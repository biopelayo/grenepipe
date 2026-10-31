#!/usr/bin/env python3
"""
Systematic discovery of public RNA-seq / DNA-seq datasets for lncRNA population studies.

Queries ENA (primary), NCBI SRA and GEO (secondary) for a given organism, pulls run- and
sample-level metadata, scores each study for suitability (stranded, rRNA-depleted / total
RNA, many accessions / genotypes, paired-end, read length), and writes ranked TSVs.

Usage:
    ./find_datasets.py --organism arabidopsis --out results/arabidopsis
    ./find_datasets.py --organism marchantia  --out results/marchantia
    ./find_datasets.py --organism chlamy      --out results/chlamy --email you@example.org

    # Only re-score an existing ENA dump (offline):
    ./find_datasets.py --organism arabidopsis --out results/arabidopsis --offline

Outputs (in --out):
    ena_runs.tsv         raw run-level ENA dump
    sra_runinfo.csv      raw SRA RunInfo (if NCBI reachable)
    studies_ranked.tsv   one row per study/BioProject, with score and evidence columns
    samples_<study>.tsv  per-study sample sheets for the top studies (ready for nf-core/rnaseq)

No third-party dependencies beyond pandas and requests.
"""

import argparse
import csv
import io
import os
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

# =================================================================================================
#     Organisms
# =================================================================================================

ORGANISMS = {
    "arabidopsis": {
        "taxid": 3702,
        "name": "Arabidopsis thaliana",
        # BioSample attribute names that identify a natural accession / genotype
        "genotype_attrs": ["ecotype", "accession", "genotype", "cultivar", "strain", "1001g_id"],
        # Known reference accession, used to flag Col-0-only studies as low value for pop. work
        "reference_genotype": r"col[- ]?0|columbia",
    },
    "marchantia": {
        "taxid": 3197,
        "name": "Marchantia polymorpha",
        "genotype_attrs": ["strain", "accession", "genotype", "isolate", "ecotype"],
        "reference_genotype": r"tak[- ]?[12]",
    },
    "chlamy": {
        "taxid": 3055,
        "name": "Chlamydomonas reinhardtii",
        "genotype_attrs": ["strain", "isolate", "genotype", "accession"],
        "reference_genotype": r"cc[- ]?(4532|503|1690|125|124|400)",
    },
}

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/search"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ENA_RUN_FIELDS = [
    "run_accession", "sample_accession", "secondary_sample_accession", "experiment_accession",
    "study_accession", "secondary_study_accession", "study_title", "experiment_title",
    "library_strategy", "library_source", "library_selection", "library_layout",
    "library_construction_protocol", "instrument_platform", "instrument_model",
    "read_count", "base_count", "nominal_length", "first_public", "last_updated",
    "sample_title", "sample_alias", "description", "scientific_name", "tax_id",
    "fastq_ftp", "fastq_md5", "sra_ftp", "submitted_ftp", "sample_description",
]

# =================================================================================================
#     Heuristics
# =================================================================================================

# Each regex is applied to the concatenation of all free-text metadata of a run.
# Positive evidence for a feature we want, negative evidence against it.
PATTERNS = {
    "ribo_depleted": re.compile(
        r"ribo[- ]?(zero|minus|depl)|rrna[- ]?depl|ribosomal rna (removal|depletion)|"
        r"total[- ]rna[- ]?seq|rrna[- ]?remov|ribo-off|ribopool|qiaseq fastselect|"
        r"nebnext rrna depletion",
        re.I,
    ),
    "polya": re.compile(r"poly[- ]?\(?a\)?[- ]?(select|enrich|tail|\+|plus)|oligo[- ]?dt|mrna[- ]?seq", re.I),
    "stranded": re.compile(
        r"strand[- ]?specific|stranded|dutp|truseq stranded|kapa stranded|"
        r"nebnext ultra (ii )?directional|smarter stranded|first[- ]strand|reverse[- ]?stranded|fr-firststrand",
        re.I,
    ),
    "small_rna": re.compile(r"small[- ]?rna|mirna|srna[- ]?seq", re.I),
    "single_cell": re.compile(r"single[- ]?(cell|nucle)|10x|scrna|snrna|drop[- ]?seq", re.I),
    "degradome": re.compile(r"degradome|pare[- ]seq|cage|ribo[- ]?seq|ribosome profiling|clip[- ]?seq|rip[- ]?seq|3'?[- ]?seq|pas[- ]?seq", re.I),
    "stress": re.compile(
        r"drought|cold|heat|salt|salinity|osmotic|aba|pathogen|infection|nitrogen|phosphate|"
        r"starvation|hypoxia|flood|uv|light|dark|circadian|temperature|stress",
        re.I,
    ),
    "tissue": re.compile(
        r"root|leaf|leaves|rosette|seedling|flower|inflorescence|seed|pollen|shoot|thallus|"
        r"gemma|gametophyte|sporophyte|whole plant",
        re.I,
    ),
}


def text_blob(row: pd.Series) -> str:
    cols = [
        "study_title", "experiment_title", "library_source", "library_selection",
        "library_construction_protocol", "sample_title", "sample_alias", "description",
        "sample_description",
    ]
    return " | ".join(str(row.get(c, "")) for c in cols if pd.notna(row.get(c, "")))


def classify_run(row: pd.Series) -> dict:
    blob = text_blob(row)
    out = {k: bool(p.search(blob)) for k, p in PATTERNS.items()}
    # ENA library_source TRANSCRIPTOMIC with selection "RANDOM"/"other" and no polyA keyword is
    # weak evidence of total RNA; "PolyA"/"cDNA" selection is evidence of polyA.
    sel = str(row.get("library_selection", "")).lower()
    if sel in ("polya",):
        out["polya"] = True
    if sel in ("random", "inverse rrna", "inverse rrna selection"):
        out["ribo_depleted"] = out["ribo_depleted"] or sel.startswith("inverse")
    out["paired"] = str(row.get("library_layout", "")).upper() == "PAIRED"
    try:
        out["read_len"] = int(float(row.get("base_count", 0)) / max(float(row.get("read_count", 1)), 1))
    except (TypeError, ValueError):
        out["read_len"] = 0
    try:
        out["reads_M"] = float(row.get("read_count", 0)) / 1e6
    except (TypeError, ValueError):
        out["reads_M"] = 0.0
    return out


# =================================================================================================
#     ENA
# =================================================================================================


def ena_search(query: str, result: str, fields: list, limit: int = 0, retries: int = 3) -> pd.DataFrame:
    params = {
        "result": result,
        "query": query,
        "fields": ",".join(fields),
        "format": "tsv",
        "limit": limit,
    }
    for attempt in range(retries):
        try:
            r = requests.get(ENA_PORTAL, params=params, timeout=300)
            r.raise_for_status()
            if not r.text.strip():
                return pd.DataFrame(columns=fields)
            return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
        except requests.RequestException as e:
            wait = 2 ** (attempt + 1)
            print(f"[ena] {e}; retry in {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError("ENA portal unreachable")


def ena_runs(taxid: int, strategies=("RNA-Seq",), extra_query: str = "") -> pd.DataFrame:
    strat = " OR ".join(f'library_strategy="{s}"' for s in strategies)
    q = f"tax_eq({taxid}) AND ({strat})"
    if extra_query:
        q += f" AND ({extra_query})"
    df = ena_search(q, "read_run", ENA_RUN_FIELDS)
    print(f"[ena] {len(df)} runs for taxid {taxid} strategies={strategies}", file=sys.stderr)
    return df


def ena_sample_attributes(sample_accessions: list, chunk: int = 200) -> pd.DataFrame:
    """Fetch BioSample attributes (ecotype/strain/etc.) via the ENA sample result."""
    fields = ["sample_accession", "secondary_sample_accession", "sample_title", "sample_alias",
              "description", "cultivar", "ecotype", "strain", "isolate", "cell_type", "tissue_type",
              "dev_stage", "environmental_sample", "collection_date", "country", "location"]
    frames = []
    accs = [a for a in dict.fromkeys(sample_accessions) if isinstance(a, str) and a]
    for i in range(0, len(accs), chunk):
        sub = accs[i:i + chunk]
        q = " OR ".join(f'sample_accession="{a}"' for a in sub)
        try:
            frames.append(ena_search(q, "sample", fields))
        except RuntimeError:
            break
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=fields)


# =================================================================================================
#     NCBI SRA / GEO (secondary; captures studies with richer GEO-side metadata)
# =================================================================================================


def sra_runinfo(organism_name: str, email: str, strategy: str = "RNA-Seq", batch: int = 5000) -> pd.DataFrame:
    term = f'"{organism_name}"[Organism] AND "{strategy}"[Strategy]'
    params = {"db": "sra", "term": term, "usehistory": "y", "retmax": 0, "email": email, "tool": "lncrna-finder"}
    r = requests.get(f"{NCBI_EUTILS}/esearch.fcgi", params=params, timeout=120)
    r.raise_for_status()
    webenv = re.search(r"<WebEnv>(\S+)</WebEnv>", r.text).group(1)
    qkey = re.search(r"<QueryKey>(\d+)</QueryKey>", r.text).group(1)
    count = int(re.search(r"<Count>(\d+)</Count>", r.text).group(1))
    print(f"[sra] {count} runs for {term}", file=sys.stderr)
    frames = []
    for start in range(0, count, batch):
        p = {"db": "sra", "query_key": qkey, "WebEnv": webenv, "rettype": "runinfo", "retmode": "text",
             "retstart": start, "retmax": batch, "email": email, "tool": "lncrna-finder"}
        rr = requests.get(f"{NCBI_EUTILS}/efetch.fcgi", params=p, timeout=600)
        rr.raise_for_status()
        frames.append(pd.read_csv(io.StringIO(rr.text), dtype=str))
        time.sleep(0.4)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# =================================================================================================
#     Scoring
# =================================================================================================


def score_studies(runs: pd.DataFrame, samples: pd.DataFrame, org: dict) -> pd.DataFrame:
    if runs.empty:
        return pd.DataFrame()
    feats = runs.apply(classify_run, axis=1, result_type="expand")
    df = pd.concat([runs, feats], axis=1)

    # Attach genotype labels from sample attributes
    if not samples.empty:
        samples = samples.drop_duplicates("sample_accession").set_index("sample_accession")
        geno_cols = [c for c in org["genotype_attrs"] if c in samples.columns]

        def genotype(acc):
            if acc not in samples.index:
                return ""
            row = samples.loc[acc]
            for c in geno_cols:
                v = row.get(c, "")
                if isinstance(v, str) and v.strip():
                    return v.strip()
            return ""

        df["genotype"] = df["sample_accession"].map(genotype)
    else:
        df["genotype"] = ""
    # Fallback: try to pull a genotype-like token from the sample title/alias
    ref_pat = re.compile(org["reference_genotype"], re.I)

    g = df.groupby("study_accession", dropna=False)
    rows = []
    for study, sub in g:
        genos = set(x for x in sub["genotype"] if x)
        n_geno = len(genos)
        n_nonref = len([x for x in genos if not ref_pat.search(x)])
        frac = lambda c: float(sub[c].mean()) if c in sub else 0.0
        rec = {
            "study_accession": study,
            "secondary_study_accession": sub["secondary_study_accession"].iloc[0],
            "study_title": sub["study_title"].iloc[0],
            "first_public": sub["first_public"].min(),
            "n_runs": len(sub),
            "n_samples": sub["sample_accession"].nunique(),
            "n_genotypes": n_geno,
            "n_nonref_genotypes": n_nonref,
            "frac_stranded_evidence": frac("stranded"),
            "frac_ribo_depleted_evidence": frac("ribo_depleted"),
            "frac_polya_evidence": frac("polya"),
            "frac_paired": frac("paired"),
            "median_read_len": float(sub["read_len"].median()),
            "median_reads_M": float(sub["reads_M"].median()),
            "has_stress_terms": bool(sub["stress"].any()),
            "has_tissue_terms": bool(sub["tissue"].any()),
            "exclude_small_rna": bool(sub["small_rna"].mean() > 0.5),
            "exclude_single_cell": bool(sub["single_cell"].mean() > 0.5),
            "exclude_other_assay": bool(sub["degradome"].mean() > 0.5),
            "platforms": ";".join(sorted(set(sub["instrument_model"].dropna()))),
            "example_genotypes": ";".join(sorted(genos)[:8]),
        }
        # Score: what makes a study useful for lncRNA x population work.
        s = 0.0
        s += 3.0 * rec["frac_stranded_evidence"]
        s += 3.0 * rec["frac_ribo_depleted_evidence"]
        s -= 1.5 * rec["frac_polya_evidence"] * (1 - rec["frac_ribo_depleted_evidence"])
        s += 1.0 * rec["frac_paired"]
        s += 1.0 if rec["median_read_len"] >= 75 else 0.0
        s += 0.5 if rec["median_reads_M"] >= 20 else 0.0
        s += min(rec["n_nonref_genotypes"], 50) / 10.0  # up to +5 for >=50 accessions
        s += 0.5 if rec["n_samples"] >= 24 else 0.0
        s += 0.5 if rec["has_stress_terms"] else 0.0
        if rec["exclude_small_rna"] or rec["exclude_single_cell"] or rec["exclude_other_assay"]:
            s = -10.0
        rec["score"] = round(s, 2)
        rows.append(rec)
    out = pd.DataFrame(rows).sort_values("score", ascending=False)
    return out


def write_sample_sheets(runs: pd.DataFrame, ranked: pd.DataFrame, outdir: Path, top: int = 10):
    """nf-core/rnaseq-style sample sheets (sample,fastq_1,fastq_2,strandedness) for top studies."""
    feats = runs.apply(classify_run, axis=1, result_type="expand")
    df = pd.concat([runs, feats], axis=1)
    for study in ranked.head(top)["study_accession"]:
        sub = df[df["study_accession"] == study]
        recs = []
        for _, r in sub.iterrows():
            ftp = str(r.get("fastq_ftp", "") or "")
            files = ["ftp://" + f for f in ftp.split(";") if f]
            if not files:
                continue
            recs.append({
                "sample": r["run_accession"],
                "fastq_1": files[0],
                "fastq_2": files[1] if len(files) > 1 else "",
                "strandedness": "auto",
                "biosample": r["sample_accession"],
                "genotype": r.get("genotype", ""),
                "sample_title": r.get("sample_title", ""),
            })
        if recs:
            pd.DataFrame(recs).to_csv(outdir / f"samples_{study}.csv", index=False)


# =================================================================================================
#     Main
# =================================================================================================


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", choices=ORGANISMS.keys(), required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""), help="for NCBI E-utilities")
    ap.add_argument("--offline", action="store_true", help="re-score existing dumps without network")
    ap.add_argument("--with-dna", action="store_true", help="also dump WGS runs (for genotype matching)")
    ap.add_argument("--skip-sra", action="store_true")
    args = ap.parse_args()

    org = ORGANISMS[args.organism]
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    runs_path = outdir / "ena_runs.tsv"
    samples_path = outdir / "ena_samples.tsv"

    if args.offline:
        runs = pd.read_csv(runs_path, sep="\t", dtype=str)
        samples = pd.read_csv(samples_path, sep="\t", dtype=str) if samples_path.exists() else pd.DataFrame()
    else:
        runs = ena_runs(org["taxid"])
        runs.to_csv(runs_path, sep="\t", index=False)
        samples = ena_sample_attributes(runs["sample_accession"].tolist())
        samples.to_csv(samples_path, sep="\t", index=False)
        if args.with_dna:
            dna = ena_runs(org["taxid"], strategies=("WGS", "WGA"))
            dna.to_csv(outdir / "ena_wgs_runs.tsv", sep="\t", index=False)
        if not args.skip_sra and args.email:
            try:
                sra = sra_runinfo(org["name"], args.email)
                sra.to_csv(outdir / "sra_runinfo.csv", index=False)
            except Exception as e:  # NCBI is flaky; ENA already covers the same runs
                print(f"[sra] skipped: {e}", file=sys.stderr)

    ranked = score_studies(runs, samples, org)
    ranked.to_csv(outdir / "studies_ranked.tsv", sep="\t", index=False)
    runs_with_geno = runs.copy()
    if not samples.empty:
        gcols = [c for c in org["genotype_attrs"] if c in samples.columns]
        s = samples.drop_duplicates("sample_accession").set_index("sample_accession")
        runs_with_geno["genotype"] = runs_with_geno["sample_accession"].map(
            lambda a: next((s.loc[a][c] for c in gcols if a in s.index and isinstance(s.loc[a].get(c), str) and s.loc[a][c].strip()), "")
        )
    write_sample_sheets(runs_with_geno, ranked, outdir)

    print(f"\nTop studies for {org['name']} (score, n_genotypes, stranded, ribo-depleted, title):")
    for _, r in ranked.head(25).iterrows():
        print(f"  {r['score']:5.1f}  {r['n_nonref_genotypes']:4d}  {r['frac_stranded_evidence']:.2f}  "
              f"{r['frac_ribo_depleted_evidence']:.2f}  {r['study_accession']}  {str(r['study_title'])[:90]}")


if __name__ == "__main__":
    main()
