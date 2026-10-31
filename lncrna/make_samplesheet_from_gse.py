#!/usr/bin/env python3
"""
Build a workflow sample sheet from one or more GEO series (GSE) via NCBI + ENA.

GEO series -> SRA study (SRP) is resolved through NCBI E-utilities; run-level FASTQ links come
from the ENA Portal API (ENA mirrors every SRA run); sample attributes from the ENA BioSample XML.
Dataset-specific knowledge (library strandedness, how accession ids and batches are encoded in
sample titles) lives in PRESETS below, derived from the papers' methods sections.

Usage:
    ./make_samplesheet_from_gse.py --gse GSE224760 --out config/samples.GSE224760.csv --email you@x.org
    ./make_samplesheet_from_gse.py --gse GSE226691 --out config/samples.GSE226691.csv --email you@x.org
    ./make_samplesheet_from_gse.py --gse GSE80744  --out config/samples.GSE80744.csv  --email you@x.org
    # subset for a pilot (first N accessions, all their samples):
    ./make_samplesheet_from_gse.py --gse GSE224760 --out config/samples.pilot.csv --email you@x.org --max-accessions 2

Columns written: sample,fastq_1,fastq_2,accession,accession_id,tissue,condition,study,batch,strandedness,
                 biosample,run,sample_title
ALWAYS review the CSV by hand before running the workflow.
"""

import argparse
import io
import re
import sys
import time

import pandas as pd
import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA = "https://www.ebi.ac.uk/ena/portal/api/search"

# -------------------------------------------------------------------------------------------------
# Dataset presets (from the papers' methods; see docs/02_repositorio_y_viabilidad_local.md)
# -------------------------------------------------------------------------------------------------
PRESETS = {
    # Kornienko et al. 2024, Plant Cell: TruSeq Stranded mRNA (dUTP) -> reverse; PE 125; title acc.<id>.rep<n>
    "GSE224760": dict(strandedness="reverse", tissue="rosette14", condition="control",
                      title_re=r"acc\.?(?P<acc>\d{3,5})\.?rep(?P<batch>\d+)", batch_prefix="rep"),
    # Kornienko et al. 2024: seedling, 9-leaf rosette, flower, pollen; PE 150; stranded
    "GSE226691": dict(strandedness="reverse", tissue=None, condition="control",
                      title_re=r"(?P<acc>\d{3,5})", batch_prefix=""),
    # Kawakatsu et al. 2016, Cell (1001 Epigenomes): TruSeq RNA, NOT stranded, SE 100; title/ecotype has 1001G id
    "GSE80744": dict(strandedness="unstranded", tissue="rosette", condition="control",
                     title_re=r"(?P<acc>\d{3,5})", batch_prefix=""),
}
TISSUE_WORDS = {"seedling": "seedling", "9-leaf": "rosette9", "9 leaf": "rosette9", "rosette": "rosette",
                "flower": "flower", "pollen": "pollen", "leaf": "rosette", "leaves": "rosette"}


def gse_to_srp(gse: str, email: str) -> list:
    r = requests.get(f"{EUTILS}/esearch.fcgi", params={"db": "gds", "term": f"{gse}[ACCN] AND gse[ETYP]",
                                                      "retmode": "json", "email": email}, timeout=60)
    r.raise_for_status()
    ids = r.json()["esearchresult"]["idlist"]
    srps = set()
    for uid in ids:
        s = requests.get(f"{EUTILS}/esummary.fcgi", params={"db": "gds", "id": uid, "retmode": "json",
                                                           "email": email}, timeout=60)
        s.raise_for_status()
        doc = s.json()["result"][uid]
        for rel in doc.get("extrelations", []):
            if rel.get("relationtype") == "SRA":
                srps.add(rel["targetobject"])
        time.sleep(0.34)
    return sorted(srps)


def ena_runs_for_study(srp: str) -> pd.DataFrame:
    fields = ["run_accession", "sample_accession", "study_accession", "secondary_study_accession",
              "library_layout", "library_strategy", "library_source", "fastq_ftp", "sample_title",
              "sample_alias", "read_count", "base_count", "instrument_model", "experiment_title"]
    r = requests.get(ENA, params={"result": "read_run", "query": f'secondary_study_accession="{srp}"',
                                  "fields": ",".join(fields), "format": "tsv", "limit": 0}, timeout=300)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str) if r.text.strip() else pd.DataFrame(columns=fields)


def ena_sample_xml_attrs(sample_acc: str) -> dict:
    r = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample_acc}", timeout=60)
    if r.status_code != 200:
        return {}
    attrs = {}
    m = re.search(r"<TITLE>(.*?)</TITLE>", r.text, re.S)
    if m:
        attrs["__title"] = m.group(1).strip()
    for m in re.finditer(r"<SAMPLE_ATTRIBUTE>\s*<TAG>(.*?)</TAG>\s*<VALUE>(.*?)</VALUE>", r.text, re.S):
        attrs[m.group(1).strip().lower()] = m.group(2).strip()
    return attrs


def pick(attrs: dict, keys) -> str:
    for k in keys:
        for ak, av in attrs.items():
            if k in ak and av and av.lower() not in ("na", "n/a", "none", "missing", "not applicable"):
                return av
    return ""


def guess_tissue(text: str, default: str) -> str:
    t = (text or "").lower()
    for w, lab in TISSUE_WORDS.items():
        if w in t:
            return lab
    return default or ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gse", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--email", required=True)
    ap.add_argument("--max-accessions", type=int, default=0, help="pilot: keep only the first N accessions")
    args = ap.parse_args()

    rows = []
    for gse in args.gse:
        preset = PRESETS.get(gse, dict(strandedness="auto", tissue=None, condition="control",
                                       title_re=r"(?P<acc>\d{3,5})", batch_prefix=""))
        srps = gse_to_srp(gse, args.email)
        print(f"[{gse}] SRA studies: {srps}", file=sys.stderr)
        for srp in srps:
            runs = ena_runs_for_study(srp)
            runs = runs[runs["library_strategy"] == "RNA-Seq"]
            print(f"[{gse}/{srp}] {len(runs)} RNA-seq runs", file=sys.stderr)
            cache = {}
            for _, r in runs.iterrows():
                sacc = r["sample_accession"]
                if sacc not in cache:
                    cache[sacc] = ena_sample_xml_attrs(sacc)
                    time.sleep(0.2)
                a = cache[sacc]
                ftp = [f"ftp://{f}" for f in str(r.get("fastq_ftp") or "").split(";") if f]
                if not ftp:
                    print(f"  ! no FASTQ on ENA for {r['run_accession']} ({sacc}); use fasterq-dump", file=sys.stderr)
                    continue
                title = a.get("__title") or str(r.get("sample_title") or "") or str(r.get("sample_alias") or "")
                m = re.search(preset["title_re"], title)
                acc_id = m.group("acc") if m else pick(a, ["1001", "accession id", "ecotype id", "cs number"])
                if not acc_id:
                    ec = pick(a, ["ecotype", "accession", "genotype", "cultivar"])
                    acc_id = ec if re.fullmatch(r"\d{3,5}", ec or "") else ""
                acc_name = pick(a, ["ecotype", "accession", "genotype", "cultivar"]) or acc_id or sacc
                batch = (preset["batch_prefix"] + m.group("batch")) if (m and "batch" in m.groupdict() and m.group("batch")) else ""
                tissue = preset["tissue"] or guess_tissue(title + " " + pick(a, ["tissue", "organ", "dev", "stage"]), "")
                rows.append({
                    "sample": re.sub(r"[^A-Za-z0-9_.-]", "_", f"{gse}_{acc_id or acc_name}_{tissue}_{batch or r['run_accession']}"),
                    "fastq_1": ftp[0],
                    "fastq_2": ftp[1] if len(ftp) > 1 else "",
                    "accession": acc_name,
                    "accession_id": acc_id,
                    "tissue": tissue,
                    "condition": pick(a, ["treatment", "condition", "stress"]) or preset["condition"],
                    "study": gse,
                    "batch": batch or gse,
                    "strandedness": preset["strandedness"],
                    "biosample": sacc,
                    "run": r["run_accession"],
                    "sample_title": title,
                    "read_count": r.get("read_count", ""),
                    "instrument": r.get("instrument_model", ""),
                })
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("no samples found")
    if args.max_accessions:
        keep = list(dict.fromkeys(df["accession_id"]))[: args.max_accessions]
        df = df[df["accession_id"].isin(keep)]
    # duplicate sample names (same accession/tissue/batch sequenced in two runs) get the run appended
    dup = df["sample"].duplicated(keep=False)
    df.loc[dup, "sample"] = df.loc[dup, "sample"] + "_" + df.loc[dup, "run"]
    df.to_csv(args.out, index=False)
    n_acc = df["accession_id"].replace("", pd.NA).nunique()
    print(f"wrote {len(df)} samples, {n_acc} accessions with 1001G id, tissues={sorted(df['tissue'].unique())}, "
          f"batches={sorted(df['batch'].unique())[:6]}... -> {args.out}", file=sys.stderr)
    if (df["accession_id"] == "").any():
        print(f"  ! {(df['accession_id'] == '').sum()} samples without accession_id: fill by hand (eQTL needs it)", file=sys.stderr)


if __name__ == "__main__":
    main()
