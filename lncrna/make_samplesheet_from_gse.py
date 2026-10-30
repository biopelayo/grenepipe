#!/usr/bin/env python3
"""
Build a workflow sample sheet from one or more GEO series (GSE) via ENA.

GEO series -> SRA study (SRP) is resolved through NCBI E-utilities; run-level FASTQ links and
sample attributes come from the ENA Portal API (ENA mirrors every SRA run).

Usage:
    ./make_samplesheet_from_gse.py --gse GSE224760 GSE226691 --out config/samples.csv \
        --email you@example.org [--tissue-from "tissue"] [--condition-from "treatment"]

Columns written: sample,fastq_1,fastq_2,accession,accession_id,tissue,condition,study,strandedness
`accession` and `accession_id` are extracted from BioSample attributes (ecotype / accession /
genotype / 1001g id); review the result by hand before running the workflow.
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
              "library_layout", "library_strategy", "fastq_ftp", "sample_title", "sample_alias",
              "read_count", "instrument_model"]
    r = requests.get(ENA, params={"result": "read_run", "query": f'secondary_study_accession="{srp}"',
                                  "fields": ",".join(fields), "format": "tsv", "limit": 0}, timeout=300)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str) if r.text.strip() else pd.DataFrame(columns=fields)


def ena_sample_xml_attrs(sample_acc: str) -> dict:
    """Full BioSample attributes via ENA browser XML (richer than the portal fields)."""
    r = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample_acc}", timeout=60)
    if r.status_code != 200:
        return {}
    attrs = {}
    for m in re.finditer(r"<SAMPLE_ATTRIBUTE>\s*<TAG>(.*?)</TAG>\s*<VALUE>(.*?)</VALUE>", r.text, re.S):
        attrs[m.group(1).strip().lower()] = m.group(2).strip()
    return attrs


def pick(attrs: dict, keys) -> str:
    for k in keys:
        for ak, av in attrs.items():
            if k in ak and av and av.lower() not in ("na", "n/a", "none", "missing"):
                return av
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gse", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--email", required=True)
    ap.add_argument("--tissue-from", default="tissue")
    ap.add_argument("--condition-from", default="treatment")
    ap.add_argument("--rna-only", action="store_true", default=True)
    args = ap.parse_args()

    rows = []
    for gse in args.gse:
        srps = gse_to_srp(gse, args.email)
        print(f"[{gse}] SRA studies: {srps}", file=sys.stderr)
        for srp in srps:
            runs = ena_runs_for_study(srp)
            if args.rna_only:
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
                    continue
                acc = pick(a, ["ecotype", "accession", "genotype", "cultivar", "strain"])
                acc_id = pick(a, ["1001", "accession id", "cs number", "ecotype id"])
                if not acc_id and re.fullmatch(r"\d{3,5}", acc or ""):
                    acc_id = acc
                sample = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{acc or sacc}_{r['run_accession']}")
                rows.append({
                    "sample": sample,
                    "fastq_1": ftp[0],
                    "fastq_2": ftp[1] if len(ftp) > 1 else "",
                    "accession": acc,
                    "accession_id": acc_id,
                    "tissue": pick(a, [args.tissue_from, "tissue", "dev", "stage", "organ"]),
                    "condition": pick(a, [args.condition_from, "treatment", "condition", "stress"]) or "control",
                    "study": gse,
                    "strandedness": "auto",
                    "biosample": sacc,
                    "sample_title": r.get("sample_title", ""),
                })
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(f"wrote {len(df)} samples to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
