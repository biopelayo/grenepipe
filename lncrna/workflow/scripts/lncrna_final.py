"""Final lncRNA set: candidates that are non-coding (CPC2 [+FEELnc] [+no Pfam]), expressed
(TPM >= min_tpm in >= min_samples), classified by locus type and TE overlap. Also builds the
combined coding + lncRNA GTF used for quantification, merging optional external lncRNA GTFs."""
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from gtf_utils import read_gtf, transcripts_from_gtf  # noqa: E402

sys.stderr = open(snakemake.log[0], "w")
I, O, P = snakemake.input, snakemake.output, snakemake.params

cand = pd.read_csv(I.cand_tsv, sep="\t", dtype={"transcript_id": str}).set_index("transcript_id")
tx = transcripts_from_gtf(I.cand_gtf)

# ---- coding potential ---------------------------------------------------------------------------
cpc = pd.read_csv(I.cpc2, sep="\t", dtype=str)
cpc.columns = [c.lstrip("#") for c in cpc.columns]
id_col = "ID" if "ID" in cpc.columns else cpc.columns[0]
lab_col = "label" if "label" in cpc.columns else cpc.columns[-1]
noncoding = set(cpc.loc[cpc[lab_col].str.lower() == "noncoding", id_col])
cand["cpc2_noncoding"] = cand.index.isin(noncoding)

if P.use_feelnc:
    feel = {a.get("transcript_id") for _, a in read_gtf(I.feelnc)}
    cand["feelnc_noncoding"] = cand.index.isin(feel)
else:
    cand["feelnc_noncoding"] = True

if P.use_pfam:
    hits = set()
    with open(I.pfam) as f:
        for line in f:
            if line.startswith("#"):
                continue
            q = line.split()[3]
            hits.add(q.rsplit("_frame=", 1)[0])
    cand["pfam_hit"] = cand.index.isin(hits)
else:
    cand["pfam_hit"] = False

# ---- expression support -------------------------------------------------------------------------
tpm = pd.read_csv(I.tpm, sep="\t", index_col=0)
tpm = tpm.reindex(cand.index).fillna(0.0)
cand["n_samples_expressed"] = (tpm >= P.min_tpm).sum(axis=1)
cand["max_tpm"] = tpm.max(axis=1)
cand["expressed"] = cand["n_samples_expressed"] >= P.min_samples

# ---- TE overlap (bedtools) ----------------------------------------------------------------------
bed = "lncrna/candidates.exons.bed"
with open(bed, "w") as f:
    for tid, t in tx.items():
        for s, e in t["exons"]:
            f.write(f"{t['chrom']}\t{s - 1}\t{e}\t{tid}\t.\t{t['strand']}\n")
res = subprocess.run(["bedtools", "intersect", "-a", bed, "-b", I.te, "-wo"], capture_output=True, text=True, check=True)
te_ov = {}
for line in res.stdout.splitlines():
    p = line.split("\t")
    tid, fam, bp = p[3], p[10], int(p[-1])
    te_ov.setdefault(tid, {"bp": 0, "fams": set()})
    te_ov[tid]["bp"] += bp
    te_ov[tid]["fams"].add(fam)
cand["te_overlap_bp"] = [te_ov.get(t, {}).get("bp", 0) for t in cand.index]
cand["te_overlap_frac"] = cand["te_overlap_bp"] / cand["length"]
cand["te_families"] = [";".join(sorted(te_ov[t]["fams"])) if t in te_ov else "" for t in cand.index]
cand["te_like"] = cand["te_overlap_frac"] >= 0.5

# ---- classification -----------------------------------------------------------------------------
cand["lncrna_type"] = cand["class_code"].map({"u": "lincRNA", "x": "antisense", "i": "intronic"})
cand["pass"] = cand["cpc2_noncoding"] & cand["feelnc_noncoding"] & ~cand["pfam_hit"] & cand["expressed"]
final = cand[cand["pass"]].copy()
final["source"] = "novel_assembly"

# ---- write final lncRNA GTF + combined quantification GTF ----------------------------------------
final_ids = set(final.index)
with open(O.gtf, "w") as out:
    for tid in final_ids:
        out.write("\n".join(tx[tid]["lines"]) + "\n")

extra_rows = []
with open(O.combined, "w") as out:
    # reference annotation (all biotypes; coding and known lncRNAs are separated later via gene lists)
    with open(I.ref) as f:
        for line in f:
            if not line.startswith("#"):
                out.write(line)
    with open(O.gtf) as f:
        out.write(f.read())
    # optional external lncRNA GTFs (e.g. Kornienko 2024). Gene/transcript ids are prefixed to avoid clashes.
    extras = I.extra if isinstance(I.extra, (list, tuple)) else [I.extra]
    for k, path in enumerate([e for e in extras if e]):
        prefix = f"EXT{k + 1}_"
        etx = transcripts_from_gtf(path)
        for tid, t in etx.items():
            if t["chrom"] not in set(cand["chrom"]):
                continue
            for line in t["lines"]:
                p = line.split("\t")
                p[8] = p[8].replace(f'transcript_id "{tid}"', f'transcript_id "{prefix}{tid}"') \
                           .replace(f'gene_id "{t["gene_id"]}"', f'gene_id "{prefix}{t["gene_id"]}"')
                out.write("\t".join(p) + "\n")
            extra_rows.append({"transcript_id": prefix + tid, "gene_id": prefix + t["gene_id"], "chrom": t["chrom"],
                               "start": t["start"], "end": t["end"], "strand": t["strand"], "length": t["length"],
                               "n_exons": len(t["exons"]), "class_code": "", "lncrna_type": "external",
                               "source": os.path.basename(path), "pass": True})

classification = pd.concat([final.reset_index(), pd.DataFrame(extra_rows)], ignore_index=True)
classification.to_csv(O.tsv, sep="\t", index=False)

print(f"candidates: {len(cand)}", file=sys.stderr)
print(f"  CPC2 noncoding: {cand['cpc2_noncoding'].sum()}", file=sys.stderr)
print(f"  expressed (TPM>={P.min_tpm} in >={P.min_samples}): {cand['expressed'].sum()}", file=sys.stderr)
print(f"  final novel lncRNAs: {len(final)} "
      f"({(final['lncrna_type'] == 'lincRNA').sum()} lincRNA, {(final['lncrna_type'] == 'antisense').sum()} antisense, "
      f"{(final['lncrna_type'] == 'intronic').sum()} intronic; {final['te_like'].sum()} TE-like)", file=sys.stderr)
print(f"  external lncRNA transcripts added: {len(extra_rows)}", file=sys.stderr)
