"""Select lncRNA candidate transcripts from the gffcompare-annotated merged assembly."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from gtf_utils import transcripts_from_gtf, write_fasta_from_exons  # noqa: E402

sys.stderr = open(snakemake.log[0], "w")
P = snakemake.params

tmap = pd.read_csv(snakemake.input.tmap, sep="\t", dtype=str)
tmap = tmap.set_index("qry_id")
tx = transcripts_from_gtf(snakemake.input.gtf)

keep = {}
rows = []
for tid, t in tx.items():
    if t["chrom"] not in P.chromosomes:
        continue
    cc = t["attrs"].get("class_code") or (tmap.loc[tid, "class_code"] if tid in tmap.index else "")
    if cc not in P.class_codes:
        continue
    if t["length"] < P.min_length:
        continue
    if len(t["exons"]) == 1 and t["length"] < P.min_length_monoexonic:
        continue
    if cc == "u" and len(t["exons"]) < P.min_exons_intergenic:
        continue
    keep[tid] = t
    rows.append({"transcript_id": tid, "gene_id": t["gene_id"], "chrom": t["chrom"], "start": t["start"],
                 "end": t["end"], "strand": t["strand"], "length": t["length"], "n_exons": len(t["exons"]),
                 "class_code": cc, "ref_gene": t["attrs"].get("cmp_ref_gene", t["attrs"].get("ref_gene_id", ""))})

with open(snakemake.output.gtf, "w") as out:
    for t in keep.values():
        out.write("\n".join(t["lines"]) + "\n")
pd.DataFrame(rows).to_csv(snakemake.output.tsv, sep="\t", index=False)
write_fasta_from_exons(keep, snakemake.input.fa, snakemake.output.fa)
print(f"{len(tx)} assembled transcripts, {len(keep)} candidates "
      f"(class codes {P.class_codes}, length >= {P.min_length}, mono-exonic >= {P.min_length_monoexonic})",
      file=sys.stderr)
