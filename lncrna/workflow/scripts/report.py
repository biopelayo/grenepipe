"""Assemble a Markdown summary of the main results."""
import pandas as pd

I = snakemake.input
out = []
out.append("# lncRNA population analysis: summary\n")

cl = pd.read_csv(I.classification, sep="\t")
out.append("## lncRNA annotation\n")
out.append(f"- Final lncRNA transcripts: **{len(cl)}** "
           f"({(cl['source'] == 'novel_assembly').sum()} novel from assembly, "
           f"{(cl['source'] != 'novel_assembly').sum()} from external annotations)")
if "lncrna_type" in cl:
    out.append("- By type: " + ", ".join(f"{k}: {v}" for k, v in cl["lncrna_type"].value_counts().items()))
if "te_like" in cl:
    out.append(f"- TE-like (>= 50% exonic overlap with TAIR10 TEs): {int(cl['te_like'].fillna(False).sum())}")

fs = pd.read_csv(I.filter_stats, sep="\t")
out.append("\n## Expression matrices\n")
out.append(fs.to_markdown(index=False))

cs = pd.read_csv(I.consensus, sep="\t")
out.append("\n## Consensus clustering (lncRNA expression)\n")
out.append(cs.to_markdown(index=False))
cmp = pd.read_csv(I.compare, sep="\t")
out.append("\n### lncRNA vs coding clustering agreement\n")
out.append(cmp.to_markdown(index=False))

sig = pd.read_csv(I.signature, sep="\t")
perf = pd.read_csv(I.performance, sep="\t")
out.append("\n## lncRNA signature (elastic-net, nested CV)\n")
out.append(f"- Signature size: **{len(sig)}** lncRNAs")
out.append(perf.to_markdown(index=False))
out.append("\nTop signature members:\n")
out.append(sig.head(25).to_markdown(index=False))

if I.get("gxe"):
    g = pd.read_csv(I.gxe, sep="\t")
    out.append("\n## Genotype x condition (DESeq2 LRT)\n")
    out.append(f"- Features with significant G x E (padj < 0.05): **{int((g['padj'] < 0.05).sum())}** "
               f"of which lncRNAs: {int(((g['padj'] < 0.05) & (g['is_lncrna'] == True)).sum())}")

if I.get("egenes"):
    e = pd.read_csv(I.egenes, sep="\t")
    out.append("\n## cis-eQTL (MatrixEQTL)\n")
    out.append(f"- lncRNA eGenes at FDR {snakemake.config['analysis']['eqtl']['fdr']}: **{len(e)}**")

with open(snakemake.output[0], "w") as f:
    f.write("\n".join(out) + "\n")
