"""Merge salmon quant.sf TPM columns into one matrix (transcripts x samples)."""
import pandas as pd

cols = {}
for sample, path in zip(snakemake.params.samples, snakemake.input):
    q = pd.read_csv(path, sep="\t", index_col=0)
    cols[sample] = q["TPM"]
pd.DataFrame(cols).to_csv(snakemake.output[0], sep="\t")
