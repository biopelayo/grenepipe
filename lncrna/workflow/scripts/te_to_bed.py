"""TAIR10_Transposable_Elements.txt -> BED (Ensembl chromosome names 1..5)."""
import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

df = pd.read_csv(snakemake.input[0], sep="\t", dtype=str)
# Columns: Transposon_Name orientation_is_5prime Transposon_min_Start Transposon_max_End Transposon_Family Transposon_Super_Family
df.columns = [c.strip() for c in df.columns]
name = df["Transposon_Name"]
chrom = name.str.extract(r"AT(\d)TE", expand=False)
start = df["Transposon_min_Start"].astype(int) - 1
end = df["Transposon_max_End"].astype(int)
strand = df["orientation_is_5prime"].map({"true": "+", "false": "-", "True": "+", "False": "-"}).fillna(".")
bed = pd.DataFrame({
    "chrom": chrom, "start": start, "end": end, "name": name,
    "family": df["Transposon_Super_Family"].fillna("."), "strand": strand,
}).dropna(subset=["chrom"])
bed.sort_values(["chrom", "start"]).to_csv(snakemake.output[0], sep="\t", header=False, index=False)
print(f"{len(bed)} TEs written", file=sys.stderr)
