#!/usr/bin/env python3
"""
Simulate a small dataset with known structure to test the downstream analysis layer
(normalize -> stratify -> signature -> de_gxe -> eqtl -> report) without sequencing data.

Writes, under <outdir>:
  matrices/{counts,tpm}_{lncrna,coding,all}.tsv, matrices/metadata.tsv
  lncrna/lncrna_classification.tsv, lncrna/lncrna.gtf
  genotypes/genotypes_012.tsv, genotypes/variant_positions.tsv
  config/samples.csv, external_labels.csv

Ground truth: 3 accession groups driving a subset of lncRNAs; a drought condition with
genotype x condition interaction in a subset of features; cis-eQTLs for a subset of lncRNAs.
EVERYTHING HERE IS SYNTHETIC. It validates code paths and output formats only.
"""
import argparse
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-accessions", type=int, default=48)
    ap.add_argument("--n-lnc", type=int, default=1500)
    ap.add_argument("--n-coding", type=int, default=4000)
    ap.add_argument("--n-snps", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = a.outdir
    for d in ["matrices", "lncrna", "genotypes", "config"]:
        os.makedirs(os.path.join(out, d), exist_ok=True)

    nA = a.n_accessions
    acc_ids = [str(7000 + i) for i in range(nA)]
    acc_names = [f"Acc{i:03d}" for i in range(nA)]
    group = rng.integers(0, 3, nA)                       # 3 hidden accession groups
    conditions = ["control", "drought"]
    samples = []
    for i in range(nA):
        for c in conditions:
            for r in range(2):                             # 2 replicates
                samples.append({"sample": f"{acc_names[i]}_{c}_r{r + 1}", "accession": acc_names[i],
                                "accession_id": acc_ids[i], "tissue": "rosette", "condition": c,
                                "study": "SIM1" if i % 2 == 0 else "SIM2", "strandedness": "reverse",
                                "fastq_1": "na", "fastq_2": ""})
    meta = pd.DataFrame(samples)
    nS = len(meta)

    # ---- genotypes: 5 chromosomes, SNPs with group-correlated allele frequencies ---------------
    chrom_len = {"1": 30_000_000, "2": 19_000_000, "3": 23_000_000, "4": 18_000_000, "5": 27_000_000}
    chroms = rng.choice(list(chrom_len), a.n_snps, p=np.array(list(chrom_len.values())) / sum(chrom_len.values()))
    pos = np.array([rng.integers(1, chrom_len[c]) for c in chroms])
    order = np.lexsort((pos, chroms))
    chroms, pos = chroms[order], pos[order]
    snp_ids = [f"{c}_{p}" for c, p in zip(chroms, pos)]
    base_af = rng.beta(0.8, 0.8, a.n_snps).clip(0.05, 0.95)
    G = np.zeros((a.n_snps, nA), dtype=int)
    for g in range(3):
        af = (base_af + rng.normal(0, 0.15, a.n_snps)).clip(0.02, 0.98)
        idx = np.where(group == g)[0]
        G[:, idx] = (rng.random((a.n_snps, len(idx))) < af[:, None]).astype(int) * 2  # inbred lines: 0/2
    gt = pd.DataFrame(G, index=snp_ids, columns=acc_ids)
    gt.index.name = "snp"
    gt.to_csv(os.path.join(out, "genotypes/genotypes_012.tsv"), sep="\t")
    pd.DataFrame({"snp": snp_ids, "chr": chroms, "pos": pos}).to_csv(
        os.path.join(out, "genotypes/variant_positions.tsv"), sep="\t", index=False)

    # ---- features ------------------------------------------------------------------------------
    def gene_coords(n, prefix):
        c = rng.choice(list(chrom_len), n)
        s = np.array([rng.integers(1000, chrom_len[x] - 10000) for x in c])
        L = rng.integers(300, 5000, n)
        return pd.DataFrame({"gene_id": [f"{prefix}{i:05d}" for i in range(n)], "chrom": c, "start": s,
                             "end": s + L, "strand": rng.choice(["+", "-"], n)})
    lnc = gene_coords(a.n_lnc, "LNC")
    cod = gene_coords(a.n_coding, "AT9G")

    # accession-level log-expression
    def simulate(features, n_group_driven, n_gxe, n_eqtl, sd_group=1.2):
        n = len(features)
        mu = rng.normal(4, 1.5, n)
        acc_eff = rng.normal(0, 0.3, (n, nA))
        driven = rng.choice(n, n_group_driven, replace=False)
        for f in driven:
            shift = rng.normal(0, sd_group, 3)
            acc_eff[f] += shift[group]
        cond_eff = np.zeros((n, nA))
        de = rng.choice(n, n // 10, replace=False)
        cond_eff[de] += rng.normal(0, 1.0, (len(de), 1))
        gxe = rng.choice(n, n_gxe, replace=False)
        for f in gxe:
            cond_eff[f] += rng.normal(0, 1.0, 3)[group]
        eqtl_snp = {}
        for f in rng.choice(n, n_eqtl, replace=False):
            fc = features.loc[f, "chrom"]; fs = features.loc[f, "start"]
            cands = np.where((chroms == fc) & (np.abs(pos - fs) < 50_000))[0]
            if len(cands) == 0:
                continue
            s = rng.choice(cands)
            acc_eff[f] += 0.6 * (G[s] / 2.0 - G[s].mean() / 2.0) * rng.choice([-1, 1])
            eqtl_snp[features.loc[f, "gene_id"]] = snp_ids[s]
        logexp = np.zeros((n, nS))
        for j, row in meta.iterrows():
            i = acc_ids.index(row["accession_id"])
            e = mu + acc_eff[:, i] + (cond_eff[:, i] if row["condition"] == "drought" else 0)
            e += 0.15 * (row["study"] == "SIM2") + rng.normal(0, 0.25, n)
            logexp[:, j] = e
        counts = rng.poisson(np.exp(logexp)).astype(int)
        tpm = counts / counts.sum(axis=0, keepdims=True) * 1e6
        return counts, tpm, {"group_driven": features.loc[driven, "gene_id"].tolist(),
                             "gxe": features.loc[gxe, "gene_id"].tolist(), "eqtl": eqtl_snp}

    c_lnc, t_lnc, truth_lnc = simulate(lnc, 300, 80, 150)
    c_cod, t_cod, truth_cod = simulate(cod, 200, 150, 0, sd_group=0.5)

    def write(m, feats, path):
        df = pd.DataFrame(m, index=feats["gene_id"], columns=meta["sample"])
        df.index.name = "feature"
        df.to_csv(path, sep="\t")
    write(c_lnc, lnc, f"{out}/matrices/counts_lncrna.tsv"); write(t_lnc, lnc, f"{out}/matrices/tpm_lncrna.tsv")
    write(c_cod, cod, f"{out}/matrices/counts_coding.tsv"); write(t_cod, cod, f"{out}/matrices/tpm_coding.tsv")
    write(np.vstack([c_lnc, c_cod]), pd.concat([lnc, cod]), f"{out}/matrices/counts_all.tsv")
    write(np.vstack([t_lnc, t_cod]), pd.concat([lnc, cod]), f"{out}/matrices/tpm_all.tsv")
    meta[["sample", "accession", "accession_id", "tissue", "condition", "study"]].to_csv(
        f"{out}/matrices/metadata.tsv", sep="\t", index=False)
    meta.to_csv(f"{out}/config/samples.csv", index=False)

    cls = lnc.copy()
    cls["transcript_id"] = cls["gene_id"] + ".1"
    cls["length"] = cls["end"] - cls["start"]
    cls["n_exons"] = rng.integers(1, 5, len(cls))
    cls["class_code"] = rng.choice(["u", "x", "i"], len(cls), p=[0.6, 0.3, 0.1])
    cls["lncrna_type"] = cls["class_code"].map({"u": "lincRNA", "x": "antisense", "i": "intronic"})
    cls["te_like"] = rng.random(len(cls)) < 0.3
    cls["source"] = "novel_assembly"
    cls["pass"] = True
    cls.to_csv(f"{out}/lncrna/lncrna_classification.tsv", sep="\t", index=False)
    with open(f"{out}/lncrna/lncrna.gtf", "w") as f:
        for _, r in lnc.iterrows():
            attrs = f'gene_id "{r.gene_id}"; transcript_id "{r.gene_id}.1";'
            f.write(f"{r.chrom}\tsim\ttranscript\t{r.start}\t{r.end}\t.\t{r.strand}\t.\t{attrs}\n")
            f.write(f"{r.chrom}\tsim\texon\t{r.start}\t{r.end}\t.\t{r.strand}\t.\t{attrs}\n")
    pd.DataFrame({"accession_id": acc_ids, "label": [f"group{g}" for g in group]}).to_csv(
        f"{out}/external_labels.csv", index=False)
    pd.DataFrame({"gene_id": list(truth_lnc["eqtl"]), "true_snp": list(truth_lnc["eqtl"].values())}).to_csv(
        f"{out}/truth_eqtl.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": truth_lnc["gxe"]}).to_csv(f"{out}/truth_gxe_lnc.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": truth_lnc["group_driven"]}).to_csv(f"{out}/truth_group_lnc.tsv", sep="\t", index=False)
    print(f"simulated {nA} accessions, {nS} samples, {a.n_lnc} lncRNAs, {a.n_coding} coding, {a.n_snps} SNPs -> {out}")


if __name__ == "__main__":
    main()
