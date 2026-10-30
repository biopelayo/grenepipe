import json, yaml
cfg = yaml.safe_load(open("/home/user/grenepipe/lncrna/config/config.yaml"))
A = cfg["analysis"]
specs = {
 "normalize": dict(input=dict(counts_lnc="matrices/counts_lncrna.tsv", tpm_lnc="matrices/tpm_lncrna.tsv", counts_cod="matrices/counts_coding.tsv", tpm_cod="matrices/tpm_coding.tsv", metadata="matrices/metadata.tsv"),
   output=dict(vst_lnc="results/normalized/vst_lncrna.tsv", vst_cod="results/normalized/vst_coding.tsv", vst_lnc_bc="results/normalized/vst_lncrna_batchcorrected.tsv", vst_cod_bc="results/normalized/vst_coding_batchcorrected.tsv", pca="results/normalized/pca.pdf", stats="results/normalized/filter_stats.tsv"),
   params=dict(expressed_tpm=A["expressed_tpm"], expressed_frac=A["expressed_frac"]), log=["logs/normalize.log"], threads=4),
 "stratify": dict(input=dict(vst_lnc="results/normalized/vst_lncrna_batchcorrected.tsv", vst_cod="results/normalized/vst_coding_batchcorrected.tsv", metadata="matrices/metadata.tsv", labels="external_labels.csv"),
   output=dict(summary="results/stratification/consensus_summary.tsv", clusters="results/stratification/clusters.tsv", pdf="results/stratification/consensus_plots.pdf", compare="results/stratification/lncrna_vs_coding.tsv"),
   params=dict(top_variable=A["top_variable"], consensus=dict(A["consensus"], reps=300)), log=["logs/stratify.log"], threads=4),
 "signature": dict(input=dict(vst_lnc="results/normalized/vst_lncrna_batchcorrected.tsv", metadata="matrices/metadata.tsv", clusters="results/stratification/clusters.tsv", classification="lncrna/lncrna_classification.tsv"),
   output=dict(signature="results/signature/signature_lncrnas.tsv", performance="results/signature/nested_cv_performance.tsv", pdf="results/signature/signature_heatmap.pdf"),
   params=dict(sig=A["signature"], seed=A["consensus"]["seed"]), log=["logs/signature.log"], threads=4),
 "de_gxe": dict(input=dict(counts="matrices/counts_all.tsv", metadata="matrices/metadata.tsv", classification="lncrna/lncrna_classification.tsv"),
   output=dict(gxe="results/de/gxe_results.tsv", condition="results/de/condition_results.tsv", pdf="results/de/de_plots.pdf"), params={}, log=["logs/de_gxe.log"], threads=4),
 "eqtl": dict(input=dict(vst="results/normalized/vst_lncrna.tsv", metadata="matrices/metadata.tsv", gt="genotypes/genotypes_012.tsv", pos="genotypes/variant_positions.tsv", gtf="lncrna/lncrna.gtf"),
   output=dict(cis="results/eqtl/cis_eqtl.tsv", egenes="results/eqtl/egenes.tsv", pdf="results/eqtl/eqtl_plots.pdf"), params=dict(eqtl=A["eqtl"]), log=["logs/eqtl.log"], threads=4),
}
for k, v in specs.items(): json.dump(v, open(f"spec_{k}.json", "w"), indent=1)
print("specs written")
