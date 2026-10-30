# =================================================================================================
#     Downstream analysis: stratification, signature, G x E, cis-eQTL, report
# =================================================================================================

A = config["analysis"]


rule normalize:
    input:
        counts_lnc="matrices/counts_lncrna.tsv",
        tpm_lnc="matrices/tpm_lncrna.tsv",
        counts_cod="matrices/counts_coding.tsv",
        tpm_cod="matrices/tpm_coding.tsv",
        metadata="matrices/metadata.tsv",
    output:
        vst_lnc="results/normalized/vst_lncrna.tsv",
        vst_cod="results/normalized/vst_coding.tsv",
        vst_lnc_bc="results/normalized/vst_lncrna_batchcorrected.tsv",
        vst_cod_bc="results/normalized/vst_coding_batchcorrected.tsv",
        pca="results/normalized/pca.pdf",
        stats="results/normalized/filter_stats.tsv",
    params:
        expressed_tpm=A["expressed_tpm"],
        expressed_frac=A["expressed_frac"],
    log:
        "logs/analysis/normalize.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/normalize.R"


rule stratify:
    input:
        vst_lnc="results/normalized/vst_lncrna_batchcorrected.tsv",
        vst_cod="results/normalized/vst_coding_batchcorrected.tsv",
        metadata="matrices/metadata.tsv",
        labels=A["external_labels"] if A["external_labels"] else [],
    output:
        summary="results/stratification/consensus_summary.tsv",
        clusters="results/stratification/clusters.tsv",
        pdf="results/stratification/consensus_plots.pdf",
        compare="results/stratification/lncrna_vs_coding.tsv",
    params:
        top_variable=A["top_variable"],
        consensus=A["consensus"],
    log:
        "logs/analysis/stratify.log",
    threads: 8
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/stratify.R"


rule signature:
    input:
        vst_lnc="results/normalized/vst_lncrna_batchcorrected.tsv",
        metadata="matrices/metadata.tsv",
        clusters="results/stratification/clusters.tsv",
        classification="lncrna/lncrna_classification.tsv",
    output:
        signature="results/signature/signature_lncrnas.tsv",
        performance="results/signature/nested_cv_performance.tsv",
        pdf="results/signature/signature_heatmap.pdf",
    params:
        sig=A["signature"],
        seed=A["consensus"]["seed"],
    log:
        "logs/analysis/signature.log",
    threads: 8
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/signature.R"


rule de_gxe:
    input:
        counts="matrices/counts_all.tsv",
        metadata="matrices/metadata.tsv",
        classification="lncrna/lncrna_classification.tsv",
    output:
        gxe="results/de/gxe_results.tsv",
        condition="results/de/condition_results.tsv",
        pdf="results/de/de_plots.pdf",
    log:
        "logs/analysis/de_gxe.log",
    threads: 8
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/de_gxe.R"


rule eqtl:
    input:
        vst="results/normalized/vst_lncrna.tsv",
        metadata="matrices/metadata.tsv",
        gt="genotypes/genotypes_012.tsv",
        pos="genotypes/variant_positions.tsv",
        gtf="lncrna/lncrna.gtf",
    output:
        cis="results/eqtl/cis_eqtl.tsv",
        egenes="results/eqtl/egenes.tsv",
        pdf="results/eqtl/eqtl_plots.pdf",
    params:
        eqtl=A["eqtl"],
    log:
        "logs/analysis/eqtl.log",
    threads: 8
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/eqtl.R"


rule report:
    input:
        classification="lncrna/lncrna_classification.tsv",
        filter_stats="results/normalized/filter_stats.tsv",
        consensus="results/stratification/consensus_summary.tsv",
        compare="results/stratification/lncrna_vs_coding.tsv",
        signature="results/signature/signature_lncrnas.tsv",
        performance="results/signature/nested_cv_performance.tsv",
        gxe="results/de/gxe_results.tsv" if HAS_CONDITION else [],
        egenes="results/eqtl/egenes.tsv" if HAS_GENOTYPES else [],
    output:
        "results/report.md",
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/report.py"
