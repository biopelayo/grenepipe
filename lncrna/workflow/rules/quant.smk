# =================================================================================================
#     Quantification (salmon) and expression matrices
# =================================================================================================
#
# Two passes: (1) candidates only, to get expression support for the lncRNA filter;
#             (2) final coding + lncRNA transcriptome with genome decoys, for the analysis matrices.


rule salmon_index_candidates:
    input:
        "lncrna/candidates.fa",
    output:
        directory("quant/candidates_index"),
    log:
        "logs/quant/salmon_index_candidates.log",
    threads: 8
    conda:
        "../envs/quant.yaml"
    shell:
        "salmon index -t {input} -i {output} -p {threads} 2> {log}"


rule salmon_quant_candidates:
    input:
        unpack(clean_reads),
        index="quant/candidates_index",
    output:
        "quant/candidates/{sample}/quant.sf",
    params:
        outdir="quant/candidates/{sample}",
        reads=lambda w, input: f"-1 {input.r1} -2 {input.r2}" if is_paired(w.sample) else f"-r {input.r1}",
    log:
        "logs/quant/candidates/{sample}.log",
    threads: 8
    conda:
        "../envs/quant.yaml"
    shell:
        "salmon quant -i {input.index} -l A -p {threads} {params.reads} -o {params.outdir} 2> {log}"


rule candidates_tpm:
    input:
        expand("quant/candidates/{sample}/quant.sf", sample=SAMPLE_NAMES),
    output:
        "quant/candidates_tpm.tsv",
    params:
        samples=SAMPLE_NAMES,
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/merge_quant.py"


rule final_transcriptome:
    input:
        gtf="lncrna/coding_plus_lncrna.gtf",
        fa=REF,
    output:
        tx="quant/final_transcripts.fa",
        gentrome="quant/gentrome.fa",
        decoys="quant/decoys.txt",
        tx2gene="quant/tx2gene.tsv",
    log:
        "logs/quant/final_transcriptome.log",
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        (
        gffread -w {output.tx} -g {input.fa} {input.gtf}
        grep '>' {input.fa} | sed 's/>//; s/ .*//' > {output.decoys}
        cat {output.tx} {input.fa} > {output.gentrome}
        awk -F'\\t' '$3=="transcript" || $3=="exon"' {input.gtf} | \
          grep -oP 'transcript_id "\\K[^"]+|gene_id "\\K[^"]+' | paste - - | \
          awk -v OFS='\\t' '{{print $2,$1}}' | sort -u > {output.tx2gene}
        ) 2> {log}
        """


rule salmon_index_final:
    input:
        gentrome="quant/gentrome.fa",
        decoys="quant/decoys.txt",
    output:
        directory("quant/final_index"),
    log:
        "logs/quant/salmon_index_final.log",
    threads: 16
    conda:
        "../envs/quant.yaml"
    shell:
        "salmon index -t {input.gentrome} -d {input.decoys} -i {output} -p {threads} 2> {log}"


rule salmon_quant_final:
    input:
        unpack(clean_reads),
        index="quant/final_index",
    output:
        "quant/final/{sample}/quant.sf",
    params:
        outdir="quant/final/{sample}",
        reads=lambda w, input: f"-1 {input.r1} -2 {input.r2}" if is_paired(w.sample) else f"-r {input.r1}",
    log:
        "logs/quant/final/{sample}.log",
    threads: 8
    conda:
        "../envs/quant.yaml"
    shell:
        """
        salmon quant -i {input.index} -l A -p {threads} {params.reads} --validateMappings \
            --gcBias --seqBias -o {params.outdir} 2> {log}
        """


rule matrices:
    input:
        quants=expand("quant/final/{sample}/quant.sf", sample=SAMPLE_NAMES),
        tx2gene="quant/tx2gene.tsv",
        classification="lncrna/lncrna_classification.tsv",
        coding="resources/annotation.coding_genes.txt",
        samples=config["samples"],
    output:
        counts_lnc="matrices/counts_lncrna.tsv",
        tpm_lnc="matrices/tpm_lncrna.tsv",
        counts_cod="matrices/counts_coding.tsv",
        tpm_cod="matrices/tpm_coding.tsv",
        counts_all="matrices/counts_all.tsv",
        tpm_all="matrices/tpm_all.tsv",
        metadata="matrices/metadata.tsv",
    log:
        "logs/quant/matrices.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/tximport_matrices.R"


rule multiqc:
    input:
        expand("results/qc/fastp/{sample}.json", sample=SAMPLE_NAMES),
        expand("align/{sample}.Log.final.out", sample=SAMPLE_NAMES),
        expand("quant/final/{sample}/quant.sf", sample=SAMPLE_NAMES),
        expand("results/qc/sortmerna/{sample}.log", sample=SAMPLE_NAMES) if config["params"]["sortmerna"] else [],
    output:
        "results/qc/multiqc.html",
    log:
        "logs/qc/multiqc.log",
    conda:
        "../envs/utils.yaml"
    shell:
        "multiqc -f -o results/qc -n multiqc.html results/qc align quant/final 2> {log}"
