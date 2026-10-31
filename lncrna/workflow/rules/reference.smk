# =================================================================================================
#     Reference genome, annotation, transposable elements
# =================================================================================================

ENS = config["reference"]
ENS_BASE = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-{rel}".format(rel=ENS["ensembl_release"])
SPECIES = ENS["species"]
SPECIES_CAP = SPECIES.capitalize()


rule download_genome:
    output:
        temp("resources/genome.raw.fa.gz"),
    params:
        url=f"{ENS_BASE}/fasta/{SPECIES}/dna/{SPECIES_CAP}.{ENS['assembly']}.dna.toplevel.fa.gz",
    log:
        "logs/reference/download_genome.log",
    conda:
        "../envs/utils.yaml"
    shell:
        "curl -L --retry 5 -o {output} {params.url} 2> {log}"


rule prepare_genome:
    input:
        "resources/genome.raw.fa.gz",
    output:
        fa=REF,
        fai=REF + ".fai",
        sizes="resources/genome.sizes",
    params:
        chroms=" ".join(ENS["chromosomes"]),
    log:
        "logs/reference/prepare_genome.log",
    conda:
        "../envs/utils.yaml"
    shell:
        """
        (
        zcat {input} | seqkit grep -p {params.chroms} | seqkit seq -w 60 > {output.fa}
        samtools faidx {output.fa}
        cut -f1,2 {output.fai} > {output.sizes}
        ) 2> {log}
        """


rule download_annotation:
    output:
        temp("resources/annotation.raw.gtf.gz"),
    params:
        url=f"{ENS_BASE}/gtf/{SPECIES}/{SPECIES_CAP}.{ENS['assembly']}.{ENS['ensembl_release']}.gtf.gz",
    log:
        "logs/reference/download_annotation.log",
    conda:
        "../envs/utils.yaml"
    shell:
        "curl -L --retry 5 -o {output} {params.url} 2> {log}"


rule prepare_annotation:
    input:
        gtf="resources/annotation.raw.gtf.gz",
        fai=REF + ".fai",
    output:
        gtf=GTF,
        coding="resources/annotation.coding_genes.txt",
        known_lnc="resources/annotation.known_lncrna_genes.txt",
    log:
        "logs/reference/prepare_annotation.log",
    conda:
        "../envs/utils.yaml"
    shell:
        """
        (
        zcat {input.gtf} | grep -v '^#' | awk -v OFS='\\t' 'NR==FNR{{ok[$1]=1; next}} ($1 in ok)' {input.fai} - > {output.gtf}
        grep -P '\\tgene\\t' {output.gtf} | grep 'gene_biotype "protein_coding"' | grep -oP 'gene_id "\\K[^"]+' | sort -u > {output.coding}
        grep -P '\\tgene\\t' {output.gtf} | grep -E 'gene_biotype "(lncRNA|lnc_RNA|antisense|lincRNA|ncRNA|other_RNA)"' | grep -oP 'gene_id "\\K[^"]+' | sort -u > {output.known_lnc} || true
        ) 2> {log}
        """


rule reference_transcripts:
    input:
        fa=REF,
        gtf=GTF,
    output:
        TX,
    log:
        "logs/reference/transcripts.log",
    conda:
        "../envs/assembly.yaml"
    shell:
        "gffread -w {output} -g {input.fa} {input.gtf} 2> {log}"


rule download_te:
    output:
        temp("resources/te.raw.txt"),
    params:
        url=ENS["te_annotation_url"],
    log:
        "logs/reference/download_te.log",
    conda:
        "../envs/utils.yaml"
    shell:
        "curl -L --retry 5 -o {output} '{params.url}' 2> {log}"


rule te_bed:
    input:
        "resources/te.raw.txt",
    output:
        TE,
    log:
        "logs/reference/te_bed.log",
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/te_to_bed.py"


rule sortmerna_db:
    output:
        directory("resources/sortmerna_db"),
    params:
        url=config["params"]["sortmerna_db_url"],
    log:
        "logs/reference/sortmerna_db.log",
    conda:
        "../envs/utils.yaml"
    shell:
        """
        (mkdir -p {output} && curl -L --retry 5 {params.url} | tar -xz -C {output}) 2> {log}
        """
