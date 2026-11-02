# =================================================================================================
#     STAR alignment (2-pass)
# =================================================================================================


rule star_index:
    input:
        fa=REF,
        gtf=GTF,
    output:
        directory("resources/star_index"),
    params:
        overhang=config["params"]["star_sjdb_overhang"],
    log:
        "logs/align/star_index.log",
    threads: 16
    conda:
        "../envs/align.yaml"
    shell:
        """
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} \
             --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.overhang} \
             --genomeSAindexNbases 12 2> {log}
        """


rule star_align:
    input:
        unpack(clean_reads),
        index="resources/star_index",
    output:
        bam="align/{sample}.Aligned.sortedByCoord.out.bam",
        sj="align/{sample}.SJ.out.tab",
        log="align/{sample}.Log.final.out",
    params:
        prefix="align/{sample}.",
        extra=config["params"]["star_extra"],
        reads=lambda w, input: f"{input.r1} {input.r2}" if is_paired(w.sample) else input.r1,
    log:
        "logs/align/star/{sample}.log",
    threads: 16
    resources:
        mem_mb=16000,
    conda:
        "../envs/align.yaml"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {params.reads} \
             --readFilesCommand zcat --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate --twopassMode Basic \
             --outSAMattributes NH HI AS nM NM MD XS --quantMode GeneCounts {params.extra} 2> {log}
        """


rule bam_index:
    input:
        "align/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        "align/{sample}.Aligned.sortedByCoord.out.bam.bai",
    conda:
        "../envs/utils.yaml"
    shell:
        "samtools index {input}"
