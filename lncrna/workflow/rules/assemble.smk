# =================================================================================================
#     Transcript assembly (StringTie) and comparison with the reference (gffcompare)
# =================================================================================================


rule stringtie_assemble:
    input:
        bam="align/{sample}.Aligned.sortedByCoord.out.bam",
        bai="align/{sample}.Aligned.sortedByCoord.out.bam.bai",
        gtf=GTF,
        strand="results/strand/{sample}.txt",
    output:
        "assembly/{sample}.gtf",
    params:
        extra=config["params"]["stringtie_assemble"],
    log:
        "logs/assembly/stringtie/{sample}.log",
    threads: 8
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        case "$(cat {input.strand})" in
            reverse) STRAND="--rf" ;;
            forward) STRAND="--fr" ;;
            *)       STRAND="" ;;
        esac
        stringtie {input.bam} -G {input.gtf} $STRAND {params.extra} -p {threads} -o {output} 2> {log}
        """


rule stringtie_merge:
    input:
        gtfs=expand("assembly/{sample}.gtf", sample=SAMPLE_NAMES),
        ref=GTF,
    output:
        "assembly/merged.gtf",
    params:
        extra=config["params"]["stringtie_merge"],
    log:
        "logs/assembly/stringtie_merge.log",
    threads: 8
    conda:
        "../envs/assembly.yaml"
    shell:
        "stringtie --merge -G {input.ref} {params.extra} -p {threads} -o {output} {input.gtfs} 2> {log}"


rule gffcompare:
    input:
        merged="assembly/merged.gtf",
        ref=GTF,
    output:
        annotated="assembly/gffcmp.annotated.gtf",
        tmap="assembly/gffcmp.merged.gtf.tmap",
    params:
        prefix="assembly/gffcmp",
    log:
        "logs/assembly/gffcompare.log",
    conda:
        "../envs/assembly.yaml"
    shell:
        "gffcompare -r {input.ref} -o {params.prefix} {input.merged} 2> {log}"
