# =================================================================================================
#     Read retrieval, trimming, rRNA removal, strandedness
# =================================================================================================


def fastq_inputs(wildcards):
    """Local files are used directly; URLs are fetched by rule fetch_fastq."""
    s = SAMPLES.loc[wildcards.sample]
    files = [s["fastq_1"]] + ([s["fastq_2"]] if s["fastq_2"] else [])
    out = []
    for i, f in enumerate(files, start=1):
        if f.startswith(("ftp://", "http://", "https://")):
            out.append(f"reads/raw/{wildcards.sample}_{i}.fastq.gz")
        else:
            out.append(f)
    return out


rule fetch_fastq:
    output:
        temp("reads/raw/{sample}_{read}.fastq.gz"),
    params:
        url=lambda w: SAMPLES.loc[w.sample, f"fastq_{w.read}"],
    log:
        "logs/reads/fetch/{sample}_{read}.log",
    conda:
        "../envs/utils.yaml"
    retries: 3
    shell:
        "curl -L --retry 5 --retry-delay 10 -o {output} '{params.url}' 2> {log}"


rule fastp:
    input:
        fastq_inputs,
    output:
        r1="reads/trimmed/{sample}_1.fastq.gz",
        r2="reads/trimmed/{sample}_2.fastq.gz",
        html="results/qc/fastp/{sample}.html",
        json="results/qc/fastp/{sample}.json",
    params:
        extra=config["params"]["fastp"],
        io=lambda w, input, output: (
            f"-i {input[0]} -I {input[1]} -o {output.r1} -O {output.r2}"
            if is_paired(w.sample)
            else f"-i {input[0]} -o {output.r1}"
        ),
    log:
        "logs/reads/fastp/{sample}.log",
    threads: 8
    conda:
        "../envs/reads.yaml"
    shell:
        """
        fastp -w {threads} {params.extra} {params.io} -h {output.html} -j {output.json} 2> {log}
        touch {output.r2}
        """


rule sortmerna:
    input:
        r1="reads/trimmed/{sample}_1.fastq.gz",
        r2="reads/trimmed/{sample}_2.fastq.gz",
        db="resources/sortmerna_db",
    output:
        r1="reads/clean/{sample}_1.fastq.gz",
        r2="reads/clean/{sample}_2.fastq.gz",
        log="results/qc/sortmerna/{sample}.log",
    params:
        paired=lambda w: is_paired(w.sample),
        workdir="reads/clean/{sample}_sortmerna",
    log:
        "logs/reads/sortmerna/{sample}.log",
    threads: 8
    conda:
        "../envs/reads.yaml"
    shell:
        """
        rm -rf {params.workdir}
        REFS=$(ls {input.db}/smr_v4.3_default_db.fasta 2>/dev/null || ls {input.db}/*.fasta | head -1)
        if [ "{params.paired}" = "True" ]; then
            sortmerna --ref $REFS --reads {input.r1} --reads {input.r2} --workdir {params.workdir} \
                --threads {threads} --fastx --paired_in --out2 --other {params.workdir}/other \
                --aligned {params.workdir}/rrna 2> {log}
            mv {params.workdir}/other_fwd.f*q.gz {output.r1}
            mv {params.workdir}/other_rev.f*q.gz {output.r2}
        else
            sortmerna --ref $REFS --reads {input.r1} --workdir {params.workdir} \
                --threads {threads} --fastx --other {params.workdir}/other \
                --aligned {params.workdir}/rrna 2> {log}
            mv {params.workdir}/other.f*q.gz {output.r1}
            touch {output.r2}
        fi
        cp {params.workdir}/rrna.log {output.log}
        rm -rf {params.workdir}
        """


def clean_reads(wildcards):
    base = "reads/clean" if config["params"]["sortmerna"] else "reads/trimmed"
    return {"r1": f"{base}/{wildcards.sample}_1.fastq.gz", "r2": f"{base}/{wildcards.sample}_2.fastq.gz"}


rule salmon_index_reference:
    input:
        TX,
    output:
        directory("resources/salmon_reference_index"),
    log:
        "logs/reads/salmon_index_reference.log",
    threads: 8
    conda:
        "../envs/quant.yaml"
    shell:
        "salmon index -t {input} -i {output} -p {threads} 2> {log}"


rule infer_strandedness:
    """salmon on a read subsample against the reference transcriptome; lib_format_counts.json gives
    the library type (ISR/ISF/IU ...), which we translate for StringTie and reporting."""
    input:
        unpack(clean_reads),
        index="resources/salmon_reference_index",
    output:
        "results/strand/{sample}.txt",
    params:
        n=config["params"]["strand_subsample"] * 4,
        paired=lambda w: is_paired(w.sample),
        declared=lambda w: SAMPLES.loc[w.sample, "strandedness"],
        outdir="results/strand/{sample}_salmon",
    log:
        "logs/reads/strand/{sample}.log",
    threads: 4
    conda:
        "../envs/quant.yaml"
    shell:
        """
        if [ "{params.declared}" != "auto" ]; then
            echo "{params.declared}" > {output}
        else
            if [ "{params.paired}" = "True" ]; then
                salmon quant -i {input.index} -l A -p {threads} \
                    -1 <(zcat {input.r1} | head -n {params.n}) -2 <(zcat {input.r2} | head -n {params.n}) \
                    -o {params.outdir} 2> {log}
            else
                salmon quant -i {input.index} -l A -p {threads} \
                    -r <(zcat {input.r1} | head -n {params.n}) -o {params.outdir} 2> {log}
            fi
            python3 -c '
import json,sys
j=json.load(open("{params.outdir}/lib_format_counts.json"))
t=j["expected_format"]
m={{"ISR":"reverse","SR":"reverse","ISF":"forward","SF":"forward","IU":"unstranded","U":"unstranded"}}
print(m.get(t,"unstranded"))' > {output}
            rm -rf {params.outdir}
        fi
        """
