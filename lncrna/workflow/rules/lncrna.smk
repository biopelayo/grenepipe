# =================================================================================================
#     lncRNA candidate selection, coding potential, final annotation
# =================================================================================================

LF = config["lncrna_filter"]


rule lncrna_candidates:
    """Novel transcripts by gffcompare class code (u/x/i), length and exon filters, on nuclear
    chromosomes; writes candidate GTF + FASTA."""
    input:
        gtf="assembly/gffcmp.annotated.gtf",
        tmap="assembly/gffcmp.merged.gtf.tmap",
        fa=REF,
    output:
        gtf="lncrna/candidates.gtf",
        fa="lncrna/candidates.fa",
        tsv="lncrna/candidates.tsv",
    params:
        min_length=LF["min_length"],
        class_codes=LF["class_codes"],
        min_exons_intergenic=LF["min_exons_intergenic"],
        chromosomes=config["reference"]["lncrna_chromosomes"],
    log:
        "logs/lncrna/candidates.log",
    conda:
        "../envs/assembly.yaml"
    script:
        "../scripts/lncrna_candidates.py"


rule cpc2:
    input:
        "lncrna/candidates.fa",
    output:
        "lncrna/cpc2.txt",
    log:
        "logs/lncrna/cpc2.log",
    conda:
        "../envs/cpc2.yaml"
    shell:
        "CPC2.py -i {input} -o {output} 2> {log}"


rule feelnc:
    input:
        cand="lncrna/candidates.gtf",
        ref=GTF,
        fa=REF,
    output:
        "lncrna/feelnc/candidates.gtf.lncRNA.gtf",
    params:
        outdir="lncrna/feelnc",
    log:
        "logs/lncrna/feelnc.log",
    threads: 8
    conda:
        "../envs/feelnc.yaml"
    shell:
        """
        FEELnc_codpot.pl -i {input.cand} -a {input.ref} -g {input.fa} --mode=shuffle \
            --outdir {params.outdir} --outname candidates.gtf -p {threads} 2> {log}
        """


rule pfam_db:
    output:
        "resources/Pfam-A.hmm",
    params:
        url=LF["pfam_url"],
    log:
        "logs/lncrna/pfam_db.log",
    conda:
        "../envs/utils.yaml"
    shell:
        "(curl -L --retry 5 {params.url} | gunzip -c > {output} && hmmpress {output}) 2> {log}"


rule pfam_scan:
    input:
        fa="lncrna/candidates.fa",
        db="resources/Pfam-A.hmm",
    output:
        "lncrna/pfam.domtblout",
    log:
        "logs/lncrna/pfam.log",
    threads: 16
    conda:
        "../envs/utils.yaml"
    shell:
        """
        (seqkit translate -f 6 -w 0 {input.fa} | sed 's/\\*//g' > lncrna/candidates.6frame.faa
         hmmscan --cpu {threads} -E 1e-5 --domtblout {output} {input.db} lncrna/candidates.6frame.faa > /dev/null) 2> {log}
        """


def lncrna_evidence(wildcards):
    inp = {"cand_gtf": "lncrna/candidates.gtf", "cand_tsv": "lncrna/candidates.tsv", "cpc2": "lncrna/cpc2.txt",
           "te": TE, "tpm": "quant/candidates_tpm.tsv"}
    if LF["use_feelnc"]:
        inp["feelnc"] = "lncrna/feelnc/candidates.gtf.lncRNA.gtf"
    if LF["use_pfam"]:
        inp["pfam"] = "lncrna/pfam.domtblout"
    return inp


rule lncrna_final:
    """Combine coding-potential, Pfam, expression support and TE overlap into the final lncRNA GTF,
    plus classification table (lincRNA / antisense / intronic; TE-overlapping; novel vs known)."""
    input:
        unpack(lncrna_evidence),
        ref=GTF,
        extra=config["reference"]["extra_lncrna_gtf"],
    output:
        gtf="lncrna/lncrna.gtf",
        tsv="lncrna/lncrna_classification.tsv",
        combined="lncrna/coding_plus_lncrna.gtf",
    params:
        min_tpm=LF["min_tpm"],
        min_samples=LF["min_samples"],
        use_feelnc=LF["use_feelnc"],
        use_pfam=LF["use_pfam"],
    log:
        "logs/lncrna/final.log",
    conda:
        "../envs/assembly.yaml"
    script:
        "../scripts/lncrna_final.py"
