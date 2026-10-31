# =================================================================================================
#     1001 Genomes genotypes for the accessions in the sample sheet
# =================================================================================================

G = config["genotypes"]


rule download_vcf:
    output:
        G["vcf"],
    params:
        url=G["vcf_url"],
    log:
        "logs/genotypes/download_vcf.log",
    conda:
        "../envs/utils.yaml"
    shell:
        "curl -L --retry 5 -o {output} '{params.url}' 2> {log}"


rule accession_list:
    input:
        config["samples"],
    output:
        "genotypes/accessions.txt",
    run:
        ids = sorted(set(x for x in SAMPLES["accession_id"] if x))
        with open(output[0], "w") as f:
            f.write("\n".join(ids) + "\n")


rule genotype_matrix:
    """Biallelic SNPs, MAF/missingness filtered, restricted to the accessions in the sample sheet;
    012 matrix (samples x variants) for MatrixEQTL plus variant positions."""
    input:
        vcf=G["vcf"],
        acc="genotypes/accessions.txt",
    output:
        vcf="genotypes/subset.vcf.gz",
        gt="genotypes/genotypes_012.tsv",
        pos="genotypes/variant_positions.tsv",
    params:
        maf=G["maf"],
        max_missing=G["max_missing"],
    log:
        "logs/genotypes/genotype_matrix.log",
    threads: 8
    conda:
        "../envs/genotypes.yaml"
    shell:
        """
        (
        bcftools view --threads {threads} -S {input.acc} --force-samples -m2 -M2 -v snps {input.vcf} | \
          bcftools view --threads {threads} -q {params.maf}:minor -i 'F_MISSING<={params.max_missing}' \
          -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' {output.vcf} | \
          awk -v OFS='\\t' 'BEGIN{{print "snp","chr","pos"}} {{id=($3=="."?$1"_"$2:$3); print id,$1,$2}}' > {output.pos}
        (printf "snp"; bcftools query -l {output.vcf} | tr '\\n' '\\t' | sed 's/\\t$/\\n/';
         bcftools query -f '%CHROM\\t%POS\\t%ID[\\t%GT]\\n' {output.vcf} | \
           awk -v OFS='\\t' '{{id=($3=="."?$1"_"$2:$3); printf id; for(i=4;i<=NF;i++){{g=$i; if(g=="0/0"||g=="0|0")v=0; else if(g=="1/1"||g=="1|1")v=2; else if(g=="./."||g==".|."||g==".")v="NA"; else v=1; printf "\\t%s",v}}; printf "\\n"}}') > {output.gt}
        ) 2> {log}
        """
