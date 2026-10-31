#!/usr/bin/env bash
# =================================================================================================
#  Step-by-step local driver for the lncRNA population workflow on the real datasets.
#
#  Usage:  ./run_local.sh <step> [GSE] [RUNDIR]
#    step 0  : create the conda env with snakemake + mamba (once)
#    step 1  : build the sample sheet for GSE (default GSE224760); review it by hand afterwards
#    step 2  : download/prepare references only (genome, annotation, TEs, 1001G VCF)
#    step 3  : pilot: first 2 accessions of GSE224760, up to the final lncRNA GTF
#    step 4  : full run of the dataset
#    step 5  : report only (after step 4)
#    dry     : dry-run of the full DAG for the dataset (no execution)
#
#  Requires: conda/mamba, ~2 TB disk for the big cohorts, internet access to ENA/NCBI/Ensembl.
#  Set EMAIL to your address (NCBI E-utilities policy).
# =================================================================================================
set -euo pipefail
STEP="${1:?step required}"
GSE="${2:-GSE224760}"
RUN="${3:-$PWD/runs/$GSE}"
EMAIL="${EMAIL:-you@example.org}"
CORES="${CORES:-16}"
HERE="$(cd "$(dirname "$0")" && pwd)"
CFG="$HERE/config/config.$GSE.yaml"
SNK="snakemake --snakefile $HERE/workflow/Snakefile --use-conda --conda-frontend mamba --cores $CORES --rerun-incomplete --keep-going --printshellcmds"

mkdir -p "$RUN/config"

case "$STEP" in
  0)
    conda create -y -n lncrna -c conda-forge -c bioconda "snakemake>=8.15" "mamba<2" pandas requests
    echo ">> conda activate lncrna"
    ;;
  1)
    python3 "$HERE/make_samplesheet_from_gse.py" --gse "$GSE" --out "$HERE/config/samples.$GSE.csv" --email "$EMAIL"
    echo ">> review $HERE/config/samples.$GSE.csv: accession_id numeric, batch filled, tissue correct"
    ;;
  2)
    cd "$RUN" && $SNK --configfile "$CFG" resources/te.bed resources/transcripts.fa resources/star_index \
        resources/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
    ;;
  3)
    python3 "$HERE/make_samplesheet_from_gse.py" --gse "$GSE" --out "$RUN/config/samples.pilot.csv" --email "$EMAIL" --max-accessions 2
    sed "s#samples: .*#samples: \"$RUN/config/samples.pilot.csv\"#" "$CFG" > "$RUN/config/config.pilot.yaml"
    cd "$RUN" && $SNK --configfile "$RUN/config/config.pilot.yaml" lncrna/lncrna.gtf lncrna/lncrna_classification.tsv
    echo ">> inspect logs/lncrna/final.log and lncrna/lncrna_classification.tsv"
    ;;
  4)
    cd "$RUN" && $SNK --configfile "$CFG"
    ;;
  5)
    cd "$RUN" && $SNK --configfile "$CFG" results/report.md results/qc/multiqc.html
    echo ">> results/report.md, results/qc/multiqc.html, results/*/**.pdf"
    ;;
  dry)
    cd "$RUN" && $SNK --configfile "$CFG" -n --quiet rules
    ;;
  *)
    echo "unknown step $STEP"; exit 1 ;;
esac
