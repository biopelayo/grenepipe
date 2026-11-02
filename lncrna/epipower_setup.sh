#!/usr/bin/env bash
# =================================================================================================
#  One-shot bootstrap for EpiPower (Linux workstation/server): clone repo, create env, check
#  network/disk, build the GSE224760 sample sheet, prepare references, run the 2-accession pilot.
#
#  Usage (on EpiPower, after VPN):
#     export EMAIL=tu@correo            # NCBI E-utilities needs an email
#     export WORK=/data/lncrna          # >= 2 TB recommended
#     export CORES=32                   # physical cores to use
#     bash <(curl -sL https://raw.githubusercontent.com/biopelayo/grenepipe/pelamovic/tender-ptolemy-8dvjns/lncrna/epipower_setup.sh)
#  or clone first and run:  bash lncrna/epipower_setup.sh
#
#  Each phase is idempotent; re-run the script to continue where it stopped.
# =================================================================================================
set -euo pipefail
EMAIL="${EMAIL:?set EMAIL=you@example.org}"
WORK="${WORK:-$HOME/lncrna-work}"
CORES="${CORES:-$(nproc)}"
BRANCH="pelamovic/tender-ptolemy-8dvjns"
REPO="https://github.com/biopelayo/grenepipe.git"

log() { printf "\n\033[1;32m[%s] %s\033[0m\n" "$(date +%H:%M:%S)" "$*"; }
fail() { printf "\n\033[1;31m!! %s\033[0m\n" "$*"; exit 1; }

# ---- 0. host checks -----------------------------------------------------------------------------
log "0. host checks"
mkdir -p "$WORK"
AVAIL_GB=$(df -BG "$WORK" | awk 'NR==2{gsub("G","",$4); print $4}')
MEM_GB=$(awk '/MemTotal/{printf "%d", $2/1024/1024}' /proc/meminfo)
echo "   work dir : $WORK  (free ${AVAIL_GB} GB)"
echo "   cores    : $(nproc) available, using $CORES"
echo "   memory   : ${MEM_GB} GB"
[ "$AVAIL_GB" -ge 300 ] || echo "   ! <300 GB free: enough for the pilot, NOT for the full GSE224760 (~350 GB fastq + ~300 GB work)"
[ "$MEM_GB" -ge 32 ] || echo "   ! <32 GB RAM: STAR needs ~10 GB per job; lower CORES or set --resources mem_mb"
for h in https://www.ebi.ac.uk https://eutils.ncbi.nlm.nih.gov https://ftp.ebi.ac.uk https://1001genomes.org https://www.arabidopsis.org https://conda.anaconda.org; do
  code=$(curl -sS -m 15 -o /dev/null -w "%{http_code}" "$h" || echo 000)
  printf "   %-40s %s\n" "$h" "$code"
  case "$code" in 000) fail "no network to $h (VPN/proxy?)";; esac
done
command -v conda >/dev/null || command -v mamba >/dev/null || command -v micromamba >/dev/null || \
  fail "conda/mamba not found. Install Miniforge: https://github.com/conda-forge/miniforge (then re-run)"

# ---- 1. repo -----------------------------------------------------------------------------------
log "1. repository"
cd "$WORK"
if [ ! -d grenepipe ]; then git clone --branch "$BRANCH" "$REPO"; fi
cd grenepipe && git fetch -q origin "$BRANCH" && git checkout -q "$BRANCH" && git pull -q --ff-only origin "$BRANCH"
echo "   $(git log --oneline -1)"
LNC="$WORK/grenepipe/lncrna"

# ---- 2. conda env --------------------------------------------------------------------------------
log "2. conda env 'lncrna' (snakemake + mamba)"
eval "$(conda shell.bash hook 2>/dev/null || true)"
if ! conda env list | grep -qE '^lncrna\s'; then
  conda create -y -n lncrna -c conda-forge -c bioconda "snakemake>=8.15,<9" "mamba<2" pandas requests pyyaml
fi
conda activate lncrna
snakemake --version

# ---- 3. sample sheet (GSE224760) ------------------------------------------------------------------
log "3. sample sheet GSE224760"
SHEET="$LNC/config/samples.GSE224760.csv"
if [ ! -s "$SHEET" ]; then
  python3 "$LNC/make_samplesheet_from_gse.py" --gse GSE224760 --out "$SHEET" --email "$EMAIL"
fi
echo "   rows: $(($(wc -l < "$SHEET")-1))   accessions: $(cut -d, -f5 "$SHEET" | tail -n +2 | sort -u | wc -l)   batches: $(cut -d, -f9 "$SHEET" | tail -n +2 | sort -u | tr '\n' ' ')"
echo "   >> REVIEW $SHEET before step 5 (accession_id numeric, batch rep1..rep4, strandedness reverse)"

# ---- 4. references + conda envs of the workflow ------------------------------------------------------
log "4. references and tool environments (one-off, ~20-40 min)"
RUN="$WORK/runs/GSE224760"; mkdir -p "$RUN/config"
CFG="$LNC/config/config.GSE224760.yaml"
cd "$RUN"
snakemake --snakefile "$LNC/workflow/Snakefile" --configfile "$CFG" --use-conda --conda-frontend mamba \
  --conda-prefix "$WORK/conda-envs" --cores "$CORES" --rerun-incomplete \
  resources/te.bed resources/transcripts.fa resources/star_index resources/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
snakemake --snakefile "$LNC/workflow/Snakefile" --configfile "$CFG" --use-conda --conda-frontend mamba \
  --conda-prefix "$WORK/conda-envs" --cores 1 --conda-create-envs-only

# ---- 5. pilot: 2 accessions, up to the lncRNA GTF -------------------------------------------------------
log "5. pilot (2 accessions, all their batches) up to lncrna/lncrna.gtf"
python3 "$LNC/make_samplesheet_from_gse.py" --gse GSE224760 --out "$RUN/config/samples.pilot.csv" --email "$EMAIL" --max-accessions 2
sed "s#samples: .*#samples: \"$RUN/config/samples.pilot.csv\"#" "$CFG" > "$RUN/config/config.pilot.yaml"
snakemake --snakefile "$LNC/workflow/Snakefile" --configfile "$RUN/config/config.pilot.yaml" --use-conda --conda-frontend mamba \
  --conda-prefix "$WORK/conda-envs" --cores "$CORES" --rerun-incomplete --keep-going --printshellcmds \
  lncrna/lncrna.gtf lncrna/lncrna_classification.tsv 2>&1 | tee "$RUN/pilot.log"
echo
tail -8 "$RUN/logs/lncrna/final.log" || true
log "pilot done. Full run:  cd $RUN && $LNC/run_local.sh 4 GSE224760 $RUN   (inside tmux)"
