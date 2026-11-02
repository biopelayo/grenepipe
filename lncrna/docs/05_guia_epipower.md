# Guía paso a paso: ejecutar todo en EpiPower

Guía para correr el proyecto completo en EpiPower (Linux) tras conectar por VPN. Todo lo que
se puede automatizar está en `epipower_setup.sh`; aquí va el porqué de cada paso, qué revisar
y qué hacer si algo falla.

## 0. Antes de empezar (5 min)

Conéctate a la VPN y entra en EpiPower por SSH. Comprueba recursos y red:

```bash
nproc; free -g; df -h /data 2>/dev/null || df -h $HOME
for h in https://www.ebi.ac.uk https://eutils.ncbi.nlm.nih.gov https://ftp.ebi.ac.uk https://1001genomes.org https://www.arabidopsis.org https://conda.anaconda.org; do printf "%-40s " $h; curl -sS -m 15 -o /dev/null -w "%{http_code}\n" $h; done
```

Necesitas: ≥ 16 núcleos, ≥ 32 GB RAM (STAR usa ~10 GB por job), **≥ 1 TB libre** para
GSE224760 completo (350 GB fastq + ~300 GB de BAM/intermedios), conda o mamba instalado
(si no: Miniforge, `https://github.com/conda-forge/miniforge`), y acceso HTTPS a los seis
hosts de arriba (todos deben devolver 200 o 301/302, nunca 000).

Decide tres variables y expórtalas en cada sesión (o ponlas en `~/.bashrc`):

```bash
export EMAIL=tu@correo        # NCBI lo exige
export WORK=/data/lncrna      # disco grande
export CORES=32               # núcleos físicos que quieras dedicar
```

## 1. Bootstrap automático (30-60 min, la mayoría descarga)

```bash
curl -sL https://raw.githubusercontent.com/biopelayo/grenepipe/pelamovic/tender-ptolemy-8dvjns/lncrna/epipower_setup.sh -o epipower_setup.sh
bash epipower_setup.sh
```

Hace, en orden e idempotente (si falla, corrige y relanza; continúa donde se quedó):

1. comprueba disco, RAM y red;
2. clona `biopelayo/grenepipe` rama `pelamovic/tender-ptolemy-8dvjns` en `$WORK/grenepipe`;
3. crea el env conda `lncrna` (snakemake 8 + mamba);
4. construye `config/samples.GSE224760.csv` desde GEO/ENA;
5. descarga y prepara referencias (TAIR10, Araport11, TEs, VCF 1001G, índice STAR) y crea
   los 9 entornos conda de herramientas en `$WORK/conda-envs`;
6. lanza el **piloto**: 2 accesiones con todos sus lotes hasta `lncrna/lncrna.gtf`.

**Punto de revisión obligatorio** entre 4 y 5: abre la sample sheet.

```bash
column -s, -t $WORK/grenepipe/lncrna/config/samples.GSE224760.csv | less -S
```

Comprueba: ~80-100 filas; `accession_id` numérico en todas (p. ej. 6909, 9888); `batch`
= rep1..rep4; `strandedness` = reverse; `fastq_2` no vacío (son PE). Si alguna fila tiene
`accession_id` vacío, rellénala a mano con el id del título GSM.

## 2. Interpretar el piloto (lo que decide si seguimos)

Al terminar, mira el embudo:

```bash
RUN=$WORK/runs/GSE224760
tail -8 $RUN/logs/lncrna/final.log
grep -E "Uniquely mapped reads %|% of reads mapped to multiple" $RUN/align/*.Log.final.out
cat $RUN/results/strand/*.txt | sort | uniq -c
```

Esperado con 2 accesiones (≈ 6 muestras): mapeo único > 85 %; strand `reverse` en todas;
varios cientos a ~2.000 candidatos; 200-800 lncRNAs finales (lincRNA > antisense > intrónico
en novedad; con la anotación de Kornienko fusionada, antisense dominará). Si CPC2, StringTie
o salmon fallan, el error está en `$RUN/logs/<fase>/<muestra>.log`; los fallos más
probables en la primera ejecución son nombres de ejecutable (`CPC2.py`) y versiones de conda,
no lógica. Corrige en `workflow/rules/*.smk` o `workflow/envs/*.yaml`, haz commit, y relanza
el mismo comando: Snakemake continúa desde lo hecho.

## 3. GTF publicado de Kornienko 2024 (manual, 10 min)

Baja el fichero de anotación del suplemento de GSE224760 (GEO → "Supplementary file"), y
adapta nombres de cromosomas a Ensembl si vienen como `Chr1`:

```bash
mkdir -p $WORK/refs && cd $WORK/refs
# descargar el .gtf(.gz) del suplemento de GSE224760 aquí, luego:
zcat -f Kornienko*.gtf* | sed -E 's/^Chr([1-5])\t/\1\t/; s/^ChrM\t/Mt\t/; s/^ChrC\t/Pt\t/' > Kornienko2024_lncRNA.ensembl.gtf
cut -f1 Kornienko2024_lncRNA.ensembl.gtf | sort | uniq -c
```

Edita `lncrna/config/config.GSE224760.yaml` → `extra_lncrna_gtf: ["/data/lncrna/refs/Kornienko2024_lncRNA.ensembl.gtf"]`.

## 4. Corrida completa de GSE224760 (1-3 días)

Siempre dentro de `tmux` (o `screen`) para que sobreviva a la desconexión de la VPN:

```bash
tmux new -s lnc
conda activate lncrna
export EMAIL=... WORK=/data/lncrna CORES=32
cd $WORK/runs/GSE224760
$WORK/grenepipe/lncrna/run_local.sh 4 GSE224760 $WORK/runs/GSE224760 2>&1 | tee full.log
```

Si `run_local.sh` no admite `--conda-prefix`, usa el comando directo:

```bash
snakemake --snakefile $WORK/grenepipe/lncrna/workflow/Snakefile \
  --configfile $WORK/grenepipe/lncrna/config/config.GSE224760.yaml \
  --use-conda --conda-frontend mamba --conda-prefix $WORK/conda-envs \
  --cores $CORES --rerun-incomplete --keep-going --printshellcmds 2>&1 | tee full.log
```

Desacoplar de tmux: `Ctrl-b d`; volver: `tmux attach -t lnc`. Progreso:

```bash
grep -E "steps \([0-9.]+%\) done" full.log | tail -1
ls align/*.bam 2>/dev/null | wc -l     # muestras alineadas
du -sh reads align quant 2>/dev/null
```

Limitar memoria si STAR satura la RAM: añade `--resources mem_mb=$((MEM_GB*900))` al comando.

Ahorrar disco: cuando `align/` esté completo y `assembly/merged.gtf` exista, los fastq
descargados en `reads/raw` ya se borraron (son `temp`); `reads/trimmed` se puede borrar
después de `quant/final/*/quant.sf`.

## 5. Resultados (qué mirar y en qué orden)

```bash
cat results/report.md
xdg-open results/qc/multiqc.html 2>/dev/null || echo "copia results/qc/multiqc.html a tu portátil con scp"
cat results/stratification/consensus_summary.tsv
cat results/stratification/lncrna_vs_coding.tsv
head results/signature/signature_lncrnas.tsv; cat results/signature/nested_cv_performance.tsv
head results/eqtl/egenes.tsv
```

Figuras (PDF): `results/normalized/pca.pdf`, `results/stratification/consensus_plots.pdf`,
`results/signature/signature_heatmap.pdf`, `results/eqtl/eqtl_plots.pdf`. Para verlas en tu
portátil: `scp -r epipower:$WORK/runs/GSE224760/results ./results_GSE224760`.

Para comparar clusters con grupos genéticos, baja la tabla de accesiones de 1001genomes.org,
construye `config/accession_metadata.csv` (plantilla en `config/accession_metadata.template.csv`)
y ponla en `analysis.external_labels` del config; relanza solo `results/stratification/clusters.tsv`.

## 6. Después: GSE226691 y GSE80744

```bash
# tejidos (PE 150, stranded); usa el GTF del run anterior como extra_lncrna_gtf
$WORK/grenepipe/lncrna/run_local.sh 1 GSE226691 && $WORK/grenepipe/lncrna/run_local.sh 4 GSE226691 $WORK/runs/GSE226691
# validación (SE 100, no stranded): ~1 TB; solo lincRNA/intrónico; eQTL con >400 accesiones
$WORK/grenepipe/lncrna/run_local.sh 1 GSE80744  && $WORK/grenepipe/lncrna/run_local.sh 4 GSE80744  $WORK/runs/GSE80744
```

Para GSE80744 lanza por lotes si el disco no llega: filtra la sample sheet en trozos de ~100
muestras, corre hasta `quant/final/*/quant.sf`, borra `align/` y sigue con el siguiente trozo;
al final une las `quant.sf` y ejecuta solo las reglas de análisis.

## 7. Problemas típicos

| Síntoma | Causa | Solución |
|---|---|---|
| `curl: (7)` / `000` en el chequeo de red | VPN caída o proxy corporativo | reconectar VPN; si hay proxy, `export https_proxy=...` |
| `CreateCondaEnvironmentException` | mamba 2 o canales | el env usa `mamba<2`; si persiste, `--conda-frontend conda` |
| `fetch_fastq` falla 3 veces | ENA sin fastq para ese run | `conda install sra-tools` y `fasterq-dump --split-files SRRxxxx`; apunta `fastq_1/2` a los ficheros locales |
| STAR `EXITING because of FATAL ERROR: not enough memory` | RAM | `--resources mem_mb=...` o menos `CORES` |
| `KeyError: 'batch'` o columna ausente | sample sheet hecha a mano | añade columna `batch` (puede ser igual a `study`) |
| `CPC2.py: command not found` | nombre del binario en bioconda | `ls $WORK/conda-envs/*/bin | grep -i cpc` y ajusta `rules/lncrna.smk` |
| eQTL: `fewer than 20 accessions with both expression and genotypes` | `accession_id` vacío o no coincide con el VCF | ids numéricos del 1001G en la sample sheet |
| R: `vst` / `GenomeInfoDbData` | paquetes de datos de Bioconductor | con red abierta se instalan solos; si no, `R -e 'BiocManager::install("GenomeInfoDbData")'` en el env `r` |

## 8. Si quieres que lo ejecute yo en EpiPower

Esta sesión no puede abrir VPN ni SSH. La vía es **Remote Control** desde EpiPower:

```bash
cd $WORK/grenepipe && claude remote-control
```

(o en la app de escritorio de Claude: Settings → Claude Code → Remote Control → añadir la
carpeta `$WORK/grenepipe`). En cuanto aparezca, me conecto a esa carpeta y lanzo
`epipower_setup.sh` y el resto desde aquí, con tus permisos.
