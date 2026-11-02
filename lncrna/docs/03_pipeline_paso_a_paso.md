# Pipeline desglosado paso a paso

Cada paso: qué hace, qué regla/script lo implementa, entradas, salidas, parámetros que
importan, cómo comprobar que ha ido bien, y cómo ejecutarlo aislado con Snakemake
(`snakemake ... <target>`). Todos los comandos se lanzan desde `lncrna/` con

```bash
SNK="snakemake --use-conda --conda-frontend mamba --cores 16 --directory $RUN --configfile $CFG"
```

donde `$RUN` es la carpeta de análisis y `$CFG` el config del dataset.

## Fase A. Preparación

### A1. Sample sheet desde GEO

- Script: `make_samplesheet_from_gse.py --gse GSE224760 --out config/samples.GSE224760.csv --email tu@correo`
- Hace: GSE → SRP (NCBI E-utilities) → runs y FTP de fastq (ENA) → atributos de BioSample
  (ENA XML) → aplica el **preset** del GSE (strandedness, batch desde el título, tejido).
- Salida: CSV con `sample, fastq_1, fastq_2, accession, accession_id, tissue, condition, study, batch, strandedness`.
- Comprobar: número de filas ≈ muestras esperadas; `accession_id` numérico en todas; `batch`
  rellenado; ningún `fastq_1` vacío. Abre el CSV y revisa a mano los 5 primeros y 5 últimos.

### A2. Referencias (automático)

- Reglas: `download_genome`, `prepare_genome`, `download_annotation`, `prepare_annotation`,
  `reference_transcripts`, `download_te`, `te_bed`, `download_vcf`.
- Hace: baja TAIR10 y Araport11 (Ensembl Plants 59), filtra a `1-5,Mt,Pt`, indexa, extrae
  listas de genes codificantes y lncRNAs ya anotados, convierte la tabla de TEs de TAIR a BED,
  baja el VCF del 1001G.
- Salidas: `resources/genome.fa(.fai)`, `resources/annotation.gtf`, `resources/transcripts.fa`,
  `resources/te.bed`, `resources/1001genomes_*.vcf.gz`.
- Target aislado: `$SNK resources/te.bed resources/transcripts.fa`.
- Comprobar: `grep -c '>' resources/genome.fa` = 7; `wc -l resources/te.bed` ≈ 31.000.

### A3. GTF externo de lncRNAs (manual, opcional pero recomendado)

- Bajar el suplemento de GSE224760 (fichero de anotación de Kornienko), renombrar cromosomas
  a Ensembl si hace falta: `sed -E 's/^Chr([1-5])\t/\1\t/; s/^ChrM\t/Mt\t/; s/^ChrC\t/Pt\t/'`.
- Poner la ruta en `reference.extra_lncrna_gtf`.

## Fase B. Lecturas

### B1. Descarga

- Regla: `fetch_fastq` (curl desde ENA, 3 reintentos). Alternativa: `fasterq-dump`.
- Salida temporal: `reads/raw/<sample>_{1,2}.fastq.gz` (se borra tras fastp).

### B2. fastp

- Regla: `fastp`. Parámetros: `--detect_adapter_for_pe --cut_right ... --length_required 36`.
- Salidas: `reads/trimmed/<sample>_{1,2}.fastq.gz`, `results/qc/fastp/<sample>.{html,json}`.
- Comprobar: en el JSON, `after_filtering.total_reads` > 80 % de `before_filtering`.

### B3. SortMeRNA (solo total RNA)

- Regla: `sortmerna`; desactivado en los configs reales (`params.sortmerna: false`), porque
  las librerías son polyA.

### B4. Strandedness por muestra

- Regla: `infer_strandedness`. Si la sample sheet dice `reverse`/`forward`/`unstranded` se
  respeta; si dice `auto`, salmon sobre 500 k lecturas contra el transcriptoma de referencia
  y se lee `expected_format` (ISR → reverse, ISF → forward, IU → unstranded).
- Salida: `results/strand/<sample>.txt`.
- Comprobar: para GSE224760 todas deben decir `reverse`; si alguna dice otra cosa, esa
  muestra tiene un problema de librería y hay que mirarla.

## Fase C. Alineamiento y ensamblaje

### C1. Índice STAR

- Regla: `star_index`; `--sjdbOverhang` = longitud − 1 (124 en GSE224760).

### C2. STAR 2-pass

- Regla: `star_align`; `--twopassMode Basic`, `--quantMode GeneCounts`, parámetros del paper.
- Salidas: `align/<sample>.Aligned.sortedByCoord.out.bam(.bai)`, `SJ.out.tab`, `Log.final.out`.
- Comprobar: `Uniquely mapped reads %` > 85 % en Arabidopsis; multi-mapped < 10 %.

### C3. StringTie por muestra

- Regla: `stringtie_assemble`; `--rf` si reverse, guiado por Araport11; `-c 2 -m 150 -j 2.5 -a 15`.
- Salida: `assembly/<sample>.gtf`.

### C4. Merge y comparación

- Reglas: `stringtie_merge` (`-m 200 -c 1 -F 0.5 -T 1`), `gffcompare`.
- Salidas: `assembly/merged.gtf`, `assembly/gffcmp.annotated.gtf`, `assembly/gffcmp.merged.gtf.tmap`.
- Comprobar: `assembly/gffcmp.stats`: sensibilidad a nivel de locus > 90 % frente a Araport11.
- Nota: Kornienko hizo merge jerárquico por tejido; con un solo tejido por config es equivalente.
  Si juntas GSE224760 + GSE226691 en una sola corrida, el merge único sigue siendo válido
  pero más permisivo; alternativa: correr cada tejido por separado y fusionar los GTFs finales.

## Fase D. Anotación de lncRNAs

### D1. Candidatos

- Regla: `lncrna_candidates` → `scripts/lncrna_candidates.py`.
- Filtros: class code en `class_codes`; cromosomas nucleares; longitud ≥ `min_length` (200)
  y mono-exónicos ≥ `min_length_monoexonic` (400); intergénicos con ≥ `min_exons_intergenic`.
- Salidas: `lncrna/candidates.{gtf,fa,tsv}`.

### D2. Soporte de expresión (salmon rápido)

- Reglas: `salmon_index_candidates`, `salmon_quant_candidates`, `candidates_tpm`.
- Salida: `quant/candidates_tpm.tsv` (transcritos × muestras).

### D3. Potencial codificante y TEs

- Reglas: `cpc2` (siempre), `feelnc` (si `use_feelnc`), `pfam_db` + `pfam_scan` (si `use_pfam`).
- Regla final: `lncrna_final` → `scripts/lncrna_final.py`: pasa si CPC2 noncoding ∧ FEELnc
  noncoding ∧ sin Pfam ∧ TPM ≥ 0.5 en ≥ 2 muestras. Calcula solapamiento exónico con TEs
  (`te_like` si ≥ 50 %), clasifica lincRNA/antisense/intrónico, fusiona GTFs externos con
  prefijo `EXT1_`.
- Salidas: `lncrna/lncrna.gtf`, `lncrna/lncrna_classification.tsv`, `lncrna/coding_plus_lncrna.gtf`.
- Comprobar: `logs/lncrna/final.log` imprime el embudo (candidatos → CPC2 → expresados →
  finales por tipo). Orden de magnitud esperado en 28 accesiones de roseta: 2.000-6.000
  lncRNAs nuevos; Kornienko con 499 accesiones y 5 tejidos llegó a 11.295 loci.

## Fase E. Cuantificación final

- Reglas: `final_transcriptome` (gffread + decoys + tx2gene), `salmon_index_final`,
  `salmon_quant_final` (`--validateMappings --gcBias --seqBias`), `matrices`
  (`scripts/tximport_matrices.R`).
- Salidas: `matrices/{counts,tpm}_{lncrna,coding,all}.tsv`, `matrices/metadata.tsv`.
- Comprobar: tasa de mapeo de salmon > 85 %; suma de TPM por muestra = 1e6.

## Fase F. Análisis

### F1. Normalización

- `scripts/normalize.R`: filtro de expresión (TPM ≥ 1 en ≥ 10 % de muestras), VST, corrección
  de **batch** (columna `batch`, si no `study`) y tejido con limma solo para análisis no
  supervisados, PCA.
- Salidas: `results/normalized/vst_*_batchcorrected.tsv`, `pca.pdf`, `filter_stats.tsv`.
- Comprobar en el PCA: las réplicas (lotes) de una accesión deben caer juntas tras corregir.

### F2. Estratificación

- `scripts/stratify.R`: media por accesión, 2.000 features más variables,
  ConsensusClusterPlus (PAM, Pearson, 1.000 remuestreos, 80 %), PAC por k, silhouette,
  mismo procedimiento con codificantes, ARI lncRNA vs coding y vs etiquetas externas.
- Salidas: `results/stratification/{consensus_summary,clusters,lncrna_vs_coding}.tsv`, `consensus_plots.pdf`.
- Interpretar: PAC mínimo marca k; si la matriz de consenso es limpia (bloques nítidos) y el
  ARI con los grupos de admixture es alto, los lncRNAs recapitulan la estructura genética;
  si el ARI con admixture es bajo pero el PAC es bueno, hay estructura de expresión no
  explicada por ancestría (lo interesante).

### F3. Firma

- `scripts/signature.R`: elastic-net (α = 0.5) multinomial, CV anidada 5×5, firma final.
- Salidas: `results/signature/{signature_lncrnas,nested_cv_performance}.tsv`, heatmap.
- Validación externa: cuantificar GSE80744 con el mismo GTF y aplicar la firma (script a
  añadir; con 28 accesiones de descubrimiento la firma será inestable y la cohorte de 461
  es donde se decide si existe de verdad).

### F4. Genotipo × condición

- `scripts/de_gxe.R`; solo corre si hay > 1 nivel en `condition`. En GSE224760/GSE226691 no
  hay tratamiento, así que no se ejecuta. Dataset candidato para esta rama: GSE71188 (20
  accesiones, diurno).

### F5. cis-eQTL

- Reglas: `accession_list`, `genotype_matrix` (bcftools: bialélicos, MAF ≥ 0.05, missing
  ≤ 10 %, matriz 0/1/2), `eqtl` (`scripts/eqtl.R`: MatrixEQTL, ±100 kb, 10 PCs ocultos,
  Bonferroni intra-gen + BH).
- Requiere ≥ 20 accesiones con `accession_id`. Con 28 la potencia es baja (solo efectos
  grandes); la cohorte de 461 de GSE80744 es la que da potencia real.
- Salidas: `results/eqtl/{cis_eqtl,egenes}.tsv`, `eqtl_plots.pdf`.

### F6. Informe

- `scripts/report.py` → `results/report.md`; MultiQC → `results/qc/multiqc.html`.

## Orden recomendado de ejecución real

1. `run_local.sh 1` → entorno y sample sheet de GSE224760; revisar CSV.
2. `run_local.sh 2` → referencias (sin lecturas).
3. `run_local.sh 3` → prueba con 4 muestras (2 accesiones × 2 lotes) hasta `lncrna/lncrna.gtf`.
4. `run_local.sh 4` → corrida completa de GSE224760.
5. Repetir con GSE226691; fusionar GTFs finales.
6. GSE80744 solo cuantificación + validación de firma + eQTL con 461 accesiones.
