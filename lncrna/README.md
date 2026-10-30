# lncRNA × genoma × expresión en poblaciones de plantas

Proyecto completo para identificar lncRNAs a nivel de población, cuantificarlos por accesión y
usarlos para **estratificar accesiones** (y para análisis caso/control × genotipo) en
*Arabidopsis thaliana*, con extensión a *Chlamydomonas* y *Marchantia* donde los datos lo
permitan. La lógica es la de las firmas de lncRNAs que estratifican pacientes en oncología
(consensus clustering sobre expresión de lncRNAs, firma por regresión penalizada, validación
en cohorte independiente), aplicada a genotipos naturales con genoma conocido.

Este directorio es independiente del pipeline grenepipe (variant calling de ADN). grenepipe
sigue siendo útil para la rama de ADN si se resecuencian accesiones nuevas.

## Contenido

```
lncrna/
├── README.md                     este documento
├── seed_datasets.tsv             datasets clave ya identificados (con accessions verificados)
├── find_datasets.py              búsqueda sistemática en ENA/SRA, ranking de estudios
├── make_samplesheet_from_gse.py  GSE -> sample sheet del workflow (vía NCBI + ENA)
├── config/
│   ├── config.yaml               parámetros del workflow
│   └── samples.example.csv       formato de la sample sheet
└── workflow/                     Snakemake 8: 9 módulos de reglas, 9 envs conda, 12 scripts
    ├── Snakefile
    ├── rules/      reference, reads, align, assemble, lncrna, quant, genotypes, analysis
    ├── envs/       utils, reads, align, assembly, cpc2, feelnc, quant, genotypes, r
    └── scripts/    Python (anotación) y R (matrices, clustering, firma, G×E, eQTL)
```

## 1. Diseño

**Hipótesis.** La expresión de lncRNAs varía mucho más entre accesiones que la de genes
codificantes y está dominada por silenciamiento tipo transposón (Kornienko et al. 2024).
Esa variación tiene estructura (geografía, clima, grupos genéticos) explotable para
estratificar accesiones, y parte de ella es cis-genética (eQTL) o metilación-dependiente.

**Tres capas por accesión:**

| Capa | Fuente | Uso |
|---|---|---|
| Genoma | 1001 Genomes v3.1 VCF (1135 accesiones) | cis-eQTL de lncRNAs, presencia/ausencia, ASE |
| Expresión | Kornienko 2024 (GSE224760, GSE226691; 28 accesiones × 5 tejidos, stranded) + 1001 Epigenomes (GSE80744 / SRP074107, ~1000 transcriptomas) | anotación, matriz, clustering, validación externa |
| Epigenoma | 1001 Epigenomes (GSE43857) + Kornienko BS-seq (GSE226560) | covariable y mecanismo |

**Dos resultados publicables:**

1. **Estratificación no supervisada** de accesiones por lncRNAs: consensus clustering con
   PAC para elegir k, comparación con el clustering por genes codificantes (ARI) y con
   etiquetas externas (grupos de admixture del 1001G, clima de origen de AraClimate,
   fenotipos de AraPheno), firma mínima de lncRNAs por elastic-net con CV anidada, validación
   en la cohorte 1001 Epigenomes.
2. **Caso/control × genotipo**: en datasets con tratamiento (sequía, frío, diurno) en varias
   accesiones, lncRNAs con interacción genotipo × tratamiento (DESeq2 LRT) y eQTL.

## 2. Búsqueda sistemática de datos

```bash
pip install pandas requests
./find_datasets.py --organism arabidopsis --out results/arabidopsis --with-dna --email tu@correo
./find_datasets.py --organism chlamy      --out results/chlamy --email tu@correo
./find_datasets.py --organism marchantia  --out results/marchantia --email tu@correo
```

Baja todos los runs RNA-seq del taxón en ENA, los atributos de BioSample (ecotype/strain),
clasifica cada run (stranded, rRNA-depleted vs polyA, paired, longitud, profundidad,
estrés/tejido; excluye sRNA, single-cell, degradome/CAGE/Ribo-seq), puntúa por estudio
(premia stranded + ribo-depleted + paired + ≥75 nt + ≥20 M + nº de accesiones no-referencia)
y escribe `studies_ranked.tsv` más una sample sheet por estudio. El flujo es ranking
automático → revisión manual del top 30 → confirmación de strandedness con datos reales (el
workflow la infiere) → tabla curada.

Para los datasets ya identificados (ver `seed_datasets.tsv`):

```bash
./make_samplesheet_from_gse.py --gse GSE224760 GSE226691 --out config/samples.csv --email tu@correo
# revisar a mano: accession, accession_id (id numérico 1001G), tissue, condition
```

**Criterios de inclusión fijados a priori**: Illumina stranded, ≥75 nt, ≥20 M lecturas;
rRNA-depleted preferido (polyA aceptado como covariable); ≥20 accesiones con genotipo público
o diseño caso/control con ≥3 réplicas; metadatos de tejido/estadio/condición; excluidos
mutantes como único contraste, single-cell, small RNA, estudios sin réplicas.

## 3. Workflow

```bash
conda create -n lncrna -c conda-forge -c bioconda snakemake=8 mamba
conda activate lncrna
cd lncrna
snakemake --use-conda --conda-frontend mamba --cores 32 \
    --directory /ruta/analisis --configfile config/config.yaml
# Verificar el plan antes: añadir -n
```

Pasos (todos con envs conda fijados, descargas automáticas de genoma TAIR10 / Araport11
de Ensembl Plants, TEs de TAIR, VCF del 1001G, base de SortMeRNA):

1. **Lecturas**: descarga desde ENA si son URLs, fastp, SortMeRNA (opcional, para total RNA).
2. **Strandedness** inferida por muestra con salmon sobre el transcriptoma de referencia
   (o declarada en la sample sheet).
3. **STAR 2-pass** con índice de referencia.
4. **StringTie** por muestra (guiado por Araport11, flag de strand por muestra),
   `stringtie --merge`, `gffcompare`.
5. **Candidatos lncRNA**: class codes u/x/i, ≥200 nt, cromosomas nucleares; cuantificación
   rápida con salmon para soporte de expresión; **CPC2** (y opcionalmente FEELnc y Pfam);
   solapamiento con TEs (bedtools) → `lncrna/lncrna.gtf` + tabla de clasificación
   (lincRNA / antisense / intrónico, TE-like, novel). Se pueden fusionar GTFs externos
   (p. ej. el de Kornienko 2024) vía `reference.extra_lncrna_gtf`.
6. **Cuantificación final** con salmon (decoys genómicos, GC/seq bias) sobre codificantes +
   lncRNAs; `tximport` → matrices de counts/TPM separadas (lncRNA, coding, all) + metadatos.
7. **Normalización**: filtro de expresión, VST, corrección de batch/tejido (limma) solo para
   análisis no supervisados, PCA.
8. **Estratificación** (`stratify.R`): media por accesión, top-2000 variables,
   ConsensusClusterPlus (PAM, Pearson, 1000 remuestreos), PAC, silhouette, k óptimo,
   mismo procedimiento con genes codificantes, ARI entre ambos y contra etiquetas externas.
9. **Firma** (`signature.R`): glmnet elastic-net (bi/multinomial) con CV anidada 5×5,
   accuracy por fold, firma final con coeficientes y tipo de lncRNA, heatmap.
10. **G × E** (`de_gxe.R`, si hay >1 condición): DESeq2, efecto de condición y LRT de
    `accession:condition`.
11. **cis-eQTL** (`eqtl.R`, si ≥20 accesiones con id 1001G): bcftools subset/filtrado
    (MAF, missingness) → matriz 012; MatrixEQTL ±100 kb con 10 PCs ocultos; eGenes con
    Bonferroni intra-gen + BH.
12. **MultiQC** y `results/report.md` con el resumen de todo.

Estado de validación: el DAG completo (271 jobs con 30 muestras de prueba, todas las ramas
activas) se construye correctamente en dry-run con Snakemake 8.15; los scripts Python
compilan; los scripts R no se han podido ejecutar en este entorno (sin R ni datos). Es
esperable algún ajuste menor en la primera ejecución real.

## 4. Estado por organismo

- **Arabidopsis**: viable y publicable con datos públicos. Core: Kornienko 2024 (discovery)
  + 1001 Epigenomes (validación) + 1001 Genomes (genética).
- **Chlamydomonas**: genoma v6.1 y resecuenciación de aislados naturales existen, pero el
  RNA-seq público es casi todo de cepas de laboratorio. El workflow sirve tal cual cambiando
  referencias en `config.yaml` (fasta/GTF de Phytozome; desactivar eQTL/TE); estratificar
  cepas naturales requeriría generar RNA-seq propio.
- **Marchantia**: anotación de lncRNAs en Tak-1 viable; población no viable con datos públicos.

## 5. Pendiente (requiere red fuera de este entorno)

- Ejecutar `find_datasets.py` y `make_samplesheet_from_gse.py`.
- Bajar el GTF de lncRNAs de Kornienko 2024 (suplemento de GSE224760) y ponerlo en
  `reference.extra_lncrna_gtf`.
- Confirmar tipo de librería de GSE80744 (polyA vs total) en los registros GSM.
- Primera ejecución real sobre GSE224760 (28 accesiones, rosette) y ajustar.
