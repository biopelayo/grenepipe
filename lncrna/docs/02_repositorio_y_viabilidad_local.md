# Qué hace este repositorio y cómo ejecutarlo en local

## 1. Dos cosas distintas en el mismo repo

### 1a. `grenepipe` (raíz del repo): variant calling de ADN

Fork sin modificar de `moiexpositoalonsolab/grenepipe` v0.15.0 (Czech & Exposito-Alonso 2022).
Workflow de Snakemake 8.15.2 que va de fastq de ADN a VCF filtrado + MultiQC, con elección de
herramienta por etapa (6 trimmers, 4 mapeadores, 3 callers, 3 filtros, anotación snpEff/VEP,
HAF-pipe para pool-seq, mapDamage para ADN antiguo). Verificado en esta sesión con dry-runs
de 9 configuraciones: todas construyen el DAG.

**No sirve para RNA-seq ni para lncRNAs.** No tiene alineador splice-aware, ni ensamblaje de
transcritos, ni cuantificación. Su papel en este proyecto es nulo mientras usemos el VCF
público del 1001 Genomes. Solo haría falta si resecuenciaras accesiones nuevas (p. ej.
accesiones sin genotipo en el 1001G) para obtener sus SNPs.

Requisitos para ejecutarlo en local, por si lo necesitas:

- conda/mamba con `workflow/envs/grenepipe.yaml` (Snakemake 8.15.2, mamba < 2).
- Correr **dentro** del env conda (lee `CONDA_DEFAULT_ENV`; fuera revienta).
- Rutas absolutas en `config.yaml` y `samples.tsv` (`example/prepare.sh` lo hace).
- Comando: `snakemake --use-conda --conda-frontend mamba --cores N --directory <carpeta>`.

### 1b. `lncrna/`: el proyecto de lncRNAs (nuevo, en el PR #1)

Workflow independiente de Snakemake 8 que implementa el estándar descrito en
`01_estado_del_arte.md` más la capa de análisis poblacional:

```
fastq ──fastp──(SortMeRNA)──► strand por muestra ──► STAR 2-pass ──► StringTie ──► merge ──► gffcompare
     ──► candidatos (u/x/i, ≥200 nt) ──► salmon rápido (soporte) ──► CPC2 (+FEELnc,+Pfam) + TE ──► lncrna.gtf
     ──► salmon final (coding+lncRNA, decoys) ──► tximport ──► matrices
     ──► VST/batch ──► consensus clustering + PAC ──► firma elastic-net ──► DESeq2 G×E ──► MatrixEQTL ──► report
```

Estado de verificación:

| Capa | Verificado | Cómo |
|---|---|---|
| DAG completo | sí | dry-run 271 jobs, 30 muestras sintéticas, todas las ramas |
| Scripts Python | sintaxis | `py_compile`; `find_datasets.py` probado con datos ENA sintéticos |
| Scripts R | **ejecutados** | entorno `r.yaml` real; datos simulados con verdad conocida; recuperan clusters (ARI 1.0), firma, eQTL (57/63 verdaderos), G×E (67/80) |
| Alineamiento/ensamblaje/cuantificación | **no ejecutados** | sin acceso de red a ENA/GEO/Ensembl en este entorno |

Es decir: la parte que más riesgo tiene de fallar en la primera ejecución real es la de
herramientas externas (STAR, StringTie, CPC2, salmon), no por lógica sino por detalles de
versiones de conda, nombres de ficheros de salida y formato de GTF. Espera 1-3 iteraciones.

## 2. Viabilidad en local con los datos reales

### Datos (metadatos leídos del texto completo de los papers)

| Dataset | GEO | Librería | Diseño | Tamaño estimado |
|---|---|---|---|---|
| Kornienko roseta 14 hojas | GSE224760 | TruSeq Stranded mRNA, PE 125, HiSeq 2500 | 28 accesiones × 2-4 réplicas (lotes de cultivo separados semanas) ≈ 80-100 muestras; títulos `acc.<id1001G>.rep<n>` | ~25-40 M pares/muestra → 2-4 GB/muestra → **200-350 GB** fastq.gz |
| Kornienko tejidos | GSE226691 | TruSeq Stranded mRNA, PE 150 | plántula y roseta 9 hojas (27 acc.), flor y polen (25 acc.); ~104 muestras | **300-400 GB** |
| 1001 Epigenomes roseta | GSE80744 / SRP074107 | TruSeq RNA (polyA, **no stranded**), SE 100 | 461 accesiones usadas por Kornienko (728 muestras en GEO; 727 a 22 °C) | ~1-2 GB/muestra → **0.7-1.4 TB** |
| 1001 Genomes VCF | 1001genomes.org v3.1 | SNP + indel cortos | 1.135 accesiones | 1.3 GB comprimido |

### Recursos necesarios (orientativo, TAIR10 = 135 Mb)

| Paso | RAM | CPU | Tiempo por muestra (PE 125, 30 M pares, 8 hilos) |
|---|---|---|---|
| fastp | 2 GB | 8 | 5 min |
| STAR índice TAIR10 | 8 GB | 16 | 10 min (una vez) |
| STAR 2-pass | 8-10 GB | 8-16 | 10-20 min |
| StringTie | 2 GB | 8 | 5 min |
| salmon | 4 GB | 8 | 5 min |
| CPC2 / FEELnc | 4 GB | 8 | 30-90 min total |
| R (clustering, eQTL) | 8-16 GB | 4-8 | minutos a 1 h |
| DESeq2 LRT G×E con 28 niveles | 8 GB | 8 | 30-60 min total |

Con una estación de **16 núcleos / 64 GB / 2 TB**: GSE224760 completo en ~2 días de cómputo
(descarga aparte), GSE226691 en ~3 días, GSE80744 en 1-2 semanas (y necesita 2 TB solo para
fastq; se puede procesar por lotes y borrar).

**Recomendación**: empezar con GSE224760 (descubrimiento), después GSE226691 (tejidos),
y GSE80744 solo como cohorte de validación cuantificando con salmon directamente contra el
GTF final (sin STAR ni ensamblaje), lo que reduce su coste a ~1/5.

### Qué hay que ajustar en los ficheros (y ya está ajustado en `config/config.GSE*.yaml`)

| Parámetro | Por qué | Valor |
|---|---|---|
| `strandedness` | TruSeq Stranded mRNA es dUTP → lectura 1 antisentido | `reverse` en GSE224760/GSE226691; `unstranded` en GSE80744 |
| `params.sortmerna` | librerías polyA, no total RNA | `false` (ahorra 30 % de tiempo) |
| `params.star_sjdb_overhang` | longitud de lectura − 1 | 124 / 149 / 99 |
| `params.star_extra` | parámetros del paper | `--alignIntronMax 6000 --alignMatesGapMax 6000 --outFilterIntronMotifs RemoveNoncanonical --outFilterMultimapNmax 10` |
| `params.stringtie_assemble` | parámetros del paper | `-c 2 -m 150 -j 2.5 -a 15` |
| `lncrna_filter.class_codes` | antisense solo fiable con stranded | `[u,x,i]` stranded; `[u,i]` para GSE80744 |
| `lncrna_filter.min_exons_intergenic` + longitud mono-exón | ruido en mono-exónicos | mono-exónicos ≥ 400 nt (nueva opción `min_length_monoexonic`) |
| `batch` (columna nueva) | réplicas de GSE224760 son lotes de cultivo | `rep1..rep4` desde el título GSM |
| `accession_id` | eQTL necesita el id numérico 1001G | parseado del título `acc.9888.rep4` → `9888` |
| `reference.extra_lncrna_gtf` | fusionar la anotación publicada | GTF del suplemento de GSE224760 (descargar a mano) |
| `genotypes.vcf` | 1001G v3.1 | descarga automática |
| `analysis.external_labels` | comparar clusters con grupos genéticos | CSV `accession_id,label` de admixture groups (1001genomes.org, tabla de accesiones) |

### Riesgos conocidos antes de la primera ejecución

1. **Nombres de cromosomas**: Ensembl usa `1..5, Mt, Pt`; el GTF de Kornienko puede usar
   `Chr1..Chr5, ChrM, ChrC` (TAIR). Hay que renombrar con `sed` antes de ponerlo en
   `extra_lncrna_gtf`. El VCF del 1001G usa `1..5`.
2. **Nombres de muestras del VCF**: numéricos (`6909`). La sample sheet debe llevar el mismo
   número en `accession_id`.
3. **CPC2 en bioconda**: el ejecutable se llama `CPC2.py`; si la versión cambia el nombre,
   ajustar `rules/lncrna.smk`.
4. **FEELnc**: lento y con dependencias Perl; déjalo en `false` en la primera pasada.
5. **Descarga de SRA**: los enlaces ENA (`ftp.sra.ebi.ac.uk`) son los más fiables; si ENA no
   tiene el fastq, usar `fasterq-dump` de sra-tools (ver `run_local.sh`).
6. **Un fastq corrupto** para la descarga; el workflow reintenta 3 veces y después para.
   Reanudar con el mismo comando.
