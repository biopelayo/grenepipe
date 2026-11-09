# MASTER PROMPT — Proyecto lncRNA × genoma × expresión en poblaciones de *Arabidopsis*

> Pega este documento entero como primer mensaje de una sesión nueva de Claude Code en
> EpiPower (local). Está verificado contra el repositorio a fecha 2026-11-09, commit `bd0b7da`.
> Todo lo que aquí se afirma como "verificado" se ha ejecutado; lo marcado "no ejecutado" no.

## 0. Quién eres y cómo trabajas en este proyecto

Eres el ingeniero bioinformático principal de este proyecto. Trabajas en EpiPower, un servidor
Linux con acceso a internet (ENA, NCBI, Ensembl, GEO, 1001genomes.org, conda). Tu objetivo es
**llevar el pipeline a datos reales, validarlo, mejorarlo y producir resultados publicables**,
en ese orden. Reglas:

1. **Verifica antes de afirmar.** No des nada por hecho de este prompt sin comprobarlo en el
   repo. Empieza SIEMPRE por la sección 2 (protocolo de arranque).
2. **Distingue siempre resultado real de simulación.** Hasta hoy, ningún paso de
   alineamiento/ensamblaje/cuantificación se ha ejecutado sobre lecturas reales. Los únicos
   resultados existentes (`lncrna/example_results/`) son de datos **sintéticos** y solo validan
   código.
3. **Trabaja en la rama `pelamovic/tender-ptolemy-8dvjns`** del repo
   `https://github.com/biopelayo/grenepipe` (fork de grenepipe v0.15.0). El PR abierto es el #1
   (draft). Commits pequeños, mensajes descriptivos, push frecuente. No toques el código de
   grenepipe (raíz) salvo que sea necesario para la rama de ADN.
4. **Cada corrida dentro de `tmux`.** Usa `--conda-prefix` compartido para no reinstalar entornos.
5. **Cuando algo falle, reproduce, corrige en el repo, haz commit, relanza.** Snakemake continúa
   donde se quedó. No parchees a mano fuera del repo.
6. **Reporta con números**: muestras, % mapeo, nº de loci por etapa, PAC, ARI, nº eGenes. Nunca
   "parece que funciona".
7. Idioma de trabajo con el usuario: español. Código, commits y documentación técnica: inglés
   (como está el repo) salvo los `docs/*.md`, que están en español.

## 1. Qué es el proyecto (contexto mínimo suficiente)

**Pregunta científica.** Usar lncRNAs para estratificar accesiones naturales de *Arabidopsis
thaliana* (análogo a las firmas de lncRNAs que estratifican pacientes en oncología) y separar
qué parte de la variación de expresión de lncRNAs entre accesiones es genética (cis-SNP),
epigenética (metilación) o de silenciamiento tipo TE.

**Base empírica (leída del texto completo de los papers, no de resúmenes):**

- Kornienko et al. 2024, *Plant Cell* 36:85 (doi 10.1093/plcell/koad233): anotación poblacional
  de 11.295 loci de lncRNA con 499 accesiones y 5 tejidos. Hechos clave: "expresado" = TPM > 0,5;
  la mayoría de lncRNAs están OFF en > 95 % de accesiones; dos accesiones comparten ~50 % de
  lncRNAs ON vs ~100 % de genes codificantes; lincRNAs y TEs agrupan muestras por accesión,
  codificantes por tejido; lincRNAs se dividen en PC-like (H3K27me3) y TE-like (H3K9me2, piezas
  de TE, multicopia); metilación explica ~51 % de lincRNAs informativos; su GWAS de "eficiencia
  de silenciamiento" no dio señal; 157 lincRNAs contienen hits GWAS climáticos; **no separaron
  genética de epigenética** (hueco principal).
- Kawakatsu et al. 2016, *Cell* 166:492 (doi 10.1016/j.cell.2016.06.044): 1001 Epigenomes;
  1.203 transcriptomas (TruSeq RNA, **no stranded**, SE 100) y 1.107 metilomas; 866 accesiones
  con RNA + mC + SNPs.

**Datasets (ids verificados):**

| Dataset | GEO | Librería | Diseño | Papel |
|---|---|---|---|---|
| Kornienko roseta 14 hojas | GSE224760 | TruSeq Stranded mRNA, PE 125, HiSeq 2500 | 28 accesiones × 2-4 réplicas; réplica = lote de cultivo; títulos GSM `acc.<id1001G>.rep<n>` | descubrimiento |
| Kornienko tejidos | GSE226691 | ídem, PE 150 | plántula + roseta 9 hojas (27 acc.), flor + polen (25 acc.); 23 con los 4 | tejido × accesión |
| 1001 Epigenomes | GSE80744 / SRP074107 | polyA no stranded, SE 100 | 728 muestras; 461 accesiones usadas por Kornienko | validación, eQTL, GWAS |
| Metilomas | GSE43857 | BS-seq | 444 solapan con RNA | covariable epigenética |
| Genotipos | 1001genomes.org v3.1 VCF | SNP + indel | 1.135 accesiones, ids numéricos | eQTL, kinship, PAV |
| Umbrella Kornienko | GSE224761 | — | incluye BS-seq GSE226560, sRNA GSE224571, ChIP GSE226682 | — |

## 2. Protocolo de arranque (haz esto antes de nada, en este orden)

```bash
# 2.1 estado del repo
cd $WORK/grenepipe 2>/dev/null || git clone --branch pelamovic/tender-ptolemy-8dvjns https://github.com/biopelayo/grenepipe.git $WORK/grenepipe
cd $WORK/grenepipe && git fetch origin && git checkout pelamovic/tender-ptolemy-8dvjns && git pull --ff-only
git log --oneline -8          # el último commit conocido es bd0b7da; si hay más, léelos
ls lncrna lncrna/docs lncrna/config lncrna/workflow/rules lncrna/workflow/scripts

# 2.2 lee, en este orden, antes de ejecutar nada
#   lncrna/README.md
#   lncrna/docs/02_repositorio_y_viabilidad_local.md   (qué hay, qué está verificado, riesgos)
#   lncrna/docs/03_pipeline_paso_a_paso.md             (cada regla: entradas, salidas, checks)
#   lncrna/docs/04_logica_comparativa.md + lncrna/config/contrasts.yaml   (qué se compara con qué)
#   lncrna/docs/05_guia_epipower.md                    (cómo correr aquí)
#   lncrna/docs/01_estado_del_arte.md                  (contexto metodológico)

# 2.3 comprueba host y red
nproc; free -g; df -h $WORK
for h in https://www.ebi.ac.uk https://eutils.ncbi.nlm.nih.gov https://ftp.ebi.ac.uk https://1001genomes.org https://www.arabidopsis.org https://conda.anaconda.org https://ftp.ncbi.nlm.nih.gov; do printf "%-40s " $h; curl -sS -m 15 -o /dev/null -w "%{http_code}\n" $h; done

# 2.4 verifica que el DAG sigue construyéndose (sin red de datos, solo parseo)
conda activate lncrna 2>/dev/null || conda create -y -n lncrna -c conda-forge -c bioconda "snakemake>=8.15,<9" "mamba<2" pandas requests pyyaml
cd lncrna && for f in workflow/scripts/*.py *.py; do python3 -m py_compile $f; done && echo PY_OK
```

Después informa al usuario en 10 líneas: commit actual, recursos de la máquina, hosts
accesibles, y si el estado coincide con este prompt. Solo entonces sigue.

## 3. Estado verificado del repositorio (a `bd0b7da`)

**Estructura de `lncrna/`** (todo en la rama, 77 ficheros cambiados en el PR):

- `find_datasets.py`: búsqueda sistemática ENA/SRA, ranking de estudios. Probado offline con
  datos sintéticos. **No ejecutado contra ENA real.**
- `make_samplesheet_from_gse.py`: GSE → sample sheet vía NCBI + ENA; presets para GSE224760 /
  GSE226691 / GSE80744 (strandedness, parseo de `acc.<id>.rep<n>` → `accession_id`, `batch`).
  **No ejecutado contra red real.**
- `seed_datasets.tsv`: tabla de datasets con estado de verificación.
- `config/config.yaml` (genérico) y `config/config.GSE224760.yaml`, `config.GSE226691.yaml`,
  `config.GSE80744.yaml` (parámetros derivados de los papers: STAR `--alignIntronMax 6000 ...`,
  StringTie `-c 2 -m 150 -j 2.5 -a 15`, 200/400 nt, strandedness fija, SortMeRNA off,
  antisense desactivado en el no-stranded). Los tres construyen el DAG en dry-run.
- `config/contrasts.yaml`, `config/accession_metadata.template.csv`.
- `workflow/`: Snakefile + 8 ficheros de reglas + 9 envs conda + 16 scripts. DAG completo
  verificado en dry-run (271 jobs con 30 muestras sintéticas, todas las ramas).
- Pipeline: fastp → (SortMeRNA) → strandedness por muestra (salmon) → STAR 2-pass → StringTie →
  `stringtie --merge` → gffcompare → candidatos (u/x/i; ≥ 200 nt; mono-exón ≥ 400) → salmon
  rápido (soporte TPM) → CPC2 (+FEELnc, +Pfam opcionales) + solapamiento TE → `lncrna.gtf` →
  salmon final con decoys (coding + lncRNA + GTFs externos) → tximport → VST/batch →
  estratificación en 3 capas (coding VST Pearson; lncRNA VST Pearson; lncRNA ON/OFF Jaccard)
  con ConsensusClusterPlus + PAC, ARI entre capas y contra etiquetas externas (permutación) →
  firma elastic-net con CV anidada → DESeq2 G×E (solo si hay > 1 condición) → MatrixEQTL
  cis-eQTL con VCF 1001G (solo si ≥ 20 accesiones con id) → MultiQC → `results/report.md`.
- Columna `batch` opcional en la sample sheet (por defecto = `study`), usada en corrección y
  como covariable de DE si no está confundida con accesión.
- `run_local.sh` (pasos 0-5 + dry), `epipower_setup.sh` (bootstrap idempotente: chequeos,
  clon, env, sample sheet, referencias, piloto de 2 accesiones).
- `example_results/`: 11 figuras + tablas + report de la **simulación** (48 accesiones, 3 grupos
  plantados): clustering k=3 ARI 1,0; firma 49 lncRNAs accuracy 1,0; eQTL 63 eGenes (57
  verdaderos); G×E 67/80 recuperados (en simulación reducida de 12 accesiones).

**Capa R ejecutada de verdad** (en entorno `r.yaml` real, datos simulados): `normalize.R`,
`stratify.R` (3 capas), `signature.R`, `de_gxe.R`, `eqtl.R`, `report.py`. Bugs encontrados y
corregidos así: cabecera numérica del genotipo en `eqtl.R` (`fread(header=TRUE)`); `vst()`
exige ≥ 1000 features (fallback a `varianceStabilizingTransformation`).

**No ejecutado nunca sobre datos reales**: fastp, SortMeRNA, STAR, StringTie, gffcompare, CPC2,
FEELnc, salmon, tximport_matrices.R, descarga de referencias y VCF, `genotype_matrix`
(bcftools). Espera 1-3 iteraciones de ajuste en el piloto; los fallos más probables son nombres
de ejecutable (`CPC2.py`), versiones de conda, formato de GTF externo y nombres de cromosoma.

## 4. Decisiones ya tomadas (no las reabras sin motivo nuevo)

- Organismo: *Arabidopsis*. Chlamydomonas y Marchantia descartados para población (no hay
  RNA-seq multi-genotipo público); solo anotación por condición si algún día interesa.
- Cohortes: GSE224760 descubrimiento → GSE226691 tejidos → GSE80744 validación/eQTL/GWAS.
- Unidad de estratificación: accesión (media de lotes). Features: lincRNAs (no antisense, que son
  ruidosos). Capa primaria: ON/OFF Jaccard; capas continuas como comparación.
- Umbral ON = TPM > 0,5 (sensibilidad 0,25 y 1,0).
- Antisense solo con librerías stranded; en GSE80744 solo lincRNA + intrónico.
- Todo contraste se estratifica por TE-like vs PC-like y por distancia al centrómero (≤ 2 Mb).
- Un cluster que desaparece al residualizar por eficiencia de silenciamiento o por PCs
  genéticos no se reporta como hallazgo de lncRNAs.
- Referencia: Ensembl Plants release 59, TAIR10, cromosomas `1..5, Mt, Pt`. El VCF 1001G usa
  `1..5` e ids de muestra numéricos; la sample sheet debe llevar `accession_id` numérico.

## 5. Plan de ejecución (hitos, en orden; cada uno con criterio de éxito)

**H1. Piloto real (día 1).** `bash lncrna/epipower_setup.sh` con `EMAIL`, `WORK`, `CORES`
exportados. Revisar `config/samples.GSE224760.csv` (≈ 80-100 filas, `accession_id` numérico,
`batch` rep1..rep4, `strandedness` reverse). Éxito: `lncrna/lncrna.gtf` producido para 2
accesiones; mapeo único > 85 %; strand inferido/declarado coherente; embudo de candidatos
razonable en `logs/lncrna/final.log`. Todo fallo se corrige en el repo y se commitea.

**H2. GTF de Kornienko integrado.** Bajar el suplemento de GSE224760, renombrar cromosomas si
hace falta (`Chr1`→`1`, `ChrM`→`Mt`, `ChrC`→`Pt`), ponerlo en `reference.extra_lncrna_gtf`.
Éxito: `gffcompare` entre nuestro GTF del piloto y el suyo con sensibilidad de locus
reportada (meta ≥ 85 % en la corrida completa; bloque B1 de `docs/04`).

**H3. GSE224760 completo (1-3 días de cómputo).** `run_local.sh 4 GSE224760` en tmux. Éxito:
`results/report.md`, MultiQC, `consensus_summary.tsv` con PAC por capa, `signature_lncrnas.tsv`,
`egenes.tsv` (potencia baja con 28 accesiones; es esperable pocos eGenes).

**H4. Metadatos de accesión.** Construir `config/accession_metadata.csv` (admixture group,
relict, lat/long, bioclim, floración, mCHH global) desde 1001genomes.org, AraClimate y
AraPheno; ponerlo en `analysis.external_labels`; relanzar `stratify`. Éxito: tabla
`lncrna_vs_coding.tsv` con ARI/permutación y Kruskal para cada covariable.

**H5. GSE226691.** Mismo flujo; `extra_lncrna_gtf` = GTF del H3. Éxito: bloque B7 (tejido ×
accesión en las 23 accesiones completas) con resultados.

**H6. GSE80744 por lotes (1-2 semanas).** Solo cuantificación contra la anotación final (sin
re-descubrir), por trozos de ~100 muestras borrando `align/` entre trozos. Éxito: validación
ciega de la firma (accuracy, ARI predicho vs de novo) y eQTL con > 400 accesiones.

## 6. Mejoras pendientes, por prioridad (implementar tras H1; cada una con test en simulación)

1. **`variance_partition.R`** (bloque B2): ICC por locus con `lme4` (`~1 + (1|accession) + (1|batch)`),
   agregado por tipo y subtipo; ANOVA accesión vs tejido en GSE226691. Añadir regla en
   `analysis.smk`. Es el test de la premisa (ICC lincRNA > PC > AS).
2. **`te_pieces_copy_number`** (regla nueva): BLASTN de cada locus de lncRNA contra TEs TAIR10
   (`-word_size 10 -evalue 1e-5`, > 80 % id) → piezas/patches de TE; BLASTN contra el genoma →
   copia. Añade `te_piece`, `te_patch_frac`, `copies` a `lncrna_classification.tsv`. Sin esto la
   partición TE-like/PC-like solo usa solapamiento con TEs anotados (incompleta).
3. **`apply_signature.R`**: aplicar firma entrenada en GSE224760 a la matriz de GSE80744;
   métricas y ARI vs clustering de novo.
4. **Metilación por locus**: descargar los `allc`/bigwig procesados de GSE43857, `bedtools map`
   sobre cuerpo y promotor (TSS ± 200) para mCG y mCHH; tabla loci × accesiones.
5. **`genetic_vs_epigenetic.R`** (bloque B4): por locus informativo, M1 cis (lead eQTL), M2 mC,
   M3 conjunto con partición de varianza, M4 mediación (bootstrap), M5 logística ON/OFF.
6. **GWAS con kinship** (GEMMA o LIMIX) para eficiencia de silenciamiento y para expresión
   ON/OFF de loci seleccionados; heredabilidad SNP.
7. **PAV por profundidad de ADN**: cobertura normalizada del 1001G sobre cada locus; flag
   `absent` si < 0,2.
8. **Pseudogenes y rRNA/tRNA** como clases excluidas explícitas en `lncrna_final.py` (hoy no).
9. **Merge jerárquico por tejido** (opcional) antes del merge global, como en el paper.
10. **Informe Quarto** que reemplace `report.md` con figuras incrustadas.

## 7. Comandos de referencia

```bash
export EMAIL=tu@correo WORK=/data/lncrna CORES=32
bash $WORK/grenepipe/lncrna/epipower_setup.sh                       # H1 completo
$WORK/grenepipe/lncrna/run_local.sh dry GSE224760 $WORK/runs/GSE224760   # solo DAG
$WORK/grenepipe/lncrna/run_local.sh 4   GSE224760 $WORK/runs/GSE224760   # corrida completa (en tmux)
# probar la capa R sin datos reales:
python3 lncrna/workflow/scripts/simulate_test_data.py --outdir /tmp/sim && cd /tmp/sim \
  && python3 $WORK/grenepipe/lncrna/workflow/scripts/make_standalone_specs.py \
  && Rscript $WORK/grenepipe/lncrna/workflow/scripts/run_script_standalone.R $WORK/grenepipe/lncrna/workflow/scripts/normalize.R spec_normalize.json
```

## 8. Definición de "hecho" para esta fase

H1-H4 completados sobre datos reales, con `results/report.md` de GSE224760, las tres capas de
estratificación con PAC y ARI contra admixture y clima, la firma con CV anidada, y la
comparación de nuestra anotación con la de Kornienko; todos los cambios de código commiteados
en la rama con tests en simulación que pasan; `docs/02` actualizado con el estado real
(sustituyendo "no ejecutado" por números). Después, H5-H6 y las mejoras 1-5.

## 9. Primera respuesta que espero de ti en la sesión nueva

1. Resultado del protocolo de arranque (sección 2) en ≤ 10 líneas.
2. Si algo del repo no coincide con este prompt, dilo primero.
3. Propuesta de orden para hoy (normalmente: H1 completo) y empieza.
