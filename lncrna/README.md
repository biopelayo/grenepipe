# lncRNA × genoma × expresión en poblaciones de plantas

Proyecto para identificar, cuantificar y usar lncRNAs como marcadores para **estratificar
accesiones** (o casos/controles) en *Arabidopsis thaliana*, con extensión exploratoria a
*Chlamydomonas reinhardtii* y *Marchantia polymorpha*. La lógica es la misma que la de las
firmas de lncRNAs usadas para estratificar pacientes en oncología (consensus clustering sobre
expresión de lncRNAs, firmas por regresión penalizada, validación en cohortes independientes),
pero con genotipos de población en lugar de tumores.

Este directorio es independiente de grenepipe (que es un pipeline de variant calling de ADN y
no sirve para RNA-seq). grenepipe sí puede usarse en la rama de ADN del proyecto.

## 1. Pregunta y diseño

**Hipótesis central.** La expresión de lncRNAs varía mucho más entre accesiones que la de genes
codificantes (Kornienko et al. 2024 muestran que la mayoría de lncRNAs están silenciados y su
expresión es extremadamente variable, en gran parte por silenciamiento tipo transposón). Esa
variabilidad contiene estructura (geográfica, climática, de historia demográfica, de estrategia
de vida) que se puede usar para **estratificar accesiones**, y lncRNAs concretos responden a
estrés de forma dependiente del genotipo (caso/control × accesión).

**Tres capas de datos, una por accesión:**

| Capa | Fuente | Para qué |
|---|---|---|
| Genoma (SNP/indel/SV) | 1001 Genomes v3.1 VCF, 1135 accesiones | eQTL de lncRNAs, presencia/ausencia de loci, ASE, GWAS de "expresión de lncRNA" como fenotipo |
| Expresión (RNA-seq stranded) | Kornienko 2024 (499 accesiones, 5 tejidos), 1001 Epigenomes (~1000 transcriptomas) | Anotación de lncRNAs de población, matriz de expresión, clustering |
| Epigenoma (mC) | 1001 Epigenomes, 444 accesiones solapadas | Explica ~50 % de la variación de lincRNAs; covariable y mecanismo |

**Dos análisis publicables con los mismos datos:**

1. **Estratificación no supervisada.** Consensus clustering (ConsensusClusterPlus / NMF) sobre
   la matriz de expresión de lncRNAs, estabilidad por remuestreo, comparación con la estructura
   genética (admixture groups del 1001G), clima de origen (WorldClim en AraClimate) y fenotipos
   (AraPheno). Pregunta: ¿los lncRNAs separan accesiones mejor o distinto que los genes
   codificantes, y qué los explica (cis-SNP, metilación, TE)?
2. **Caso/control × genotipo.** Dataset de estrés en muchas accesiones (sequía leve en ~100
   accesiones, Clauw et al.; diurno/frío en 20 accesiones GSE71188) para lncRNAs con
   interacción genotipo × tratamiento, y firma de lncRNAs predictiva de la respuesta.

## 2. Sistema de búsqueda de datos

`find_datasets.py` hace la búsqueda **sistemática y reproducible**, no manual:

1. Descarga todos los runs RNA-seq de ENA para el taxón (ENA Portal API, campo a campo).
2. Descarga atributos de BioSample (ecotype / strain / accession) para contar genotipos.
3. Clasifica cada run con heurísticas sobre el texto de metadatos: stranded, rRNA-depleted
   vs polyA, paired, longitud de lectura, profundidad, términos de estrés/tejido; excluye
   small RNA, single-cell, degradome, CAGE, Ribo-seq.
4. Agrupa por estudio y asigna una puntuación que premia stranded + ribo-depleted + paired +
   ≥75 nt + ≥20 M lecturas + número de accesiones no-referencia + diseño con estrés.
5. Escribe `studies_ranked.tsv` y una sample sheet lista para nf-core/rnaseq por estudio.
6. Opcionalmente vuelca los runs WGS (`--with-dna`) para cruzar genotipos RNA/DNA por BioSample.

```bash
pip install pandas requests
./find_datasets.py --organism arabidopsis --out results/arabidopsis --with-dna --email tu@correo
./find_datasets.py --organism chlamy      --out results/chlamy --email tu@correo
./find_datasets.py --organism marchantia  --out results/marchantia --email tu@correo
```

Las heurísticas son imperfectas: la mitad de los submitters no describen el protocolo.
El flujo es: ranking automático → revisión manual del top 30 → confirmación de strandedness
con los datos reales (nf-core/rnaseq lo infiere con salmon/RSeQC) → tabla final curada.
`seed_datasets.tsv` contiene los datasets ya identificados en la revisión bibliográfica, con
columna `status` que indica cuáles tienen el accession verificado y cuáles hay que confirmar
leyendo el Data Availability del paper.

## 3. Criterios de inclusión (fijados antes de mirar resultados)

- RNA-seq Illumina, **stranded**, ≥ 75 nt, paired preferido, ≥ 20 M lecturas por muestra.
- rRNA-depleted / total RNA preferido (captura lncRNAs no poliadenilados y antisense);
  polyA aceptado como capa secundaria y marcado como covariable.
- ≥ 20 accesiones no-referencia con genotipo público (1001G o resecuenciado en el mismo
  estudio) o diseño caso/control con ≥ 3 réplicas biológicas por grupo.
- Metadatos suficientes para asignar tejido, estadio y condición a cada muestra.
- Excluidos: mutantes como único contraste, single-cell, small RNA, estudios sin réplicas.

## 4. Pipeline de análisis (resumen)

1. **Preprocesado y cuantificación**: nf-core/rnaseq (fastp → SortMeRNA para total RNA →
   STAR 2-pass → salmon), con GTF = Araport11 + anotación de lncRNAs de Kornienko 2024.
2. **Anotación de lncRNAs de novo por accesión**: StringTie2 guiado por genoma, merge con
   TACO, filtro de potencial codificante (CPC2 + FEELnc + PLEK, Pfam/hmmscan), longitud ≥ 200,
   exclusión de solapamiento sense con genes codificantes, clasificación lincRNA / antisense /
   intrónico, etiquetado de solapamiento con TEs (TAIR10 TE annotation). Mapeo **por accesión
   contra pseudogenomas** (SNPs del 1001G inyectados con `g2gtools` o `vcf2fasta`) para evitar
   sesgo de referencia en la cuantificación y en la ASE.
3. **Matriz y normalización**: counts por locus, TMM/VST, corrección de batch (ComBat-seq o
   covariables de estudio/tejido), filtrado de lncRNAs expresados en ≥ 10 % de accesiones.
4. **Estratificación**: ConsensusClusterPlus (PAM/hclust, 1000 remuestreos, k 2–10, PAC para
   elegir k), NMF como método alternativo, comparación con PCA de SNPs y con los grupos
   genéticos del 1001G mediante ARI; firma mínima por LASSO/elastic-net multinomial con
   validación cruzada anidada y validación externa en 1001 Epigenomes.
5. **Genética**: cis-eQTL (MatrixEQTL / tensorQTL, ±100 kb, PEER factors), GWAS de expresión
   (GEMMA/LIMIX con kinship del 1001G), ASE en heterocigotos raros, asociación con metilación
   (mC en el locus) y con solapamiento con TE.
6. **Caso/control**: DESeq2/limma-voom con modelo `~ accesión + tratamiento + accesión:tratamiento`,
   o modelos mixtos con kinship; lncRNAs con G×E significativo; co-expresión (WGCNA) para
   asignar función por "guilt by association".
7. **Validación y reproducibilidad**: Snakemake/Nextflow para todo, conda/containers
   fijados, semillas, informe Quarto, datos derivados en Zenodo.

## 5. Estado por organismo

- **Arabidopsis**: viable y de nivel publicable con datos públicos. Dependencia crítica:
  conseguir los accessions de Kornienko 2024 y de 1001 Epigenomes (bloqueados en esta sesión
  por red, ver abajo).
- **Chlamydomonas**: genoma v6 excelente y resecuenciación de decenas de aislados naturales
  (Flowers 2015, Craig 2019), pero el RNA-seq público es casi todo de cepas de laboratorio.
  Estratificación de cepas naturales requeriría generar RNA-seq propio; la anotación de
  lncRNAs condición-dependiente sí es viable.
- **Marchantia**: anotación de lncRNAs viable (Tak-1), población no viable con datos públicos.

## 6. Pendiente

- Ejecutar `find_datasets.py` desde una máquina con acceso a ENA/NCBI (en esta sesión el proxy
  bloquea ebi.ac.uk y ncbi.nlm.nih.gov).
- Rellenar `status` de `seed_datasets.tsv` leyendo los Data Availability de Kornienko 2024,
  Kawakatsu 2016 y Clauw 2016.
- Decidir tejido principal (rosette/seedling es el que tiene más accesiones).
