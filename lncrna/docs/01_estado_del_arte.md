# Estado del arte: bioinformática de lncRNAs con DNA-seq y RNA-seq (2024-2026)

Evaluación paso a paso de qué se considera hoy práctica correcta, qué herramientas dominan,
dónde están los límites, y cómo encaja todo esto con el proyecto en *Arabidopsis*.
Fuentes principales: Kornienko et al. 2024 *Plant Cell* ([DOI](https://doi.org/10.1093/plcell/koad233), texto completo vía PubMed/PMC),
Kawakatsu et al. 2016 *Cell* ([DOI](https://doi.org/10.1016/j.cell.2016.06.044), texto completo vía PubMed/PMC),
Plant-LncPipe 2024 ([Hortic Res](https://academic.oup.com/hr/article/11/4/uhae041/7604008)),
LncRAnalyzer 2025 ([Plant J](https://onlinelibrary.wiley.com/doi/10.1111/tpj.70509)),
PLncFire 2025 ([bioRxiv](https://www.biorxiv.org/content/10.1101/2025.10.13.682253v1.full.pdf)),
pangenoma de 32 accesiones ([Nat Commun 2023](https://www.nature.com/articles/s41467-023-42029-4)),
variación de TEs en 1.114 accesiones ([PMC10980334](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10980334/)).

## Paso 0. Qué es un lncRNA operativamente

Transcrito ≥ 200 nt sin potencial codificante detectable, poliadenilado o no, nuclear o
citoplásmico. La definición es **operacional y depende del protocolo de librería**: una
librería polyA solo ve lncRNAs poliadenilados; una librería total-RNA con depleción de rRNA ve
también los no-polyA (muchos antisense, enhancer-RNAs, transcritos de TEs). Sin librería
*stranded* no se puede distinguir antisense de sense y la clase "antisense" es ruido.

Consecuencia directa para nosotros: los datos de Kornienko 2024 son **TruSeq Stranded mRNA**
(polyA, stranded, PE 125/150 nt, HiSeq 2500) y los de 1001 Epigenomes son **TruSeq RNA**
(polyA, **no stranded**, SE 100 nt). Así que: (a) el universo de lncRNAs al alcance es el
poliadenilado; (b) antisense solo es fiable en el dataset de Kornienko; (c) en 1001 Epigenomes
solo deben cuantificarse lincRNAs e intrónicos.

## Paso 1. Datos de partida: qué pide hoy una anotación poblacional seria

1. **Muchos genotipos, pocas réplicas** descubre más loci que pocos genotipos con muchas
   réplicas. Kornienko 2024 lo demuestra con curvas de saturación: 461 accesiones en roseta
   descubren muchos más loci que 12 muestras × 14 réplicas de Col-0 con el mismo número de
   lecturas. Los lncRNAs son accesión-específicos; muestrear accesiones es muestrear loci.
2. **Varios tejidos**: 23 accesiones × 4 tejidos (plántula, roseta, flor, polen); el polen y
   la flor aportan loci exclusivos.
3. **Stranded** obligatorio para antisense; **PE ≥ 100 nt** para ensamblar isoformas.
4. **Genotipos de la misma población**: en Arabidopsis el VCF del 1001 Genomes (1.135
   accesiones) permite eQTL, presencia/ausencia y análisis alelo-específico sin secuenciar ADN.
5. **Epigenoma emparejado**: en Arabidopsis existe para 444 accesiones (metilación) y explica
   ~51 % de la variación de lincRNAs y ~22 % de antisense (Kornienko 2024). Ningún otro
   sistema vegetal tiene esto.

## Paso 2. Preprocesado y alineamiento

Consenso actual: fastp (adaptadores, calidad) → alineador splice-aware **STAR** (2-pass para
splice junctions nuevas) o HISAT2 → BAM ordenado. Para total-RNA se añade SortMeRNA para
quitar rRNA residual. Kornienko usó STAR con `--alignIntronMax 6000 --alignMatesGapMax 6000
--outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonical`; los intrones de
Arabidopsis son cortos y limitar `alignIntronMax` evita quimeras.

Punto no resuelto: **sesgo de referencia**. Mapear 500 accesiones contra Col-0 pierde lecturas
en regiones divergentes y en loci ausentes del genoma de referencia, justo donde abundan
TEs y lncRNAs. Las soluciones en uso son (a) pseudogenomas por accesión sustituyendo SNPs del
1001G (lo hizo 1001 Epigenomes para metilomas), (b) pangenomas de referencia con grafos. Para
lncRNAs aún no hay un estándar; la práctica de 2024 sigue siendo Col-0 + aceptar el sesgo, y
declararlo.

## Paso 3. Ensamblaje de transcritos y definición de loci

**StringTie2** guiado por anotación es el estándar (Kornienko: `-c 2 -m 150 -j 2.5 -a 15`,
v2.1.5); Scallop2 es la alternativa. Después, **merge jerárquico**: primero por tejido/tipo de
dato, después entre tejidos (Kornienko usó Cuffmerge; hoy `stringtie --merge` o TACO hacen lo
mismo y están mantenidos). Luego **gffcompare** contra la anotación de referencia para asignar
class codes: `u` (intergénico), `x` (antisense exónico), `i` (intrónico), `o`/`e`/`j` (solapan
en sentido con genes conocidos: se descartan o se tratan como isoformas).

Filtros que hoy se consideran mínimos:

- Longitud ≥ 200 nt (multi-exón) y, cada vez más, **≥ 400 nt para mono-exónicos** (Kornienko),
  porque los mono-exónicos cortos son mayoritariamente ruido o fragmentos de TEs.
- Separación explícita de **genes TE, fragmentos de TE (> 60 % solapamiento exónico en la
  misma hebra) y pseudogenes** antes de llamar lncRNA. Es el paso que la mayoría de pipelines
  genéricos omiten y en Arabidopsis es decisivo: una gran parte de los "lncRNAs" nuevos son
  transcritos de TEs silenciados de forma variable.
- Soporte de expresión: TPM ≥ 0.5–1 en ≥ 2 muestras o en ≥ 1 accesión con réplicas.

## Paso 4. Potencial codificante

No hay una herramienta ganadora; el consenso de 2024-2025 es **combinar 2-3 predictores y
usar modelos reentrenados en plantas**. Plant-LncPipe (2024) muestra que LncFinder-plant y
CPAT-plant reentrenados superan a los modelos originales; LncRAnalyzer (2025) encuentra que
CPC2 y FEELnc son los más precisos en sorgo; PLncFire (2025, preprint) usa consenso CPC2 +
FEELnc + Plant-LncPipe y reporta que FEELnc recupera 85 % de lncRNAs anotados frente a 44 %
de CPC2 usado solo. RNAsamba advierte en su propio artículo de transferencia pobre entre
especies. A esto se añade **Pfam/hmmscan en 6 marcos** para quitar ORFs con dominios y, en
trabajos cuidadosos, exclusión de sORFs con evidencia Ribo-seq.

Decisión para el proyecto: CPC2 como filtro base (rápido, sin entrenamiento), FEELnc en modo
shuffle como segundo voto, Pfam opcional. Es lo que implementa el workflow.

## Paso 5. Cuantificación

**salmon** con decoys genómicos (o kallisto) sobre un transcriptoma combinado codificantes +
lncRNAs es hoy preferible a featureCounts: corrige sesgos GC/secuencia, maneja multi-mapeo
(crítico en lncRNAs derivados de TEs) y devuelve TPM y counts comparables. Kornienko usó
featureCounts exónico; la diferencia práctica es pequeña para loci bien separados, pero
salmon reparte mejor las lecturas ambiguas. La agregación a gen se hace con tximport
(`lengthScaledTPM`).

## Paso 6. Variación de expresión entre genotipos

Lo que se sabe: los lncRNAs tienen **mucha mayor varianza inter-accesión** que los genes
codificantes, la mayoría están silenciados en la mayoría de accesiones, y esa variación está
asociada a marcas de silenciamiento tipo TE (mCHH, H3K9me2, siRNAs de 24 nt). Es exactamente la
estructura que hace falta para estratificar.

Métodos que se consideran correctos:

- **Normalización**: VST o rlog de DESeq2; corrección de lote con covariables explícitas
  (limma `removeBatchEffect`) solo para análisis no supervisados; nunca en el modelo de DE.
- **Clustering no supervisado**: ConsensusClusterPlus (remuestreo 80 %, 1000 iteraciones) con
  PAC para elegir k, NMF como alternativa; comparar con etiquetas externas por ARI. Es el
  esquema que la oncología ha estandarizado (subtipos por lncRNAs en TCGA) y se traslada
  sin cambios.
- **Firmas**: elastic-net multinomial con validación cruzada anidada y validación en cohorte
  independiente. En oncología el paso siguiente es Cox; aquí el equivalente es asociar
  clusters con clima de origen, floración, grupos genéticos.

## Paso 7. Integración con el genoma (la parte DNA-seq)

1. **cis-eQTL** sobre expresión de lncRNAs con SNPs del 1001G: MatrixEQTL/tensorQTL, ventana
   ±100 kb, factores ocultos (PEER o PCs), kinship si se hace GWAS genómico (GEMMA/LIMIX,
   que es lo que usaron en 1001 Epigenomes). El antecedente más claro en plantas es maíz, donde
   dominan los trans-eQTL para lncRNAs; en Arabidopsis esto no se ha hecho sistemáticamente
   para la anotación de 2024. **Hueco publicable.**
2. **Presencia/ausencia de loci de lncRNA** entre accesiones: requiere ensamblajes. Hoy hay 32
   accesiones con HiFi (Nat Commun 2023) y matrices PAV de TEs; con cobertura de lecturas
   cortas del 1001G se llama PAV con ~66 % de verdaderos positivos (estudio de TEs en 1.114
   accesiones). Se puede hacer un primer análisis con profundidad normalizada de lecturas de ADN
   sobre cada locus de lncRNA.
3. **Expresión alelo-específica** en heterocigotos: casi irrelevante en Arabidopsis
   (autógama, líneas puras).
4. **Metilación y small RNAs** del mismo tejido como variables explicativas (ya hecho por
   Kornienko a nivel descriptivo; falta modelar conjuntamente SNP + mC + TE).

## Paso 8. Validación funcional computacional

Conservación (sintenia y secuencia en Brassicaceae), co-expresión con genes codificantes
(WGCNA, "guilt by association"), solapamiento con picos de GWAS fenotípico del AraGWAS Catalog,
estructura (RNAfold) y mecanismos (cis-NAT con el gen sentido). Ninguno es prueba; todos
aportan priorización.

## Paso 9. Fronteras y lo que todavía no es estándar

- **Long reads** (Iso-Seq / ONT): resuelven isoformas completas y suben la proporción de
  novedad (~30 % de isoformas nuevas en nueve especies vegetales); no hay atlas long-read de
  lncRNAs en poblaciones de Arabidopsis. Añadir Iso-Seq de 5-10 accesiones sería un
  diferenciador claro.
- **Single-cell / single-nucleus**: prácticamente sin trabajos de lncRNAs en plantas.
- **Pangenomas de transcritos**: anotar lncRNAs sobre cada ensamblaje de accesión y proyectar
  a coordenadas comunes. Falta herramienta madura.
- **Chlamydomonas y Marchantia**: tienen genoma excelente (CC-4532 v6.1; Tak-1 v6.1) pero no
  hay paneles de RNA-seq multi-genotipo; estratificación poblacional inviable sin generar datos.

## Resumen operativo

| Componente | Estándar actual | Lo que hace el workflow |
|---|---|---|
| Lecturas | fastp, STAR 2-pass, SortMeRNA si total RNA | sí |
| Strand | inferir por muestra | sí (salmon sobre transcriptoma de referencia) |
| Ensamblaje | StringTie2 guiado, merge jerárquico, gffcompare | StringTie2 + merge + gffcompare (merge en un paso) |
| Filtros | ≥ 200/400 nt, TE/pseudogén explícitos, soporte expresión | ≥ 200 nt, TE por solapamiento exónico, soporte TPM; **falta** 400 nt mono-exón y pseudogenes |
| Codificante | consenso CPC2 + FEELnc (+Pfam) | sí, FEELnc y Pfam opcionales |
| Cuantificación | salmon con decoys + tximport | sí |
| Estratificación | consensus clustering + PAC + firma penalizada | sí |
| Genoma | cis-eQTL con 1001G, kinship | MatrixEQTL con PCs; **falta** kinship/GWAS |
| Epigenoma | mC y sRNA como covariables | no implementado |
