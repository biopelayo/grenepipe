# Lógica comparativa de los análisis

Qué se compara con qué, por qué, con qué test, controlando qué, y qué resultado decide el
paso siguiente. Construido sobre lo que Kornienko et al. 2024 ya demostraron
([DOI](https://doi.org/10.1093/plcell/koad233)) y sobre lo que dejaron abierto, con los
genotipos y metadatos del 1001 Genomes / 1001 Epigenomes (Kawakatsu et al. 2016,
[DOI](https://doi.org/10.1016/j.cell.2016.06.044)). La versión legible por máquina está en
`config/contrasts.yaml`.

## 0. Hechos de partida que condicionan el diseño

| Hecho (Kornienko 2024) | Consecuencia para nosotros |
|---|---|
| "Expresado" = TPM > 0,5; la mayoría de lncRNAs están OFF en > 95 % de accesiones; ~50 % de PC genes están ON en todas | La variable primaria de un lncRNA entre accesiones es **binaria (ON/OFF)**, no continua. Hay que analizar la matriz de presencia además de la de TPM |
| Dos accesiones comparten ~50 % de lncRNAs expresados vs ~100 % de PC | La distancia entre accesiones por lncRNAs debe ser **Jaccard sobre ON/OFF**; Pearson sobre VST infraestima la estructura |
| lincRNAs y TEs agrupan muestras por **accesión**; PC genes por **tejido** | Clustering de accesiones: usar lincRNAs; clustering de tejidos: PC. Multi-tejido (GSE226691) separa las dos señales |
| Variación intra-accesión (lotes) < inter-accesión para todo salvo AS lncRNAs | AS lncRNAs son ruidosos: la estratificación se hace con **lincRNAs**; AS entra solo en análisis cis con su gen pareja |
| lincRNAs = PC-like (K27, 1 copia, sin TE) + TE-like (K9, piezas de TE, multicopia, mC, sRNA 24 nt) | **Todo contraste se estratifica por TE-like vs PC-like**: son dos poblaciones de loci con biología distinta |
| Nº de lincRNAs TE-like ON por accesión ∝ nº de TE genes ON (r alto); varía 3× entre accesiones; GWAS sin señal clara | Definir el fenotipo **"eficiencia de silenciamiento"** por accesión y tratarlo como covariable y como rasgo (GWAS con kinship, relict vs non-relict, clima) |
| Metilación explica ~51 % de lincRNAs informativos; **no se separó genético de epigenético** | Modelo conjunto SNP-cis + mC + TE-content → partición de varianza. Hueco principal |
| 157 lincRNAs contienen top hits GWAS de 65 fenotipos, sobre todo climáticos | Enriquecimiento de clusters/firma en esos loci; asociación de expresión con clima de origen |
| 14 AS lncRNAs anticorrelacionados con su PC pareja | Contraste cis AS–PC entre accesiones para generalizarlo |
| Tejidos: flores y polen expresan muchos más lincRNAs y TEs; 23 accesiones tienen los 4 tejidos | Diseño emparejado accesión × tejido en las 23 completas |

## 1. Unidades, estratos y cohortes

**Unidad de observación**: accesión (media de réplicas/lotes). **Unidad de inferencia**: locus.

**Estratos de loci** (columnas de `lncrna_classification.tsv`, ampliadas):

- tipo: lincRNA / antisense / intrónico / codificante (control) / TE gene (control)
- subtipo lincRNA: TE-like (te_overlap_frac ≥ 0,5 o pieza de TE por BLAST) vs PC-like
- copia: single vs multicopy (BLASTN del locus contra TAIR10, > 80 % id y longitud)
- origen: novel / Kornienko 2024 / Araport11
- frecuencia ON en la cohorte: privado (1 acc.), raro (< 5 %), intermedio (5-50 %), común (> 50 %)

**Estratos de accesiones** (tabla `config/accession_metadata.csv`, a construir desde la tabla de
accesiones de 1001genomes.org, AraClimate y AraPheno):

- grupo de admixture del 1001G (9 grupos + relict); relict vs non-relict
- clima de origen (bioclim de WorldClim: temperatura media anual, precipitación, estacionalidad)
- latitud / longitud
- fenotipos clave de AraPheno: tiempo de floración 10 °C/16 °C, metilación global (mCHH) de 1001 Epigenomes
- eficiencia de silenciamiento: nº de TE genes ON y nº de lincRNAs TE-like ON en roseta (calculado)

**Cohortes y papel**:

| Cohorte | Papel | Restricciones |
|---|---|---|
| GSE224760 (28 acc. × 2-4 lotes, roseta, stranded) | descubrimiento de anotación; estimación de varianza intra vs inter accesión; **no** para GWAS (n = 28) | batch = lote |
| GSE226691 (25-27 acc. × 4 tejidos, stranded) | efecto tejido; diseño emparejado 23 acc. × 4 tejidos | sin réplicas |
| GSE80744 (461-727 acc., roseta, no stranded, SE 100) | validación de clusters y firma; eQTL y GWAS con potencia; solo lincRNA + intrónico | batch desconocido: usar factores ocultos |
| Metilomas GSE43857 (444 solapan) | covariable epigenética por locus | mismo tejido |
| VCF 1001G v3.1 | genotipos para todo | ids numéricos |

## 2. Bloques de análisis, en orden de dependencia

### B1. Reproducción y extensión de la anotación (GSE224760 → GSE226691)

- Comparación: nuestra anotación vs GTF publicado de Kornienko (`gffcompare`): sensibilidad y
  precisión a nivel de locus; novedad por tipo.
- Decisión: si recuperamos < 70 % de sus loci en roseta con los mismos 28 genotipos, revisar
  parámetros antes de seguir. Si > 85 %, fusionar y continuar.

### B2. Estructura de la variación (descriptivo pero cuantificado)

| Contraste | Métrica | Test | Control |
|---|---|---|---|
| inter- vs intra-accesión por tipo de locus | CV entre medias de accesión vs CV entre lotes | Mann-Whitney por tipo; ICC por locus con modelo mixto `~1 + (1|accesión) + (1|lote)` | solo loci ON en ≥ 2 accesiones |
| TE-like vs PC-like lincRNA | frecuencia ON, CV, ICC | Mann-Whitney | distancia al centrómero (≤ 2 Mb vs > 2 Mb) como estrato |
| accesión vs tejido como factor dominante | % varianza explicada por accesión y por tejido en GSE226691 (23 acc. × 4 tejidos) | ANOVA por locus, agregado por tipo | PC genes como referencia |

Resultado esperado si el paper se sostiene: ICC(lincRNA) > ICC(PC) > ICC(AS); accesión > tejido
para lincRNA TE-like y tejido > accesión para PC.

### B3. Estratificación de accesiones (el núcleo)

Tres matrices de distancia, mismas accesiones, mismo procedimiento (ConsensusClusterPlus +
PAC):

1. **PC-VST / Pearson** (línea base: debería reproducir la estructura genética)
2. **lincRNA-VST / Pearson**
3. **lincRNA ON-OFF / Jaccard** (la distancia que corresponde a la biología descrita)

y por separado para TE-like y PC-like.

| Comparación | Métrica | Pregunta |
|---|---|---|
| clusters(2) y clusters(3) vs clusters(1) | ARI | ¿los lncRNAs aportan estructura distinta de la codificante? |
| clusters vs admixture group | ARI, test de permutación (1000) | ¿recapitulan ancestría? |
| clusters vs clima de origen | ANOVA / Kruskal de bioclim por cluster; regresión multinomial cluster ~ bioclim + PCs genéticos | ¿hay señal climática independiente de ancestría? |
| clusters vs eficiencia de silenciamiento | Kruskal | ¿un cluster es "silenciamiento fallido" global? Si sí, es un confusor: repetir el clustering residualizando por esa covariable |
| clusters vs mCHH global | Kruskal | ídem para epigenoma |

Regla: un cluster que desaparece al residualizar por eficiencia de silenciamiento o por los
10 primeros PCs genéticos no es un hallazgo de lncRNAs; se reporta como tal.

Firma: elastic-net sobre la matriz que dé el mejor PAC (binaria o continua), CV anidada en
GSE224760, **validación ciega en GSE80744** (misma anotación, solo lincRNA/intrónico,
cuantificación salmon no stranded). Métricas: accuracy, kappa, y ARI entre clusters
predichos y clusters recalculados de novo en la cohorte de validación.

### B4. Genética vs epigenética de la expresión de cada lncRNA (el hueco principal)

Para cada lincRNA informativo (ON en ≥ 5 % y OFF en ≥ 5 % de las 444 accesiones con RNA + mC):

| Modelo | Variables | Salida |
|---|---|---|
| M1 cis-genético | SNPs ± 100 kb (lead cis-eQTL, MAF ≥ 5 %) | β, p, FDR; eGene sí/no |
| M2 epigenético | mCG y mCHH del cuerpo y promotor del locus (444 acc.) | R² |
| M3 conjunto | lead SNP + mCG + mCHH + PCs genéticos + eficiencia de silenciamiento | partición de varianza (R² parcial de cada término) |
| M4 mediación | SNP → mC → expresión | ¿la metilación media el efecto del SNP? (ACME, bootstrap) |

Contrastes agregados: proporción de eGenes y R² por estrato (TE-like vs PC-like, single vs
multicopy, raro vs común); comparación con PC genes como referencia (el 1001 Epigenomes ya
tiene eQTL de PC con LIMIX). Hipótesis: TE-like → dominan mC y trans (silenciamiento);
PC-like → dominan cis-SNP.

Para ON/OFF, la variable es binaria: regresión logística con kinship (GEMMA/SAIGE-style) o
test de Fisher por SNP lead; reportar también "presencia/ausencia del locus" estimada con
profundidad de lecturas de ADN del 1001G normalizada (si el locus falta en una accesión, no es
silenciamiento, es deleción: estrato aparte).

### B5. Eficiencia de silenciamiento como rasgo de accesión

Fenotipo: nº de lincRNAs TE-like ON y nº de TE genes ON (roseta, 461 acc.), residualizados por
profundidad. Análisis: (a) heredabilidad SNP (GCTA/LIMIX); (b) GWAS con kinship (ellos no
vieron señal clara: confirmar o refinar con el fenotipo mejor definido y covariables); (c)
asociación con clima y relict; (d) correlación con mCHH global de 1001 Epigenomes.

### B6. Antisense cis (AS lncRNA ↔ PC pareja)

Para cada par AS–PC con ambos ON en ≥ 10 % de accesiones: correlación de Spearman entre
accesiones (GSE224760 y GSE226691, stranded); distribución de ρ vs pares aleatorios; lista de
pares con ρ < −0,5 (generaliza los 14 del paper); contraste de si el SNP lead del AS es también
eQTL del PC (co-localización simple por LD, r² > 0,8).

### B7. Tejido × accesión (GSE226691, 23 acc. completas)

Modelo por locus `log TPM ~ accesión + tejido + accesión:tejido`; proporción de loci con
interacción por tipo; conjuntos ON específicos de polen/flor y su solapamiento entre
accesiones (Jaccard); ¿la estratificación por lincRNAs es la misma en los 4 tejidos? (ARI
entre clusterings por tejido).

### B8. Priorización funcional

Loci de la firma y eGenes cruzados con: los 157 lincRNAs con hits GWAS climáticos del paper,
AraGWAS Catalog (hits dentro del locus ± 5 kb), conservación en Brassicaceae (sintenia),
co-expresión con PC (WGCNA en 461 acc.), y categoría K9/K27 cuando haya ChIP (14 acc.).

## 3. Confusores y cómo se controlan

| Confusor | Dónde muerde | Control |
|---|---|---|
| Lote de cultivo (GSE224760) | varianza intra-accesión | columna `batch`; modelo mixto; limma solo en no supervisado |
| Profundidad y % mapeo | nº de loci ON por accesión | residualizar por log(lecturas) y % uniquely mapped |
| Librería no stranded (GSE80744) | antisense | excluir AS; solo lincRNA/intrónico |
| Sesgo de referencia (Col-0) | loci ausentes o divergentes en accesiones lejanas; relicts | flag de PAV con ADN; sensibilidad excluyendo relicts |
| Estructura poblacional | eQTL, GWAS, clusters vs clima | kinship / PCs genéticos; permutaciones estratificadas por admixture group |
| Distancia al centrómero | TE-like vs PC-like | estrato ≤ 2 Mb / > 2 Mb |
| Multimapeo (multicopia) | cuantificación de lincRNAs TE-like | salmon reparte; reportar con y sin loci multicopia |
| Umbral ON/OFF (0,5 TPM) | todo lo binario | sensibilidad con 0,25 y 1,0 |

## 4. Criterios de decisión al final de cada bloque

- B1: recuperación ≥ 85 % de la anotación publicada → seguir; si no, ajustar.
- B2: si ICC(lincRNA) no supera ICC(PC), la premisa de estratificación falla; parar y reportar.
- B3: PAC < 0,1 y k estable en remuestreos de accesiones (jackknife de 10 %) → cluster
  robusto; ARI con admixture > 0,6 → "recapitula ancestría" (resultado menor); ARI con
  admixture < 0,3 pero asociación climática significativa tras controlar PCs genéticos →
  resultado principal.
- B3 validación: firma con accuracy > 0,8 y ARI predicho-vs-de-novo > 0,5 en GSE80744 → firma
  válida.
- B4: si > 30 % de lincRNAs informativos tienen cis-eQTL a FDR 5 % en 444 acc. → genética
  sustancial; partición de varianza por estrato es el resultado central del bloque.

## 5. Qué falta implementar para cubrir esta lógica

| Elemento | Estado |
|---|---|
| Matriz ON/OFF + clustering Jaccard + comparación con continuo | **añadido** a `stratify.R` (`binary` set) |
| Tabla de metadatos de accesión y su uso en `stratify.R` | plantilla `config/accession_metadata.template.csv`; `external_labels` acepta cualquier columna |
| ICC / modelo mixto intra vs inter accesión | pendiente (`variance_partition.R`) |
| TE pieces por BLAST y copia por BLASTN | pendiente (regla `te_pieces_copy_number`); hoy solo solapamiento exónico con TEs anotados |
| Metilación por locus (444 acc.) | pendiente: descarga de GSE43857 procesado + `bedtools map` |
| Modelo conjunto y mediación | pendiente (`genetic_vs_epigenetic.R`) |
| GWAS con kinship (GEMMA) | pendiente |
| PAV por profundidad de ADN | pendiente; opcionalmente grenepipe para BAMs |
| Validación de firma en cohorte externa | pendiente (`apply_signature.R`) |
