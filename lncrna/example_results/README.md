# Example results: SYNTHETIC validation run of the analysis layer

**Todo lo que hay en este directorio proviene de datos simulados** con
`workflow/scripts/simulate_test_data.py`. No es un resultado biológico. Sirve para:

1. Verificar que los scripts R del workflow (`normalize.R`, `stratify.R`, `signature.R`,
   `de_gxe.R`, `eqtl.R`) se ejecutan de principio a fin y producen los ficheros esperados.
2. Comprobar que recuperan una estructura conocida (validación positiva del método).
3. Mostrar qué figuras y tablas genera el proyecto antes de tener datos reales.

## Simulación

48 accesiones en 3 grupos ocultos, 2 condiciones (control / drought) × 2 réplicas = 192
muestras, 2 "estudios" como batch, 1.500 lncRNAs y 4.000 genes codificantes, 20.000 SNPs
en 5 cromosomas con frecuencias alélicas correlacionadas con el grupo. Verdad simulada:
300 lncRNAs dirigidos por grupo, 80 lncRNAs con interacción genotipo × condición, 150 lncRNAs
con cis-eQTL (SNP a < 50 kb).

## Resultados de la validación

| Paso | Resultado | Verdad |
|---|---|---|
| Consensus clustering lncRNA | k = 3 por PAC (PAC = 0), silhouette 0.68 | 3 grupos |
| Clustering lncRNA vs coding | ARI = 1.00 | mismos grupos |
| Clustering vs etiquetas externas | ARI = 1.00 | `external_labels.csv` |
| Firma elastic-net (CV anidada 5×5) | 49 lncRNAs, accuracy 1.00 en todos los folds | las 49 están entre los 300 lncRNAs dirigidos por grupo |
| cis-eQTL (MatrixEQTL, FDR 0.05) | 63 eGenes | 57 verdaderos (lead SNP correcto en 55), 6 falsos positivos (~10 % FDR empírico) |
| G × E (DESeq2 LRT), simulación reducida de 12 accesiones (`*_small`) | 203 features con interacción (padj < 0.05), 77 lncRNAs | 80 lncRNAs con G×E simulada, 67 recuperados |
| Efecto de condición (DESeq2), simulación reducida | 290 features DE | 10 % de features con efecto simulado |

La simulación completa de 48 accesiones hace el LRT de interacción (48 niveles × condición)
demasiado lento para este entorno de 4 núcleos; en un servidor normal es cuestión de minutos.

Las figuras están en `figures/` (PCA, PAC por k, matrices de consenso, heatmap de la firma,
distribución de p y tamaños de efecto de eQTL) y las tablas en `tables/`, incluidos los
ficheros de verdad simulada (`truth_*.tsv`) para reproducir la comparación.

## Reproducir

```bash
python3 workflow/scripts/simulate_test_data.py --outdir /tmp/simrun
# generar los spec JSON como en el README principal y ejecutar cada script con
Rscript workflow/scripts/run_script_standalone.R workflow/scripts/normalize.R spec_normalize.json
```
