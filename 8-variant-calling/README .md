# Variant Calling w R 
## Variant calling w RStudio to proces identyfikacji wariantów genetycznych, takich jak SNP, InDel, CNV czy SV, w danych sekwencjonowania DNA. Służy do badania mutacji genetycznych, analiz populacyjnych oraz poszukiwania biomarkerów czy celów terapeutycznych. W RStudio analizę przeprowadza się za pomocą pakietów Bioconductor, takich jak VariantAnnotation czy Rsamtools, pozwalających na przetwarzanie plików w formatach FASTQ, BAM, CRAM i VCF. Proces obejmuje mapowanie odczytów do genomu referencyjnego, identyfikację wariantów, ich filtrację, annotację oraz wizualizację wyników. Poniżej znajduje się kod, umożliwiający analizę wariantów w RStudio.
---

## Zadanie 1: Instalacja i załadowanie niezbędnych pakietów

Zainstaluj (jeżeli jest to konieczne) i załaduj następujące pakiety z `Bioconductor`: `VariantTools`, `Rsamtools`, `GenomicRanges`, `GenomicFeatures`, `VariantAnnotation`, `BiocParallel`.

**Instrukcje:**

```{R}
# Instalacja menedżera pakietów Bioconductor (jeśli nie jest zainstalowany)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

```{r Instalacja pakietów}
BiocManager::install(c("VariantTools", "Rsamtools", "GenomicRanges", "GenomicFeatures", "VariantAnnotation", "BiocParallel"))
```
