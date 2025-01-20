# Cała analiza opisana poniżej służy do przetwarzania, filtrowania, anotacji i eksploracji wariantów genetycznych w celu lepszego zrozumienia ich znaczenia biologicznego. Umożliwia identyfikację kluczowych cech wariantów, takich jak jakość danych (np. kolumna QUAL), lokalizacja w genomie (regiony kodujące, UTR, międzygenowe) oraz potencjalny wpływ na funkcje genów.
# Krok 1: instalacja i załadowanie niezbędnych pakietów: 
```{r}
install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "GenomicRanges", "AnnotationHub"))
library(VariantAnnotation)
library(GenomicRanges)
library(AnnotationHub)
```
