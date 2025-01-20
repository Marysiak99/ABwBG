# Cała analiza opisana poniżej służy do przetwarzania, filtrowania, anotacji i eksploracji wariantów genetycznych w celu lepszego zrozumienia ich znaczenia biologicznego. Umożliwia identyfikację kluczowych cech wariantów, takich jak jakość danych (np. kolumna QUAL), lokalizacja w genomie (regiony kodujące, UTR, międzygenowe) oraz potencjalny wpływ na funkcje genów.
## Krok 1: instalacja i załadowanie niezbędnych pakietów: 
```{r}
install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "GenomicRanges", "AnnotationHub"))
library(VariantAnnotation)
library(GenomicRanges)
library(AnnotationHub)
```
## Krok 2: Wczytanie i eksploracja danych:
### Znajduję ścieżkę do przykładowego pliku VCF z pakietu VariantAnnotation:
```{r}
vcf_file <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
```
"C:/Users/marys/AppData/Local/R/win-library/4.4/VariantAnnotation/extdata/chr22.vcf.gz"
### Wczytuję plik VCF do obiektu typu VCF
```{r}
vcf_data <- readVcf(vcf_file, "hg19")
```
### Wyświetlam podstawowe informacje o obiekcie VCF
```{r}
vcf_
```
### Pobieram nagłówki i metadane
```{r}
vcf_header <- header(vcf_data)

# Wyświetl nagłówki
vcf_header
```
```{r}
### Pola INFO
info(vcf_header)

### Pola FORMAT
geno(vcf_header)
```
```{r}
### Liczba wariantów w pliku VCF
num_variants <- nrow(vcf_data)
num_variants
