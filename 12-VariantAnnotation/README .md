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
```
### Wyświetlam nagłówki
```{r}
vcf_header
```
### Pola INFO
```{r}
info(vcf_header)
```
### Pola FORMAT
```{r}
geno(vcf_header)
```
### Liczba wariantów w pliku VCF
```{r}
num_variants <- nrow(vcf_data)
num_variants
```
## Krok 3: Filtracja i analiza jakości
```{r}
# Pobieram kolumnę QUAL i usuwam NA
qual_values <- fixed(vcf_data)$QUAL
qual_values_clean <- qual_values[!is.na(qual_values)]

# Filtracja wariantów o jakości QUAL >= 30, pomijając NA
high_quality_variants <- vcf_data[!is.na(qual_values) & qual_values >= 30, ]

# Porównuję liczby wariantów przed i po filtracji
num_variants_before <- nrow(vcf_data)
num_variants_after <- nrow(high_quality_variants)

cat("Liczba wariantów przed filtracją:", num_variants_before, "\n")
cat("Liczba wariantów po filtracji (QUAL >= 30):", num_variants_after, "\n")
```
Liczba wariantów przed filtracją: 10376 
Liczba wariantów po filtracji (QUAL >= 30): 10364 
