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
### Liczba wariantów: 10376 
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
### Liczba wariantów przed filtracją: 10376 
### Liczba wariantów po filtracji (QUAL >= 30): 10364 
# Krok 4: Anotacja wariantów
```{r}
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
```
```{r}
### Sprawdzenie dostępnych kolumn w INFO - ładuję wymagane pakiety
```{r}
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
```
### Przekształcam obiekt VCF na GRanges
```{r}
granges_vcf <- as(vcf_data, "GRanges")
```
### Przygotowanie obiektu TxDb
```{r}
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
```
### Anotacja wariantów
```{r}
annotated_variants <- locateVariants(granges_vcf, txdb)
```
### Sprawdzenie wyników
```{r}
head(annotated_variants)
```
### Pobieranie regionów 5'UTR i 3'UTR
```{r}
utr5 <- fiveUTRsByTranscript(txdb)
utr3 <- threeUTRsByTranscript(txdb)

utr_regions <- c(utr5, utr3)

utr5_variants <- findOverlaps(granges_vcf, utr5)
utr3_variants <- findOverlaps(granges_vcf, utr3)

cat("Liczba wariantów w 5'UTR:", length(utr5_variants), "\n")
cat("Liczba wariantów w 3'UTR:", length(utr3_variants), "\n")
```
# Podsumowanie 
### W analizie wariantów genetycznych na chromosomie 22 przed filtracją znajdowało się 10376 wariantów, z czego po zastosowaniu filtru jakościowego (QUAL ≥ 30) liczba ta zmniejszyła się do 10364. Po wyodrębnieniu wariantów z regionów kodujących, porównano ich liczbę z wariantami znajdującymi się w regionach międzygenowych. W analizie uwzględniono także warianty z regionów 5'UTR i 3'UTR, które następnie zostały porównane pod kątem liczebności. Możliwości rozszerzenia analizy obejmują ocenę funkcji genów, powiązanie wariantów z chorobami, analizę częstotliwości wariantów w różnych populacjach, oraz zastosowanie narzędzi bioinformatycznych, takich jak SIFT czy PolyPhen, do oceny biologicznego wpływu tych wariantów.
