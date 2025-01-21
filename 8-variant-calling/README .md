# Variant Calling w R 
## Variant calling w RStudio to proces identyfikacji wariantów genetycznych, takich jak SNP, InDel, CNV czy SV, w danych sekwencjonowania DNA. Służy do badania mutacji genetycznych, analiz populacyjnych oraz poszukiwania biomarkerów czy celów terapeutycznych. W RStudio analizę przeprowadza się za pomocą pakietów Bioconductor, takich jak VariantAnnotation czy Rsamtools, pozwalających na przetwarzanie plików w formatach FASTQ, BAM, CRAM i VCF. Proces obejmuje mapowanie odczytów do genomu referencyjnego, identyfikację wariantów, ich filtrację, annotację oraz wizualizację wyników. Poniżej znajduje się kod, umożliwiający analizę wariantów w RStudio.
---

## Instalacja i załadowanie niezbędnych pakietów
```{R}
# Instalacja menedżera pakietów Bioconductor (jeśli nie jest zainstalowany)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
```{r Instalacja pakietów}
BiocManager::install(c("VariantTools", "Rsamtools", "GenomicRanges", "GenomicFeatures", "VariantAnnotation", "BiocParallel"))
```
```{r Załadowanie pakietów}
library(VariantTools)
library(Rsamtools)
library(GenomicRanges)
library(GenomicFeatures)
library(VariantAnnotation)
library(BiocParallel)
```
## Zapoznanie się z pakietami do wykrywania wariantów
### Funkcja ?? w R służy do wyszukiwania dokumentacji i pomocy dotyczącej określonych terminów lub pakietów. Wyszukuje zarówno w nazwach funkcji, jak i w opisie dokumentacji dostępnej w załadowanych pakietach oraz zainstalowanych w systemie. Na przykład:
1. Wyświetlam pomoc dla pakietu `VariantTools`:
```{R}
??VariantTools
```
2. Wprowadzenie do pakietu:

```{R}
vignette("VariantTools")
```
## Konfiguracja środowiska pracy
### Ustawiam katalog roboczy i sprawdzam dostępność danych.
```{R}
setwd("C:/Users/marys/Dropbox/Mój komputer (LAPTOP-6CIOMO42)/Desktop/alignment")
```
### Sprawdzam, czy pliki są dostępne:

```{R}
list.files()
```
### Należy upewnić się, że w katalogu znajdują się pliki niezbędne do analizy. W tym przypadku będą to:
- Plik BAM z odczytami (`aligned_sample.BAM`)
- Genom referencyjny w formacie FASTA (`ecoli_reference.fasta`)
  
##  Wczytanie danych

### Wczytuję plik BAM i genom referencyjny. Indeksuję plik FASTA.
```{R}
bamfile <- "C:/Users/marys/Dropbox/Mój komputer (LAPTOP-6CIOMO42)/Desktop/alignment/aligned_sample.BAM"
bam <- BamFile(bamfile)
```
```{R}
ref_genome <- "C:/Users/marys/Dropbox/Mój komputer (LAPTOP-6CIOMO42)/Desktop/alignment/ecoli_reference.fasta"
fa <- FaFile(ref_genome)
```
### Sortuję plik BAM według współrzędnych
```{r}
input_bam <- "C:/Users/marys/Dropbox/Mój komputer (LAPTOP-6CIOMO42)/Desktop/alignment/aligned_sample.BAM"
output_bam <- "C:/Users/marys/Dropbox/Mój komputer (LAPTOP-6CIOMO42)/Desktop/alignment/sorted_aligned_sample.BAM"

sortBam(file = input_bam, destination = output_bam, overwrite = TRUE)

# Definiuję przesortowany plik: 
sorted_bam <- "C:/Users/marys/Dropbox/Mój komputer (LAPTOP-6CIOMO42)/Desktop/alignment/sorted_aligned_sample.BAM.bam"    
```
### Indeksuję plik FASTA, jeśli indeks nie istnieje, oraz przesortowany plik BAM:

```{R}
indexFa(ref_genome)
indexBam(sorted_bam)
```
## Podstawowa kontrola jakości danych sekwencyjnych przed wykrywaniem wariantów
### Sprawdzam nagłówek pliku BAM:
```{R}
scanBamHeader(bam)
```
### Sprawdzam podstawowe statystyki pliku BAM:
```{R}
idxstats <- idxstatsBam(sorted_bam)
print(idxstats)
```
### Obliczam i wizualizuję pokrycie genomu:
```{R}
coverage_data <- coverage(sorted_bam)
summary(coverage_data[[1]]) # dla genomów prokariota
plot(coverage_data[[1]], main="Pokrycie genomu dla sekwencji U00096.3", ylab="Pokrycie", xlab="Pozycja w genomie") 
```
## Wykrywanie wariantów za pomocą funkcji `callVariants()`.

```{R}
pileup_param <- PileupParam(
    distinguish_strands = FALSE,
    distinguish_nucleotides = TRUE,
    min_base_quality = 20
)

pile <- pileup(sorted_bam, scanBamParam = ScanBamParam(), pileupParam = pileup_param)

```

### Konwertuję dane pileup do ramki danych z uzgodnieniem nazw sekwencji
```{r}
library(dplyr)

pile_df <- as.data.frame(pile)
pile_df <- pile_df %>%
    mutate(seqnames = as.character(seqnames)) %>%
    mutate(seqnames = ifelse(seqnames == "U00096.3", "NC_000913.3", seqnames))

```

### Grupuję dane według pozycji
```{r}
variant_candidates <- pile_df %>%
    group_by(seqnames, pos) %>%
    summarise(
        total = sum(count),
        A = sum(count[nucleotide == "A"]),
        C = sum(count[nucleotide == "C"]),
        G = sum(count[nucleotide == "G"]),
        T = sum(count[nucleotide == "T"]),
        .groups = 'drop'
    ) %>%
    mutate(
        ref = as.character(getSeq(fa, GRanges(seqnames, IRanges(pos, pos))))
    ) %>%
    rowwise() %>%
    mutate(
        # Obliczanie alternatywnych alleli
        alt_alleles = list(setdiff(c("A", "C", "G", "T"), ref)),
        # Liczenie odczytów dla referencyjnego i alternatywnych alleli
        ref_count = sum(c_across(c("A", "C", "G", "T"))[ref]),
        alt_count = sum(c_across(c("A", "C", "G", "T"))[alt_alleles])
    ) %>%
    ungroup() %>%
    # Filtracja na podstawie minimalnej liczby odczytów dla wariantu
    filter(alt_count >= 5) %>%
    # Opcjonalne filtrowanie na podstawie proporcji
    filter((alt_count / total) >= 0.2)

```

```{r}
# Przykład wyświetlenia wariantów:
head(variant_candidates)
```
## Filtracja i eksportowanie wyników do pliku
```{r}
# Filtracja wariantów na podstawie jakości i głębokości pokrycia
filtered_variants <- variant_candidates %>%
    filter(total >= 10, alt_count / total >= 0.2, alt_count >= 5)

# Wyświetlenie liczby wariantów przed i po filtrowaniu
cat("Liczba wariantów przed filtrowaniem:", nrow(variant_candidates), "\n")
cat("Liczba wariantów po filtrowaniu:", nrow(filtered_variants), "\n")

# Konwersja do data.frame dla eksportu
df_variants <- as.data.frame(filtered_variants)

# Eksport do pliku CSV
write.csv(df_variants, "ścieżka/do/pliku/wyniki_wariantow.csv", row.names = FALSE)
```
---

## Podsumowanie i dyskusja
### Dokładność wykrywania wariantów genetycznych zależy od wielu czynników, takich jak jakość danych sekwencyjnych, głębokość pokrycia oraz parametry użyte podczas procesu variant calling. Istotne jest również unikanie błędów w mapowaniu odczytów do genomu referencyjnego, które mogą znacząco wpłynąć na wiarygodność wyników. Filtrowanie wariantów odgrywa kluczową rolę w eliminowaniu fałszywie pozytywnych wyników, co pozwala skupić się na biologicznie istotnych wariantach i poprawia jakość analiz. Źródła błędów w procesie variant calling mogą obejmować artefakty PCR, błędy sekwencjonowania, zanieczyszczenia próbek czy niewłaściwe mapowanie odczytów. Dlatego kluczowe jest zastosowanie odpowiednich metod filtracji oraz dokładne sprawdzanie wyników. Po wykryciu wariantów kolejne kroki obejmują ich annotację, analizę funkcjonalną oraz – w kontekście medycznym – interpretację kliniczną. 
