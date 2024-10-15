## Zadanie 2: Instalacja i użycie pakietów

**Cel:** Nauka instalowania i używania pakietów w R

### Instrukcje i rozwiązanie:

1.  **Zainstaluj pakiet `ggplot2`.**
  
  ``` r
install.packages("ggplot2")
```

2.  **Załaduj pakiet.**
  
  ``` r
library(ggplot2)
```

3.  **Sprawdź dokumentację pakietu.**
  
  ``` r
?ggplot2
```

4.  **Zainstaluj i załaduj** dodatkowy pakiet, który Cię zainteresuje. Listę pakietów w repozytorium CRAN wraz z opisami znajdziesz tutaj: <https://cran.r-project.org/web/packages/available_packages_by_name.html>
  
  ``` r
install.packages("plotly")
library(plotly)

# Pakiet 'plotly' służy do tworzenia interaktywnych wykresów.
# Zainteresował mnie, ponieważ umożliwia tworzenie dynamicznych wizualizacji danych.
```
install.packages("geno2proteo")
library(geno2proteo)
#pakiet ułatwia szukanie sekwencji DNA oraz białek 
#nie udało się zainstalować, ponieważ pakiet "GenomicRanges" jest dostępny w innej wersji R

