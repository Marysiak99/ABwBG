# Proste wykresy eksploracyjne
## Wykres violin przedstawia rozkład danych numerycznych, łącząc elementy wykresu pudełkowego (boxplot) i estymacji gęstości rozkładu (kde plot). Pokazuje zarówno podsumowanie statystyczne, jak mediana i kwartyle, jak i kształt rozkładu danych w postaci symetrycznej krzywej. Dzięki temu umożliwia łatwe porównanie rozkładów między różnymi grupami oraz identyfikację asymetrii lub wielomodalności. Wykres violin jest użyteczny, gdy chcemy uzyskać pełniejszy obraz danych niż ten, który oferuje klasyczny boxplot. Jednak jego interpretacja może być trudniejsza w przypadku małych próbek lub dużego szumu w danych. 
```{r}
# install.packages("ggplot2") 
library(ggplot2)

# Przykładowy zbiór danych: iris
data(iris)

## 1a. Boxplot
ggplot(iris, aes(x = Species, y = Sepal.Length, fill=Species)) +
  geom_boxplot() +
  labs(title="Boxplot - Długość działki kielicha (Sepal.Length) wg gatunków irysa")

## 1b. Histogram
ggplot(iris, aes(x = Sepal.Width)) +
  geom_histogram(binwidth = 0.2, fill="lightblue", color="black") +
  labs(title="Histogram - Szerokość działki kielicha (Sepal.Width)")

## 1c. Scatter plot
ggplot(iris, aes(x=Sepal.Length, y=Petal.Length, color=Species)) +
  geom_point() +
  labs(title="Scatter plot - Iris")

## 1d. Violin + Boxplot (hybryda)
ggplot(iris, aes(x=Species, y=Sepal.Width, fill=Species)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.1, color="black", outlier.shape=NA) +
  labs(title="Violin + Boxplot - Iris")
```



# Stacked Bar Plot (skumulowane słupki)

```{r}
# Dane przykładowe
df_bar <- data.frame(
  Sample = rep(c("S1","S2","S3"), each=3),
  Category = rep(c("A","B","C"), times=3),
  Count = c(10, 5, 15, 8, 12, 6, 20, 10, 5)
)

ggplot(df_bar, aes(x=Sample, y=Count, fill=Category)) +
  geom_bar(stat="identity") +
  labs(title="Skumulowany wykres słupkowy")
```

# Waffle Plot
## Waffle plot to wykres służący do wizualizacji udziałów procentowych lub proporcji różnych kategorii w zbiorze danych. Składa się z siatki kwadratów, gdzie każdy kwadrat reprezentuje określoną liczbę jednostek lub procent. Dzięki swojej strukturze graficznej jest bardziej atrakcyjny wizualnie niż tradycyjny wykres słupkowy czy kołowy, a jednocześnie czytelnie przedstawia względne wielkości kategorii. Waffle plot jest szczególnie przydatny do prezentowania udziałów w prosty i intuicyjny sposób, ale bywa mniej efektywny przy dużej liczbie kategorii. 

```{r}
# install.packages("waffle")  # w razie potrzeby
library(waffle)

parts <- c(`Cat A (50)`=50, `Cat B (30)`=30, `Cat C (20)`=20)

waffle(parts, rows=5, 
       title="Przykładowy Waffle Plot",
       legend_pos = "bottom")
```


# Time Series Plot (analiza czasowa)
## Przedstawia zmiany wartości danych w czasie, umożliwiając identyfikację trendów, cykli oraz nagłych zmian. W kontekście badań genów wykres ten jest niezwykle użyteczny do analizy ekspresji genów w różnych punktach czasowych, np. podczas procesu rozwoju organizmu lub w odpowiedzi na bodźce. Pozwala śledzić dynamikę działania genów, wykrywać wzorce czasowe oraz porównywać ekspresję między różnymi genami lub warunkami eksperymentalnymi. 

```{r}
# Dane przykładowe (np. ekspresja genu w 5 punktach czasowych, dla 2 genów)
df_time <- data.frame(
  Time = rep(1:5, times=2),
  Expression = c(1,2,3,2.5,4, 2,3.5,3,4,5),
  Gene = rep(c("GeneA","GeneB"), each=5)
)

ggplot(df_time, aes(x=Time, y=Expression, color=Gene)) +
  geom_line() +
  geom_point() +
  labs(title="Analiza czasowa ekspresji genów")
```


# Waterfall Plot
## wykres służący do wizualizacji skumulowanych zmian wartości w kolejnych krokach lub obserwacjach. Składa się z serii słupków, które mogą być dodatnie lub ujemne, a każdy z nich zaczyna się tam, gdzie kończy poprzedni. Dzięki temu pokazuje, jak poszczególne składowe przyczyniają się do całkowitej zmiany wartości. W badaniach biomedycznych, w tym genetycznych, waterfall plot jest często wykorzystywany do przedstawiania reakcji pacjentów na terapię, np. zmian wielkości guzów nowotworowych w odpowiedzi na leczenie. Umożliwia łatwą identyfikację pacjentów z pozytywną, negatywną lub neutralną reakcją na terapię, co pomaga w ocenie skuteczności leczenia.

## 5a. Klasyczny Waterfall (prosta wersja)

```{r}
# Dane sztuczne: zmiana wielkości guza w % (hipotetyczny przykład)
set.seed(123)
df_wf <- data.frame(
  Pacjent = paste0("P", 1:20),
  Zmiana = sample(seq(-100, 100, by=10), 20)
)
# Sortujemy wg wartości
df_wf <- df_wf[order(df_wf$Zmiana), ]

# Prosty "waterfall" z ggplot2 (barplot z uporządkowanymi słupkami)
df_wf$Pacjent <- factor(df_wf$Pacjent, levels=df_wf$Pacjent)  # kolejność

ggplot(df_wf, aes(x=Pacjent, y=Zmiana, fill=Zmiana>0)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=c("red","forestgreen"), name="Zmiana dodatnia?") +
  labs(title="Klasyczny Waterfall Plot - Zmiana wielkości",
       x="Pacjent", y="Zmiana (%)")
```

## 5b. Waterfall w analizach mutacji (pakiet maftools)

```{r}
BiocManager::install("maftools")
library(maftools)
# Przykładowe dane w pakiecie maftools: "tcga_laml.maf.gz"
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf)
oncoplot(maf = laml, top=10)  # typowy "onco waterfall" w stylu mutacji
```


# Volcano Plot
## wykres używany do wizualizacji wyników analiz statystycznych, najczęściej w badaniach różnicowej ekspresji genów. Przedstawia zależność między istotnością statystyczną (p-value) a wielkością efektu (fold change) w formie punktów. Oś X reprezentuje logarytmicznie przekształcony fold change, a oś Y – ujemny logarytm p-value, co pozwala uwypuklić istotnie zmienione geny. Wykres ten jest szczególnie przydatny do szybkiej identyfikacji genów, które wykazują znaczące zmiany w ekspresji oraz istotne statystycznie różnice między badanymi warunkami, np. w porównaniu zdrowych i chorych próbek. Geny, które są zarówno istotne statystycznie, jak i mają duży fold change, pojawiają się jako punkty położone wysoko i daleko od środka wykresu (tzw. „wulkany”).

## Metoda 1: base R (przykład minimalny)

```{r}
# Przykładowe dane
set.seed(123)
df_volcano <- data.frame(
  gene = paste0("Gene", 1:100),
  log2FC = rnorm(100, 0, 1),
  pval = runif(100, 0, 0.05)
)
df_volcano$negLogP <- -log10(df_volcano$pval)

plot(df_volcano$log2FC, df_volcano$negLogP,
     pch=20, col="grey50",
     xlab="log2 Fold Change", ylab="-log10(p-value)",
     main="Volcano Plot (base R)")

abline(h=-log10(0.05), col="red", lty=2)
abline(v=c(-1, 1), col="blue", lty=2)
```

### Metoda 2: pakiet EnhancedVolcano

```{r}
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# Dane przykładowe - minimalnie
EnhancedVolcano(df_volcano,
  lab = df_volcano$gene,
  x = 'log2FC',
  y = 'pval',
  pCutoff = 0.01,
  FCcutoff = 1,
  title = 'Przykładowy volcano plot',
  subtitle = 'DE analysis',
  legendPosition = "right")
```

# Heatmap
## Heatmap (mapa cieplna) to wykres, który przedstawia wartości danych w formie kolorowej siatki, gdzie intensywność koloru reprezentuje wielkość danej wartości. Jest często stosowany w badaniach genetycznych do wizualizacji dużych zbiorów danych, takich jak poziomy ekspresji genów w różnych warunkach lub próbkach.W heatmapach każdy wiersz może reprezentować gen, a każda kolumna warunek eksperymentalny lub próbkę. Dzięki zastosowaniu skali kolorów badacze mogą szybko zidentyfikować wzorce, np. grupy genów, które wykazują podobną ekspresję, lub różnice w ekspresji między różnymi stanami biologicznymi. Heatmapy są szczególnie użyteczne w połączeniu z analizą klastrowania, co pozwala grupować geny lub próbki o podobnych właściwościach.
```{r}
install.packages("pheatmap")
library(pheatmap)

# Tworzymy losową macierz 10 genów x 5 próbek
set.seed(123)
mat <- matrix(rnorm(50), nrow=10, ncol=5)
rownames(mat) <- paste0("Gene", 1:10)
colnames(mat) <- paste0("Sample", 1:5)

pheatmap(mat, 
         scale="row", 
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         main="Heatmap - 10 genów x 5 próbek")
```


# Wykresy redukcji wymiarów (PCA, t-SNE)
## to wizualizacje, które przedstawiają dane wielowymiarowe w dwóch lub trzech wymiarach, ułatwiając ich analizę i interpretację. Stosuje się je, gdy dane zawierają wiele cech (np. ekspresje tysięcy genów), a chcemy zobaczyć ich ogólną strukturę, znaleźć podobieństwa między próbkami lub zidentyfikować grupy o podobnych właściwościach.

## 8a. PCA

```{r}
data(iris)
pca_result <- prcomp(iris[,1:4], center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Species = iris$Species
)

ggplot(pca_df, aes(x=PC1, y=PC2, color=Species)) +
  geom_point() +
  labs(title="PCA - Iris")
```

## 8b. t-SNE

```{r}
# install.packages("Rtsne")
library(Rtsne)

# Usuwamy duplikaty względem kolumn 1:4 (Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
iris_nodup <- iris[!duplicated(iris[,1:4]), ]

# Teraz wywołujemy Rtsne
tsne_out <- Rtsne(iris_nodup[,1:4], pca=FALSE, perplexity=20, max_iter=500)

library(ggplot2)

# Tworzymy data.frame z wynikami t-SNE
tsne_df <- data.frame(
  X = tsne_out$Y[,1],
  Y = tsne_out$Y[,2],
  Species = iris_nodup$Species  # bo usunęliśmy te same wiersze
)

ggplot(tsne_df, aes(x=X, y=Y, color=Species)) +
  geom_point() +
  labs(title="t-SNE - Iris (bez duplikatów)")

```


# Manhattan Plot
## Używany w genetyce, szczególnie w badaniach asocjacyjnych (GWAS), do wizualizacji wyników poszukiwań związków między zmiennymi genotypowymi a cechami fenotypowymi. Jego celem jest pokazanie statystycznej istotności związku między poszczególnymi markerami genetycznymi (np. SNP) a cechami, takimi jak choroby.

```{r}
# install.packages("qqman")
library(qqman)

# Generujemy przykład: 500 SNP w 5 chromosomach
set.seed(123)
chrom <- rep(1:5, each=100)
bp <- rep(1:100, times=5)
pval <- runif(500, min=0, max=0.1)
df_gwas <- data.frame(CHR=chrom, BP=bp, P=pval, SNP=paste0("rs",1:500))

manhattan(df_gwas,
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-3),
          main="Przykładowy Manhattan Plot")
```



```{r}
library(ggplot2)
library(dplyr)
# Generowanie danych
set.seed(123)
chrom <- rep(1:5, each=100)
bp <- rep(1:100, times=5)
pval <- runif(500, min=0, max=0.1)
df_gwas <- data.frame(CHR=chrom, BP=bp, P=pval, SNP=paste0("rs",1:500))

# Logarytm dziesiętny z p-value
df_gwas <- df_gwas %>%
  mutate(logP = -log10(P))

# Tworzenie wykresu Manhattan z użyciem słupków
ggplot(df_gwas, aes(x=BP, y=logP, fill=factor(CHR))) +
  geom_bar(stat="identity", width=0.8) +  # Słupki
  facet_wrap(~ CHR, scales="free_x", nrow=1) +  # Oddzielne panele dla chromosomów
  geom_hline(yintercept = -log10(5e-8), color="red", linetype="dashed") +  # Linia genome-wide
  geom_hline(yintercept = -log10(1e-3), color="blue", linetype="dashed") +  # Linia suggestive
  labs(title = "Przykładowy Manhattan Plot - Słupki",
       x = "Pozycja na chromosomie",
       y = "-log10(p-value)") +
  theme_minimal()

```

# Venn Diagrams i UpSet Plots
## Graficzne przedstawienie relacji między różnymi zbiorami danych. Wykorzystują one nakładające się okręgi, gdzie każdy okrąg reprezentuje jeden zbiór, a części wspólne okręgów pokazują elementy, które należą do więcej niż jednego zbioru. W kontekście badań genetycznych wykresy Venn'a są szczególnie użyteczne do wizualizacji wspólnych i różniących się genów lub innych markerów genetycznych między różnymi grupami, takimi jak chorzy vs zdrowi, różne grupy czasowe lub różne warunki eksperymentalne.
## diagram Venna

```{r}
# install.packages("VennDiagram")
library(VennDiagram)
library(grid)  # do grid.draw

setA <- paste0("Gene", 1:10)
setB <- paste0("Gene", 6:15)

venn <- venn.diagram(
  x = list(A=setA, B=setB),
  filename = NULL,
  fill = c("skyblue", "pink"),
  alpha = c(0.5, 0.5),
  cex = 2,
  cat.cex = 2
)
grid.newpage()
grid.draw(venn)
```

## UpSet Plot

```{r}
# install.packages("UpSetR")
library(UpSetR)

listInput <- list(
  SetA = setA,
  SetB = setB,
  SetC = paste0("Gene", 8:12)
)

upset(fromList(listInput), 
      order.by = "freq", 
      main.bar.color = "steelblue",
      sets.bar.color = "tomato")
```


# Pathway and Annotation Plots (np. KEGG)

```{r}
BiocManager::install("pathview")
library(pathview)

# Przykładowe sztuczne dane (Entrez ID i logFC)
genelist <- c("1840"=2, "4609"=-1, "7124"=1)  # Entrez ID
# KEGG pathway (np. "hsa04110" - Cell cycle dla Homo sapiens)
# Zwróć uwagę, że w realnych analizach używasz prawidłowych ID
pv.out <- pathview(gene.data = genelist,
                   pathway.id = "hsa04110",
                   species = "hsa",
                   out.suffix="example")
```


# Drzewo filogenetyczne
## Wykres, który przedstawia ewolucyjny stosunek między różnymi organizmami, gatunkami, lub sekwencjami genetycznymi. Drzewa te pokazują, jak poszczególne jednostki są ze sobą powiązane przez wspólnych przodków, na podstawie podobieństw i różnic w ich cechach, takich jak sekwencje DNA, białek czy morfologia.

```{r}
# install.packages("ape")
library(ape)

tree <- rtree(10)  # losowe drzewo z 10 taksonami
plot(tree, main="Random Phylogenetic Tree")
```

# Synteny Plots (np. genoPlotR)
## Służą do wizualizacji układu genów na chromosomach różnych organizmów lub różnych regionów tego samego genomu, w celu porównania ich rozmieszczenia i identyfikowania zachowanych układów genów (tzw. syntenii). Syntenia odnosi się do sytuacji, w której geny o podobnych funkcjach lub sekwencjach są utrzymywane w tych samych pozycjach na chromosomach u różnych organizmów, co może wskazywać na ich wspólne pochodzenie lub ewolucyjne konserwowanie tych układów.

```{r}
install.packages("genoPlotR")
library(genoPlotR)
```

```{r}
# Przykładowe dane w pakiecie
data("barto", package="genoPlotR")

plot_gene_map(dna_segs = barto$dna_segs,
              comparisons = barto$comparisons,
              main = "Synteny plot - Bartonella genomes (genoPlotR)")
```


# Circos Plots
## Służą do wizualizacji złożonych relacji i zależności między danymi w obrębie pojedynczego genomu lub między różnymi genomami. Są szeroko stosowane w bioinformatyce, zwłaszcza do analizy danych genomowych, gdzie pomagają przedstawiać skomplikowane wzorce, takie jak syntenia, pozycje genów, zmiany strukturalne w genomach czy interakcje między różnymi regionami DNA.
```{r}
install.packages("circlize")
library(circlize)
```

Przygotwoanie zakresów sektorów
```{r}
library(dplyr)
library(circlize)

bed <- data.frame(
  chr   = c("chr1","chr1","chr2","chr2"),
  start = c(1, 50, 1, 50),
  end   = c(25, 75, 25, 75),
  value = c(10, 20, 5, 15)
)

# Grupujemy, żeby wyliczyć minimalny start i maksymalny end dla każdego chr
chr_ranges <- bed %>%
  group_by(chr) %>%
  summarise(
    min_start = min(start),
    max_end   = max(end)
  )

```


```{r}
install.packages("circlize")
library(circlize)

circos.clear()  # czyścimy stan przed nową inicjalizacją

circos.initialize(
  factors = chr_ranges$chr, 
  xlim    = cbind(chr_ranges$min_start, chr_ranges$max_end)
)

circos.trackPlotRegion(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    # Odczytujemy informację o sektorze
    sector.name = CELL_META$sector.index
    # Rysujemy napis na środku sektora
    circos.text(
      CELL_META$xcenter,
      0.5,
      sector.name,
      facing = "bending.inside"
    )
  }
)

for(i in seq_len(nrow(bed))) {
  # Wyciągamy chrom, start, end
  chr   = bed$chr[i]
  start = bed$start[i]
  end   = bed$end[i]
  val   = bed$value[i]

  # Rysujemy prostokąt w sektorze `chr`
  # "Wysokość" prostokąta zrobimy, np. od 0 do val/20 (tak, by coś było widać)
  circos.rect(
    xleft       = start, 
    ybottom     = 0, 
    xright      = end, 
    ytop        = val/20, 
    sector.index= chr,
    col         = "skyblue", 
    border      = "black"
  )
}

circos.clear()
```


# Ideograms (np. karyoploteR)
## Wizualne przedstawienie chromosomów, które są używane do ukazywania ich struktury i organizacji. W kontekście genomiki, ideogramy to schematyczne diagramy, które przedstawiają chromosomy w postaci prostych, zazwyczaj prostokątnych lub cylindrycznych bloków, z podziałami na różne regiony. Mogą one zawierać dodatkowe informacje o lokalizacji genów, markerów genetycznych, mutacjach, a także innych istotnych elementach w obrębie chromosomu.

```{r}
BiocManager::install("karyoploteR")
library(karyoploteR)

kp <- plotKaryotype(genome="hg19")  # lub inny genom
# Przykładowo zaznaczmy region na chr1
region <- toGRanges(data.frame(chr="chr1", start=1e6, end=2e6))
kpRect(kp, data=region, y0=0, y1=1, col="red", border=NA)
```

