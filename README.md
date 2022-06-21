# DAMMA
Distillation of Animal-associated Microorganisms' Metabolic Annotations

### Install DAMMA
DAMMA can be directly installed from this Github repository using install_github() function.
```
install.packages("devtools") #only if devtools is not installed
library(devtools)
install_github("anttonalberdi/DAMMA")
```

### Input data
DAMMA only requires an annotations table containing MAG identifiers and annotations. Such tables are usually quite large, so using data.table is recommended for smooth processing of the data.

```
#Load annotations
library(data.table)
#annotations_file="/mydir/annotations.tsv"
#annotations <- fread(annotations_file)
```

DAMMA contains an example annotation table that can be loaded along with the functions table. These data are used in this documentation for showcasing DAMMA scripts.

```
#Load DAMMA library
library(DAMMA)

#Load DAMMA support data
data(damma_data)

#Visualise example annotations
head(annotations_example)

#Visualise functions table
head(functions_table)
```

### Run distillation
The damma() function requires specifying in which column(s) to find MAG identifiers and annotation data.
- magcol: number of column containing MAG identifies.
- keggcol: index(es) of column(s) containing KO codes (eg: K00169).
- eccol: index(es) of column(s) containing EC codes (eg: EC:3.2.4.15).
- pepcol: index(es) of column(s) containing peptidase codes (eg: C03H.

```
#Using example data
distilled_table <- damma(annotations_example,functions_table,magcol=2,keggcol=9,eccol=c(10,19),pepcol=12)
```

### Aggregrate raw distillates into 73 compounds

The raw fullness data can be aggregated to the compound level using the damma_compounds() function. These data are useful to obtain high-resolution functional information for statistical analyses that enable a large number of features.
```
distilled_table_compounds <- damma_compounds(distilled_table,functions_table)
```

### Convert compounds into a binary table

DAMMA offers the possibility to convert the fullness table into a binary format by setting a cutoff threshold.
```
distilled_table_compounds_bin <- damma_bin(distilled_table_compounds,threshold=0.9)
```

This conversion can account for the completeness of MAGs, and decrease the threshold when the completeness of the MAG is lower than 100%.
```
completeness_table <- cbind(genome=c("bin_m1.cct123","bin_m1.mtb106","bin_m1.mtb2","bin_m1.mxb107_sub","bin_m1.vmb35","bin_m1.vmb46","bin_m9.vmb60"),completeness=c(100,98,85.8,94.5,97,100,70))
distilled_table_compounds_bin <- damma_bin(distilled_table_compounds,threshold=0.9,completeness=completeness_table)
```

### Aggregrate compounds into 10 functions

The compounds data can be aggregated to the main functional levels using the damma_functions() function. These data are useful to obtain overall functional information for statistical analyses that required a reduced number of features.
```
distilled_table_functions <- damma_functions(distilled_table_compounds,functions_table,normalise=FALSE)
distilled_table_functions_bin <- aggregate_functions(distilled_table_compounds_bin,functions_table,normalise=TRUE)

```

### Compound-level heatmap

Functional data can be visualised as a heatmap using ggplot2.
```
library(data.table)
library(ggplot2)
library(RColorBrewer)

#Prepare input table
compounds_table_df <- melt(distilled_table_compounds)
colnames(compounds_table_df) <- c("MAGs","Compounds","Fullness")
compounds_table_df2 <- merge(compounds_table_df,functions_table,by.x="Compounds",by.y="Compound")
compounds_table_df2$Function <- as.factor(compounds_table_df2$Function)
compounds_table_df2$Function <- factor(compounds_table_df2$Function, levels=c("Dietary carbohydrate degradation","Dietary lipid degradation","Protein degradation","Mucin degradation","SCFA production","Organic anion production","Secondary bile acid production","Amino acid production","Amino acid derivative production","Vitamin production"))

#Plot heatmap
ggplot(compounds_table_df2, aes(x=MAGs, y=Compounds, fill=Fullness, group=Function))+
  geom_tile(colour="white", size=0.1)+
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  #scale_fill_gradientn(limits = c(0,1), colours = rev(c("#781a25","#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da","#f1faf0","#f4f4f4")))+
  scale_fill_gradientn(limits = c(0,1), colours=brewer.pal(7, "YlGnBu"))+
  facet_grid(Function ~ ., scales = "free", space = "free")+
  theme_grey(base_size=8)+
  theme(strip.text.y = element_text(angle = 0),axis.text.x=element_blank())
```

![Comparison of fullness and binary outcome of compound metabolism.](images/compound_fullness-binary.jpg)

### Function-level heatmap

Functional data can be visualised as a heatmap using ggplot2.
```
library(data.table)
library(ggplot2)
library(RColorBrewer)

#Prepare input table
functions_table_df <- melt(distilled_table_functions)
colnames(functions_table_df) <- c("MAGs","Functions","Index")
functions_table_df$Function <- as.factor(functions_table_df$Function)
functions_table_df$Function <- factor(functions_table_df$Function, levels=c("Dietary carbohydrate degradation","Dietary lipid degradation","Protein degradation","Mucin degradation","SCFA production","Organic anion production","Secondary bile acid production","Amino acid production","Amino acid derivative production","Vitamin production"))

#Plot heatmap
ggplot(functions_table_df, aes(x=MAGs, y=Functions, fill=Index))+
  geom_tile(colour="white", size=0.1)+
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  #scale_fill_gradientn(limits = c(0,1), colours = rev(c("#781a25","#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da","#f1faf0","#f4f4f4")))+
  scale_fill_gradientn(limits = c(0,1), colours=brewer.pal(7, "YlGnBu"))+
  facet_grid(Function ~ ., scales = "free", space = "free")+
  theme_grey(base_size=8)+
  theme(strip.text.y = element_text(angle = 0),axis.text.x=element_blank())
```
