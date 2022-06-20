# DAMMA
Distillation of Animal-associated Microorganisms' Metabolic Annotations

### Install DAMMA
```
install.packages("devtools") #only if devtools is not installed
library(devtools)
install_github("anttonalberdi/DAMMA")
```

### Run distillation
- **Required files / tables**: A) single table containing all MAG annotations; B) functions table provided with DAMMA.
- **Arguments**: magcol=[number of column containing MAG identifies], keggcol=[index(es) of column(s) containing KO codes (eg: K00169)], eccol=[index(es) of column(s) containing EC codes (eg: EC:3.2.4.15)], pepcol=[index(es) of column(s) containing peptidase codes (eg: C03H]

```
#Load annotations
library(data.table)
#annotations_file="/mydir/annotations.tsv"
#annotations <- fread(annotations_file)

#Load example annotations
data(damma_data)
head(annotations_example)

#Load functions
data(damma_data)
head(functions_table)

#Run distillation
distilled_table <- damma(annotations_example,functions_table,magcol=2,keggcol=9,eccol=c(10,19),pepcol=12)
```

### Aggregrate compounds
```
compounds_table <- aggregate_compounds(distilled_table,functions_table)
```

### Compound-level heatmap
```
library(ggplot2)
library(RColorBrewer)

compounds_table_df <- melt(compounds_table)
colnames(compounds_table_df) <- c("MAGs","Compounds","Fullness")
compounds_table_df2 <- merge(compounds_table_df,functions_table,by.x="Compounds",by.y="Compound")

ggplot(compounds_table_df2, aes(x=MAGs, y=Compounds, fill=Fullness, group=Function))+
  geom_tile(colour="white", size=0.1)+
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  #scale_fill_gradientn( colours = rev(c("#781a25","#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da","#f1faf0","#f4f4f4")))+
  scale_fill_gradientn(colours=brewer.pal(7, "YlGnBu"))+
  facet_grid(Function ~ ., scales = "free", space = "free")+
  theme_grey(base_size=8)+
  theme(strip.text.y = element_text(angle = 0),axis.text.x=element_blank())
```
