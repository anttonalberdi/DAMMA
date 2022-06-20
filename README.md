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
annotations_file="/Users/anttonalberdi/annotations.tsv"
functions_file="/Users/anttonalberdi/github/DAMMA/functions/DAMMA_functions_20220619.tsv"

annotations <- fread(annotations_file)
functions <- read.table(functions_file,header=TRUE,row.names=1,sep="\t")
distilled_table <- fullness_matrix(annotations,functions,magcol=2,keggcol=9,eccol=c(10,19),pepcol=12)
```

### Transform fullness table
```
#Raw
heatmap(distilled_table, Colv=NA, scale="none")

#Binary
threshold=0.9
distilled_table_binary <- distilled_table
distilled_table_binary[distilled_table_binary >= threshold] <- 1
distilled_table_binary[distilled_table_binary < threshold] <- 0

#Binary corrected
threshold=0.9
completeness <- read.csv("/Users/anttonalberdi/genomeInfo.csv",row.names=1)
completeness <- completeness[rownames(completeness) %in% rownames(distilled_table),]
threshold_corrected <- completeness[,1]/100 * threshold

distilled_table_binary_corrected <- distilled_table
for(r in c(1:nrow(distilled_table_binary_corrected))){
  row <- distilled_table_binary_corrected[r, ]
  row[row >= threshold_corrected[r]] <- 1
  row[row < threshold_corrected[r]] <- 0
  distilled_table_binary_corrected[r, ] <- row
}
```

### Compound-level heatmap
```
library(ggplot2)
library(RColorBrewer)

distilled_table2=distilled_table
distilled_table2=distilled_table_binary_corrected

distilled_table2[is.na(distilled_table2)] <- 0
distilled_df <- melt(distilled_table2)
colnames(distilled_df) <- c("MAGs","Functions","Fullness")
distilled_df2 <- merge(distilled_df,functions,by.x="Functions",by.y="row.names")

ggplot(distilled_df2, aes(x=Functions, y=MAGs, fill=Fullness, group=Function))+
  geom_tile(colour="white", size=0.1)+
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  #scale_fill_gradientn( colours = rev(c("#781a25","#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da","#f1faf0","#f4f4f4")))+
  scale_fill_gradientn(colours=brewer.pal(7, "YlGnBu"))+
  facet_grid(~ Function, scales = "free", space = "free")
  theme_grey(base_size=8)
```

### Function-level heatmap
```
library(ggplot2)
library(RColorBrewer)

distilled_table2=distilled_table
distilled_table2=distilled_table_binary_corrected

distilled_table2[is.na(distilled_table2)] <- 0
distilled_df <- melt(distilled_table2)
colnames(distilled_df) <- c("MAGs","Functions","Fullness")
distilled_df2 <- merge(distilled_df,functions,by.x="Functions",by.y="row.names")

#Aggregate compounds
distilled_df3 <- aggregate(distilled_df2$Fullness, by=list(distilled_df2$MAGs, distilled_df2$Compound),FUN=sum)
colnames(distilled_df3) <- c("MAGs","Compounds","Fullness")
distilled_df3$Fullness[distilled_df3$Fullness > 1] <- 1
distilled_df3 <- merge(distilled_df3,unique(functions[,c(1,3)]),by.x="Compounds",by.y="Compound")
distilled_df3$Function <- as.factor(distilled_df3$Function)
distilled_df3$Function <- factor(distilled_df3$Function, levels=c("Dietary carbohydrate degradation","Dietary lipid degradation","Protein degradation","Mucin degradation","SCFA production","Organic anion production","Secondary bile acid production","Amino acid production","Amino acid derivative production","Vitamin production"))

ggplot(distilled_df3, aes(x=MAGs, y=Compounds, fill=Fullness, group=Function))+
  geom_tile(colour="white", size=0.1)+
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  #scale_fill_gradientn( colours = rev(c("#781a25","#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da","#f1faf0","#f4f4f4")))+
  scale_fill_gradientn(colours=brewer.pal(7, "YlGnBu"))+
  facet_grid(Function ~ ., scales = "free", space = "free")+
  theme_grey(base_size=8)+
  theme(strip.text.y = element_text(angle = 0),axis.text.x=element_blank())

```
