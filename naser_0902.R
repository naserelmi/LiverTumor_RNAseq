

library(openxlsx)
df <- read.xlsx("Expression_Profile.GRCh38.gene.xlsx",
                sheet = 1,
                rowNames = TRUE,
                startRow = 1,
                detectDates = FALSE,
                skipEmptyRows = TRUE)

#find index of repetitive gene symbols
repetitive_indices <- which(duplicated(df$Gene_Symbol) | duplicated(df$Gene_Symbol, fromLast = TRUE))
print(repetitive_indices)
df_nr<-df[-c(repetitive_indices),]


#Select raw data
df_nr_raw<-df_nr[,1:21]


#select nessecarry colums (Gene Symbol+ raw count)

data_raw<-cbind(df_nr_raw$Gene_Symbol,df_nr_raw[,9:dim(df_nr_raw)[2]])
colnames(data_raw)[1]<-'Gene_Symbol'
raw_data<-data_raw[,-1]
rownames(raw_data)<-data_raw[,1]

countdata<-raw_data
# Eliminate the genes(rows) with low count based on CPM 
library(edgeR)
myCPM <- cpm(countdata)
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1], xlab="CPM", ylab="Raw Count", main=colnames(myCPM)[1], 
     ylim=c(0,50), xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.25)
abline(v=1)
abline(h=15)
# We set CPM>1 to get at least 15 reads for each gene in at least 3 samples
isexpr <- rowSums(myCPM>0.20) >= 5
summary(isexpr)
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[isexpr,]
dim(countdata)
dim(counts.keep)

df <-na.omit(counts.keep)
ex <- data.matrix(df)
head(ex)
logcounts <- log(ex+1,2)
# Box and whisker plot for unnormalized sample gene expression 
boxplot(logcounts,
        ylab= 'Expression',
        varwidth=TRUE,
        notch=TRUE,
        border=TRUE,outline = FALSE)
library(reshape)
melted_data <- melt(logcounts)
colnames(melted_data)<-c('Gene_Symbol','Sample_Name','Expression_Value')
melted_data$`Sample Name` <- gsub("_Read_Count", "", melted_data$`Sample Name`)
library(ggplot2)
ggplot(melted_data, aes(x = `Sample_Name`, y = `Expression_Value`, fill = Sample_Name)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Sample", y = "Expression", title = "Box Plot of Gene/Transcript Expression")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(
    axis.text.x = element_text(size = 12),  # Adjust x-axis label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis label font size
    axis.title = element_text(size = 14)    # Adjust axis title font size
  )+coord_cartesian(ylim = c(0, 17))  # Set the y-axis limits

#normalize with DEseq2
library(DESeq2)

# Define the conditions for each sample
meta_data <- data.frame(condition = c('Tumor',rep('Adjacent_A',8),rep('Adjacent_R',4)))
meta_data$condition <- factor(meta_data$condition)

dds <- DESeqDataSetFromMatrix(countData = ex,
                              colData = meta_data,
                              design = ~ condition)
dds <- DESeq(dds)

dds <- estimateSizeFactors(dds)
dds <- rlog(dds) #use rlog!
normalized_transformed_data_rlog <- assay(dds)
melted_normalized_data <- melt(normalized_transformed_data_rlog)
colnames(melted_normalized_data)<-c('Gene_Symbol','Sample_Name','Expression_Value')
melted_normalized_data$Sample_Name <- gsub("_Read_Count", "", melted_normalized_data$Sample_Name)

ggplot(melted_normalized_data, aes(x = Sample_Name, y = Expression_Value, fill = Sample_Name)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Sample", y = "Expression", title = "Box Plot of Normalized Gene/Transcript Expression")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(
    axis.text.x = element_text(size = 12),  # Adjust x-axis label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis label font size
    axis.title = element_text(size = 14)    # Adjust axis title font size
  )+coord_cartesian(ylim = c(0, 17))  # Set the y-axis limits


#normalize and differential gene expression using edgeR exacttest






















#load again raw data
df_edger <- read.xlsx("Expression_Profile.GRCh38.gene.xlsx",
                sheet = 1,
                rowNames = TRUE,
                startRow = 1,
                detectDates = FALSE,
                skipEmptyRows = TRUE)

df_edger<-df_edger[,1:21]
df_edgee<-cbind(df_edger$Gene_Symbol,df_edger$gene_biotype,df_edger$HGNC,df_edger$MIM,df_edger[,9:dim(df_edger)[2]])
colnames(df_edgee)[1:4]=c('Gene_Symbol','Gene_biotype','HGCN','MIM')
colnames(df_edgee) <- gsub("_Read_Count", "", colnames(df_edgee))

y <- DGEList(counts=df_edgee[,5:dim(df_edgee)[2]], genes=df_edgee[,1:4])
keep <- filterByExpr(y)
table(keep)
keep
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y)
y$samples

o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Gene_Symbol)
y <- y[!d,]
nrow(y)

#calculate dispersion
y<-normLibSizes(y)
y$samples

#design matrix for comparison
tissue<-factor(colnames(y))
tissue_cluster<-factor(c('c0','c0','c1','c1','c2','c3','c3','c3','c3','c3','c1','c2','c2'))
data.frame(Sample=colnames(y),tissue,tissue_cluster)
design<-model.matrix(~0+tissue_cluster)
rownames(design)<-colnames(y)
install.packages('statmod')
library(statmod)
library(limma)
library(edgeR)
#dispersion estimation
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

plotBCV(y)

design_tissue<-model.matrix(~0+tissue,data=y$samples)


y$samples$group<-tissue
levels(y$samples$group)
#pair=c('a','b') means b versus a

#get entrez id
BiocManager::install("org.Hs.eg.db")  # Install the org.Hs.eg.db package
library(org.Hs.eg.db)
library(dplyr)
library(AnnotationDbi)
library(biomaRt)

# Vector of HGNC gene symbols
hgnc_symbols <- y$genes$Gene_Symbol # Add your gene symbols here
# Convert HGNC symbols to Entrez Gene IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys = hgnc_symbols, keytype = "SYMBOL", column = "ENTREZID")
unique_entrez_ids <- sapply(entrez_ids, function(x) x[1])
result_df <- data.frame(GeneSymbol = hgnc_symbols, EntrezID = unique_entrez_ids)
y$gene$entrezID<-result_df$EntrezID
y$genes$entrezID<-result_df$EntrezID



#BiocManager::install("GO.db")
#BiocManager::install("clusterProfiler")
library(GO.db)
library(clusterProfiler)

####get output for differnt samples
###test1: TA2mm vs CENTER
et_A2_C <- exactTest(y, pair = c('CENTER','TA2mm'),dispersion=y$common.dispersion)
et_A2_C_df<-as.data.frame(topTags(et_A2_C,n=dim(et_A2_C)[1]))
up_A2_C <- et_A2_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_A2_C <- et_A2_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A2_C,p.value = 0.01,lfc=3))
plotMD(et_A2_C)
abline(h=c(-3, 3), col="blue")
go_up_A2_C <- goana(as.vector(up_A2_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A2_C<-data.frame(topGO(go_up_A2_C, ontology = "BP",p.value = 0.05))
#####-----------------------------------------------------------
##test2 TA4mm vs Center

et_A4_C <- exactTest(y, pair = c('CENTER','TA4mm'),dispersion=y$common.dispersion)
et_A4_C_df<-as.data.frame(topTags(et_A4_C,n=dim(et_A4_C)[1]))
up_A4_C <- et_A4_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_A4_C <- et_A4_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A4_C,p.value = 0.01,lfc=3))
plotMD(et_A4_C)
abline(h=c(-3, 3), col="blue")
go_up_A4_C <- goana(as.vector(up_A4_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A4_C<-data.frame(topGO(go_up_A4_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test3 TA6mm vs Center

et_A6_C <- exactTest(y, pair = c('CENTER','TA6mm'),dispersion=y$common.dispersion)
et_A6_C_df<-as.data.frame(topTags(et_A6_C,n=dim(et_A6_C)[1]))
up_A6_C <- et_A6_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_A6_C <- et_A6_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A6_C,p.value = 0.01,lfc=3))
plotMD(et_A6_C)
abline(h=c(-3, 3), col="blue")
go_up_A6_C <- goana(as.vector(up_A6_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A6_C<-data.frame(topGO(go_up_A6_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test4 TA10mm vs Center

et_A10_C <- exactTest(y, pair = c('CENTER','TA10mm'),dispersion=y$common.dispersion)
et_A10_C_df<-as.data.frame(topTags(et_A10_C,n=dim(et_A10_C)[1]))
up_A10_C <- et_A10_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_A10_C <- et_A10_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A10_C,p.value = 0.01,lfc=3))
plotMD(et_A10_C)
abline(h=c(-3, 3), col="blue")
go_up_A10_C <- goana(as.vector(up_A10_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A10_C<-data.frame(topGO(go_up_A10_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test5 TA14mm vs Center

et_A14_C <- exactTest(y, pair = c('CENTER','TA14mm'),dispersion=y$common.dispersion)
et_A14_C_df<-as.data.frame(topTags(et_A14_C,n=dim(et_A14_C)[1]))
up_A14_C <- et_A14_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_A14_C <- et_A14_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A14_C,p.value = 0.01,lfc=3))
plotMD(et_A14_C)
abline(h=c(-3, 3), col="blue")
go_up_A14_C <- goana(as.vector(up_A14_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A14_C<-data.frame(topGO(go_up_A14_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test6 TA20mm vs Center

et_A20_C <- exactTest(y, pair = c('CENTER','TA20mm'),dispersion=y$common.dispersion)
et_A20_C_df<-as.data.frame(topTags(et_A20_C,n=dim(et_A20_C)[1]))
up_A20_C <- et_A20_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_A20_C <- et_A20_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A20_C,p.value = 0.01,lfc=3))
plotMD(et_A20_C)
abline(h=c(-3, 3), col="blue")
go_up_A20_C <- goana(as.vector(up_A20_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A20_C<-data.frame(topGO(go_up_A20_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test7 TA24mm vs Center

et_A24_C <- exactTest(y, pair = c('CENTER','TA24mm'),dispersion=y$common.dispersion)
et_A24_C_df<-as.data.frame(topTags(et_A24_C,n=dim(et_A24_C)[1]))
up_A24_C <- et_A24_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_A24_C <- et_A24_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A24_C,p.value = 0.01,lfc=3))
plotMD(et_A24_C)
abline(h=c(-3, 3), col="blue")
go_up_A24_C <- goana(as.vector(up_A24_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A24_C<-data.frame(topGO(go_up_A24_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test7 TA30mm vs Center

et_A30_C <- exactTest(y, pair = c('CENTER','TA30mm'),dispersion=y$common.dispersion)
et_A30_C_df<-as.data.frame(topTags(et_A30_C,n=dim(et_A30_C)[1]))
up_A30_C <- et_A30_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_A30_C <- et_A30_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A30_C,p.value = 0.01,lfc=3))
plotMD(et_A30_C)
abline(h=c(-3, 3), col="blue")
go_up_A30_C <- goana(as.vector(up_A30_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A30_C<-data.frame(topGO(go_up_A30_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test8 TAR2mm vs Center

et_AR2_C <- exactTest(y, pair = c('CENTER','TAR2mm'),dispersion=y$common.dispersion)
et_AR2_C_df<-as.data.frame(topTags(et_AR2_C,n=dim(et_AR2_C)[1]))
up_AR2_C <- et_AR2_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_AR2_C <- et_AR2_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_AR2_C,p.value = 0.01,lfc=3))
plotMD(et_AR2_C)
abline(h=c(-3, 3), col="blue")
go_up_AR2_C <- goana(as.vector(up_AR2_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_AR2_C<-data.frame(topGO(go_up_AR2_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test9 TAR4mm vs Center

et_AR4_C <- exactTest(y, pair = c('CENTER','TAR4mm'),dispersion=y$common.dispersion)
et_AR4_C_df<-as.data.frame(topTags(et_AR4_C,n=dim(et_AR4_C)[1]))
up_AR4_C <- et_AR4_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_AR4_C <- et_AR4_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_AR4_C,p.value = 0.01,lfc=3))
plotMD(et_AR4_C)
abline(h=c(-3, 3), col="blue")
go_up_AR4_C <- goana(as.vector(up_AR4_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_AR4_C<-data.frame(topGO(go_up_AR4_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test10 TAR6mm vs Center

et_AR6_C <- exactTest(y, pair = c('CENTER','TAR6mm'),dispersion=y$common.dispersion)
et_AR6_C_df<-as.data.frame(topTags(et_AR6_C,n=dim(et_AR6_C)[1]))
up_AR6_C <- et_AR6_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_AR6_C <- et_AR6_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_AR6_C,p.value = 0.01,lfc=3))
plotMD(et_AR6_C)
abline(h=c(-3, 3), col="blue")
go_up_AR6_C <- goana(as.vector(up_AR6_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_AR6_C<-data.frame(topGO(go_up_AR6_C, ontology = "BP",p.value = 0.05))

#####-----------------------------------------------------------
##test11 TAR10mm vs Center
et_AR10_C <- exactTest(y, pair = c('CENTER','TAR10mm'),dispersion=y$common.dispersion)
et_AR10_C_df<-as.data.frame(topTags(et_AR10_C,n=dim(et_AR10_C)[1]))
up_AR10_C <- et_AR10_C_df %>% filter(FDR<0.01 & (logFC) >3)
down_AR10_C <- et_AR10_C_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_AR10_C,p.value = 0.01,lfc=3))
plotMD(et_AR10_C)
abline(h=c(-3, 3), col="blue")
go_up_AR10_C <- goana(as.vector(up_AR10_C$entrezID), species="Hs",geneid = entrezID)
go_BP_up_AR10_C<-data.frame(topGO(go_up_AR10_C, ontology = "BP",p.value = 0.05))

#########################################################################
####comparison against A30 (the most normal sample)
#test1 Center vs A30

et_CENTER_A30 <- exactTest(y, pair = c('TA30mm','CENTER'),dispersion=y$common.dispersion)
et_CENTER_A30_df<-as.data.frame(topTags(et_CENTER_A30,n=dim(et_CENTER_A30)[1]))
up_CENTER_A30 <- et_CENTER_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_CENTER_A30 <- et_CENTER_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_CENTER_A30,p.value = 0.01,lfc=3))
plotMD(et_CENTER_A30)
abline(h=c(-3, 3), col="blue")
go_up_CENTER_A30 <- goana(as.vector(up_CENTER_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_CENTER_A30<-data.frame(topGO(go_up_CENTER_A30, ontology = "BP",p.value = 0.05))

####comparison against A30 (the most normal sample)
#test1 A2 vs A30

et_A2_A30 <- exactTest(y, pair = c('TA30mm','TA2mm'),dispersion=y$common.dispersion)
et_A2_A30_df<-as.data.frame(topTags(et_A2_A30,n=dim(et_A2_A30)[1]))
up_A2_A30 <- et_A2_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_A2_A30 <- et_A2_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A2_A30,p.value = 0.01,lfc=3))
plotMD(et_A2_A30)
abline(h=c(-3, 3), col="blue")
go_up_A2_A30 <- goana(as.vector(up_A2_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A2_A30<-data.frame(topGO(go_up_A2_A30, ontology = "BP",p.value = 0.05))

#test1 A4 vs A30

et_A4_A30 <- exactTest(y, pair = c('TA30mm','TA4mm'),dispersion=y$common.dispersion)
et_A4_A30_df<-as.data.frame(topTags(et_A4_A30,n=dim(et_A4_A30)[1]))
up_A4_A30 <- et_A4_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_A4_A30 <- et_A4_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A4_A30,p.value = 0.01,lfc=3))
plotMD(et_A4_A30)
abline(h=c(-3, 3), col="blue")
go_up_A4_A30 <- goana(as.vector(up_A4_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A4_A30<-data.frame(topGO(go_up_A4_A30, ontology = "BP",p.value = 0.05))

#test1 A6 vs A30

et_A6_A30 <- exactTest(y, pair = c('TA30mm','TA6mm'),dispersion=y$common.dispersion)
et_A6_A30_df<-as.data.frame(topTags(et_A6_A30,n=dim(et_A6_A30)[1]))
up_A6_A30 <- et_A6_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_A6_A30 <- et_A6_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A6_A30,p.value = 0.01,lfc=3))
plotMD(et_A6_A30)
abline(h=c(-3, 3), col="blue")
go_up_A6_A30 <- goana(as.vector(up_A6_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A6_A30<-data.frame(topGO(go_up_A6_A30, ontology = "BP",p.value = 0.05))


#test1 A10 vs A30

et_A10_A30 <- exactTest(y, pair = c('TA30mm','TA10mm'),dispersion=y$common.dispersion)
et_A10_A30_df<-as.data.frame(topTags(et_A10_A30,n=dim(et_A10_A30)[1]))
up_A10_A30 <- et_A10_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_A10_A30 <- et_A10_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A10_A30,p.value = 0.01,lfc=3))
plotMD(et_A10_A30)
abline(h=c(-3, 3), col="blue")
go_up_A10_A30 <- goana(as.vector(up_A10_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A10_A30<-data.frame(topGO(go_up_A10_A30, ontology = "BP",p.value = 0.05))

#test1 A14 vs A30

et_A14_A30 <- exactTest(y, pair = c('TA30mm','TA14mm'),dispersion=y$common.dispersion)
et_A14_A30_df<-as.data.frame(topTags(et_A14_A30,n=dim(et_A14_A30)[1]))
up_A14_A30 <- et_A14_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_A14_A30 <- et_A14_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A14_A30,p.value = 0.01,lfc=3))
plotMD(et_A14_A30)
abline(h=c(-3, 3), col="blue")
go_up_A14_A30 <- goana(as.vector(up_A14_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A14_A30<-data.frame(topGO(go_up_A14_A30, ontology = "BP",p.value = 0.05))

#test1 A20 vs A30

et_A20_A30 <- exactTest(y, pair = c('TA30mm','TA20mm'),dispersion=y$common.dispersion)
et_A20_A30_df<-as.data.frame(topTags(et_A20_A30,n=dim(et_A20_A30)[1]))
up_A20_A30 <- et_A20_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_A20_A30 <- et_A20_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A20_A30,p.value = 0.01,lfc=3))
plotMD(et_A20_A30)
abline(h=c(-3, 3), col="blue")
go_up_A20_A30 <- goana(as.vector(up_A20_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A20_A30<-data.frame(topGO(go_up_A20_A30, ontology = "BP",p.value = 0.05))


#test1 A24 vs A30

et_A24_A30 <- exactTest(y, pair = c('TA30mm','TA24mm'),dispersion=y$common.dispersion)
et_A24_A30_df<-as.data.frame(topTags(et_A24_A30,n=dim(et_A24_A30)[1]))
up_A24_A30 <- et_A24_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_A24_A30 <- et_A24_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_A24_A30,p.value = 0.01,lfc=3))
plotMD(et_A24_A30)
abline(h=c(-3, 3), col="blue")
go_up_A24_A30 <- goana(as.vector(up_A24_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_A24_A30<-data.frame(topGO(go_up_A24_A30, ontology = "BP",p.value = 0.05))

#test1 AR2 vs A30

et_AR2_A30 <- exactTest(y, pair = c('TA30mm','TAR2mm'),dispersion=y$common.dispersion)
et_AR2_A30_df<-as.data.frame(topTags(et_AR2_A30,n=dim(et_AR2_A30)[1]))
up_AR2_A30 <- et_AR2_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_AR2_A30 <- et_AR2_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_AR2_A30,p.value = 0.01,lfc=3))
plotMD(et_AR2_A30)
abline(h=c(-3, 3), col="blue")
go_up_AR2_A30 <- goana(as.vector(up_AR2_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_AR2_A30<-data.frame(topGO(go_up_AR2_A30, ontology = "BP",p.value = 0.05))

#test1 AR4 vs A30

et_AR4_A30 <- exactTest(y, pair = c('TA30mm','TAR4mm'),dispersion=y$common.dispersion)
et_AR4_A30_df<-as.data.frame(topTags(et_AR4_A30,n=dim(et_AR4_A30)[1]))
up_AR4_A30 <- et_AR4_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_AR4_A30 <- et_AR4_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_AR4_A30,p.value = 0.01,lfc=3))
plotMD(et_AR4_A30)
abline(h=c(-3, 3), col="blue")
go_up_AR4_A30 <- goana(as.vector(up_AR4_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_AR4_A30<-data.frame(topGO(go_up_AR4_A30, ontology = "BP",p.value = 0.05))

#test1 AR6 vs A30

et_AR6_A30 <- exactTest(y, pair = c('TA30mm','TAR6mm'),dispersion=y$common.dispersion)
et_AR6_A30_df<-as.data.frame(topTags(et_AR6_A30,n=dim(et_AR6_A30)[1]))
up_AR6_A30 <- et_AR6_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_AR6_A30 <- et_AR6_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_AR6_A30,p.value = 0.01,lfc=3))
plotMD(et_AR6_A30)
abline(h=c(-3, 3), col="blue")
go_up_AR6_A30 <- goana(as.vector(up_AR6_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_AR6_A30<-data.frame(topGO(go_up_AR6_A30, ontology = "BP",p.value = 0.05))

#test1 AR10 vs A30

et_AR10_A30 <- exactTest(y, pair = c('TA30mm','TAR10mm'),dispersion=y$common.dispersion)
et_AR10_A30_df<-as.data.frame(topTags(et_AR10_A30,n=dim(et_AR10_A30)[1]))
up_AR10_A30 <- et_AR10_A30_df %>% filter(FDR<0.01 & (logFC) >3)
down_AR10_A30 <- et_AR10_A30_df %>% filter(FDR<0.01 & logFC < -3)
summary(decideTests(et_AR10_A30,p.value = 0.01,lfc=3))
plotMD(et_AR10_A30)
abline(h=c(-3, 3), col="blue")
go_up_AR10_A30 <- goana(as.vector(up_AR10_A30$entrezID), species="Hs",geneid = entrezID)
go_BP_up_AR10_A30<-data.frame(topGO(go_up_AR10_A30, ontology = "BP",p.value = 0.05))

save.image(file = "my_workspace_20230829.RData")
load('C:/Users/User/Desktop/Dr_samanian/my_workspace_20230829.RData')



#this does not work for samples with one replicate
#rownames(design)<-colnames(y)
# "CENTER"  "TA2mm"   "TA4mm"   "TA6mm"   "TA10mm"  "TA14mm"  "TA20mm"  "TA24mm"  "TA30mm"  "TAR2mm"  "TAR4mm"  "TAR6mm"  "TAR10mm"
#fit <- glmQLFit(y, design)
#fit<-glmQLFit(y, design=design, dispersion=y$common.dispersion)

DEG_list <- list()
DEG_list_df<-list()
up_list<-list()
down_list<-list()
dim1=length(seq(1,4,0.25))
dim2=length(colnames(y))
dim3=length(colnames(y))
dim1_names <- as.character(seq(1,4,0.25))
dim2_names <- colnames(y)
dim3_names <- colnames(y)
# Assign names to the dimensions using the dimnames attribute
matrix_up <- array(NA, dim = c(dim1, dim2, dim3))
dimnames(matrix_up) <- list(dim1_names, dim2_names, dim3_names)
matrix_down <- array(NA, dim = c(dim1, dim2, dim3))
dimnames(matrix_down) <- list(dim1_names, dim2_names, dim3_names)
matrix_sum <- array(NA, dim = c(dim1, dim2, dim3))
dimnames(matrix_sum) <- list(dim1_names, dim2_names, dim3_names)
lfc_count=1
for (lfc in seq(1,4,0.25)){
  x_item_count=1
  for (x_item in colnames(y)){
    y_item_count=1
    for (y_item in colnames(y)){
      
      file_et_name<-paste0(x_item,'_',y_item)
      
      if (is.null(DEG_list[[as.character(lfc)]])) {
        DEG_list[[as.character(lfc)]] <- list()
        DEG_list_df[[as.character(lfc)]] <- list()
        up_list[[as.character(lfc)]] <- list()
        down_list[[as.character(lfc)]] <- list()
        
      }
      if (is.null(DEG_list[[as.character(lfc)]][[file_et_name]])) {
        DEG_list[[as.character(lfc)]][[file_et_name]] <- list()
        DEG_list_df[[as.character(lfc)]][[file_et_name]] <- list()
        up_list[[as.character(lfc)]][[file_et_name]] <- list()
        down_list[[as.character(lfc)]][[file_et_name]] <- list()
        
      }
      
      # Perform exactTest and store the result in the list
      DEG_list[[as.character(lfc)]][[file_et_name]] <- exactTest(y, pair = c(x_item,y_item),dispersion=y$common.dispersion)
      DEG_list_df[[as.character(lfc)]][[file_et_name]]<-as.data.frame(topTags( DEG_list[[as.character(lfc)]][[file_et_name]],n=dim( DEG_list[[as.character(lfc)]][[file_et_name]])[1]))
      up_file <- DEG_list_df[[as.character(lfc)]][[file_et_name]] %>% filter(FDR<0.01 & (logFC) >3)
      down_file <- DEG_list_df[[as.character(lfc)]][[file_et_name]] %>% filter(FDR<0.01 & logFC < -3)
      up_list[[as.character(lfc)]][[file_et_name]]<-up_file
      down_list[[as.character(lfc)]][[file_et_name]]<-down_file
      deg_summary<-summary(decideTests(DEG_list[[as.character(lfc)]][[file_et_name]],p.value = 0.01,lfc=lfc))
      num_up<-deg_summary[1]
      num_down<-deg_summary[3]
      #plotMD(DEG_list[[as.character(lfc)]][[file_et_name]])
      #abline(h=c(-lfc, lfc), col="blue")
      
      matrix_up[lfc_count, x_item_count, y_item_count] <- num_up
      matrix_down[lfc_count, x_item_count, y_item_count] <- num_down
      matrix_sum[lfc_count, x_item_count, y_item_count] <- num_up+num_down
      
      y_item_count<-y_item_count+1
      

      
      
    }
    x_item_count<-x_item_count+1
  }
  lfc_count<-lfc_count+1
}




##create a datframe for plotting

col_names <- c("Num_up", "Num_down", "Num_sum",'distance','distance_class','lfc_value')
# Create an empty data frame with predefined column names
data_frame <- data.frame(matrix(nrow =0 , ncol = length(col_names)))
colnames(data_frame) <- col_names
lfc_val=seq(1,4,0.25)
distance=c(0,2,4,6,10,14,20,24,30,2,4,6,10)
distane_cl<-colnames(y)
for (lfc in 1:length(seq(1,4,0.25))){
  
  # Create an empty data frame with predefined column names
  temp_df <- data.frame(matrix(nrow =(dim(matrix_sum)[1])^2 , ncol = length(col_names)))
  colnames(temp_df) <- col_names
  temp_df$Num_up<-matrix_up[lfc,1,]
  temp_df$Num_down<-matrix_down[lfc,1,]
  temp_df$Num_sum<-matrix_sum[lfc,1,]
  temp_df$distance<-distance
  temp_df$distance_class<-distane_cl
  temp_df$lfc_value<-rep(lfc_val[lfc],length(seq(1,4,0.25)))
  data_frame<-rbind(data_frame,temp_df)
  temp_df<-data.frame()
  
}

data_frame_center<-data_frame

color_palette <- colorRampPalette(c("darkred", "lightblue"))

# Generate category_colors for 13 categories
num_categories <- length(seq(1,4,0.25))
category_colors <- color_palette(num_categories)


library(ggplot2)
library(ggplotify)
library(plotly)
#up
boxplot_up_center<-ggplot(data_frame_center, aes(x = distance_class, y = Num_up, fill = distance_class)) +
  geom_boxplot() +
  scale_fill_manual(values = category_colors) +
  labs(title = "Box Plot of Up regulated genes in different distances against the CENTER", x = "Sample Name", y = "Number of UP DEGs") +
  scale_x_discrete(limits = unique(data_frame_center$distance_class))+guides(fill = "none")+  # Remove the legend for fill
# Maintain the original order
  theme_minimal()
print(boxplot_up_center)
coord_fixed()
ggplotly(boxplot_up_center)

#down
boxplot_down_center<-ggplot(data_frame_center, aes(x = distance_class, y = Num_down, fill = distance_class)) +
  geom_boxplot() +
  scale_fill_manual(values = category_colors) +
  labs(title = "Box Plot of Down regulated genes in different distances against the CENTER", x = "Sample Name", y = "Number of Down DEGs") +
  scale_x_discrete(limits = unique(data_frame_center$distance_class))+guides(fill = "none")+  # Remove the legend for fill
  # Maintain the original order
  theme_minimal()
print(boxplot_down_center)
coord_fixed()
ggplotly(boxplot_down_center)

#sum

boxplot_sum_center<-ggplot(data_frame_center, aes(x = distance_class, y = Num_sum, fill = distance_class)) +
  geom_boxplot() +
  scale_fill_manual(values = category_colors) +
  labs(title = "Box Plot of Number of DEG genes in different distances against the CENTER", x = "Sample Name", y = "Number of DEGs") +
  scale_x_discrete(limits = unique(data_frame_center$distance_class))+guides(fill = "none")+  # Remove the legend for fill
  # Maintain the original order
  theme_minimal()
print(boxplot_sum_center)
coord_fixed()
ggplotly(boxplot_sum_center)

library(patchwork)

plot1 <- ggplot(data = data.frame(x, y1), aes(x, y1)) + geom_line() + ggtitle("Plot 1")
plot2 <- ggplot(data = data.frame(x, y2), aes(x, y2)) + geom_line() + ggtitle("Plot 2")
plot3 <- ggplot(data = data.frame(x, y3), aes(x, y3)) + geom_line() + ggtitle("Plot 3")

data_filtered=data_frame_center[data_frame_center$lfc_value==1,]
plot_lfc1_sum<-ggplot(data_filtered, aes(x = data_filtered$distance_class, y = data_filtered$Num_sum)) +
  geom_point(size = 3) +
  labs(title = "Dot Plot, logfc>1", x = "Sample Name", y = "Number of DEGs") +
  # Maintain the original order
  theme_minimal()+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_x_discrete(limits = unique(data_filtered$distance_class))

plot_lfc2_sum<-ggplot(data_frame_center[data_frame_center$lfc_value=='2',], aes(x = distance_class, y = Num_sum)) +
  geom_point(size = 3) +
  labs(title = "Dot Plot, logfc>2", x = "Sample Name", y = "Number of DEGs") +
  # Maintain the original order
  theme_minimal()+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_x_discrete(limits = unique(data_filtered$distance_class))


plot_lfc3_sum<-ggplot(data_frame_center[data_frame_center$lfc_value=='3',], aes(x = distance_class, y = Num_sum)) +
  geom_point(size = 3) +
  labs(title = "Dot Plot, logfc>3", x = "Sample Name", y = "Number of DEGs") +
  # Maintain the original order
  theme_minimal()+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_x_discrete(limits = unique(data_filtered$distance_class))


# Arrange the plots in a subplot
subplot_plots <- plot_lfc1_sum | plot_lfc2_sum | plot_lfc3_sum

print(subplot_plots)
coord_fixed()
ggplotly(plot_lfc2_sum)




# Now normalize the data for clustering and correlation analysis
normalized_expression_edgeR <- cpm(y,log = TRUE)

#pca plot on normalized data







#yy<-DGEList(counts = normalized_expression_edgeR,genes=y$genes)
#yy <- normLibSizes(yy)
#yy$samples
#yy<-calcNormFactors(yy)
#normalized_expression_edgeR2 <- cpm(yy,log = TRUE)
plotMDS(y)
#plotMDS(yy)


pca_result <- prcomp(t(normalized_expression_edgeR))

# Create a data frame with PCA results
pca_df <- data.frame(PC1 = pca_result$x[, 1],
                     PC2 = pca_result$x[, 2],
                     Group = colnames(normalized_expression_edgeR))



# Create a PCA plot using ggplot2
pca_plot <- ggplot(data = pca_df, aes(x = PC1, y = PC2, color = Group, label = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Group), nudge_x = 8, nudge_y = -8, size = 3) +  # Add labels
  labs(title = "PCA Plot",
       x = "PC1",
       y = "PC2",
       color = "Group") +
  theme_minimal()
# Print the PCA plot
print(pca_plot)
coord_fixed()
ggplotly(pca_plot)

library(Rtsne)

perplexity_value <- 3  # Adjust as needed

# Perform t-SNE
set.seed(123)  # For reproducibility
tsne_result <- Rtsne(t(normalized_expression_edgeR),perplexity = perplexity_value,theta=0.01)

# Create a data frame with t-SNE results
tsne_df <- data.frame(X = tsne_result$Y[, 1],
                      Y = tsne_result$Y[, 2],
                      Group = colnames(normalized_expression_edgeR))

tsne_plot <- ggplot(data = tsne_df, aes(x = X, y = Y, color = Group, label = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Group), nudge_x = 3, nudge_y = -3, size = 3) +  # Add labels
  labs(title = "tsne Plot",
       x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2",
       color = "Group") +
  theme_minimal()
# Print the tsne plot
print(tsne_plot)
coord_fixed()
ggplotly(tsne_plot)

# Create a data frame with umap results

library(umap)
# Perform UMAP
# Determine an appropriate number of neighbors (e.g., 5-50)
num_neighbors <- 5  # Adjust as needed
custom.settings = umap.defaults
custom.settings$n_neighbors = num_neighbors
umap_result <- umap(t(normalized_expression_edgeR),config=custom.settings)


# Create a data frame with UMAP results
umap_df <- data.frame(X = umap_result$layout[, 1],
                      Y = umap_result$layout[, 2],
                      Group = colnames(normalized_expression_edgeR))


umap_plot <- ggplot(data = umap_df, aes(x = X, y = Y, color = Group, label = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Group), nudge_x = 0.01, nudge_y = -.2, size = 3) +  # Add labels
  labs(title = "tsne Plot",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2",
       color = "Group") +
  theme_minimal()
# Print the tsne plot
print(umap_plot)
coord_fixed()
ggplotly(umap_plot)

save.image('C:/Users/User/Desktop/Dr_samanian/my_data_2023_8_30.RData')
setwd('C:/Users/User/Desktop/Dr_samanian')
load('Dr_samanian_20230817.RData')



#heatmap
library(pheatmap)
# Calculate variance across genes
gene_variance <- apply(normalized_expression_edgeR, 1, var)

# Select the top 100 variable genes
top_variable_genes <- names(sort(gene_variance, decreasing = TRUE))[1:100]
heatmap_data <- normalized_expression_edgeR[top_variable_genes, ]

# Calculate distance matrix for clustering
dist_matrix <- dist(heatmap_data, method = "euclidean")

# Perform hierarchical clustering
hc_rows <- hclust(dist_matrix, method = "complete")
hc_cols <- hclust(dist(t(heatmap_data)), method = "complete")

# Order the rows and columns based on the clustering
ordered_rows <- heatmap_data[order.dendrogram(as.dendrogram(hc_rows)), ]
ordered_cols <- heatmap_data[, order.dendrogram(as.dendrogram(hc_cols))]
# Plot the heatmap

pheatmap(ordered_rows, 
         clustering_distance_rows = dist_matrix,
         clustering_distance_cols = dist(t(ordered_cols)),
         col = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "none", # or "row" if you want to scale rows
         fontsize = 8,   # adjust font size
         cellwidth = 20, # adjust cell width
         cellheight = 4, # adjust cell height
         main = "Gene expression Heatmap plot",
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_col = 10) # adjust column annotation font size)

plot(hc_cols)
plot(hc_rows)


#correlation matrix
# Calculate the correlation matrix
correlation_matrix <- cor(normalized_expression_edgeR,method = 'kendall')
dim(correlation_matrix)

pheatmap(correlation_matrix, 
         col = colorRampPalette(c("blue", "white", "red"))(50),
         symm=TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none", # or "row" if you want to scale rows
         fontsize = 10,   # adjust font size
         cellwidth = 30, # adjust cell width
         cellheight = 30, # adjust cell height
         main = "Samples Correlation",
         show_rownames = TRUE,
         show_colnames = TRUE,
         xlab='Samples',
         ylab='Samples',
         margins = c(10, 10),
         cexRow = 0.8, cexCol = 0.8,
         fontsize_col = 10) # adjust column annotation font size)



#merge degs all vs center
library(clusterProfiler)
library(org.Hs.eg.db)

deg_Center_TA30<-rbind(up_list$'3'$CENTER_TA30mm,down_list$'3'$CENTER_TA30mm)
hgnc_oncorenrich<-as.data.frame(deg_Center_TA30$Gene_Symbol)
write.table(deg_Center_TA30,'C:/Users/User/Desktop/Dr_samanian/deg_Center_TA30.txt',row.names = F)

BP_CENTER_TA30 <- enrichGO(gene= deg_Center_TA30$entrezID,
                ont = "BP",
                OrgDb = org.Hs.eg.db,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable = TRUE)
BP_CENTER_TA30_df<-as.data.frame(BP_CENTER_TA30)
write.table(BP_CENTER_TA30_df,'C:/Users/User/Desktop/Dr_samanian/BP_CENTER_TA30.txt',row.names = F)

MF_CENTER_TA30 <- enrichGO(gene= deg_Center_TA30$entrezID,
                           ont = "MF",
                           OrgDb = org.Hs.eg.db,
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           readable = TRUE)
MF_CENTER_TA30_df<-as.data.frame(MF_CENTER_TA30)
write.table(MF_CENTER_TA30_df,'C:/Users/User/Desktop/Dr_samanian/MF_CENTER_TA30.txt',row.names = F)


#create genelist for GSEA (GeneID,Foldchange, sorted)

## feature 1: numeric vector
geneList = deg_Center_TA30$logFC

## feature 2: named vector
names(geneList) = as.character(deg_Center_TA30$entrezID)

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)

GSEA_BP_Center_TA30 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 20,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)


goplot(GSEA_BP_Center_TA30)
GSEA_BP_results<-GSEA_BP_Center_TA30@result
write.csv(GSEA_BP_results,'C:/Users/User/Desktop/Dr_samanian/GSEA_BP_Center_TA30.csv',row.names = F)

GSEA_BP_Neg<-GSEA_BP_Center_TA30
GSEA_BP_Pos<-GSEA_BP_Center_TA30

GSEA_BP_Neg@result<-GSEA_BP_Neg@result[GSEA_BP_Neg@result$enrichmentScore<0,]
GSEA_BP_Pos@result<-GSEA_BP_Pos@result[GSEA_BP_Pos@result$enrichmentScore>0,]

p1 <- dotplot(GSEA_BP_Pos, showCategory = 25, font.size=12)
p2 <- dotplot(GSEA_BP_Neg, showCategory = 25, font.size=12)


cowplot::plot_grid(p1,p2,labels=c('          Positively enriched BP (TA30mm-Tumor)','          Negatively enriched BP (TA30mm-Tumor)'))




#KEGG Pathway Enrichment

GSEA_KEGG_Center_TA30 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 20,
               pvalueCutoff = 0.05,
               verbose      = FALSE)


GSEA_KEGG_results<-GSEA_KEGG_Center_TA30@result
write.csv(GSEA_KEGG_results,'C:/Users/User/Desktop/Dr_samanian/GSEA_KEGG_Center_TA30.csv',row.names = F)

p1 <- dotplot(GSEA_KEGG_Center_TA30, showCategory = 25, font.size=12)

cowplot::plot_grid(p1,labels=c('        Enriched KEGG Pathways (TA30mm-Tumor)'))

library("pathview")
hsa05207 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa05207",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))


######CELL MARKER Enrichment

cell_marker_data <- vroom::vroom('http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx')


###disease gene
library(DOSE)
GSEA_DIONT_Center_TA30 <- gseDO(geneList,
           minGSSize     = 20,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)
GSEA_DIONT_results<-GSEA_DIONT_Center_TA30@result
write.csv(GSEA_DIONT_results,'C:/Users/User/Desktop/Dr_samanian/GSEA_DIONT_Center_TA30.csv',row.names = F)

GSEA_DIONT_Neg<-GSEA_DIONT_Center_TA30
GSEA_DIONT_Pos<-GSEA_DIONT_Center_TA30

GSEA_DIONT_Neg@result<-GSEA_DIONT_Neg@result[GSEA_DIONT_Neg@result$qvalue<0.05 & GSEA_DIONT_Neg@result$enrichmentScore<0,]
GSEA_DIONT_Pos@result<-GSEA_DIONT_Pos@result[GSEA_DIONT_Pos@result$qvalue<0.05 & GSEA_DIONT_Pos@result$enrichmentScore>0,]

p1 <- dotplot(GSEA_DIONT_Pos, showCategory = 15, font.size=12)
p2 <- dotplot(GSEA_DIONT_Neg, showCategory = 15, font.size=12)


cowplot::plot_grid(p1,p2,labels=c('Positively enriched Disease Ontology Terms (TA30mm vs Tumor)','Netgatively enriched Disease Ontology Terms (TA30mm vs Tumor)'))



#read data from oncoEnrichR results
oncoenrich_DEG_CENTER_A30_Association <- read.xlsx("C:/Users/User/Desktop/Dr_samanian/Galaxy3-[DEG_CENTER_TA30mm_-_xlsx].xlsx",
                sheet = 'CANCER_ASSOCIATION',
                rowNames = F,
                startRow = 1,
                detectDates = FALSE,
                skipEmptyRows = TRUE)
onco_DEG_Driver<-oncoenrich_DEG_CENTER_A30_Association[oncoenrich_DEG_CENTER_A30_Association$cancer_driver=='TRUE',]

onco_DEG__Hallmark <- read.xlsx("C:/Users/User/Desktop/Dr_samanian/Galaxy3-[DEG_CENTER_TA30mm_-_xlsx].xlsx",
                                                   sheet = 'CANCER_HALLMARK',
                                                   rowNames = F,
                                                   startRow = 1,
                                                   detectDates = FALSE,
                                                   skipEmptyRows = TRUE)
onco_DEG__tissue <- read.xlsx("C:/Users/User/Desktop/Dr_samanian/Galaxy3-[DEG_CENTER_TA30mm_-_xlsx].xlsx",
                                                   sheet = 'CELL_TISSUE',
                                                   rowNames = F,
                                                   startRow = 1,
                                                   detectDates = FALSE,
                                                   skipEmptyRows = TRUE)

onco_DEG__Hallmark_filt<-onco_DEG__Hallmark[onco_DEG__Hallmark$symbol %in% onco_DEG_Driver$symbol,]
onco_DEG__tissue_filt<-onco_DEG__tissue[onco_DEG__tissue$symbol %in% onco_DEG_Driver$symbol,]
onco_DEG_tissue_bulk<-onco_DEG__tissue_filt[onco_DEG__tissue_filt$category=='tissue',]
onco_DEG_tissue_sc<-onco_DEG__tissue_filt[onco_DEG__tissue_filt$category=='scRNA',]
onco_DEG_Driver$Tissue<-onco_DEG_tissue_bulk$tissue_or_celltype
onco_DEG_Driver$scRNA<-onco_DEG_tissue_sc$tissue_or_celltype

write.csv(onco_DEG_Driver,'C:/Users/User/Desktop/Dr_samanian/onco_DEG_Driver_Center_TA30.csv',row.names = F)
write.csv(onco_DEG__Hallmark_filt,'C:/Users/User/Desktop/Dr_samanian/onco_DEG_Hallmark.csv',row.names = F)

refrence_DEG<-as.data.frame(onco_DEG_Driver$symbol)


DEG_raw_Center_TA2mm<-rbind(up_list$'3'$CENTER_TA2mm,down_list$'3'$CENTER_TA2mm)
DEG_raw_Center_TA4mm<-rbind(up_list$'3'$CENTER_TA4mm,down_list$'3'$CENTER_TA4mm)
DEG_raw_Center_TA6mm<-rbind(up_list$'3'$CENTER_TA6mm,down_list$'3'$CENTER_TA6mm)
DEG_raw_Center_TA10mm<-rbind(up_list$'3'$CENTER_TA10mm,down_list$'3'$CENTER_TA10mm)
DEG_raw_Center_TA14mm<-rbind(up_list$'3'$CENTER_TA14mm,down_list$'3'$CENTER_TA14mm)
DEG_raw_Center_TA20mm<-rbind(up_list$'3'$CENTER_TA20mm,down_list$'3'$CENTER_TA20mm)
DEG_raw_Center_TA24mm<-rbind(up_list$'3'$CENTER_TA24mm,down_list$'3'$CENTER_TA24mm)
DEG_raw_Center_TA30mm<-rbind(up_list$'3'$CENTER_TA30mm,down_list$'3'$CENTER_TA30mm)
DEG_raw_Center_TAR2mm<-rbind(up_list$'3'$CENTER_TAR2mm,down_list$'3'$CENTER_TAR2mm)
DEG_raw_Center_TAR4mm<-rbind(up_list$'3'$CENTER_TAR4mm,down_list$'3'$CENTER_TAR4mm)
DEG_raw_Center_TAR6mm<-rbind(up_list$'3'$CENTER_TAR6mm,down_list$'3'$CENTER_TAR6mm)
DEG_raw_Center_TAR10mm<-rbind(up_list$'3'$CENTER_TAR10mm,down_list$'3'$CENTER_TAR10mm)

DEG_raw_CENTER_list<-list(TA2mm=rbind(up_list$'3'$CENTER_TA2mm,down_list$'3'$CENTER_TA2mm),
                          TA4mm=rbind(up_list$'3'$CENTER_TA4mm,down_list$'3'$CENTER_TA4mm),
                          TA6mm=rbind(up_list$'3'$CENTER_TA6mm,down_list$'3'$CENTER_TA6mm),
                          TA10mm=rbind(up_list$'3'$CENTER_TA10mm,down_list$'3'$CENTER_TA10mm),
                          TA14mm=rbind(up_list$'3'$CENTER_TA14mm,down_list$'3'$CENTER_TA14mm),
                          TA20mm=rbind(up_list$'3'$CENTER_TA20mm,down_list$'3'$CENTER_TA20mm),
                          TA24mm=rbind(up_list$'3'$CENTER_TA24mm,down_list$'3'$CENTER_TA24mm),
                          TA30mm=rbind(up_list$'3'$CENTER_TA30mm,down_list$'3'$CENTER_TA30mm),
                          TAR2mm=rbind(up_list$'3'$CENTER_TAR2mm,down_list$'3'$CENTER_TAR2mm),
                          TAR4mm=rbind(up_list$'3'$CENTER_TAR4mm,down_list$'3'$CENTER_TAR4mm),
                          TAR6mm=rbind(up_list$'3'$CENTER_TAR6mm,down_list$'3'$CENTER_TAR6mm),
                          TAR10mm=rbind(up_list$'3'$CENTER_TAR10mm,down_list$'3'$CENTER_TAR10mm))


temp<-refrence_DEG
colnames(temp)<-'HGNC'
for (sample in colnames(y)){
  if (sample!='CENTER')
  {
    temp<-merge(temp,DEG_raw_CENTER_list[[sample]][,c('Gene_Symbol','logFC')],by.x='HGNC',by.y='Gene_Symbol',all.x=T)
  }
  
 
}

colnames(temp)<-c('HGNC',colnames(y)[2:13])
temp[is.na(temp)] <- 0
temp_round <- round(temp[,-1], digits = 0)
temp_round$HGNC<-temp$HGNC

write.csv(temp,'C:/Users/User/Desktop/Dr_samanian/DriverGene_DEG_vs_Center.csv',row.names = F)
write.csv(temp_round,'C:/Users/User/Desktop/Dr_samanian/DriverGene_DEG_vs_Center_round.csv',row.names = F)
#bubble chart
library(ggplot2)
library(reshape2)
dev.off()  # Close the current graphics device
dev.new()  # Open a new graphics device
pcm = melt(temp_round, id = c("HGNC"))
pcm$HGNC <- factor(pcm$HGNC,levels=unique(pcm$HGNC))
colnames(pcm)<-c('HGNC','variable','LogFC')

bubble<-ggplot(pcm, aes(x = variable, y = HGNC, size = abs(LogFC), fill = LogFC)) +
  geom_point(shape = 21, color = "black", alpha = 0.7) +
  scale_size_continuous(range = c(3, 10)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Sample Name", y = "Gene symbol (HGNC)", title = "Bubble Chart of Driver genes LogFC in tumor adjacent samples versus Tumor, abs(logfc>3)") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 12, face = "bold",angle = 90, vjust = 0.5),  # Adjust x-axis label font size and weight
    axis.text.y = element_text(size = 12, face = "bold")   # Adjust y-axis label font size and weight
  )
print(bubble)
coord_fixed()
ggplotly(bubble)

save.image('C:/Users/User/Desktop/Dr_samanian/my_data_0903.RData')


