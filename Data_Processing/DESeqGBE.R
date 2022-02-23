library(DESeq2)
library(tximport)
library(readr)
library("RColorBrewer")
library(gplots)
library(stringr)
library("pheatmap")


#Salmon counts are found at /2/scratch/arteen/GBE_rnaseq/salmon_counts

#Import
#To the directory of counts
setwd("/Users/arteen/Desktop/School/Projects/GBE")

#I renamed the folder containing all the counts (salmon_counts) to quants
quant_files <- file.path("quants", list.files("quants"), "quant.sf")
file.exists(quant_files)

samples <- c("ORE_sd1_r1_GTGGCC_L002_quant",
  "ORE_sd1_r2_TTAGGC_L003_quant",
  "ORE_sdE3_r1_GTGGCC_L004_quant",
  "ORE_sdE3_r2_TGACCA_L005_quant",
  "ORE_sdETX4_r1_ATCACG_L003_quant",
  "ORE_sdETX4_r2_GTTTCG_L004_quant",
  "ORE_vg1_r1_TGACCA_L007_quant",
  "ORE_vg1_r2_TTAGGC_L001_quant",
  "ORE_vg213_r1_AGTCAA_L006_quant",
  "ORE_vg213_r2_CGATGT_L007_quant",
  "ORE_vg2a33_r1_CGATGT_L005_quant",
  "ORE_vg2a33_r2_AGTTCC_L006_quant",
  "ORE_w_r1_ATCACG_L001_quant",
  "ORE_w_r2_GTTTCG_L002_quant",
  "OREf_SAMm_sd1_CTTGTA_L007_quant",
  "OREf_SAMm_sdE3_ATTCCT_L002_quant",
  "OREf_SAMm_sdETX4_GGCTAC_L001_quant",
  "OREf_SAMm_vg1_CTTGTA_L005_quant",
  "OREf_SAMm_vg213_ATTCCT_L004_quant",
  "OREf_SAMm_vg2a33_GGCTAC_L003_quant",
  "OREf_SAMm_w_GTCCGC_L006_quant",
  "SAM_sd1_r1_CGTACG_L004_quant",
  "SAM_sd1_r2_GCCAAT_L005_quant",
  "SAM_sdE3_r1_ATGTCA_L006_quant",
  "SAM_sdE3_r2_GCCAAT_L007_quant",
  "SAM_sdETX4_r1_ACAGTG_L005_quant",
  "SAM_sdETX4_r2_CCGTCC_L006_quant",
  "SAM_vg1_r1_CGTACG_L002_quant",
  "SAM_vg1_r2_GATCAG_L003_quant",
  "SAM_vg213_r1_ACTTGA_L001_quant",
  "SAM_vg213_r2_GAGTGG_L002_quant",
  "SAM_vg2a33_r1_ACAGTG_L007_quant",
  "SAM_vg2a33_r2_GATCAG_L001_quant",
  "SAM_w_r1_ACTTGA_L003_quant",
  "SAM_w_r2_GAGTGG_L004_quant",
  "SAMf_OREm_sd1_GTGAAA_L006_quant",
  "SAMf_OREm_sdE3_TAGCTT_L001_quant",
  "SAMf_OREm_sdETX4_CAGATC_L007_quant",
  "SAMf_OREm_vg1_ACTGAT_L004_quant",
  "SAMf_OREm_vg213_TAGCTT_L003_quant",
  "SAMf_OREm_vg2a33_ACTGAT_L002_quant",
  "SAMf_OREm_w_CAGATC_L005_quant")


names(quant_files) <- samples

library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(x = txdb, keys = k, "GENEID", "TXNAME")

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)


strain <- c(rep("ORE",14),rep("OREf_SAMm",7),rep("SAM",14),rep("SAMf_OREm",7))
strain <- as.factor(strain)

mutation <- c("sd1","sd1","sdE3","sdE3","sdETX4","sdETX4","vg1",
  "vg1","vg213","vg213","vg2a33","vg2a33","w","w",
  "sd1","sdE3","sdETX4","vg1","vg213","vg2a33",
  "w","sd1","sd1","sdE3","sdE3","sdETX4","sdETX4",
  "vg1","vg1","vg213","vg213","vg2a33","vg2a33",
  "w","w","sd1","sdE3","sdETX4","vg1","vg213",
  "vg2a33","w")
mutation <- as.factor(mutation)
levels(mutation) <- c ("w","sd1","vg2a33","sdETX4","sdE3","vg213","vg1")

lane <- str_sub(samples, - 10, - 1) 
lane <- str_sub(lane,1,4)
lane <- as.factor(lane)

rna.design <- data.frame(sample=samples, 
                         file=quant_files,
                         strain=strain,
                         mutation=mutation,
                         lane=lane)

#not loading a model, just using ~lane to look at dispersion and PCA
#pretty sure you can just add the real/another model here if youd like
load.model <- formula(~lane)

all.data <- DESeqDataSetFromTximport(txi,rna.design,design=load.model)


#Filtering out rows with 1 or 0 counts total
keep <- rowSums(counts(all.data)) > 1
all.data <- all.data[keep,]
nrow(all.data)


#creating DESeq object
all.data <- DESeq(all.data)

plotDispEsts(all.data,CV=T)

for_pca <- vst(all.data,blind=TRUE)

plotPCA(for_pca,intgroup=c("lane"),ntop=2000)

plotPCA(for_pca,intgroup=c("strain"),ntop=2000)

#heatmap
rlogMat <- assay(for_pca)
distsRL <- dist(t(rlogMat))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(all.data),paste(strain,mutation,lane,sep=" : "))
hc <-  hclust(distsRL)
hmcol <-  colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=rev(hmcol),margin=c(10,10))



#Fitting a model on
#I am replacing the old all.data (used previously) with this one that uses the model below
load.model <- formula(~lane + strain + mutation + strain:mutation)
all.data <- DESeqDataSetFromTximport(txi,rna.design,design=load.model)
all.data <- DESeq(all.data)

#Filtering out rows with 1 or 0 counts total
keep <- rowSums(counts(all.data)) > 1
all.data <- all.data[keep,]
nrow(all.data)

allResults <- results(all.data)

#Showing contrasts
resultsNames(all.data)

strainResult <- results(tmp, contrast = c("strain", "ORE","SAM"))
strainResult


#results - y stuff

#Subset only results w/ p-adjusted < 0.1
resSig <- subset(strainResult, padj < 0.1)

#Strongest down regulation
head(resSig[ order(resSig$log2FoldChange), ])
#Strongest up regulation
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])



#Another PCA + Heatmap because why not lol
tmpPCA <- vst(all.data,blind=TRUE)
plotPCA(tmpPCA,intgroup=c("lane"),ntop=2000)


rlogMat <- assay(tmpPCA)
distsRL <- dist(t(rlogMat))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(tmp),paste(strain,mutation,lane,sep=" : "))
hc <-  hclust(distsRL)
hmcol <-  colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=rev(hmcol),margin=c(10,10))





