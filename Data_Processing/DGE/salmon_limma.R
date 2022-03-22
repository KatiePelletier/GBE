library(DESeq2)
library(tximport)
library(readr)
library("RColorBrewer")
library(gplots)
library(stringr)
library(edgeR)
library(limma)


#Salmon counts are found at /2/scratch/arteen/GBE_rnaseq/salmon_counts

#Import
#To the directory of counts
setwd("/Users/amandan/Desktop/Dworkin/background_effects/data")

#I renamed the folder containing all the counts (salmon_counts) to quants
quant_files <- file.path("quants", list.files("quants"), "quant.sf")
quant_files <- quant_files[ -c(5) ]
file.exists(quant_files)

samples <- c("ORE_sd1_r1_GTGGCC_L002_quant",
             "ORE_sd1_r2_TTAGGC_L003_quant",
             "ORE_sdE3_r1_GTGGCC_L004_quant",
             "ORE_sdE3_r2_TGACCA_L005_quant",
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


strain <- c(rep("ORE", 13), rep("F1", 7),
            rep("SAM", 14), rep("F1", 7))

strain <- as.factor(strain)
strain <- relevel(strain, "ORE")

mutation_type <- c("sd","sd","sd","sd","sd","vg",
                   "vg","vg","vg","vg","vg","w","w",
                   "sd","sd","sd","vg","vg","vg",
                   "w","sd","sd","sd","sd","sd","sd",
                   "vg","vg","vg","vg","vg","vg",
                   "w","w","sd","sd","sd","vg","vg",
                   "vg","w")

mutation_type <- as.factor(mutation_type)
mutation_type <- relevel(mutation_type, "w")

pertubation <- c(8.571557971,8.571557971,
                 5.140398551,5.140398551,
                 6.419982317,
                 1.475362319,1.475362319,
                 3.693979933,3.693979933,
                 7.78021978,7.78021978,
                 10,10,
                 8.571557971,5.140398551,
                 6.419982317,1.475362319,
                 3.693979933,7.78021978,
                 10,
                 8.571557971,8.571557971,
                 5.140398551,5.140398551,
                 6.419982317,6.419982317,
                 1.475362319,1.475362319,
                 3.693979933,3.693979933,
                 7.78021978,7.78021978,
                 10,10,
                 8.571557971,5.140398551,
                 6.419982317,1.475362319,
                 3.693979933,7.78021978,
                 10)

pertubation <- pertubation - min(pertubation)

lane <- str_sub(samples, - 10, - 1) 
lane <- str_sub(lane, 1, 4)
lane <- as.factor(lane)



rna.design <- data.frame(sample=samples,
                         file=quant_files,
                         strain=strain,
                         mutation_type=mutation_type,
                         pertubation=pertubation,
                         lane=lane)

#Actual limma voom part
dgeList <- DGEList(txi$counts) # Create a DGEList object 

design <- model.matrix(~ strain + pertubation + pertubation:strain)


keep2 <- filterByExpr(dgeList, design, min.count = 1) # Determine which genes have sufficiently large counts to be retained ina  statistical analysis 
# Keep genes that have at least min.count reads in a worthwile number of samples

# So, Arteen originally had filterByExpr(dgeList, design) but I added min.count = 1 because his way
# was removing 4943 samples which left us with 8758 which idk seems like not a lot?
#

length(which(keep2 == FALSE))
# Removing 2846 genes

dgeList <- dgeList[keep2, ] # Remove rows that have consistently low or very low counts

dim(dgeList)
# now left with 10855 genes  

# normalize and run voom transformation
dgeList <- calcNormFactors(dgeList)

v <- voom(dgeList, design,plot=TRUE) # Transform count data to log2-counts per million (logCPM) using the 
# calculated normalization factors 


#Fitting
fit <- lmFit(v, design) # Fit the linear model
fit <- eBayes(fit) # Compute moderated t-statistics, moderated F-statistic, and log-odds of differential 
# expression by empirical Bayes moderation of the SE 

summary(decideTests(fit)) # identify which genes are significantly differentially expressed for each 
# contrast from a fit object containing p values and tests statistics 

#Pulling sex and interaction term again
pert.limma <- topTable(fit, coef = 4, sort.by = "p", n = Inf)
# Log fold-changes do look like they make more sense here now

#write.csv(pert.limma, "limma_salmon_mar22_2022.csv")
