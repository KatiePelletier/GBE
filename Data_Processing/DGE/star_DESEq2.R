
### Script to import STAR counts, run same analyses as Arteen with DESeq2

### I processed the ReadsPerGene.out file to grab the 1st and 4th columns already)

#############################################

library(DESeq2)
library(stringr)

# Create the data table to work with (not the counts, for now):
file_names <- list.files(path = "/Users/amandan/Desktop/Dworkin/background_effects/data/STAR_counts_4thcol")

file_name <- file_names

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

file_names <- str_split_fixed(file_names, pattern = "_", n = 6)

lane <- file_names[,5]
lane <- as.factor(lane)

sample_name <- paste0(1:42, "_", file_names[,1], "_", file_names[,2], "_", file_names[,3])
sample_name <- gsub("_r2|_r1", replacement = "", sample_name)

#levels(mutation) <- c ("w","sd1","vg2a33","sdETX4","sdE3","vg213","vg1")

mutation_type <- c("sd","sd","sd","sd","sd","sd","vg",
                   "vg","vg","vg","vg","vg","w","w",
                   "sd","sd","sd","vg","vg","vg",
                   "w","sd","sd","sd","sd","sd","sd",
                   "vg","vg","vg","vg","vg","vg",
                   "w","w","sd","sd","sd","vg","vg",
                   "vg","w")

mutation_type <- as.factor(mutation_type)

#Idk just by hand, its not that serious?
perturbation <- c (0.158665903, 0.158665903,
                  1.03033975, 1.03033975,
                  0.63442789,0.63442789,
                  1.716418055,1.716418055,
                  1.403581039, 1.403581039,
                  0.615181681, 0.615181681,
                  0,0,
                  0.158665903,1.03033975,
                  0.63442789,1.716418055,
                  1.403581039,0.615181681,
                  0,
                  0.158665903, 0.158665903,
                  1.03033975, 1.03033975,
                  0.63442789,0.63442789,
                  1.716418055,1.716418055,
                  1.403581039, 1.403581039,
                  0.615181681, 0.615181681,
                  0,0,
                  0.158665903,1.03033975,
                  0.63442789,1.716418055,
                  1.403581039,0.615181681,
                  0
)

sample_table <- data.frame(sample_name, file_name, strain, mutation, lane, mutation_type, perturbation)

load.model <- formula(~lane + mutation_type + perturbation + strain + perturbation:strain)
#Read in the counts:

star.counts <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                      directory = "/Users/amandan/Desktop/Dworkin/background_effects/data/STAR_counts_4thcol",
                                      design = load.model)


head(counts(star.counts)) # Note that the first 4 rows are not actually counts but metadata

#write.csv(counts(star.counts)[5:nrow(counts(star.counts)),], "star_counts.csv")

#filter out lowly mapped reads:
keep <- rowSums(counts(star.counts)) > 5
star.counts <- star.counts[keep,]
nrow(star.counts)

#remove the first 4 lines:
star.counts <- star.counts[-c(1:4),]
nrow(star.counts)

#creating DESeq object
all.data <- DESeq(star.counts)

# plot dispersion estimates
plotDispEsts(all.data,CV=T)

# exploratory PCAs
# Play around with the 2000 and 500 most variable genes... a weird guy jumps out at 500..
for_pca <- vst(all.data,blind=FALSE) #variance stabilized counts
# what would regularized log transformed look like?

plotPCA(for_pca,intgroup=c("strain"),ntop=500)
plotPCA(for_pca,intgroup=c("mutation"),ntop=500)
plotPCA(for_pca,intgroup=c("lane"),ntop=500)
plotPCA(for_pca,intgroup=c("mutation_type"),ntop=500)

plotPCA(for_pca,intgroup=c("strain"),ntop=2000)
plotPCA(for_pca,intgroup=c("mutation"),ntop=2000)
plotPCA(for_pca,intgroup=c("lane"),ntop=2000)

# Look into that sample that is the strange out liar 

z <- plotPCA(for_pca,intgroup=c("strain"),ntop=500)
z + geom_text(aes(label = sample_name)) + xlim(-30,20)
# So the outlier is sample 5 (ORE_sdETX4)

# What's causing the weird gene? Maybe somebody grabbed a bit of another tissue...
# Might have to do analyses with this included and excluded 

weirdo <- counts(all.data)[,5]
weirdo <- weirdo[order(weirdo, decreasing = TRUE)]
weirdo[1:20]
bro <- counts(all.data)[,6]
bro <- bro[order(bro, decreasing = TRUE)]
bro["FBgn0003373"]
weirdo["FBgn0003373"] #salivary gland secretion looks weird...
counts(all.data)[rownames(counts(all.data))=="FBgn0003373",] #yep...

counts(all.data)[rownames(counts(all.data))=="FBgn0003975",] #vestigial
counts(all.data)[rownames(counts(all.data))=="FBgn0003345",] #scalloped 

# Model stuff:
# Showing contrasts

# resultsNames(all.data)
# strainResult <- results(all.data, contrast = c("strain", "ORE","SAM"))
# strainResult
#write.csv(strainResult, "deseq_star_model1_results.csv")

# These gave us crazy results from the salmon counts... let's see what happens here

resultsNames(all.data)
allResults <- results(all.data)
summary(allResults)

plotMA(allResults)

plot(allResults$log2FoldChange)

###########################################3

####### Notes from meeting on March 2, 2022 #########

# Driad chandler paper stuff for wing size, and cell counts from juveniles? 
# Instead of working with gene counts

# Design matrix: 
# what we are currently doing is arbitrary dominance (unordered factor for the three ORE/SAM/het)
#1) additive (how many alleles or ore vs. sam): 0,1,2
#2) dominance deviation: 0,1,0
#3) instead of mutation as an unordered factor: could have a column for degree of perturbation (wing area or cell count)
#4) is the mutant allele scalloped or vestigial? (gene level effect)

# Vector correlations between different backgrounds/mutations...

# Are we getting the same genes that are differential expressed between the Salmon and STAR counts?
# similar expression values? especially for vg1? Ian thinks these will be very different

# Ian doesn't expect the results between DESEq2,limmavoom, edgeR to be similar?
# sooooo... check that out!

# Could do a smaller model, like comparing sdE3 to w, wrap head around that, then expand


# Here is the wing size data stuff

size <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/Combined_Final_Data_05_14.csv")
dim(size)

# why do we want to include the alleles that we don't have, will this be overfitting our model?


###### Looking at the differences in expression coutns with STAR vs. Salmon

# Load in the data

# Salmon counts from Arteen:
salmon_counts <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/salmon_counts/salmon_counts.csv")
dim(salmon_counts)

# Star counts from me:
star_counts <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/star_counts.csv")
dim(star_counts)

# First off, STAR gave us so many more genes than Salmon did.... let's think about this more?

length(which(star_counts$X %in% salmon_counts$X == TRUE))
length(which(salmon_counts$X %in% star_counts$X == TRUE))

# There are 13669 matching genes between the two... which is aaaalmost all of the star count
star_genes <-star_counts$X[which(star_counts$X %in% salmon_counts$X == FALSE)]
# sooo... these are probably genes that are splice variants
#write.csv(star_genes[1:1000], "star_unique_genes_1000.csv", row.names = FALSE)

# Checked with FlyBase and looks like a lot of pseudo genes, transfer RNAs... etc.

# So, how similar are the values that we have calculated for the 13669 genes that are shared between them?

matched_star <- star_counts[star_counts$X %in% salmon_counts$X,]
matched_salmon <- salmon_counts[salmon_counts$X %in% star_counts$X,]

# order the samples so the genes are in the same order
matched_star <- matched_star[order(matched_star$X),]
matched_salmon <- matched_salmon[order(matched_salmon$X),]

# do the column names match?
colnames(matched_star)
colnames(matched_salmon)
# names aren't the same but the samples are in the same order... so that's fine

# Can I just plot? would that work? Perhaps one sample at a time? Generate some correlations?
plot(x = matched_star[,2], y = matched_salmon[,2])
plot(x = matched_star[,3], y = matched_salmon[,3])
plot(x = matched_star[,4], y = matched_salmon[,4])
plot(x = matched_star[,5], y = matched_salmon[,5])
plot(x = matched_star[,6], y = matched_salmon[,6])
plot(x = matched_star[,7], y = matched_salmon[,7])
plot(x = matched_star[,8], y = matched_salmon[,8])
plot(x = matched_star[,9], y = matched_salmon[,9])
plot(x = matched_star[,10], y = matched_salmon[,10])
plot(x = matched_star[,11], y = matched_salmon[,11])
plot(x = matched_star[,12], y = matched_salmon[,12])
plot(x = matched_star[,13], y = matched_salmon[,13])
plot(x = matched_star[,14], y = matched_salmon[,14])
plot(x = matched_star[,15], y = matched_salmon[,15])

# So i think to interpret this plot, if it was a straight diagonal line there would be perfect correlation
# between the counts, but as you can see for most, the salmon counts are slightly higher but I think that this makes sense
# ???? Don't think that this is telling me what I want so I shall try another way... subtract the star counts from the salmon counts


sub_counts <- matched_salmon[,2:ncol(matched_salmon)] - matched_star[,2:ncol(matched_star)] 

plot(x = 1:nrow(sub_counts), y = sub_counts[,1])
max(sub_counts[,1])
matched_salmon$X[sub_counts[,1] == 229291]
matched_star$X[sub_counts[,1] == 229291]

matched_salmon[matched_salmon$X == "FBgn0284245",]
matched_star[matched_star$X == "FBgn0284245",]

matched_star[1:5,1:3]
matched_salmon[1:5,1:3]

df2 <-data.frame(c(matched_salmon[,2], matched_star[,2]), rep(1:nrow(matched_salmon), times = 2), 
                   c(rep("blue", times = nrow(matched_salmon)), rep("red", times = nrow(matched_star))))

plot(x = df2[,2], y = df2[,1], col = df2[,3])

df3 <- data.frame(c(matched_salmon[,3], matched_star[,3]), rep(1:nrow(matched_salmon), times = 2), 
                  c(rep("blue", times = nrow(matched_salmon)), rep("red", times = nrow(matched_star))))

plot(x = df3[,2], y = df3[,1], col = df3[,3])

df4 <- data.frame(c(matched_salmon[,4], matched_star[,4]), rep(1:nrow(matched_salmon), times = 2), 
                  c(rep("blue", times = nrow(matched_salmon)), rep("red", times = nrow(matched_star))))

plot(x = df4[,2], y = df4[,1], col = df4[,3])

df5 <- data.frame(c(matched_salmon[,5], matched_star[,5]), rep(1:nrow(matched_salmon), times = 2), 
                  c(rep("blue", times = nrow(matched_salmon)), rep("red", times = nrow(matched_star))))

plot(x = df5[,2], y = df5[,1], col = df5[,3])

df6 <- data.frame(c(matched_salmon[,6], matched_star[,6]), rep(1:nrow(matched_salmon), times = 2), 
                  c(rep("blue", times = nrow(matched_salmon)), rep("red", times = nrow(matched_star))))

plot(x = df6[,2], y = df6[,1], col = df6[,3])

# Salmon counts seem to be higher usually than star counts, I think that this makes sense

######################

# Compare the salmon model results to the star model results (same models, are the top genes the same?)

salmon_res <- read.csv("deseq_salmon_model1_results.csv")
star_res <- read.csv("deseq_star_model1_results.csv")

# Well, first off the two models had different genes being input so maybe that will affect things

salmon_res <- salmon_res[order(salmon_res$padj, decreasing = FALSE),]
star_res <- star_res[order(star_res$padj, decreasing = FALSE),]

head(salmon_res)
head(star_res)

# How much of the top 100 match?

length(which(salmon_res$X[1:100] %in% star_res$X[1:100] == TRUE)) #67 match! not bad

length(which(salmon_res$X[1:250] %in% star_res$X[1:250] == TRUE)) #166 matching 

# Are the not matching ones even in the salmon results?

salmon_no <- salmon_res$X[which(salmon_res$X[1:250] %in% star_res$X[1:250] == FALSE)]
star_no <- star_res$X[which(star_res$X[1:250] %in% salmon_res$X[1:250] == FALSE)]

length(which(star_no %in% salmon_res$X == FALSE)) #ehhh about 41 out of 84 not matching aren't even in the dataset

### Should do a GO term thing with this...
