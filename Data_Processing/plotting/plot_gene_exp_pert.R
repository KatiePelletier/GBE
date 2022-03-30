###################################################################
## Script to make a plot to inspect expression of different
## genes by perturbation 
###################################################################

# Libraries
library(ggplot2) # plotting
library(stringr)

### Create sample information for plotting. Skip this and read in the output
# of it below ####

# # perturbation values
# pert_vals <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/final_alleleEffects_semiQscale_noBackground.csv")
# # get rid of the alleles that we don't have which are... row 10,11,12,3,6
# p_sub <- pert_vals[-c(3,6,10,11,12),]
# 
# # Pull out a gene, the background, and perturbation
# p_sub$perturbation <- p_sub$emmean - min(p_sub$emmean)
# 
# p_sub$allele <- gsub("\\[|\\]|-|t", "", p_sub$Maternal_Allele)
# 
# 
# ### Create a matrix that has everything I need, one for STAR counts and one for SALMON counts
# ### [gene expression] [background] [allele] [perturbation]
# 
# split_names <- str_split_fixed(string = names(my_gene)[-1],
#                                pattern = "_", 
#                                n = 4)
# 
# split_names
# 
# background <- c(rep("ORE", times = 14),
#                 rep("hybrid", times = 7),
#                 rep("SAM", times = 14),
#                 rep("hybrid", times = 7))
# 
# allele <- c(split_names[,3][1:14],
#             split_names[,4][15:21],
#             split_names[,3][22:35],
#             split_names[,4][36:42])
# 
# perturbation <- c()
# 
# for(i in 1:length(allele)){
#   perturbation <- c(perturbation, p_sub[which(p_sub$allele == allele[i]),7])
# }
# 
# # Define the sample information
# sample_info <- data.frame(background = background,
#                           allele = allele,
#                           perturbation = perturbation)
# 
# write.csv(sample_info, "plot_sample_info.csv", row.names = FALSE)

#### Read in data ####

# Read in the normalized count data from STAR and SALMON
# salmon counts
salmon_counts <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/salmon_norm_counts.csv")
dim(salmon_counts)

# STAR counts
star_counts <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/star_norm_counts.csv")
dim(star_counts)

# Sample information, for creating matrices for plotting:
sample_info <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/plot_sample_info.csv")

#### Plotting Function ####

# Note that this function will use the star counts by default,
# but you can change it to the salmon_counts if you wish
plot_pert_exp <- function(gene,method=star_counts){
  gene_exp <- as.numeric(method[which(method$X == gene),][-1])
  plot_this <- cbind.data.frame(gene_exp, sample_info)
  ggplot(plot_this, aes(x = perturbation, y = gene_exp, colour = allele, shape = background)) +
    geom_point(size = 4, alpha = 0.65) + 
    labs(y = "normalized gene expression", x = "perturbation") + 
    ggtitle(gene)
}

# Dll
plot_pert_exp("FBgn0000157")

# Vg
plot_pert_exp("FBgn0003975")
plot_pert_exp("FBgn0003975", method = salmon_counts)

# Sd
plot_pert_exp("FBgn0003345")
plot_pert_exp("FBgn0003345", method = salmon_counts)

# omb/bifid
plot_pert_exp("FBgn0000179")

#optix
plot_pert_exp("FBgn0025360")
