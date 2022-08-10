##############################################################################
#### Given a gene, plot logCPM values on Y axis, perturbation on X axis, CI and fit estimates from the PredictorEffect function, and expression values for each background/allele combination
## Amanda Neves August 10, 2022
##############################################################################

# libraries

library(lme4)
library(glmmTMB)
library(emmeans)
library(parallel)
library(stringr)

# Load in data, run model
setwd("/Users/amandan/Desktop/Dworkin/background_effects/data")
logCPM_norm <- read.csv("gaussian_salmon_logcpm_aug10.csv", row.names = 1)

gene_names <- read.delim("FlyGeneDictionary.txt",
                         header = TRUE,
                         sep = "\t",
                         dec = ".") # for adding gene name to plot 

background <- c(rep("ORE", 14), rep("SAM", 14))
background <- as.factor(background)
background <- relevel(background, "ORE")

perturbation <- c(7.434782609, 7.434782609,
                  3.541666667, 3.541666667,
                  4.564102564, 4.564102564,
                  1.733333333, 1.733333333,
                  3.692307692, 3.692307692,
                  5.846153846, 5.846153846,
                  10, 10,
                  9.708333333, 9.708333333,
                  6.739130435, 6.739130435,
                  8.275862069, 8.275862069,
                  1.217391304, 1.217391304,
                  3.695652174, 3.695652174,
                  9.714285714, 9.714285714,
                  10, 10)

perturbation <- abs(perturbation - max(perturbation))


lane <- c(2,3,4,5,3,4,7,1,6,7,5,6,1,2,4,5,6,7,5,6,2,3,1,2,7,1,3,4)
lane <- as.factor(lane)

nt <- min(parallel::detectCores(),5)

# Plotting function

plot_PE_gene <- function(gene, model = c("linear", "quadratic")){
  geneID <- gene
  single_cpm <- logCPM_norm[geneID,]
  rownames(single_cpm) <- NULL
  tmb.data <- data.frame(single_cpm = as.numeric(single_cpm[1,]), background, perturbation, lane)
  
  if (model == "linear"){
    tmb.mod.outlier <- glmmTMB(single_cpm ~ background + 
                                 perturbation + 
                                 background:perturbation + 
                                 (1 | lane),
                               data = tmb.data,
                               control = glmmTMBControl(parallel = nt))
  } else {
    tmb.mod.outlier <- glmmTMB(single_cpm ~ background *
                                           poly(perturbation,2) +
                                           (1 | lane),
                                         data = tmb.data,
                                         control = glmmTMBControl(parallel = nt))
  }
  
  plot_df <- as.data.frame(predictorEffect("perturbation", tmb.mod.outlier, xlevels = 13))
  
  split_ids <- str_split(names(single_cpm), "_")
  
  base_background <- sapply(split_ids,"[[",1)
  
  base_allele <- c("sd[1]", "sd[1]", "sd[E3]", "sd[E3]", "sd[ETX4]", "sd[ETX4]", "vg[1]", "vg[1]", "vg[21-3]", "vg[21-3]", "vg[2a33]", "vg[2a33]", "wild type", "wild type", "sd[1]", "sd[1]", "sd[E3]", "sd[E3]", "sd[ETX4]", "sd[ETX4]", "vg[1]", "vg[1]", "vg[21-3]", "vg[21-3]", "vg[2a33]", "vg[2a33]", "wild type", "wild type")
  
  
  # and i need to add the perturbation of each... 
  
  base_perturbation <- perturbation
  
  point_df <- data.frame(background = base_background,
                         allele = base_allele,
                         expression = as.numeric(single_cpm[1,]),
                         perturbation = base_perturbation)
  
  plot_title <- ifelse(gene %in% gene_names$X.submitted_item,
                       gene_names$current_symbol[match(gene, gene_names$X.submitted_item)],
                       gene)
  
  ggplot() + geom_line(data = plot_df, aes(x = perturbation, y = fit, colour = background)) + geom_ribbon(data=plot_df,aes(ymin=lower,ymax=upper,fill = background, x = perturbation),alpha=0.3) + geom_point(data = point_df, size = 4, alpha = 0.65, aes(x = perturbation, y = expression, shape = allele, colour = background)) + scale_shape_manual(values=1:7) + labs(y = "Gene Expression (CPM)", x = "Perturbation") + theme_classic() + ggtitle(substitute(italic(x), list(x=plot_title)))
}

# test it out on a gene!

plot_PE_gene(gene = "FBgn0003882", model = "linear" )
plot_PE_gene(gene = "FBgn0003882", model = "quadratic" )
