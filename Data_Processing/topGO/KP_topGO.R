### Script to do GO analysis on output of differential expression analyses
# (here is limma specifically)

#Amanda wrote most of this then Katie edited. 

###########################

# libraries

library(topGO)

# Load in the data

pert.limma <- read.csv("limma_salmon_mar22_2022.csv")

str(pert.limma)
# This is the output of DGE analysis done in limma by Arteen, using salmon counts
# the model was ~ strain + pertubation + pertubation:strain
# I already have ordered the results by adjusted p-value

# This file is from Katie (previously from Sam George), it basicaly assigns
# GO terms to genes (is adjusted from the gene_association.fb from FlyBase to be 
# readable by topGO)
gene_GO <- readMappings("fly_to_GO.delim")


####Need to get the sig up and down genes (using adj p <0.05 and sign of lfc)
gene_signifcances_up <- pert.limma[pert.limma$logFC > 0,]$adj.P.Val 
names(gene_signifcances_up) <- pert.limma[pert.limma$logFC > 0,]$X
gene_signifcances_up <- gene_signifcances_up[complete.cases(gene_signifcances_up)]

gene_signifcances_down <- pert.limma[pert.limma$logFC < 0,]$adj.P.Val 
names(gene_signifcances_down) <- pert.limma[pert.limma$logFC < 0,]$X
gene_signifcances_down <- gene_signifcances_down[complete.cases(gene_signifcances_down)]

all_genes <- pert.limma$adj.P.Val
names(all_genes) <- pert.limma$X
all_genes <- all_genes[complete.cases(all_genes)]


#for filtering based on p-val. this can be changed
gene_filter <- function(allScore){
  return(allScore < 0.05)
}

#Now I will make these topGO objects
Positive_Genes <- new("topGOdata",
                      ontology = "BP", 
                      allGenes = gene_signifcances_up,  
                      geneSel = gene_filter, 
                      nodeSize = 10, 
                      annotationFun = annFUN.gene2GO, 
                      gene2GO = gene_GO) 

negitive_Genes <- new("topGOdata",
                   ontology = "BP", 
                   allGenes = gene_signifcances_down,  
                   geneSel = gene_filter, 
                   nodeSize = 10, 
                   annotationFun = annFUN.gene2GO, 
                   gene2GO = gene_GO) 

#KP to KP can I do this with all the genes later (up and down)?

all_genes <- new("topGOdata",
                      ontology = "BP", 
                      allGenes = all_genes,  
                      geneSel = gene_filter, 
                      nodeSize = 10, 
                      annotationFun = annFUN.gene2GO, 
                      gene2GO = gene_GO) 



####Now a go analysis 

#####just testing for one term (wing development)
wing_genes <- genesInTerm(Positive_Genes, "GO:0035220")[[1]]
Pos_genes <- genes(Positive_Genes)
sig_genes <- sigGenes(Positive_Genes)

wing_test <- new("classicCount", testStatistic = GOFisherTest, name = "fisher",
                 allMembers = Pos_genes, groupMembers = wing_genes,
                 sigMembers = sig_genes)

contTable(wing_test)
runTest(wing_test)

#################################################

###This is looking at all the groups (defults)
Positive_Fisher <- runTest(Positive_Genes, algorithm = "classic", statistic = "fisher")
Positive_Fisher

pos_GO <- Positive_Fisher@score

Pos.Table <- GenTable(Positive_Genes, Fisher = Positive_Fisher,    
                      topNodes = 10)
Pos.Table


###and negitive
neg_Fisher <- runTest(negitive_Genes, algorithm = "classic", statistic = "fisher")

neg_GO <- neg_Fisher@score

neg.Table <- GenTable(negitive_Genes, Fisher = neg_Fisher,    
                      topNodes = 10)

neg.Table


#and all 
all_Fisher <- runTest(all_genes, algorithm = "classic", statistic = "fisher")

all_GO <- neg_Fisher@score

all.Table <- GenTable(all_genes, Fisher = all_Fisher,    
                      topNodes = 10)

all.Table


