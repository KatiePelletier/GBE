### Script to do GO analysis on output of differential expression analyses
# (here is limma specifically)

###########################

# libraries

library(topGO)

# Load in the data

pert.limma <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/limma_salmon_mar22_2022.csv")
# This is the output of DGE analysis done in limma by Arteen, using salmon counts
# the model was ~ strain + pertubation + pertubation:strain
# I already have ordered the results by adjusted p-value

# This file is from Katie (previously from Sam George), it basicaly assigns
# GO terms to genes (is adjusted from the gene_association.fb from FlyBase to be 
# readable by topGO)
gene_GO <- readMappings("/Users/amandan/Desktop/Dworkin/background_effects/data/fly_to_GO.delim")

# Create the file which contains our "interesting" genes (p-adj < 0.05)
pert.limma$adj.P.Val < 0.05 # top 251 are < 0.05
my_genes <- pert.limma$X[1:251]

# If it's an interesting gene give a 1, if not give a 0:
interesting_genes <- factor(as.numeric(pert.limma$X %in% my_genes))

# give gene names to interesting_genes
names(interesting_genes) <- pert.limma$X 

topGOobj <- new("topGOdata", 
            ontology = "BP", 
            allGenes = interesting_genes, 
            annotationFun = annFUN.gene2GO, 
            gene2GO = gene_GO)

resultFisher <- runTest(topGOobj, algorithm = "classic", statistic = "fisher")

result_table <- GenTable(topGOobj, classic = resultFisher, ranksOf = "classic", topNodes = 20)

