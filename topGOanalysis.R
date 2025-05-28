library(topGO)
library(data.table)
library(splitstackshape)
##need to filter by to just keep the genes in the main chromossomes in case it adds to some aritficial imbalance when i diddnt allow some genes to be conseidred for differences in methylation
setwd("./methyl_analysis/goanalysis/")


length_data <- read.table("./STRINGS/src/gene_length.csv")
names(length_data) <- c('name','length')
all_genes <- length_data$name

categories <- readMappings(file = "./categories_topGO.csv")

enriched_all <- read.table("./finalcov10xwcover10icover3nsamples12_annotated_all_formated.csv",header=TRUE)



enriched_genes <- unique(enriched_all$name)

geneList <- factor(as.integer(all_genes %in% enriched_genes))
names(geneList) <- all_genes

GOdata <- new("topGOdata",ontology="BP",allGenes=geneList,annot=annFUN.gene2GO,gene2GO=categories)

test.stat <- new("classicCount",testStatistic= GOFisherTest,name="Fisher test")
resultFisher <- getSigGroups(GOdata,test.stat)


allRes <- GenTable(GOdata, Fisher = resultFisher,topNodes=100)
allRes$adjusted <- p.adjust(allRes$Fisher,method="BH")
allRes
