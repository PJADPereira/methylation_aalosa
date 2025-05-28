library(goseq)


setwd("./methyl_analysis/goanalysis/")


length_data <- read.table("./STRINGS/src/gene_length.csv")
names(length_data) <- c('name','length')

categories <- read.table("./STRINGS/src/goterms.csv")
names(categories) <- c("name","category")

enriched_all <- read.table("./differentiallymethylatedregions.csv",header=TRUE,sep=";")
enriched_genes <- unique(enriched_all$name)




all_genes <- length_data$name

goseq_vector = as.integer(all_genes %in% enriched_genes)

names(goseq_vector)<-all_genes

pwf <- nullp(goseq_vector,bias.data=length_data$length)


GO.wall = goseq(pwf,gene2cat=categories)
GO.nobias=goseq(pwf,gene2cat=categories,method="Hypergeometric")
plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.wall[,1],GO.nobias[,1]),2]), xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",xlim=c(-3,0))
abline(0,1,col=3,lty=2)
GO.nobias$adjusted <- p.adjust(GO.nobias$over_represented_pvalue, method="BH")


# GO.nobias.adjusted <- GO.nobias
# GO.nobias.adjusted$adjusted <- p.adjust(GO.nobias$over_represented_pvalue, method="BH")
# enriched.GO=p.adjust(GO.nobias$over_represented_pvalue, method="BH")
# #enriched.GO=GO.wall$category[GO.nobias$over_represented_pvalue<.05]
