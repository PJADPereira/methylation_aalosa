require(methylKit)
require(ggplot2)
require(genomation)
require(dplyr)
library(ggfortify)

## FUNCTION NECESSARY TO PREP THE FILES

makeList <- function(vec) {
  len <- length(vec)
  out <- as.list(rep(1:len), len)
  for (i in 1:len){
    out[i] <- vec[i]
  }
  out
}



setwd("./methyl_analysis/")


## LOAD THE FILES

tej_files <- "./bedGraphs/TEJ/individuals/"
cbd_files <- "./bedGraphs/CBD/individuals/"

tej_files_list <- makeList(paste0(tej_files,list.files(tej_files,pattern = "bismark.cov.gz")))
cbd_files_list <- makeList(paste0(cbd_files,list.files(cbd_files,pattern = "bismark.cov.gz")))
all_files <- c(tej_files_list,cbd_files_list)
samples <- strsplit(c(list.files(tej_files,pattern = "bismark.cov.gz"),list.files(cbd_files,pattern = "bismark.cov.gz")),"_local_bed.gz.bismark.cov.gz")
#condit <- c(rep(0,24),rep(1,24))
condit <- c(rep(0,18),rep(1,19))
indcover <- 3
myobj <- methRead(location = all_files,
                  sample.id = samples,
                  assembly="ncbi",
                  treatment= as.vector(condit),
                  context="CpG",
                  pipeline = "bismarkCoverage",
                  mincov = indcover)


### FILTER AND NORMALIZE COVERAGE

filtered_myobj <- filterByCoverage(myobj,hi.perc=99)
normalized_myobj <- normalizeCoverage(filtered_myobj)


### SETUP THRESHOLDS

chrms <-  c("NC_063189.1", "NC_063190.1", "NC_063191.1", "NC_063192.1", "NC_063193.1", "NC_063194.1", "NC_063195.1", "NC_063196.1", "NC_063197.1", "NC_063198.1", "NC_063199.1", "NC_063200.1", "NC_063201.1", "NC_063202.1", "NC_063203.1", "NC_063204.1", "NC_063205.1", "NC_063206.1", "NC_063207.1", "NC_063208.1", "NC_063209.1", "NC_063210.1", "NC_063211.1", "NC_063212.1")

wsize <- 2000
wcover <- 10
qval <-  0.05
diff <- 25
nsamples <- 10L


### LOAD THE COVARIATES - File covariate_Table.csv in github

covariates <- read.table("covariate_Table.csv",sep=";",header= T)

samples_to_order_covariates <-  as.vector(sapply(unlist(samples, use.names = FALSE),gsub,pattern="_local_bed.gz_filtered.bismark.cov.gz",replacement=""))
to_keep_covs <- subset(covariates,covariates$Sample %in% samples_to_order_covariates)
order_of_covs <- match(samples_to_order_covariates,to_keep_covs$Sample) 

to_keep_covs <- to_keep_covs[order_of_covs,]
to_keep_covs <- to_keep_covs[,c(2,3)]
to_keep_covs <- data.frame("Sex"=to_keep_covs$Sex)





### PREPARE THE WINDOWS


tiles = tileMethylCounts(normalized_myobj,win.size=wsize,step.size=(wsize/2),cov.bases = wcover,mc.cores=8)
meth=unite(tiles, destrand=FALSE,min.per.group = nsamples,mc.cores = 8)

### If necessary to rerun the analysis, save the meth object as a RDS in the file which can be loaded later bypassing the long tilling process
#saveRDS(meth,paste0("./individual_results/RDS/","wsize-",wsize,"_wcover-",wcover,"_diff-",diff,"_nsamples-",nsamples,"_indcover-",indcover,"_nwinwos-",nrow(meth),".RDS"))


### PLOT PCA
pca <- PCASamples(meth,adj.lim = c(1,0.1),obj.return = T,center = T,transpose=T)
data <- data.frame("sample"=meth@sample.ids,"treatment"=meth@treatment)
data$population <- "CBD"
data$population[which(data$treatment==0)] <- "TEJ"

## since adding this later, to ensure code integrity create a new dataframe
forplot_cov <- covariates
names(forplot_cov) <- c("sample","sex","year","pop")

data <- merge(data,forplot_cov,by="sample")

autoplot(pca, data = data, colour = 'population',shape="sex")+
  theme(
    #axis.title.x=element_text("Genome position"),
    plot.margin=unit(c(0,0,0,0), "null"),panel.spacing=unit(c(0,0,0,0), "null"),
    #axis.title.y=element_text("Coverage"),
    
    panel.background=element_blank()#,panel.border=element_blank(),panel.grid.major=element_blank(),
  )+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5),
        axis.text.x = element_text(family="sans",size=9),
        axis.text.y = element_text(family="sans",size=9),
        axis.title.x = element_text(family="sans",size=11),
        legend.text = element_text(family="sans",size=9),
        legend.title = element_text(family="sans",size=11),
  )




### CALCULATE DIFF METHYLATION WITH AND WITHOUT COVARIATES
# Year is not a good covariate because only one of the groups has a difference in year and therefore cannot be statistically estimated it's effect!

myDiff_all_cov=calculateDiffMeth(meth,covariates=to_keep_covs)
myDiff_all_nor=calculateDiffMeth(meth)
myDiff_cov <- myDiff_all_cov[which(myDiff_all_cov$chr %in% chrms),]
myDiff_nor <- myDiff_all_nor[which(myDiff_all_nor$chr %in% chrms),]
chrmcolors <- data.frame(chr=chrms,color=rep(c("black","darkgray"),times=length(chrms)/2))
diffmeth_cov <- data.frame(chr=myDiff_cov$chr,diff=myDiff_cov$meth.diff,qval=myDiff_cov$qvalue)
diffmeth_cov <- merge(diffmeth_cov,chrmcolors,by="chr",all=F)

diffmeth_nor <- data.frame(chr=myDiff_nor$chr,diff=myDiff_nor$meth.diff,qval=myDiff_nor$qvalue)
diffmeth_nor <- merge(diffmeth_nor,chrmcolors,by="chr",all=F)



dmc <- ggplot(diffmeth_cov)+
  geom_point(aes(x=1:length(diff),y=diff,colour=color))+
  labs(y = "different_methylation", x = "Genome position")+
  scale_color_identity() +
  geom_hline(yintercept =25)+
  geom_hline(yintercept=-25)+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        #axis.title.x=element_text("Genome position"),
        plot.margin=unit(c(0,0,0,0), "null"),panel.spacing=unit(c(0,0,0,0), "null"),
        #axis.title.y=element_text("Coverage"),
        legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

dmc


dmn <- ggplot(diffmeth_nor)+
  geom_point(aes(x=1:length(diff),y=diff,colour=color))+
  labs(y = "different_methylation", x = "Genome position")+
  scale_color_identity() +
  geom_hline(yintercept =25)+
  geom_hline(yintercept=-25)+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        #axis.title.x=element_text("Genome position"),
        plot.margin=unit(c(0,0,0,0), "null"),panel.spacing=unit(c(0,0,0,0), "null"),
        #axis.title.y=element_text("Coverage"),
        legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

dmn




myDiff25pall=getMethylDiff(myDiff_nor,difference=diff,qvalue=qval)
myDiff25pcov=getMethylDiff(myDiff_cov,difference=diff,qvalue=qval)


### MEASURE WHAT TYPE OF FEATURES ARE DIFFERENTIALLY METHYLATED


gene.obj <- readTranscriptFeatures("./ncbi_1.1/genomic_transcoder.bed12",remove.unusual=FALSE)
#a <- annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
#a <- annotateWithFeatures(as(myDiff25p,"GRanges"),gff)

diffAnn=annotateWithGeneParts(as(myDiff25pcov,"GRanges"),gene.obj)
associations <- getAssociationWithTSS(diffAnn)
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

myDiff25ptable <- getData(myDiff25pcov)
myDiff25ptable$linenumber <- 1:nrow(myDiff25ptable)
associations$linenumber <- associations$target.row
myDiff25ptable <- merge(myDiff25ptable,associations,by="linenumber")











## go enrichment
library(goseq)

##need to filter by to just keep the genes in the main chromossomes in case it adds to some aritficial imbalance when i diddnt allow some genes to be conseidred for differences in methylation
setwd("/methyl_analysis/goanalysis/")

### PARSE STRINGS RESULTS USING build_gene_to_go.py

length_data <- read.table("/methyl_analysis/goanalysis/STRINGS/src/gene_length.csv")
names(length_data) <- c('name','length')

categories <- read.table("/methyl_analysis/goanalysis/STRINGS/src/goterms.csv")
names(categories) <- c("name","category")

enriched_genes <- unique(sapply(strsplit(associations$feature.name,";"),"[",3))
#enriched_genes_promoters <- unique(sapply(strsplit(associationspromoter$feature.name,";"),"[",3))


all_genes <- length_data$name

goseq_vector = as.integer(all_genes %in% enriched_genes)
#goseq_promoters_vector = as.integer(all_genes %in% enriched_genes_promoters)
names(goseq_vector)<-all_genes
#names(goseq_promoters_vector) <- all_genes

pwf <- nullp(goseq_vector,bias.data=length_data$length)
#pwf_promoter <- nullp(goseq_promoters_vector,bias.data=length_data$length)


GO.wall = goseq(pwf,gene2cat=categories)
GO.nobias=goseq(pwf,gene2cat=categories,method="Hypergeometric")
plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.wall[,1],GO.nobias[,1]),2]), xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",xlim=c(-3,0))
abline(0,1,col=3,lty=2)
GO.nobias$adjusted <- p.adjust(GO.nobias$over_represented_pvalue, method="BH")
GO.wall$adjusted <- p.adjust(GO.wall$over_represented_pvalue,method="BH")

write.table(GO.wall,"wall_all.csv")
write.table(GO.nobias,"nobias_all.csv")
