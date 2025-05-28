require(methylKit)
require(ggplot2)
require(genomation)

makeList <- function(vec) {
  len <- length(vec)
  out <- as.list(rep(1:len), len)
  for (i in 1:len){
    out[i] <- vec[i]
  }
  out
}



setwd("./methyl_analysis/revisions/methyilanalysis_cpg")


tej_files <- "./filteredbedgraphs/TEJ/"
cbd_files <- "./filteredbedgraphs/CBD/"

tej_files_list <- makeList(paste0(tej_files,list.files(tej_files,pattern = "bismark.cov.gz")))
cbd_files_list <- makeList(paste0(cbd_files,list.files(cbd_files,pattern = "bismark.cov.gz")))
all_files <- c(tej_files_list,cbd_files_list)
samples <- strsplit(c(list.files(tej_files,pattern = "bismark.cov.gz"),list.files(cbd_files,pattern = "bismark.cov.gz")),"_local_bed.gz.bismark.cov.gz")
condit <- c(rep(0,24),rep(1,24))
#condit <- c(rep(0,18),rep(1,19))
indcover <- 3
myobj <- methRead(location = all_files,
                  sample.id = samples,
                  assembly="ncbi",
                  treatment= as.vector(condit),
                  context="CpG",
                  pipeline = "bismarkCoverage",
                  mincov = indcover)



filtered_myobj <- filterByCoverage(myobj,hi.perc=99)
normalized_myobj <- normalizeCoverage(filtered_myobj)

for (i in 1:48){
  coverage <- data.frame(chr = myobj@.Data[[i]]$chr, coverage= myobj@.Data[[i]]$coverage)
  print(mean(coverage$coverage))
  
}

chrms <-  c("NC_063189.1", "NC_063190.1", "NC_063191.1", "NC_063192.1", "NC_063193.1", "NC_063194.1", "NC_063195.1", "NC_063196.1", "NC_063197.1", "NC_063198.1", "NC_063199.1", "NC_063200.1", "NC_063201.1", "NC_063202.1", "NC_063203.1", "NC_063204.1", "NC_063205.1", "NC_063206.1", "NC_063207.1", "NC_063208.1", "NC_063209.1", "NC_063210.1", "NC_063211.1", "NC_063212.1")

qval <-  0.05
diff <- 25
nsamples <- 10L


meth=unite(normalized_myobj, destrand=FALSE,min.per.group = nsamples,mc.cores = 8)


PCASamples(meth,adj.lim=c(1,0.1))#,adj.lim = c(1,0.1))


myDiff_all=calculateDiffMeth(meth)
myDiff <- myDiff_all[which(myDiff_all$chr %in% chrms),]
chrmcolors <- data.frame(chr=chrms,color=rep(c("black","darkgray"),times=length(chrms)))
diffmeth <- data.frame(chr=myDiff$chr,diff=myDiff$meth.diff,qval=myDiff$qvalue)
diffmeth <- merge(diffmeth,chrmcolors,by="chr",all=TRUE)



dm <- ggplot(diffmeth)+
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

dm

myDiff25p=getMethylDiff(myDiff,difference=diff,qvalue=qval)

file_name <- paste0("cov10wsize-",wsize,"_wcover-",wcover,"_diff-",diff,"_nsamples-",nsamples,"_indcover-",indcover,"_nwinwos-",nrow(meth),".csv")
png_name  <- paste0("cov10wsize-",wsize,"_wcover-",wcover,"_diff-",diff,"_nsamples-",nsamples,"_indcover-",indcover,"_nwinwos-",nrow(meth),".png")


gene.obj <- readTranscriptFeatures("./alosa_genome/ncbi_1.1/genomic_transcoder.bed12",remove.unusual=FALSE)


diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
associations <- getAssociationWithTSS(diffAnn)
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

myDiff25ptable <- getData(myDiff25p)
myDiff25ptable$linenumber <- 1:nrow(myDiff25ptable)
associations$linenumber <- associations$target.row
myDiff25ptable <- merge(myDiff25ptable,associations,by="linenumber")

write.table(myDiff25ptable,"./finalcov10xwcover10icover3nsamples12_annotated_all.csv",row.names = F)



