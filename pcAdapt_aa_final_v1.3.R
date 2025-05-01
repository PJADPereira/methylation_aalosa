library(pcadapt)
library(data.table)
library(qvalue)
library(ggplot2)
library(gridExtra)
setwd("./methyl_analysis/")

#data <- fread("aalo_snape_20kbwindows.pcadapt",sep="\t",header = TRUE)
data <- fread("aalo_snape_20snpswindows_finalset.pcadapt",sep="\t",header = TRUE)

pops <- data.frame(data[,1])[,1]
data <- read.pcadapt(data[,2:dim(data)[2]], type = "pool")

rownames(data) <- pops 

res <- pcadapt(data,K=4,min.maf = 0)


plot(res,option="screeplot")


#plot(res, option = "qqplot")

plot(res, option = "scores", pop = rownames(data))

## only use pca2 if there are components 3 and 4, if not adjust or comment it out


hist(res$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

#plot(-log10(res$pvalues))

chrms <-  c("NC_063189.1", "NC_063190.1", "NC_063191.1", "NC_063192.1", "NC_063193.1", "NC_063194.1", "NC_063195.1", "NC_063196.1", "NC_063197.1", "NC_063198.1", "NC_063199.1", "NC_063200.1", "NC_063201.1", "NC_063202.1", "NC_063203.1", "NC_063204.1", "NC_063205.1", "NC_063206.1", "NC_063207.1", "NC_063208.1", "NC_063209.1", "NC_063210.1", "NC_063211.1", "NC_063212.1")

window_chrms <- sapply(strsplit(colnames(data), split = ":"),"[",1)
chrmcolors <- data.frame(chr=chrms,color=rep(c("black","darkgray"),times=length(chrms)/2))
start_positions <- sapply(strsplit(sapply(strsplit(colnames(data), split = ":"),"[",2),split="-"),"[",1)
stop_positions <- sapply(strsplit(sapply(strsplit(colnames(data), split = ":"),"[",2),split="-"),"[",2)

results <- data.frame(chr=window_chrms,pvalues=res$pvalues,start=start_positions,stop=stop_positions)
results <- merge(results,chrmcolors,by="chr")
results$padj <- p.adjust(results$pvalues,method="fdr")
percentile1cutoff <-  sort(results$pvalues,decreasing = FALSE)[ceiling(length(res$pvalues)*0.01)]
percentile.1cutoff <-  sort(results$pvalues,decreasing = FALSE)[ceiling(length(res$pvalues)*0.001)]
results$line <- 1:dim(results)[1]
#results$color[which(results$pvalues<percentile1cuttoff)] <- "red"


mh <- ggplot(results,aes(x=line))+
  geom_point(aes(x=line,y=-log10(pvalues),colour=color))+
  labs(y = "-log10(pvalues)", x = "Genome position")+
  geom_hline(yintercept = -log10(percentile1cutoff))+
  #geom_hline(yintercept = -log10(percentile.1cutoff))+
  scale_color_identity() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        #axis.title.x=element_text("Genome position"),
        plot.margin=unit(c(0,0,0,0), "null"),panel.spacing=unit(c(0,0,0,0), "null"),
        #axis.title.y=element_text("Coverage"),
        legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())





png("/methyl_analysis/manhattan_finalset_v1.3.png",height=6,width=12,units = "in",res=420)
plot(mh)
dev.off()
