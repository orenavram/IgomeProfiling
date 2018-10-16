library("RColorBrewer")
library("gplots")
#setwd("/groups/pupko/orenavr2/gershoni/Experiments/Exp_DP_4/analyses/1_ALL/Motifs/")
setwd("/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/1_ALL/mAbs_motifs_vs_mAbs_hits")
#setwd("/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_5/analyses/Analysis1/results")
#DataFile="merged_df_top_k_features.csv"
#setwd("/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_2/analyses/9_donors_vs_585and3662")
DataName="mAbs_motifs_merged_hits_df"
DataFile=paste(DataName, ".csv", sep="")
Data=read.table(file=DataFile,sep=",",header=T)
row.names(Data) <- Data$sample_name
Data= Data [,-(dim(Data)[2])]  # remove sample name
#Data= Data [,-(dim(Data)[2])]  # remove label
#Data= Data [,-1]  # remove index

names(Data)
ColSideColors = c(    # grouping col-variables into different bioogical conditions motifs
 rep("green", 156),    # First 160 columns: green C300
 rep("orange", 168),   # Next 252 columns: blue - C300S
 rep("blue", 169),      # Next 383 columns: yellow - H5N1
rep("black", 157),    # Next 349 columns: black - Naive
rep("red", 144),    # Next 349 columns: black - Naive
rep("yellow", 159))    # Next 349 columns: black - Naive

# ColSideColors = c(    # grouping col-variables into different bioogical conditions motifs
#   rep("green", 8),    # First 160 columns: green C300
#   rep("orange", 4),   # Next 252 columns: blue - C300S
#   rep("blue", 7),      # Next 383 columns: yellow - H5N1
#   rep("black", 4),    # Next 349 columns: black - Naive
#   rep("red", 2),    # Next 349 columns: black - Naive
#   rep("yellow", 3))    # Next 349 columns: black - Naive

NumOfMotifs=length(Data)
pdf(paste(DataName, ".pdf", sep=""), width = 18, height = 16)
#pdf("heat_map_all_significant_motifs_exp14.pdf", width = 11.69, height = 8.27)
#pdf("ferrets_vs_texas.pdf", width = 11.69, height = 8.27)
res=heatmap.2(as.matrix(log(Data+1)),col=brewer.pal(8,"Blues"),trace="none",dendrogram='row', Rowv=TRUE, Colv=FALSE,cexRow=1.1,cexCol=1.1,labCol = "",xlab=paste(NumOfMotifs,"motifs"),margins=c(3,10),key=TRUE, symkey=FALSE, density.info="none",key.xlab="Log(Number of hits+1)", ColSideColors=ColSideColors)
#res=heatmap.2(as.matrix(-log(Data)),col=brewer.pal(8,"Blues"),trace="none",dendrogram='row', Rowv=TRUE, Colv=FALSE,cexRow=1.1,cexCol=1.1,labCol = "",xlab=paste(NumOfMotifs,"motifs"),margins=c(3,10),key=TRUE, symkey=FALSE, density.info="none",key.xlab="-Log2(pvalue)", ColSideColors=ColSideColors)
#res=heatmap.2(as.matrix(log(Data+1)),col=brewer.pal(8,"Blues"),trace="none",dendrogram='row', Rowv=TRUE, Colv=FALSE,cexRow=1.1,cexCol=1.1,labCol = "",xlab=paste(NumOfMotifs,"motifs"),margins=c(3,10),key=TRUE, symkey=FALSE, density.info="none",key.xlab="Log(Number of hits)", ColSideColors=ColSideColors)
par(lend = 20)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       # legend = c('(8) mAb_c54', '(4) mAb_CR8020', '(7) mAb_c108', '(4) mAb_c585', '(2) mAb_F045', '(3) mAb_c3662'), # category labels
       legend = c('(157) mAb_c54', '(169) mAb_CR8020', '(170) mAb_c108', '(158) mAb_c585', '(145) mAb_c3662', '(160) mAb_F045'), # category labels
       col = c("green", "orange", "blue", "black", "red", "yellow"),  # color key
      lty= 1,             # line style
      lwd = 10            # line width
)
dev.off()