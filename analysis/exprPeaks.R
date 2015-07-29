############################################
#
#  To run this file, use the folowing command -> 
#   source('~/work/research/researchProjects/encode/ctcf/analysis/main.R')
#
#

dataDir <- getFullPath("data")

annot.file <- paste0(dataDir,"/peakCount_foldChange_annot.tab")
peak.file <- paste0(dataDir,"/peakCount_foldChange.tab")
expr.file <- paste0(dataDir,"/exonCosi_dsRegion.tab")
expr.psi.file <- paste0(dataDir,"/exonPsi_dsRegion.tab")

peakExpr.file <- paste0(dataDir,"/exonCosi_peakCount_foldChange.tab")
peakExpr.psi.file<- paste0(dataDir,"/exonCosi_peakCount_foldChange_PSI.tab")

getPeakColsNames_cosi <- function(){
  peak <- read.csv(peak.file, sep="\t", stringsAsFactors=FALSE)
  annot <- read.csv(annot.file, sep="\t", stringsAsFactors=FALSE)
  expr <- read.csv(expr.file, sep="\t", stringsAsFactors=FALSE)
  expr.cells <- colnames(expr)[4:21]
  annot.cells.tmp <- annot[which(annot$cell %in% expr.cells &((annot$replicate == 1))),]
  annot.drop <- c(13,16,19,23,24,25,28)
  found <- c()
  drop <- c()
  for(i in seq_along(annot.cells.tmp$cell)){
    localCell <- annot.cells.tmp$cell[i]
    if(localCell %in% found){
      drop <- c(drop,i)
    } else{
      found <- c(found,localCell)
    }
  }
  annot.cells <- annot.cells.tmp[-drop,]
  
  expr.cols <- annot.cells$cell
  peak.cols <- annot.cells$tag
  
  expr.match <- expr[c("label",expr.cols)]
  peak.count.match <- peak[c("label",paste("peakCount",peak.cols,sep="."))]
  colnames(peak.count.match) <- c("label",expr.cols)
  peak.fold.match <- peak[c("label",paste("sumFoldChange",peak.cols,sep="."))]
  colnames(peak.fold.match) <- c("label",paste(expr.cols,"foldChange",sep="."))
  
  combtmp <- merge(expr.match,peak.count.match,by="label",suffix=c(".expr",".peakCount"))
  comb <- merge(combtmp,peak.fold.match,by="label")
  write.table(comb,
              file=peakExpr.file,sep="\t",quote=FALSE,row.names=FALSE)
}
labelToView <- function(label="chr16:27709649-27709821_+"){
  chr = as.character(unlist(strsplit(label,split=":"))[1])
  label2 = as.character(unlist(strsplit(label,split=":"))[2])
  strand = as.character(unlist(strsplit(label2,split="_"))[2])
  ranges = as.character(unlist(strsplit(label2,split="_"))[1])
  startPos = as.numeric(unlist(strsplit(ranges,split="-"))[1])
  stopPos = as.numeric(unlist(strsplit(ranges,split="-"))[1])
  viewStart <- ifelse(strand=="+",startPos,startPos-1000)
  viewStop <- ifelse(strand=="+",stopPos+1000,stopPos)
  paste0(chr,":",viewStart,"-",viewStop)
  
}
getPeakColsNames_PSI <- function(){
  peak <- read.csv(peak.file, sep="\t", stringsAsFactors=FALSE)
  annot <- read.csv(annot.file, sep="\t", stringsAsFactors=FALSE)
  expr <- read.csv(expr.psi.file, sep="\t", stringsAsFactors=FALSE)
  expr.cells <- colnames(expr)[4:21]
  annot.cells.tmp <- annot[which(annot$cell %in% expr.cells &((annot$replicate == 1))),]
  annot.drop <- c(13,16,19,23,24,25,28)
  found <- c()
  drop <- c()
  for(i in seq_along(annot.cells.tmp$cell)){
    localCell <- annot.cells.tmp$cell[i]
    if(localCell %in% found){
      drop <- c(drop,i)
    } else{
      found <- c(found,localCell)
    }
    
  }
  annot.cells <- annot.cells.tmp[-drop,]
  expr.cols <- annot.cells$cell
  peak.cols <- annot.cells$tag
  
  expr.match <- expr[c("label",expr.cols)]
  peak.count.match <- peak[c("label",paste("peakCount",peak.cols,sep="."))]
  colnames(peak.count.match) <- c("label",expr.cols)
  peak.fold.match <- peak[c("label",paste("sumFoldChange",peak.cols,sep="."))]
  colnames(peak.fold.match) <- c("label",paste(expr.cols,"foldChange",sep="."))
  
  combtmp <- merge(expr.match,peak.count.match,by="label",suffix=c(".expr",".peakCount"))
  comb <- merge(combtmp,peak.fold.match,by="label")
  write.table(comb,
              file=peakExpr.psi.file,sep="\t",quote=FALSE,row.names=FALSE)
}
countNA <- function(vec){
  sum(is.na(vec))
}
countNotNa<- function(vec){
  sum(!is.na(vec))
}

sumRmNa <- function(vec){
  sum(vec,rm.na=TRUE)
}

applyFnByNa <- function(vecNA,vec2,fn){
  stopifnot( (length(vecNA) == length(vec2)))
  idx <- which(!is.na(vecNA))
  vec2good <- vec2[idx]
  do.call(fn,list(vec2good))
}

applyCorrByNa <- function(vecNA,vec2,method){
  stopifnot( (length(vecNA) == length(vec2)))
  if(sum(!is.na(vecNA)) <3){
    return(NA)
  }
  idx <- which(!is.na(vecNA))
  vecNAgood <- vecNA[idx]
  vec2good <- vec2[idx]
  cor(vecNA[idx],vec2[idx],method=method)
}

applyCorrByNaPvalue <- function(vecNA,vec2,method){
  stopifnot( (length(vecNA) == length(vec2)))
  if(sum(!is.na(vecNA)) <3){
    return(NA)
  }
  idx <- which(!is.na(vecNA))
  vecNAgood <- vecNA[idx]
  vec2good <- vec2[idx]
  cor.test(vecNA[idx],vec2[idx],method=method)$p.value
}

getNumberOfNonNaBoth <- function(vecNA,vec2,method){
  stopifnot( (length(vecNA) == length(vec2)))
  if(sum(!is.na(vecNA)) <3){
    return(NA)
  }
  # NA * NA -> NA, NA * (real number) -> NA, (real) * (real) -> (real)
  idx <- which(!is.na(vecNA * vec2))
  length(idx)
}

applyCorrByNaColsBothCols <- function(vecNA,vec2,method){
  stopifnot( (length(vecNA) == length(vec2)))
  if(sum(!is.na(vecNA)) <3){
    return(NA)
  }
  # NA * NA -> NA, NA * (real number) -> NA, (real) * (real) -> (real)
  idx <- which(!is.na(vecNA * vec2))
  if(length(idx) < 3){return(NA)}
  #dput(vecNA[idx])
  #dput(vec2[idx])
  cor(vecNA[idx],vec2[idx],method=method)
}
applyCorrByNaPvalueBothCols <- function(vecNA,vec2,method){
  stopifnot( (length(vecNA) == length(vec2)))
  if(sum(!is.na(vecNA)) <3){
    return(NA)
  }
  idx <- which(!is.na(vecNA* vec2))
  if(length(idx) < 3){return(NA)}
  #dput(vecNA[idx])
  #dput(vec2[idx])
  cor.test(vecNA[idx],vec2[idx],method=method)$p.value
}

#wilcox.test(1:3,2:4,paired=TRUE)$p.value
applyWilcoxTestByNaPvalue <- function(vecNA,vec2,method){
  stopifnot( (length(vecNA) == length(vec2)))
  if(sum(!is.na(vecNA)) <3){
    return(NA)
  }
  idx <- which(!is.na(vecNA))
  v1 <- vecNA[idx]
  v2 <- vec2[idx]
  
  v.low <- v1[which(v2 == 0)]
  v.high <- v1[which(v2 > 0)]
  if(length(v.low) < 3 || length(v.high) < 3){
    return(NA)
  }
  
  wilcox.test(v.low,v.high)$p.value
}
fracGreaterZero <- function(vec){
  sum(vec>0)/length(vec)
}

shuffle <- function(vec){
  vecClass <- class(vec)
  out <- sample(x=vec,size=length(vec),replace=FALSE)
  as(object=out,Class=vecClass)
} 

shuffleRows <- function(df,cols){
  colClasses <- as.character(apply(df[,cols],2, class))
  df[,cols] <- apply(df[,cols],2,shuffle)
  for(i in seq_along(cols)){
    df[,i] <- as(object=df[,i],Class=colClasses[i])
  }
  df
}

zeroToNa <-function(vec){
  vec[vec == 0] <- NA
  vec
}

replaceRowZeroToNa <- function(df,cols,value=0,replacement=NA){
  colClasses <- as.logical(apply(df[,cols],2, is.numeric))
  if(sum(colClasses) != length(colClasses)){
    stop(" all specified df columns must be of type numeric")
  }
  
  df[,cols] <- apply(df[,cols],2,zeroToNa)
  df
}



getCorrPeaksCosi <- function(){
  plotDir <- getFullPath("plots")
  outdir <- file.path(plotDir,"CoSI-FoldChange/")  
  if(!file.exists(outdir)){dir.create(path=outdir,recursive=TRUE)}
  
  
  pe <- read.csv(peakExpr.file, sep="\t", stringsAsFactors=FALSE)
  pe.cols <- colnames(pe)
  pe.expr.cols <- pe.cols[grep(x=pe.cols,pattern="expr")]
  pe.fold.cols <- pe.cols[grep(x=pe.cols,pattern="foldChange")]
  pe.count.cols <- pe.cols[grep(x=pe.cols,pattern="peakCount")]
  pe$peakCountSumRaw <- apply(pe[pe.count.cols],1,sum)
  pe$peakCountSum   <- apply(pe[c(pe.expr.cols,pe.count.cols)],1,
                             function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                                    x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.count.cols))],
                                                    sum))
  
  
  pe$exprCount <- apply(pe[pe.expr.cols],1,countNotNa)
  pe$coSIave   <- apply(pe[pe.expr.cols],1,function(x)mean(x,na.rm=TRUE))
  #pe$FoldChangeave   <- apply(pe[pe.fold.cols],1,function(x)mean(x,na.rm=TRUE))
  pe$FoldChangeave   <- apply(pe[c(pe.expr.cols,pe.fold.cols)],1,
                              function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                                     x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                     mean))
  
  pe$FoldChangeVar  <- apply(pe[c(pe.expr.cols,pe.fold.cols)],1,
                             function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                                    x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                    var))
  pe$coSIVar   <- apply(pe[c(pe.expr.cols,pe.expr.cols)],1,
                        function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                               x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.expr.cols))],
                                               var))
  
  
  
  #pe$FoldChangeaveGroup <- 
  # fraction of expressed exons with peaks
  pe$peakFracExpr   <- apply(pe[c(pe.expr.cols,pe.count.cols)],1,
                             function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                                    x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.count.cols))],
                                                    fracGreaterZero))
  pe$peakFracExprGroup   <- cut(x=pe$peakFracExpr,breaks=(1:length(pe.count.cols))/length(pe.count.cols))
  
  pe$pearsonCorr         <- apply(pe[c(pe.expr.cols,pe.fold.cols)],1,
                                  function(x)applyCorrByNa(x[1:length(pe.expr.cols)],
                                                           x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                           method="pearson"))
  pe$pearsonCorrPvalue   <- apply(pe[c(pe.expr.cols,pe.fold.cols)],1,
                                  function(x)applyCorrByNaPvalue(x[1:length(pe.expr.cols)],
                                                                 x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                 method="pearson"))
  
  #wilcox.test(1:3,2:4,paired=TRUE)$p.value
  
  pe$wilcoxCorrPvalue   <- apply(pe[c(pe.expr.cols,pe.count.cols)],1,
                                 function(x)applyWilcoxTestByNaPvalue(x[1:length(pe.expr.cols)],
                                                                      x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.count.cols))]))
  
  pe.nona <- pe[!is.na(pe$pearsonCorrPvalue),]
  pe.pearson <- pe.nona[order(pe.nona$pearsonCorrPvalue,decreasing=FALSE),]
  pe.pearson$alpha <- 0.05
  pe.pearson$M <- length(pe.nona$pearsonCorrPvalue)
  pe.pearson$j <- seq_along(pe.nona$pearsonCorrPvalue) 
  FDRpass <- with(pe.pearson, which(pearsonCorrPvalue < alpha*(j/M)))
  pe.pearson.fdr <- pe.pearson[FDRpass,]
  
  pw.nona <- pe[!is.na(pe$wilcoxCorrPvalue),]
  pw.wilcox <- pw.nona[order(pw.nona$wilcoxCorrPvalue,decreasing=FALSE),]
  pw.wilcox$alpha <- 0.05
  pw.wilcox$M <- length(pw.wilcox$wilcoxCorrPvalue)
  pw.wilcox$j <- seq_along(pw.wilcox$wilcoxCorrPvalue) 
  FDRpassW <- with(pw.wilcox, which(wilcoxCorrPvalue < alpha*(j/M)))
  pw.wilcox.fdr <- pw.wilcox[FDRpassW,]
  
  pe.peakShuffle <- shuffleRows(pe,pe.fold.cols)
  pe.peakShuffle$pearsonCorrPvalue   <- apply(pe.peakShuffle[c(pe.expr.cols,pe.fold.cols)],1,
                                              function(x)applyCorrByNaPvalue(x[1:length(pe.expr.cols)],
                                                                             x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                             method="pearson"))
  pe.peakShuffle$pearsonCorr   <- apply(pe.peakShuffle[c(pe.expr.cols,pe.fold.cols)],1,
                                        function(x)applyCorrByNa(x[1:length(pe.expr.cols)],
                                                                 x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                 method="pearson"))
  
  pes.nona <- pe.peakShuffle[!is.na(pe.peakShuffle$pearsonCorrPvalue),]
  pes.pearson <- pes.nona[order(pes.nona$pearsonCorrPvalue,decreasing=FALSE),]
  pes.pearson$alpha <- 0.05
  pes.pearson$M <- length(pes.pearson$pearsonCorrPvalue)
  pes.pearson$j <- seq_along(pes.pearson$pearsonCorrPvalue) 
  FDRpassShuffle <- with(pes.pearson, which(pearsonCorrPvalue < alpha*(j/M)))
  pes.pearson.fdr <- pes.pearson[FDRpassShuffle,]
  
  
  
  pe.peakShuffle.corr <- shuffleRows(pe.nona,pe.fold.cols)
  pe.peakShuffle.corr$pearsonCorrPvalue   <- apply(pe.peakShuffle.corr[c(pe.expr.cols,pe.fold.cols)],1,
                                                   function(x)applyCorrByNaPvalue(x[1:length(pe.expr.cols)],
                                                                                  x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                                  method="pearson"))
  pe.peakShuffle.corr$pearsonCorr   <- apply(pe.peakShuffle.corr[c(pe.expr.cols,pe.fold.cols)],1,
                                             function(x)applyCorrByNa(x[1:length(pe.expr.cols)],
                                                                      x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                      method="pearson"))
  
  
  
  
  pesc.nona <- pe.peakShuffle.corr[!is.na(pe.peakShuffle.corr$pearsonCorrPvalue),]
  pesc.pearson <- pesc.nona[order(pesc.nona$pearsonCorrPvalue,decreasing=FALSE),]
  pesc.pearson$alpha <- 0.05
  pesc.pearson$M <- length(pesc.pearson$pearsonCorrPvalue)
  pesc.pearson$j <- seq_along(pesc.pearson$pearsonCorrPvalue) 
  FDRpassShuffleCorr <- with(pesc.pearson, which(pearsonCorrPvalue < alpha*(j/M)))
  pesc.pearson.fdr <- pesc.pearson[FDRpassShuffleCorr,]
  
  
  pe.short <- pe[1:10000,]
  #   ggplot(pe[pe$FoldChangeave > 0,], aes(x=FoldChangeaveGroup,y=coSIave))+geom_boxplot() + 
  #     xlab("ave CTCF peaks fold change in all celltypes") + ylab("ave coSI in all cell type") + theme_bw()+
  #     ggtitle("average coSI for CTCF peaks found\nversus average fold change above bkgd")
  #   ggsave(paste0(outdir,"foldChangeGroup-avecoSI-boxplot.png"), height=5,width=10)
  
  ggplot(pe, aes(x=FoldChangeave,y=coSIave,fill=peakCountSum))+geom_point() + 
    xlab("ave CTCF peaks fold change") + ylab("ave coSI") + theme_bw()+
    ggtitle("average coSI for CTCF peaks found")
  ggsave(paste0(outdir,"foldChangeAve-avecoSI-point.png"), height=7,width=7)
  
  #wide
  
  ggplot(pe[which(pe$peakCountSum %in% c(0,15)),], aes(coSIave,fill=factor(peakCountSum)))+geom_density(alpha=I(0.5)) + 
    xlab("average coSI") + ylab("density") + theme_bw()+
    ggtitle("average coSI\nCTCF occupied in all celltypes vs. none")
  ggsave(paste0(outdir,"coSI-density-colorByPeakCount.png"), height=7,width=7)
  
  
  
  ggplot(pe, aes(x=peakCountSum,y=coSIave,fill=factor(peakCountSum)))+geom_boxplot() + 
    xlab("CTCF peaks found") + ylab("density") + theme_bw()+
    ggtitle("average coSI for CTCF peaks found")
  ggsave(paste0(outdir,"coSI-peakCountSum-colorByPeakCount-boxplot.png"), height=7,width=10)
  
  
  ggplot(pe[which(pe$peakCountSum > 0),], aes(peakCountSum))+geom_bar(binwidth=1) + 
    scale_y_log10() + xlab("CTCF peaks") + ylab("count of exons") + 
    ggtitle("number of CTCF binding events detected\ndownstream of exon\nin all cell types")
  ggsave(paste0(outdir,"coSI-peakCountSum-bar.png"), height=7,width=10)
  
  ggplot(pe[which(pe$peakCountSum > 0),], aes(x=peakFracExpr,y=peakCountSum))+geom_point() + 
    xlab("fraction of expressed exons with CTCF peaks") + ylab("number downstream CTCF peaks") + 
    ggtitle("number of CTCF binding events detected\ndownstream of exon\nin all cell types") +
    theme_bw()
  ggsave(paste0(outdir,"peakCountSum-peakFrac-point.png"), height=7,width=7)
  
  ggplot(pe, aes(x=peakFracExpr,y=coSIave))+geom_point() + 
    xlab("fraction of expressed exons with CTCF peaks") + ylab("average coSI") + theme_bw()+
    ggtitle("frac CTCF peaks vs. ave coSI")
  ggsave(paste0(outdir,"peakFracExpr-coSIave-point.png"), height=7,width=7)
  
  
  ggplot(pe, aes(x=peakFracExprGroup,y=coSIave))+geom_boxplot() + 
    xlab("fraction of expressed exons with CTCF peaks") + ylab("mean coSI") + theme_bw()+
    ggtitle("fraction of expressed exons w/ CTCF peaks \nvs. mean coSI of expressed exons")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(outdir,"peakFracExpr-coSIave-boxplot.png"), height=7,width=12)
  
  ggplot(pe[which(!is.na(pe$coSIVar) & !is.na(pe$coSIVar)),], aes(x=log10(FoldChangeVar),y=coSIVar))+geom_point() + 
    xlab("variance of CTCF peak fold change") + ylab("variance of coSI") + theme_bw()+
    ggtitle("variance in CTCF foldChange vs. variance coSI\nin expressed exons only")
  ggsave(paste0(outdir,"foldChangeVar-coSIVar-point.png"), height=7,width=7)
  
  ggplot(pe, aes(x=-log10(pearsonCorrPvalue)))+geom_density() + 
    xlab("-log10(pearsonCorr P-value)") + ylab("density") + theme_bw()+
    ggtitle("Pearson Correlation P-value Between\n coSI and foldChange\nof expressed exons")
  ggsave(paste0(outdir,"coSI-CorrPvalue-density.png"), height=7,width=7)
  
  ggplot(pe.nona, aes(x=-log10(abs(pearsonCorr))))+geom_density() + 
    xlab("-log10(abs(pearsonCorr))") + ylab("density") + theme_bw()+
    ggtitle("Pearson Correlation Between\n coSI and foldChange\nof expressed exons")
  ggsave(paste0(outdir,"coSI-LogCorr-density.png"), height=7,width=7)
  
  ggplot(pe.nona, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle("Pearson Correlation Between\n coSI and foldChange\nof expressed exons")
  ggsave(paste0(outdir,"coSI-LogCorr-density.png"), height=7,width=7)
  
  
  gplot(pe.nona, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("Pearson Correlation Coeff\nBetween coSI and foldChangeof expressed exons\n(in at least 3 cell types)\nN=",
                  dim(pe.nona)[1]))
  ggsave(paste0(outdir,"coSI-LogCorr-density.png"), height=7,width=7)
  
  ggplot(pe.pearson.fdr, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("Pearson Correlation Coeff\nBetween coSI and foldChangeof expressed exons\n(pass FDR of 0.05)\nN=",
                  dim(pe.pearson.fdr)[1]))
  ggsave(paste0(outdir,"coSI-LogCorr-FDR_pass-density.png"), height=7,width=7)
  
  
  
  
  ggplot(pes.nona, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("Shuffled downstream regions\nPearson Correlation Coeff\nBetween coSI and foldChangeof shuffled exons\n(in at least 3 cell types)\nN=",
                  dim(pes.nona)[1]))
  ggsave(paste0(outdir,"coSI-shuffle-LogCorr-gt3Expr-density.png"), height=7,width=7)
  
  ggplot(pes.pearson.fdr, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("Shuffled downstream regions\nPearson Correlation Coeff\nBetween coSI and foldChangeof shuffled exons\n(pass FDR of 0.05)\nN=",
                  dim(pes.pearson.fdr)[1]))
  ggsave(paste0(outdir,"coSI-shuffle-LogCorr-FDR_pass-density.png"), height=7,width=7)
  
  
  ggplot(pesc.nona, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("Shuffled downstream regions(regions w/ defined corr only)\nPearson Correlation Coeff\nBetween coSI and foldChangeof shuffled exons\n(in at least 3 cell types)\nN=",
                  dim(pesc.nona)[1]))
  ggsave(paste0(outdir,"coSI-shuffleCorrOnly-LogCorr-gt3Expr-density.png"), height=7,width=7)
  
  ggplot(pesc.pearson.fdr, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("Shuffled downstream regions(regions w/ defined corr only)\nPearson Correlation Coeff\nBetween coSI and foldChangeof shuffled exons\n(pass FDR of 0.05)\nN=",
                  dim(pesc.pearson.fdr)[1]))+xlim(-1,1)
  ggsave(paste0(outdir,"coSI-shuffleCorrOnly-LogCorr-FDR_pass-density.png"), height=7,width=7)
  
}


getCorrPeaksPSI <- function(){
  plotDir <- getFullPath("plots")
  outdir <- file.path(plotDir,"PSI-FoldChange/")  
  if(!file.exists(outdir)){dir.create(path=outdir,recursive=TRUE)}
  
  
  
  pe <- unique(read.csv(peakExpr.psi.file, sep="\t", stringsAsFactors=FALSE))
  pe.cols <- colnames(pe)
  pe.expr.cols <- pe.cols[grep(x=pe.cols,pattern="expr")]
  pe.fold.cols <- pe.cols[grep(x=pe.cols,pattern="foldChange")]
  pe.count.cols <- pe.cols[grep(x=pe.cols,pattern="peakCount")]
  pe$peakCountSumRaw <- apply(pe[pe.count.cols],1,sum)
  pe$peakCountSum   <- apply(pe[c(pe.expr.cols,pe.count.cols)],1,
                             function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                                    x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.count.cols))],
                                                    sum))
  
  
  pe$exprCount <- apply(pe[pe.expr.cols],1,countNotNa)
  pe$PSIave   <- apply(pe[pe.expr.cols],1,function(x)mean(x,na.rm=TRUE))
  #pe$FoldChangeave   <- apply(pe[pe.fold.cols],1,function(x)mean(x,na.rm=TRUE))
  pe$FoldChangeave   <- apply(pe[c(pe.expr.cols,pe.fold.cols)],1,
                              function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                                     x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                     mean))
  pe$FoldChangeaveGroup <- cut(x=pe$FoldChangeave,breaks=seq(from=0,to=80,by=5))
  
  pe$FoldChangeVar  <- apply(pe[c(pe.expr.cols,pe.fold.cols)],1,
                             function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                                    x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                    var))
  pe$PSIVar   <- apply(pe[c(pe.expr.cols,pe.expr.cols)],1,
                       function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                              x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.expr.cols))],
                                              var))
  
  
  
  #pe$FoldChangeaveGroup <- 
  # fraction of expressed exons with peaks
  pe$peakFracExpr   <- apply(pe[c(pe.expr.cols,pe.count.cols)],1,
                             function(x)applyFnByNa(x[1:length(pe.expr.cols)],
                                                    x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.count.cols))],
                                                    fracGreaterZero))
  pe$cellsWithPeaks = with(pe,exprCount * peakFracExpr)
  ggplot(pe, aes(cellsWithPeaks))+geom_bar()+
    ggtitle("Number of expressed exons\nwith peaks")+
    theme_bw()
  
  
  PSI.mat <- as.matrix(pe[pe.expr.cols])
  foldChange.mat <- as.matrix(pe[pe.fold.cols])
  #> dim(foldChange.mat)
  #[1] 165352     15
  fc.vec <- as.vector(foldChange.mat)
  psi.vec <- as.vector(PSI.mat)
  
  
  pe$peakFracExprGroup   <- cut(x=pe$peakFracExpr,breaks=(1:length(pe.count.cols))/length(pe.count.cols))
  
  pe$pearsonCorr   <- apply(pe[c(pe.expr.cols,pe.fold.cols)],1,
                            function(x)applyCorrByNa(x[1:length(pe.expr.cols)],
                                                     x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                     method="pearson"))
  
  pe$pearsonCorrPvalue   <- apply(pe[c(pe.expr.cols,pe.fold.cols)],1,
                                  function(x)applyCorrByNaPvalue(x[1:length(pe.expr.cols)],
                                                                 x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                 method="pearson"))
  
  #wilcox.test(1:3,2:4,paired=TRUE)$p.value
  
  pe$wilcoxCorrPvalue   <- apply(pe[c(pe.expr.cols,pe.count.cols)],1,
                                 function(x)applyWilcoxTestByNaPvalue(x[1:length(pe.expr.cols)],
                                                                      x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.count.cols))]))
  
  
  
  pe.nona <- pe[!is.na(pe$pearsonCorrPvalue),]
  pe.pearson <- pe.nona[order(pe.nona$pearsonCorrPvalue,decreasing=FALSE),]
  pe.pearson$alpha <- 0.05
  pe.pearson$M <- length(pe.nona$pearsonCorrPvalue)
  pe.pearson$j <- seq_along(pe.nona$pearsonCorrPvalue) 
  FDRpass <- with(pe.pearson, which(pearsonCorrPvalue < alpha*(j/M)))
  pe.pearson.fdr <- pe.pearson[FDRpass,]
  write.table(pe.pearson.fdr,file=paste0(outdir,"FDRpass-peakFoldChange-vs-PSI.tab"),sep="\t")
  
  
  pe.nona6 <- pe[which(pe$exprCount >= 6),]
  pe.pearson6 <- pe.nona6[order(pe.nona6$pearsonCorrPvalue,decreasing=FALSE),]
  pe.pearson6$alpha <- 0.05
  pe.pearson6$M <- length(pe.pearson6$pearsonCorrPvalue)
  pe.pearson6$j <- seq_along(pe.pearson6$pearsonCorrPvalue) 
  FDRpass6 <- with(pe.pearson6, which(pearsonCorrPvalue < alpha*(j/M)))
  pe.pearson.fdr6 <- pe.pearson6[FDRpass6,]
  
  df= data.frame(psi= as.numeric(pe.pearson.fdr6[,pe.expr.cols]),
                 foldChange = as.numeric(pe.pearson.fdr6[,pe.fold.cols]))
  
  pw.nona <- pe[!is.na(pe$wilcoxCorrPvalue),]
  pw.wilcox <- pw.nona[order(pw.nona$wilcoxCorrPvalue,decreasing=FALSE),]
  pw.wilcox$alpha <- 0.05
  pw.wilcox$M <- length(pw.wilcox$wilcoxCorrPvalue)
  pw.wilcox$j <- seq_along(pw.wilcox$wilcoxCorrPvalue) 
  FDRpassW <- with(pw.wilcox, which(wilcoxCorrPvalue < alpha*(j/M)))
  pw.wilcox.fdr <- pw.wilcox[FDRpassW,]
  
  
  
  
  
  
  pe.peakShuffle <- shuffleRows(pe,pe.fold.cols)
  pe.peakShuffle$pearsonCorrPvalue   <- apply(pe.peakShuffle[c(pe.expr.cols,pe.fold.cols)],1,
                                              function(x)applyCorrByNaPvalue(x[1:length(pe.expr.cols)],
                                                                             x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                             method="pearson"))
  pe.peakShuffle$pearsonCorr   <- apply(pe.peakShuffle[c(pe.expr.cols,pe.fold.cols)],1,
                                        function(x)applyCorrByNa(x[1:length(pe.expr.cols)],
                                                                 x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                 method="pearson"))
  
  pes.nona <- pe.peakShuffle[!is.na(pe.peakShuffle$pearsonCorrPvalue),]
  pes.pearson <- pes.nona[order(pes.nona$pearsonCorrPvalue,decreasing=FALSE),]
  pes.pearson$alpha <- 0.05
  pes.pearson$M <- length(pes.pearson$pearsonCorrPvalue)
  pes.pearson$j <- seq_along(pes.pearson$pearsonCorrPvalue) 
  FDRpassShuffle <- with(pes.pearson, which(pearsonCorrPvalue < alpha*(j/M)))
  pes.pearson.fdr <- pes.pearson[FDRpassShuffle,]
  
  
  
  pe.peakShuffle.corr <- shuffleRows(pe.nona,pe.fold.cols)
  pe.peakShuffle.corr$pearsonCorrPvalue   <- apply(pe.peakShuffle.corr[c(pe.expr.cols,pe.fold.cols)],1,
                                                   function(x)applyCorrByNaPvalue(x[1:length(pe.expr.cols)],
                                                                                  x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                                  method="pearson"))
  pe.peakShuffle.corr$pearsonCorr   <- apply(pe.peakShuffle.corr[c(pe.expr.cols,pe.fold.cols)],1,
                                             function(x)applyCorrByNa(x[1:length(pe.expr.cols)],
                                                                      x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                      method="pearson"))
  
  
  
  
  pesc.nona <- pe.peakShuffle.corr[!is.na(pe.peakShuffle.corr$pearsonCorrPvalue),]
  pesc.pearson <- pesc.nona[order(pesc.nona$pearsonCorrPvalue,decreasing=FALSE),]
  pesc.pearson$alpha <- 0.05
  pesc.pearson$M <- length(pesc.pearson$pearsonCorrPvalue)
  pesc.pearson$j <- seq_along(pesc.pearson$pearsonCorrPvalue) 
  FDRpassShuffleCorr <- with(pesc.pearson, which(pearsonCorrPvalue < alpha*(j/M)))
  pesc.pearson.fdr <- pesc.pearson[FDRpassShuffleCorr,]
  
  
  
  
  peakExprTable <- dcast(as.data.frame(group_by(pe,cellsWithPeaks,exprCount) %.% summarise(count=length(label))), formula= cellsWithPeaks ~ exprCount,value.var="count",)
  plot(peakExprTable)
  
  # sanity :  sum(p.adjust(pes.nona$pearsonCorrPvalue, method="fdr", n=length(pes.nona$pearsonCorrPvalue)) < 0.05) == dim(pes.pearson.fdr)
  
  pe.short <- pe[1:10000,]
  
  pe.peakNa <- replaceRowZeroToNa(pe,pe.fold.cols) #getNumberOfNonNaBoth
  pe.peakNa$cellTypesBothValues   <- apply(pe.peakNa[c(pe.expr.cols,pe.fold.cols)],1,
                                           function(x)getNumberOfNonNaBoth(x[1:length(pe.expr.cols)],
                                                                           x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                           method="pearson"))
  
  pe.peakNa$pearsonCorrPvalue   <- apply(pe.peakNa[c(pe.expr.cols,pe.fold.cols)],1,
                                         function(x)applyCorrByNaPvalueBothCols(x[1:length(pe.expr.cols)],
                                                                                x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                                method="pearson"))
  pe.peakNa$pearsonCorr   <- apply(pe.peakNa[c(pe.expr.cols,pe.fold.cols)],1,
                                   function(x)applyCorrByNaColsBothCols(x[1:length(pe.expr.cols)],
                                                                        x[(length(pe.expr.cols)+1):(length(pe.expr.cols)+length(pe.fold.cols))],
                                                                        method="pearson"))
  pena.nona <- pe.peakNa[!is.na(pe.peakNa$pearsonCorrPvalue),]
  pena.pearson <- pena.nona[order(pena.nona$pearsonCorrPvalue,decreasing=FALSE),]
  pena.pearson$alpha <- 0.25
  pena.pearson$M <- length(pena.pearson$pearsonCorrPvalue)
  pena.pearson$j <- seq_along(pena.pearson$pearsonCorrPvalue) 
  pena.pearson$FDRpvalueCutoff <- with(pena.pearson,alpha*(j/M))
  FDRpassShuffleCorrPeakNa <- with(pena.pearson, which(pearsonCorrPvalue < alpha*(j/M)))
  pena.pearson.fdr <- pena.pearson[FDRpassShuffleCorrPeakNa,]
  pena.pearson$cutoff_0_25 <-  with(pena.pearson, alpha*(j/M))
  pena.pearson$cutoff_0_05 <-  with(pena.pearson, 0.05*(j/M))

  
  
  ggplot(as.data.frame(group_by(pe,cellsWithPeaks,exprCount) %.% summarise(count=length(label))), aes(x=cellsWithPeaks,y=exprCount,label=count,fill=count)) + 
    geom_text(size=2) + theme_bw()+
    ggtitle("Table of counts for number of cells \nwith peaks and expressed exons(PSI)")
  ggsave(paste0(outdir,"peak-vs-expr-table.png"), height=5,width=7)
  
  
  ggplot(pe, aes(cellsWithPeaks))+geom_bar(binwidth=1)+
    ggtitle("Number of cells w/ peaks\nfor expressed")+
    theme_bw()
  ggsave(paste0(outdir,"peakDistro-boxplot.png"), height=5,width=5)
  
  
  ggplot(pe[which(pe$cellsWithPeaks != 0),], aes(cellsWithPeaks))+geom_bar(binwidth=1)+
    ggtitle(paste("Number of cells w/ peaks\nfor expressed\nN=",
                  dim(pe[which(pe$cellsWithPeaks != 0),])[1]))+
    theme_bw()
  ggsave(paste0(outdir,"peakDistro-noZero-boxplot.png"), height=5,width=5)
  
  
  ggplot(pe[pe$FoldChangeave > 0,], aes(x=FoldChangeaveGroup,y=PSIave))+geom_boxplot() + 
    xlab("ave CTCF peaks fold change in all celltypes") + ylab("ave PSI in all cell type") + theme_bw()+
    ggtitle("average PSI for CTCF peaks found\nversus average fold change above bkgd")
  ggsave(paste0(outdir,"foldChangeGroup-avePSI-boxplot.png"), height=5,width=10)
  
  ggplot(pe, aes(x=FoldChangeave,y=PSIave,fill=peakCountSum))+geom_point() + geom_density2d() +
    xlab("ave CTCF peaks fold change") + ylab("ave PSI") + theme_bw()+
    ggtitle("average PSI for CTCF peaks found")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  ggsave(paste0(outdir,"foldChangeAve-avePSI-point.png"), height=7,width=7)
  
  cellExprFc <- data.frame(foldChange=fc.vec,PSI=psi.vec)
  ggplot(cellExprFc, aes(x=foldChange,y=PSI)) + geom_density2d() + 
    xlab("CTCF peaks fold change") + ylab("ave PSI") + theme_bw()+
    ggtitle("PSI of exon vs. downstream CTCF\nFold change")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  ggplot(cellExprFc, aes(x=PSI)) + geom_density() + 
    xlab("PSI of exon") + theme_bw() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  
  ggplot(cellExprFc, aes(x=foldChange)) + geom_density() + 
    xlab("CTCF peaks fold change") + theme_bw() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  
  #wide
  
  ggplot(pe[which(pe$peakCountSum %in% c(0,15)),], aes(PSIave,fill=factor(peakCountSum)))+geom_density(alpha=I(0.5)) + 
    xlab("average PSI") + ylab("density") + theme_bw()+
    ggtitle("average PSI\nCTCF occupied in all celltypes vs. none")
  ggsave(paste0(outdir,"PSI-density-colorByPeakCount.png"), height=7,width=7)
  
  
  
  ggplot(pe, aes(x=peakCountSum,y=PSIave,fill=factor(peakCountSum)))+geom_boxplot() + 
    xlab("CTCF peaks found") + ylab("density") + theme_bw()+
    ggtitle("average PSI for CTCF peaks found")
  ggsave(paste0(outdir,"PSI-peakCountSum-colorByPeakCount-boxplot.png"), height=7,width=10)
  
  
  ggplot(pe[which(pe$peakCountSum > 0),], aes(peakCountSum))+geom_bar(binwidth=1) + 
    scale_y_log10() + xlab("CTCF peaks") + ylab("count of exons") + 
    ggtitle("number of CTCF binding events detected\ndownstream of exon\nin all cell types")
  ggsave(paste0(outdir,"PSI-peakCountSum-bar.png"), height=7,width=10)
  
  ggplot(pe[which(pe$peakCountSum > 0),], aes(x=peakFracExpr,y=peakCountSum))+geom_point() + 
    xlab("fraction of expressed exons with CTCF peaks") + ylab("number downstream CTCF peaks") + 
    ggtitle("number of CTCF binding events detected\ndownstream of exon\nin all cell types") +
    theme_bw()
  ggsave(paste0(outdir,"peakCountSum-peakFrac-point.png"), height=7,width=7)
  
  ggplot(pe, aes(x=peakFracExpr,y=PSIave))+geom_point() + 
    xlab("fraction of expressed exons with CTCF peaks") + ylab("average PSI") + theme_bw()+
    ggtitle("frac CTCF peaks vs. ave PSI")
  ggsave(paste0(outdir,"peakFracExpr-PSIave-point.png"), height=7,width=7)
  
  
  ggplot(pe, aes(x=peakFracExprGroup,y=PSIave))+geom_boxplot() + 
    xlab("fraction of expressed exons with CTCF peaks") + ylab("mean PSI") + theme_bw()+
    ggtitle("fraction of expressed exons w/ CTCF peaks \nvs. mean PSI of expressed exons")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(outdir,"peakFracExpr-PSIave-boxplot.png"), height=7,width=12)
  
  
  ggplot(pe[which(!is.na(pe$PSIVar) & !is.na(pe$PSIVar)),], aes(x=log10(FoldChangeVar),y=PSIVar))+geom_point() + 
    xlab("variance of CTCF peak fold change") + ylab("variance of PSI") + theme_bw()+
    ggtitle("variance in CTCF foldChange vs. variance PSI\nin expressed exons only")
  ggsave(paste0(outdir,"foldChangeVar-PSIVar-point.png"), height=7,width=7)
  
  ggplot(pe, aes(x=-log10(pearsonCorrPvalue)))+geom_density() + 
    xlab("-log10(pearsonCorr P-value)") + ylab("density") + theme_bw()+
    ggtitle("Pearson Correlation P-value Between\n PSI and foldChange\nof expressed exons")
  ggsave(paste0(outdir,"PSI-CorrPvalue-density.png"), height=7,width=7)
  
  ggplot(pe.nona, aes(x=-log10(abs(pearsonCorr))))+geom_density() + 
    xlab("-log10(abs(pearsonCorr))") + ylab("density") + theme_bw()+
    ggtitle("Pearson Correlation Between\n PSI and foldChange\nof expressed exons")
  ggsave(paste0(outdir,"PSI-LogCorr-density.png"), height=7,width=7)
  
  ggplot(pe.nona, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle("Pearson Correlation Between\n PSI and foldChange\nof expressed exons")
  ggsave(paste0(outdir,"PSI-LogCorr-density.png"), height=7,width=7)
  
  ggplot(pe.nona, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("Pearson Correlation Coeff\nBetween PSI and foldChangeof expressed exons\n(in at least 3 cell types)\nN=",
                  dim(pe.nona)[1]))
  ggsave(paste0(outdir,"PSI-LogCorr-density.png"), height=7,width=7)
  
  ggplot(pe.pearson.fdr, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("Pearson Correlation Coeff\nBetween PSI and foldChangeof expressed exons\n(pass FDR of 0.05)\nN=",
                  dim(pe.pearson.fdr)[1]))
  ggsave(paste0(outdir,"PSI-LogCorr-FDR_pass-density.png"), height=7,width=7)
  
  ggplot(pes.nona, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("Shuffled downstream regions\nPearson Correlation Coeff\nBetween PSI and foldChangeof shuffled exons\n(in at least 3 cell types)\nN=",
                  dim(pes.nona)[1]))
  ggsave(paste0(outdir,"PSI-shuffle-LogCorr-gt3Expr-density.png"), height=7,width=7)
  
  ggplot(pes.pearson.fdr, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("shuffled downstream regions\npearson correlation coeff\nbetween PSI and foldChangeof shuffled exons\n(pass FDR of 0.05)\nN=",
                  dim(pes.pearson.fdr)[1]))
  ggsave(paste0(outdir,"PSI-shuffle-LogCorr-FDR_pass-density.png"), height=7,width=7)
  
  
  ggplot(pesc.nona, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("shuffled downstream regions(regions w/ defined corr only)\npearson correlation coeff\nbetween PSI and foldChangeof shuffled exons\n(in at least 3 cell types)\nN=",
                  dim(pesc.nona)[1]))
  ggsave(paste0(outdir,"PSI-shuffleCorrOnly-LogCorr-gt3Expr-density.png"), height=7,width=7)
  
  ggplot(pesc.pearson.fdr, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("shuffled downstream regions(regions w/ defined corr only)\npearson correlation coeff\nbetween PSI and foldChangeof shuffled exons\n(pass FDR of 0.05)\nN=",
                  dim(pesc.pearson.fdr)[1]))+xlim(-1,1)
  ggsave(paste0(outdir,"PSI-shuffleCorrOnly-LogCorr-FDR_pass-density.png"), height=7,width=7)
  
  
  ggplot(pena.pearson, aes(x=pearsonCorr))+geom_density() + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("downstream regions(regions w/ defined correlation only)\npearson correlation coeff\nbetween PSI and foldChangeof  exons\n(if no peak in found  -> foldChange for exon in cell = NA)\nN=",
                  dim(pena.pearson)[1]))
  ggsave(paste0(outdir,"PSI-noPeakNa-LogCorr-gt3Expr-density.png"), height=7,width=7)
  
  ggplot(pena.pearson, aes(x=pearsonCorr))+geom_bar(binwidth=0.05) + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("downstream regions(regions w/ defined correlation only)\npearson correlation coeff\nbetween PSI and foldChangeof  exons\n(if no peak in found  -> foldChange for exon in cell = NA)\nN=",
                  dim(pena.pearson)[1]))
  ggsave(paste0(outdir,"PSI-noPeakNa-LogCorr-gt3Expr-bars.png"), height=7,width=7)
  
  ggplot(pena.pearson.fdr, aes(x=pearsonCorr))+geom_bar(binwidth=0.05) + 
    xlab("pearsonCorr") + ylab("density") + theme_bw()+
    ggtitle(paste("downstream regions(regions w/ defined correlation only)\npearson correlation coeff\nbetween PSI and foldChange of  exons\n(if no peak in found  -> foldChange for exon in cell = NA)\nN=",
                  dim(pena.pearson.fdr)[1]))
  ggsave(paste0(outdir,"PSI-noPeakNa-LogCorr-gt3Expr-fdrPass-bars.png"), height=7,width=7)
  
  #cellTypesBothValues
  ggplot(pena.pearson, aes(x=cellTypesBothValues))+geom_bar(binwidth=0.5) + 
    xlab("# of cells types w/ (PSI and fold change > 0)") + ylab("number of exons found") + theme_bw()+
    ggtitle(paste("downstream regions(regions w/ defined correlation only)\ncell types with exon experssion\nand MACS peak foldChange > 0\n(if no peak in cell  -> foldChange in cell = NA)\nN=",
                  dim(pena.pearson)[1]))+
    xlim(0,max(pena.pearson$cellTypesBothValues)+1)
  ggsave(paste0(outdir,"PSI-noPeakNa-cellTypesFound-bars.png"), height=7,width=7)
  
  pena.melt <- melt(pena.pearson[c("j","pearsonCorrPvalue","cutoff_0_25","cutoff_0_05")],id.var="j")
  
  ggplot(pena.pearson, aes(x=j,y=pearsonCorrPvalue)) + geom_point() +
    scale_y_log10() + scale_x_log10()+
    annotation_logticks() + theme_bw()+
    xlab("exons ordered by p-value")+
    ylab("p-value")+
    ggtitle(paste("p-value vs. FDR cutoff @ 0.25(blue) & 0.05(red)\ndownstream regions(regions w/ defined correlation only)\ncell types with exon experssion\nand MACS peak foldChange > 0\n(if no peak in cell  -> foldChange in cell = NA)\nN =",
                  dim(pena.pearson)[1]))+
    annotate("text",x=2,y=0.00004,label="SEMA-3B",color="black",size=4)+
    annotate("line",x=pena.pearson$j,y= pena.pearson$cutoff_0_05,color="red")+
    annotate("line",x=pena.pearson$j,y= pena.pearson$cutoff_0_25,color="blue")+
    annotate("point",x=pena.pearson$j[2],y= pena.pearson$pearsonCorrPvalue[2],color="yellow",alpha=I(0.5),size=I(10))
  ggsave(paste0(outdir,"PSI-noPeakNa-pValVsCutoff-line.png"), height=7,width=7)
  
}



