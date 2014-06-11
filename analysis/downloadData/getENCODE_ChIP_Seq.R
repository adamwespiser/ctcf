homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}
projectDir <- normalizePath(file.path(script.dir(), ".."))
getFullPath <- function(subpath){ file.path(projectDir, subpath) }

clusterDir <<- "/project/umw_zhiping_weng/wespisea/ctcf/data/"

jan11FreezeUrl <<- "http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/files.txt"
jan11FreezeLocal <<- getFullPath("data/integration_data_jan2011-files.txt")
jan11FreezeLocalTab <<- getFullPath("data/integration_data_jan2011-files.tab")


cshlRNAseqUrlBase <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/"
cshlRNAseqUrl <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/files.txt"
cshlRNAseqLocal <<- getFullPath("data/wgEncodeCshlLongRnaSeq-files.txt")
cshlRNAseqLocalTab <<- getFullPath("data/wgEncodeCshlLongRnaSeq-files.tab")
cshl.list <- list(base=cshlRNAseqUrlBase,url=cshlRNAseqUrl,local=cshlRNAseqLocal,localTab=cshlRNAseqLocalTab)

sydhTfbsUrl <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/files.txt"
sydhTfbsUrlBase <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/"
sydhTfbsLocal <<- getFullPath("data/wgEncodeSydhTfbs-files.txt")
sydhTfbsLocalTab <<- getFullPath("data/wgEncodeSydhTfbs-files.tab")
sydh.list <- list(base=sydhTfbsUrlBase,url=sydhTfbsUrl,local=sydhTfbsLocal,localTab=sydhTfbsLocalTab)

haibTfbsUrl <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/files.txt"
haibTfbsUrlBase <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/"
haibTfbsLocal <<- getFullPath("data/wgEncodeHaibTfbs-files.txt")
haibTfbsLocalTab <<- getFullPath("data/wgEncodeHaibTfbs-files.tab")
haib.list <- list(base=haibTfbsUrlBase,url=haibTfbsUrl,local=haibTfbsLocal,localTab=haibTfbsLocalTab)

UchicagoTfbsUrl <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUchicagoTfbs/files.txt"
UchicagoTfbsUrlBase <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUchicagoTfbs/"
UchicagoTfbsLocal <<- getFullPath("data/wgEncodeUchicagoTfbs-files.txt")
UchicagoTfbsLocalTab <<- getFullPath("data/wgEncodeUchicagoTfbs-files.tab")
Uchicago.list <- list(base=UchicagoTfbsUrlBase,url=UchicagoTfbsUrl,local=UchicagoTfbsLocal,localTab=UchicagoTfbsLocalTab)

uwTfbsUrl <<-   "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/files.txt"
uwTfbsUrlBase <<-   "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/"
uwTfbsLocal <<- getFullPath("data/wgEncodeUwTfbs-files.txt")
uwTfbsLocalTab <<- getFullPath("data/wgEncodeUwTfbs-files.tab")
uw.list <- list(base=uwTfbsUrlBase,url=uwTfbsUrl,local=uwTfbsLocal,localTab=uwTfbsLocalTab)





convertToDfWithCols <- function(str,recSep,labelSep,valRm="",cols){
  records <- unique(as.character(unlist(strsplit(str, recSep))))
  nms <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[1])))
  vals <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[2])))
  # vals <- sapply(vals, function(x)gsub(x,";",""))
  l <- as.list(vals)
  names(l) <- nms 
  df <- as.data.frame(l)
  cols.notfound <- cols[-which(cols %in% nms)]
  for ( i in cols.notfound){
    df[[ i ]] <- NA
  }
  #remove extra cols 
  df.cols <- colnames(df)
  for (j in df.cols[-which(df.cols %in% cols)]){
    df[[ j ]] <- NULL
  }
  
  df
}


joinUnevenDfs <- function(df1,df2){
  stopifnot(is.data.frame(df1),is.data.frame(df2))
  
  col1 <-   colnames(df1)
  col2 <- colnames(df2)
  all.cols <- unique(c(col1,col2))
  
  if(TRUE == all.equal(col1,col2)[1]){
    return(rbind(df1[all.cols],df2[all.cols]))
  }
  
  df1.col.need <- all.cols[-which(all.cols %in% colnames(df1) )]
  df2.col.need <- all.cols[-which(all.cols %in% colnames(df2) )]
  
  if(length(df1.col.need) > 0){
    for(i in df1.col.need){
      df1[[ i ]] <- NA
    }
  }
  
  if(length(df2.col.need) > 0){
    for(j in df2.col.need){
      df2[[ j ]] <- NA
    }
  }
  
  o.df <- rbind(df1[all.cols],df2[all.cols])
  
  colo <- colnames(o.df)
  o.df.dim1 <- dim(x= o.df )[1]
  for(k in colo){
    if(length(which(is.na(o.df[[ k ]]))) == o.df.dim1 ){
      o.df[[ k ]] <- NULL
    }
  }
  o.df
  
}




convertToList <- function(str,recSep,labelSep,valRm=""){
  records <- as.character(unlist(strsplit(str, recSep)))
  nms <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[1])))
  vals <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[2])))
  # vals <- sapply(vals, function(x)gsub(x,";",""))
  l <- as.list(vals)
  names(l) <- nms 
  l}

shortenIdVec <- function(gene_id, sep = "\\."){
  as.character(sapply(gene_id, function(x)unlist(strsplit(x,split=sep))[1]))
}




fetchEncodeFilesTxt <- function(urlAddr,localAddr,localTabAddr){
  if (!file.exists(localAddr)){
    download.file(urlAddr,dest=localAddr)
  }
  
  d <- readLines(localAddr)
  if(length(grep(d,pattern="#"))){
    d <- d[-grep(d,pattern="#")]
  }
  
  d.split <- lapply(d, function(x)strsplit(x,split="\t"))
  d.names <- sapply(d.split, function(x)unlist(x)[1])
  d.info <- sapply(d.split, function(x)unlist(x)[2])
  d.list <- lapply(d.info,function(x)convertToList(x,"; ", "=",";"))
  colAnnot <-  unique(do.call("c",lapply(d.list, function(x)names(x))))
  df <- do.call(rbind,lapply(d.info,function(x)convertToDfWithCols(x,"; ", "=",";",cols=colAnnot)))
  df$filename <- d.names
  
  exportAsTable(df=df,file=localTabAddr)
}


getCTCFchipAllFilesCellTypes <- function(){
  allExp <- read.csv(file=jan11FreezeLocalTab,sep="\t",stringsAsFactor=FALSE)
  allExp.chSeq <- allExp[which(allExp$dataType %in% c("ChipSeq", "ChIP-seq")),]
  allExp.all.celltypes <- unique(allExp.chSeq$antibody)
  ctcf.antibody <- allExp.all.celltypes[grep(allExp.all.celltypes,pattern="CTCF",ignore.case=TRUE)]
  if(grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE)){
    ctcf.antibody <- ctcf.antibody[-grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE)]
  }
  
  allExp.ctcf <- allExp[which(allExp$antibody %in% ctcf.antibody),]
  unique(allExp.ctcf$cell)
  
}


getCTCFPeakSeq <- function(){
  allExp <- read.csv(file=jan11FreezeLocalTab,sep="\t",stringsAsFactor=FALSE)
  allExp.chSeq <- allExp[which(allExp$dataType %in% c("ChipSeq", "ChIP-seq")),]
  allExp.all.celltypes <- unique(allExp.chSeq$antibody)
  ctcf.antibody <- allExp.all.celltypes[grep(allExp.all.celltypes,pattern="CTCF",ignore.case=TRUE)]
  if(grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE)){
    ctcf.antibody <- ctcf.antibody[-grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE)]
  }
  
  allExp.ctcf <- allExp[which(allExp$antibody %in% ctcf.antibody),]
  allExp[which(allExp.ctcf$softwareVersion == "PeakSeq"),]
}

getCTCFSPP <- function(){
  allExp <- read.csv(file=jan11FreezeLocalTab,sep="\t",stringsAsFactor=FALSE)
  allExp.chSeq <- allExp[which(allExp$dataType %in% c("ChipSeq", "ChIP-seq")),]
  allExp.all.celltypes <- unique(allExp.chSeq$antibody)
  ctcf.antibody <- allExp.all.celltypes[grep(allExp.all.celltypes,pattern="CTCF",ignore.case=TRUE)]
  if(grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE)){
    ctcf.antibody <- ctcf.antibody[-grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE)]
  }
  
  allExp.ctcf <- allExp[which(allExp$antibody %in% ctcf.antibody),]
  allExp[which(allExp.ctcf$softwareVersion == "SPP"),]
}

getCTCFchipAllFiles <- function(){
  allExp <- read.csv(file=jan11FreezeLocalTab,sep="\t",stringsAsFactor=FALSE)
  allExp.chSeq <- allExp[which(allExp$dataType %in% c("ChipSeq", "ChIP-seq")),]
  allExp.all.celltypes <- unique(allExp.chSeq$antibody)
  ctcf.antibody <- allExp.all.celltypes[grep(allExp.all.celltypes,pattern="CTCF",ignore.case=TRUE)]
  if(grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE)){
    ctcf.antibody <- ctcf.antibody[-grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE)]
  }
  
  allExp[which(allExp$antibody %in% ctcf.antibody),]
  
  
}


getCTCF_tfbs<- function(labList){
  chSeq <- read.csv(file=labList$localTab,sep="\t",stringsAsFactor=FALSE)
  chSeq$filename <- paste0(labList$base,chSeq$filename)
  chSeq.ab <- unique(chSeq$antibody)
  ctcf.antibody <- chSeq.ab[grep(chSeq.ab,pattern="CTCF",ignore.case=TRUE)]
  if(0 == length(ctcf.antibody)){
    return(data.frame())
  }
  if(length(grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE) > 0)){
    ctcf.antibody <- ctcf.antibody[-grep(ctcf.antibody,pattern="CTCFL",ignore.case=TRUE)]
  }
  
  chSeq[which(chSeq$antibody %in% ctcf.antibody),]
}


getInput_tfbs<- function(labList){
  chSeq <- read.csv(file=labList$localTab,sep="\t",stringsAsFactor=FALSE)
  chSeq$filename <- paste0(labList$base,chSeq$filename)
  
  chSeq[which(chSeq$setType == "input"),]
}



getTfbsUniqueCells <- function(df){
  stopifnot("cell" %in% colnames(df))
  if(dim(df)[1] == 0){
    return(c())
  }
  unique(df[["cell"]])
}

getAllTfbs <- function(){
  uw <- getCTCF_tfbs(uw.list)
  haib <- getCTCF_tfbs(haib.list)
  sydh <- getCTCF_tfbs(sydh.list)
  comb <- joinUnevenDfs(joinUnevenDfs(uw,haib),sydh) 
}

getAllTfbsInput <- function(){
  uw <- getInput_tfbs(uw.list)
  haib <- getInput_tfbs(haib.list)
  sydh <- getInput_tfbs(sydh.list)
  comb <- joinUnevenDfs(joinUnevenDfs(uw,haib),sydh) 
  
}




getCtcfWithRNASeq <- function(){
   
  inputExp.all <- getAllTfbsInput()
  inputExp.all$fileEnding <- sapply(inputExp.all$filename,getLastUrl)
  
  inputExp <- inputExp.all[which(inputExp.all$type %in% c("bai","bam")),]
  
  
  tfbs.all <- getAllTfbs()
  tfbs.all$fileEnding <- sapply(tfbs.all$filename,getLastUrl)
  
  tfbs.ab <- tfbs.all[which(tfbs.all$type %in% c("bai","bam")),]
  tfbs <- merge(tfbs.ab,inputExp, by=c("project","grant","lab","cell","type","replicate","protocol"),suffix=c("",".cntrl"))
  
  tfbs.cells <- getTfbsUniqueCells(tfbs)
   
  rnaSeq <- read.csv(file=cshlRNAseqLocalTab,sep="\t",stringsAsFactor=FALSE)
  rnaSeq <- rnaSeq[which(rnaSeq$type == "fastq"),]
  rnaSeq$fileEnding <- rnaSeq$filename
  rnaSeq$filename <- paste0(cshlRNAseqUrlBase,rnaSeq$filename)
  rnaSeq.lpa.wholeCell <- rnaSeq[which(rnaSeq$localization == "cell" & rnaSeq$rnaExtract == "longPolyA"),]
  rna.cells <- unique(rnaSeq.lpa.wholeCell$cell)
  
  cellTypes <- intersect(rna.cells,tfbs.cells)
  
  rnaSeq.peaks.df <-joinUnevenDfs(tfbs[which(tfbs$cell %in% cellTypes),],
                                  rnaSeq.lpa.wholeCell[which(rnaSeq.lpa.wholeCell$cell %in% cellTypes),])  
  rnaSeq.peaks.df <- rnaSeq.peaks.df[which(rnaSeq.peaks.df$treatment %in% c("None",NA)),]
  rnaSeq.peaks.df <- rnaSeq.peaks.df[which(rnaSeq.peaks.df$objStatus != "replaced - object is identifying wrong peaks" |
                                             is.na(rnaSeq.peaks.df$objStatus)),]
  rnaSeq.peaks.df <- rnaSeq.peaks.df[-which(rnaSeq.peaks.df$controlId.cntrl %in% c("wgEncodeEH000706","wgEncodeEH000771")),]
  
  
  
  
  
  rnaSeq.peaks.df$localFile <- paste0(clusterDir,rnaSeq.peaks.df$fileEnding)
  rnaSeq.peaks.df$downloadCmd <- paste0("wget --continue ",rnaSeq.peaks.df$filename,
                                        " -O ",rnaSeq.peaks.df$localFile)
  
  rnaSeq.Cntr.idx <- !is.na(rnaSeq.peaks.df$fileEnding.cntrl)
  rnaSeq.peaks.df$localFile.cntrl <- paste0(clusterDir,rnaSeq.peaks.df$fileEnding.cntrl)
  
  cntrlDownload <-  paste0("wget --continue ",rnaSeq.peaks.df$filename.cntrl,
                           " -O ",rnaSeq.peaks.df$localFile.cntrl)[rnaSeq.Cntr.idx]
  
  totalDownload <- c(cntrlDownload,rnaSeq.peaks.df$downloadCmd)
  
  writeLines(con="~/sandbox/downloadCTCF_bamBaiControl",
             text=paste(totalDownload,rep(c("&",""),length.out=length(totalDownload))))
  # scpFile(file.local="~/sandbox/downloadCTCF_bamBaiControl", dir.remote="~/bin/") 
  
  
  tfbs.peaks.full <- rnaSeq.peaks.df[rnaSeq.Cntr.idx,]

  tfbs.peaks <- tfbs.peaks.full[-grep(x=tfbs.peaks.full$localFile, pattern="bai"),]
  tfbs.peaks$cntrlName <- gsub(x=gsub(x=tfbs.peaks$fileEnding.cntrl,pattern="wgEncode",replacement=""),
                                    pattern=".bam",replacement="")
  tfbs.peaks$expName <- gsub(x=gsub(x=tfbs.peaks$fileEnding,pattern="wgEncode",replacement=""),
                               pattern=".bam",replacement="")
  tfbs.peaks$name <- with(tfbs.peaks,paste0(expName,"_",cntrlName))
  macs2dir <- "/project/umw_zhiping_weng/wespisea/ctcf/macs2/"
  tfbs.peaks$nameFile <- with(tfbs.peaks,paste0(macs2dir,name))
  write(file="~/sandbox/macs2_ctcf_manifest", tfbs.peaks$nameFile)
  # scpFile(file.local="~/sandbox/macs2_ctcf", dir.remote=macs2dir) 
  
  # /home/aw30w/bin/MACS-master/bin/macs2 callpeak -t ChIP.bam -c Control.bam -f BAM -g hs -n test -B -q 0.01
  # https://github.com/taoliu/MACS/
  macs2.cmd <-  paste0("python /home/aw30w/bin/MACS-master/bin/macs2 callpeak -t ",tfbs.peaks$localFile,
                           " -c ",tfbs.peaks$localFile.cntrl,
                           " -f BAM -g hs -B -q 0.01 ",
                           " -n ",tfbs.peaks$nameFile)
  write(file="~/sandbox/macs2_ctcf", macs2.cmd)
  
  
  
  
  # scpFile(file.local="~/sandbox/macs2_ctcf", dir.remote="~/bin/") 
  # perl /home/aw30w/bin/runTask.pl -f ~/bin/macs2_ctcf -c 24 -m 2024 -W 720 -Q short -t ctcfMacs2

  localDir <- "/home/wespisea/data/ctcf/MACS2/"
  
  if(TRUE == getFile){
    dir.create(localDir,recursive=TRUE)
    o1 <- paste0("scp aw30w@ghpcc06.umassrc.org:",tfbs.peaks$nameFile, "*{narrowPeak,pdf}  ",localDir)
    write(o1,file="~/sandbox/ctcfMacsFetch")
    system("chmod u+x ~/sandbox/ctcfMacsFetch")
    suppressWarnings(system("~/sandbox/ctcfMacsFetch"))
    
  }
  tfbs.peaks$narrowPeak <- with(tfbs.peaks,paste0(localDir,"/",name,"_peaks.narrowPeak"))
  

}

getLastUrl <- function(x){
  splitUrl <- as.character(unlist(strsplit(x,"/")))
  splitUrl[length(splitUrl)]
}


getCtcfPeakSeqWithRNASeq <- function(){
  
  tfbs <- getCTCFSPP()
  tfbs$fileEnding <- sapply(tfbs$filename,getLastUrl)
  tfbs.cells <- getTfbsUniqueCells(tfbs)
  
  rnaSeq <- read.csv(file=cshlRNAseqLocalTab,sep="\t",stringsAsFactor=FALSE)
  rnaSeq <- rnaSeq[which(rnaSeq$type == "fastq"),]
  rnaSeq$fileEnding <- rnaSeq$filename
  rnaSeq$filename <- paste0(cshlRNAseqUrlBase,rnaSeq$filename)
  rnaSeq.lpa.wholeCell <- rnaSeq[which(rnaSeq$localization == "cell" & rnaSeq$rnaExtract == "longPolyA"),]
  rna.cells <- unique(rnaSeq.lpa.wholeCell$cell)
  
  cellTypes <- intersect(rna.cells,tfbs.cells)
  
  rnaSeq.peaks.df <-joinUnevenDfs(tfbs[which(tfbs$cell %in% cellTypes),],
                                  rnaSeq.lpa.wholeCell[which(rnaSeq.lpa.wholeCell$cell %in% cellTypes),])  
  
  rnaSeq.peaks.df$localFile <- paste0(clusterDir,rnaSeq.peaks.df$fileEnding)
  rnaSeq.peaks.df$downloadCmd <- paste0("wget --continue ",rnaSeq.peaks.df$filename," -O ",rnaSeq.peaks.df$localFile)
  
  writeLines(con="~/sandbox/downloadCTCF",
            text=paste(rnaSeq.peaks.df$downloadCmd,rep(c("&",""),times=length(rnaSeq.peaks.df$downloadCmd))))
  
  # scpFile(file.local="~/sandbox/downloadCTCF", dir.remote="~/bin/") 
  
  
  
  #/home/wespisea/bin/bedtools2/bin/intersectBed -wo -a /home/wespisea/work/research/researchProjects/coexpr/lncNET/data/exonCosi_dsRegion.bed6 -b UwTfbsNhlfCtcfStdAlnRep1_UwTfbsNhlfInputStdAlnRep1_peaks.narrowPeak > ~/sandbox/bedBug
  
  
  
  
}



getCosiFileExons <- funciton(){
  dataDir <- getFullPath("data")
  trg <- paste0(trg,"/tx_paper.cosi.tsv")
  download.file(url="http://genome.crg.es/~dmitri/export/2014-06-02/tx_paper.cosi.tsv",
                destfile=trg)
}

getCosiFileCellTypes <- funciton(){
  dataDir <- getFullPath("data")
  trg <- paste0(dataDir,"/tx_paper.cosi.tsv")
  finalOut <- paste0(dataDir,"/exonCosi_dsRegion.tab")
  dsFinalOut <- paste0(dataDir,"/exonCosi_dsRegion.bed6")
  
  
  exonBedFile <- paste0(dataDir,"/tx_paper.exon.bed6")
  exonDownstreamBedFile <- paste0(dataDir,"/tx_paper.exonDownStream.bed6")
  exonBedFileSort <- paste0(dataDir,"/tx_paper.exon.sort.bed6")
  exonDownstreamBedFileSort <- paste0(dataDir,"/tx_paper.exonDownStream.sort.bed6")
  exonDownstreamBedFileSortNoExonOverlap <- paste0(dataDir,"/tx_paper.exonDownStream.noExOverlap.sort.bed6")
  
  # 1000 bp downstream, stop at next starting exon 
  exonDownstreamBedFileSortNoExonOverlapFirst <- paste0(dataDir,"/tx_paper.exonDownStream.noExOverlap.first.sort.bed6")
  
  
  trg.df <- read.csv(file=trg,stringsAsFactors=FALSE,sep="\t")
  cols <- colnames(trg.df)
  rw <-strsplit(rownames(trg.df),"_")
  trg.df$chr <- sapply(rw, function(x)x[1])  
  # format(moDat2, scientific=FALSE)
  trg.df$startPos <- sapply(rw, function(x)as.numeric(x[2])) 
  trg.df$stopPos <- sapply(rw, function(x)as.numeric(x[3]))
  trg.df$strand <- sapply(rw, function(x)x[4])
  trg.df$label <- with(trg.df,paste0(chr,":",format(startPos, scientific=FALSE,collapse=TRUE),"-",format(stopPos, scientific=FALSE),"_",strand,sep=""))
  trg.df$label <- gsub(trg.df$label, pattern=" ",replacement="")
  trg.df$zero <- 0
  # chr5 100000 500000 read1 50 +
  write.table(format(trg.df[c("chr","startPos","stopPos","label","zero","strand")], scientific=FALSE),file=exonBedFile,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)  
  
  trg.df$dsStart <- ifelse(trg.df$strand == "+", trg.df$stopPos, trg.df$startPos - 1000)
  trg.df$dsStop <- ifelse(trg.df$strand == "+", trg.df$stopPos + 1000, trg.df$startPos)
  write.table(format(trg.df[c("chr","dsStart","dsStop","label","zero","strand")], scientific=FALSE),
    file=exonDownstreamBedFile,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  
  cmd1 <- paste("/home/wespisea/bin/bedtools2/bin/sortBed",
         " -i ", exonBedFile,
          " > ",exonBedFileSort) 
  
  cmd2 <- paste("/home/wespisea/bin/bedtools2/bin/sortBed",
        " -i ", exonDownstreamBedFile,
        " > ", exonDownstreamBedFileSort) 
  
  
  cmd3 <- paste("/home/wespisea/bin/bedtools2/bin/subtractBed -s ",
        " -a ", exonDownstreamBedFileSort, 
        " -b ", exonBedFileSort,
        " > ",exonDownstreamBedFileSortNoExonOverlap)
  
  system(cmd1);system(cmd2);system(cmd3)
  
  ds.df <- read.csv(file=exonDownstreamBedFileSortNoExonOverlap, stringsAsFactors=FALSE,
                    sep="\t",header=FALSE)
  colnames(ds.df) <- c("chr","startPos","stopPos","label","zero","strand")
  ds.plus.df <- ds.df[which(ds.df$strand == "+"),]
  ds.minus.df <- ds.df[which(ds.df$strand == "-"),]
  
  ds.plus.first <- as.data.frame(group_by(ds.plus.df, label) %.% 
                                   mutate(minStartPos = min(startPos)) %.%
                                   filter(startPos == minStartPos))
  ds.plus.first$minStartPos <- NULL
  
  ds.minus.first <- as.data.frame(group_by(ds.minus.df, label) %.% 
                                   mutate(maxStartPos = max(startPos)) %.%
                                   filter(startPos == maxStartPos))
  ds.minus.first$maxStartPos <- NULL
  
  ds.first.comb <- rbind(ds.plus.first,ds.minus.first)
  tmp <- paste0(tempfile(),".bed6")
  write.table(ds.first.comb,
              file=tmp,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  cmd4 <- paste("/home/wespisea/bin/bedtools2/bin/sortBed",
                " -i ", tmp,
                " > ", exonDownstreamBedFileSortNoExonOverlapFirst) 
  system(cmd4)
  
  ggplot(ds.first.comb, aes(stopPos - startPos)) + geom_bar(binwith=10)
  
  rnaSeq <- read.csv(file=cshlRNAseqLocalTab,sep="\t",stringsAsFactor=FALSE)
  rnaSeq <- rnaSeq[which(rnaSeq$type == "fastq"),]
  rnaSeq$fileEnding <- rnaSeq$filename
  rnaSeq$filename <- paste0(cshlRNAseqUrlBase,rnaSeq$filename)
  rnaSeq$exp <- with(rnaSeq,paste(cell,localization,rnaExtract,replicate,sep="_"))
  rnaSeqLID <- unique(rnaSeq[c("labExpId","exp","cell")])
  
  cols.split <- strsplit(cols, split="_")
  cols.one <- sapply(cols.split, function(x)x[1])  
  cols.two <- sapply(cols.split, function(x)x[2])  
  
  cols.one.exp <- sapply(cols.one, function(x)rnaSeqLID[which(rnaSeqLID$labExpId == x),"exp"])
  cols.two.exp <- sapply(cols.two, function(x)rnaSeqLID[which(rnaSeqLID$labExpId == x),"exp"])
  
  cols.one.cell <- sapply(cols.one, function(x)rnaSeqLID[which(rnaSeqLID$labExpId == x),"cell"])
  cols.two.cell <- sapply(cols.two, function(x)rnaSeqLID[which(rnaSeqLID$labExpId == x),"cell"])
  
  newCols <- gsub(x=gsub(x=as.character(unlist(cols.one.cell)),pattern="-",replacement=""),pattern="\\+",replacement="")
  
  oldCols <- colnames(trg.df)
  oldCols[1:length(newCols)] <- newCols
  colnames(trg.df) <- oldCols
  
  
  paste(cols.one.exp,cols.two.exp,sep="--")
  
  
  comb <- merge(y=ds.first.comb,x=trg.df, by =  c("label","zero","strand"),suffixes=c("",".ds"))
  write.table(comb,file=finalOut,sep="\t")
  
  t2 <- tempfile()
  write.table(comb[c("chr.ds","dsStart","dsStop","label","zero","strand")],
              file=t2,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  cmd5 <- paste("/home/wespisea/bin/bedtools2/bin/sortBed",
                " -i ", t2,
                " > ", dsFinalOut) 
  
}






