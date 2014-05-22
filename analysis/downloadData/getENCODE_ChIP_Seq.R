homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}
projectDir <- normalizePath(file.path(script.dir(), ".."))
getFullPath <- function(subpath){ file.path(projectDir, subpath) }



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






getCtcfWithRNASeq <- function(){
   
  tfbs <- getAllTfbs()
  tfbs.cells <- getTfbsUniqueCells(tfbs)
   
  rnaSeq <- read.csv(file=cshlRNAseqLocalTab,sep="\t",stringsAsFactor=FALSE)
  rnaSeq$filename <- paste0(cshlRNAseqUrlBase,rnaSeq$filename)
  rnaSeq.lpa.wholeCell <- rnaSeq[which(rnaSeq$localization == "cell" & rnaSeq$rnaExtract == "longPolyA"),]
  rna.cells <- unique(rnaSeq.lpa.wholeCell$cell)
  
  cellTypes <- intersect(rna.cells,tfbs.cells)
    
  rnaSeq.peaks.df <-joinUnevenDfs(tfbs[which(tfbs$cell %in% cellTypes),],
                                  rnaSeq.lpa.wholeCell[which(rnaSeq.lpa.wholeCell$cell %in% cellTypes),])  
  
}


sydhTfbsLocalTab

