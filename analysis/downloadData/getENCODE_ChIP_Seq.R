homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}
projectDir <- normalizePath(file.path(script.dir(), ".."))
getFullPath <- function(subpath){ file.path(projectDir, subpath) }


jan11FreezeUrl <<- "http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/files.txt"
jan11FreezeLocal <<- getFullPath("data/integration_data_jan2011-files.txt")
jan11FreezeLocalTab <<- getFullPath("data/integration_data_jan2011-files.tab")


cshlRNAseqUrl <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/files.txt"
cshlRNAseqLocal <<- getFullPath("data/wgEncodeCshlLongRnaSeq-files.txt")
cshlRNAseqLocalTab <<- getFullPath("data/wgEncodeCshlLongRnaSeq-files.tab")

sydhTfbsUrl <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/files.txt"
sydhTfbsLocal <<- getFullPath("data/wgEncodeSydhTfbs-files.txt")
sydhTfbsLocalTab <<- getFullPath("data/wgEncodeSydhTfbs-files.tab")



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
    df[[ j ]] <- NA
  }
  
  df
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


getCtcfWithRNASeq <- function(){
  rnaSeq <- read.csv(file=cshlRNAseqLocalTab,sep="\t",stringsAsFactor=FALSE)
  
  
}


sydhTfbsLocalTab

