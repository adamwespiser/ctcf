homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}
projectDir <- normalizePath(file.path(script.dir(), ".."))
getFullPath <- function(subpath){ file.path(projectDir, subpath) }




convertToDfWithCols <- function(str,recSep,labelSep,valRm="",cols){
	  records <- as.character(unlist(strsplit(str, recSep)))
  nms <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[1])))
    vals <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[2])))
    # vals <- sapply(vals, function(x)gsub(x,";",""))
    l <- as.list(vals)
	  names(l) <- nms 
	  df <- as.data.frame(l)
	    cols.notfound <- cols[-which(cols %in% nms)]
	    for ( i in cols.notfound){
			    df[[i ]] <- NA
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








fetchEncodeDccFilesText <- function(filesTxt="~/data/wgEncodeCshlLongRnaSeqFiles.txt",
									                                    filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab"){
	  if (!file.exists(filesTxt)){
		      download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/files.txt",dest=filesTxt)
	    }
  
  d <- readLines(filesTxt)
    d.split <- lapply(d, function(x)strsplit(x,split="\t"))
    d.names <- sapply(d.split, function(x)unlist(x)[1])
	  d.info <- sapply(d.split, function(x)unlist(x)[2])
	  d.list <- lapply(d.info,function(x)convertToList(x,"; ", "=",";"))
	    colAnnot <-  unique(do.call("c",lapply(d.list, function(x)names(x))))
	    df <- do.call(rbind,lapply(d.info,function(x)convertToDfWithCols(x,"; ", "=",";",cols=colAnnot)))
		  df$filename <- d.names
		  
		  exportAsTable(df=df,file=filesTxtTab)
}


