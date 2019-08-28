#' TotalReadsGen
#' This function uses the total number of raw reads and BLAST alignment info to determine the proportion of reads from each taxa. This is used for negative control analysis.
#' @param MetaDataObject Metadata matrix to describe biomass and batches.
#' @param TotalReadsOutput Path and filename to store TotalReadsFrame. Character. NA by default.
#' @param RawDataPath Path where zipped fastq files reside. Character. NA by default.
#' @keywords totalreads
#' @import data.table
#' @export
#' @examples
#' TotalReadsGen()


TotalReadsGen <- function(MetaDataObject, TotalReadsOutput = NA, RawDataPath = NA){
  
  if(!file.exists(TotalReadsOutput)){
      sample.vec = unique(MetaDataObject$Sample) ;
    total.reads.vec = c() ;
  
    for( j in sample.vec){
      wcl = system(command = paste0("zcat ", RawDataPath, "/",j, ".R1.fastq.gz | wc -l"), intern = T) ;
      total.reads = as.numeric(as.matrix(strsplit(wcl, split = " ")[[1]][1]))/4 ;
      total.reads.vec = c(total.reads.vec, total.reads) ; 
    }
    
    read.frame = data.frame("Sample" = sample.vec, "Total.Reads" = total.reads.vec) ;
    write.table(x = read.frame, file = TotalReadsOutput, quote = F, sep = "\t", row.names = F, col.names = T) ;
  }else{
    read.frame = data.frame(read.table(TotalReadsOutput, sep="\t", header = T)) ;
  }

  return(read.frame) ;
}

