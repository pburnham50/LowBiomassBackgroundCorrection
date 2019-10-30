#' ReadProportions
#' This function uses the total number of raw reads and BLAST alignment info to determine the proportion of reads from each taxa. This is used for negative control analysis.
#' @param sample Sample name. Character.
#' @param out.path Pathway to .totalreads.tab file. Current directory by default.
#' @param raw.data.path Pathway to find zipped fastq files. Character.
#' @param tblat.path Pathway to alignment files (tblat.1 format). Character.
#' @param ReadAbundMatrix Matrix listing samples and read abundance.
#' @keywords proportions
#' @import data.table
#' @export
#' @examples
#' ReadProportions()

ReadProportions <- function(sample, out.path ="./", raw.data.path, tblat.path, ReadAbundMatrix = c()){
  if(!(file.exists(paste0(out.path,"/",sample,".totalreads.tab")))){
      
      if(is.null(ReadAbundMatrix)){
        # get total number of raw reads
        wcl = system(command = paste0("zcat ", raw.data.path, "/",sample, ".R1.fastq.gz | wc -l"), intern = T) ;
        total.reads = as.numeric(as.matrix(strsplit(wcl, split = " ")[[1]][1]))/4 ;
      }else{
        total.reads = ReadAbundMatrix[ReadAbundMatrix$Sample == sample,]$Total.Reads ;
      }

      # determine proportionate reasd from alignment data
      tblat = data.frame(fread(paste0(tblat.path,"/",sample,".tblat.1"), header = F, sep = "\t")) ;
      tblat.counts = data.frame(table(tblat[,13])) ;
      colnames(tblat.counts) = c("Tax.ID","Frequency") ;
      tblat.counts$Proportion = tblat.counts$Frequency / total.reads ;
      tblat.counts$Total = total.reads ;

      write.table(x = tblat.counts, file = paste0(out.path,"/",sample,".totalreads.tab"),
                  row.names = F, col.names = T, quote = F, sep ="\t") ;
  }else{
    tblat.counts = read.table(paste0(out.path,"/",sample,".totalreads.tab"), sep = "\t", header = T) ;
  }

  return(tblat.counts)
}
