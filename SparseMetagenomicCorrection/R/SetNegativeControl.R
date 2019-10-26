#' SetNegativeControl
#' This function is used to find the proportion of all taxa present in the negative controls.
#' @param sample.vector Vector of sample names. Character vector.
#' @param raw.data.path Pathway to find zipped fastq files. Character.
#' @param tblat.path Pathway to alignment files (tblat.1 format). Character.
#' @param table.path Pathway to read proportiona table. Character.
#' @keywords negatives
#' @import data.table
#' @export
#' @examples
#' SetNegativeControl()

SetNegativeControl <- function(sample.vector, raw.data.path, tblat.path, table.path){
  #instantiate list of all contaminant taxa
  contam.taxa.list = c()

  # generate all taxa
  for (i in sample.vector){
    tax.temp = unique(data.frame(fread(paste0(tblat.path, "/",i,".tblat.1"), header = F, sep = "\t"))[,13]) ;
    contam.taxa.list = c(contam.taxa.list, tax.temp)
  }
  contam.taxa_all.list = data.frame("Tax.ID"=unique(contam.taxa.list))


  all.contam.data = c()

  for ( j in sample.vector){
    tmp.contam.data = SparseMetagenomicCorrection::ReadProportions(sample = j,
                                                                    raw.data.path = raw.data.path,
                                                                    tblat.path = tblat.path,out.path = table.path) ;
    tmp.contam.data = merge(contam.taxa_all.list, tmp.contam.data, by = "Tax.ID",all.x=T)
    tmp.contam.data[is.na(tmp.contam.data$Frequency),]$Frequency = 0
    tmp.contam.data[is.na(tmp.contam.data$Proportion),]$Proportion = 0
    tmp.contam.data[is.na(tmp.contam.data$Total),]$Total =  tmp.contam.data[!is.na(tmp.contam.data$Total),]$Total[1]

    all.contam.data = rbind(all.contam.data,tmp.contam.data) ;

  }

  ### aggregate contaminant taxa and find the mean
  agg.data = aggregate( Proportion ~ Tax.ID, all.contam.data, FUN = mean) ;

  return(agg.data)

}
