#' FilterWithNegativeControl
#' This function is used to find the proportion of all taxa present in the negative controls.
#' @param sample Sample name. Character.
#' @param contam.set Data frame of taxa in the contaminant samples.
#' @param factor.contam Multiplier to set limits for contaminant taxa. Numeric.
#' @param filter.out.list Determines if the output is filtered out taxa (TRUE) or filtered in taxa (FLASE). Logical. TRUE by default.
#' @param raw.data.path Pathway to find zipped fastq files. Character.
#' @param tblat.path Pathway to alignment files (tblat.1 format). Character.
#' @param out.path Pathway to store total reads. Character.
#' @keywords filter
#' @import data.table
#' @export
#' @examples
#' FilterWithNegativeControl()

FilterWithNegativeControl <- function(sample, contam.set, out.path = "./",
                                      factor.contam = 10,
                                      filter.out.list = T,
                                      raw.data.path, tblat.path){

  # create read proportion dataset
  sample.set = ReadProportions(sample = sample, out.path = out.path,
                               raw.data.path = raw.data.path, tblat.path = tblat.path) ;

  # set the upper bound for all contaminant taxa
  contam.set$Upper = factor.contam*contam.set$Proportion ;

  # merge sets
  contam.set = subset.data.frame(contam.set, select = c("Tax.ID", "Upper")) ;
  merge.df = merge(sample.set, contam.set, by = "Tax.ID",all.x=T) ;
  merge.df[is.na(merge.df)] = 0 ;

  # list the sample-taxa which are under the contaminant limits
  merge.df$Drop = (merge.df$Upper > merge.df$Proportion) ;

  return(data.frame("Sample" = sample, "Tax.ID" = merge.df[merge.df$Drop == filter.out.list,]$Tax.ID))
}
