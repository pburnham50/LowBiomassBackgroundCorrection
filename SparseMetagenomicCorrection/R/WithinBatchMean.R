#' WithinBatchMean
#' This function looks for correlation between a particular taxon abundance and biomass.
#' @param AbundanceObject Subset object with abundance
#' @param AbundanceVariable Character value on which to produce correlation.
#' @param MetaDataObject Metadata object containing biomass.
#' @param BatchName Character value for Batch column. NULL by DEFAULT.
#' @param BatchID Character value for specific batch. NULL by DEFAULT.
#' @importFrom reshape2 acast
#' @keywords variance
#' @export
#' @examples
#' WithinBatchMean()
#'


WithinBatchMean = function(AbundanceObject,AbundanceVariable, MetaDataObject,BatchName = NULL,BatchID = NULL){
  # Clean frames

  if(!is.null(BatchName) & !is.null(BatchID )){
    sub.meta.frame = subset.data.frame(MetaDataObject, select =  c("Sample",BatchName)) ;
    colnames(sub.meta.frame) = c("Sample", "Batch") ;
    sub.meta.frame = subset.data.frame(sub.meta.frame, subset = (Batch == BatchID))
    samples.save = as.character(sub.meta.frame$Sample) ;
    sub.data.frame = subset.data.frame(AbundanceObject, subset = (Sample %in% samples.save), select = c("Sample","Taxon",AbundanceVariable )) ;

  }else{
    samples.save = unique(as.character(as.matrix(AbundanceObject$Sample)));
    sub.data.frame = subset.data.frame(AbundanceObject, select = c("Sample","Taxon",AbundanceVariable )) ;
    BatchName = "NoBatch"
    BatchID = "NoBatch"
  }

  colnames(sub.data.frame)[3] = "Measure" ;

  # produce abundance matrix
  abund.matrix = acast(sub.data.frame, Sample ~ Taxon, value.var = "Measure") ;
  abund.matrix[is.na(abund.matrix)] = 0

  # produce mean data frame
  df.mean = data.frame("Tax.ID"=colnames(abund.matrix),
                              "Mean" = colSums(abund.matrix))

  return(list("Mean.Batch.Abundance" = df.mean,
              "Samples" = samples.save,
              "BatchName" = BatchName,
              "BatchID" = BatchID)) ;
}
