#' WithinBatchVariation
#' This function looks for correlation between a particular taxon abundance and biomass.
#' @param AbundanceObject Subset object with abundance
#' @param AbundanceVariable Character value on which to produce correlation.
#' @param MetaDataObject Metadata object containing biomass.
#' @param BatchName Character value for Batch column. NULL by DEFAULT.
#' @param BatchID Character value for specific batch. NULL by DEFAULT.
#' @param Correlate.threshold Float value of lowerbound of correlation between taxa. 10**-4 by DEFAULT
#' @param BiomassInputAdjust Boolean value to declare adjustment by biomass.FALSE by DEFAULT
#' @param BiomassVariable string to describe column name of biomass vector in metadata frame. NULL by DEFAULT
#' @importFrom reshape2 acast
#' @keywords variance
#' @export
#' @examples
#' WithinBatchVariation()
#'


WithinBatchVariation = function(AbundanceObject,AbundanceVariable, MetaDataObject,BiomassInputAdjust = F, BiomassVariable = NULL,BatchName = NULL,BatchID = NULL,Correlate.threshold = 10**-4){
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

  if(BiomassInputAdjust == T){
    tmp.sub.meta = subset.data.frame(MetaDataObject, select = c("Sample",BiomassVariable))
    tmp.merge.frame = merge(sub.data.frame, tmp.sub.meta) ;
    tmp.merge.frame$Measure = tmp.merge.frame$Measure * tmp.merge.frame[[BiomassVariable]] ;
    sub.data.frame = tmp.merge.frame[,1:3]
  }

  if(nrow(sub.data.frame)>0){
      # produce abundance matrix
  abund.matrix = acast(sub.data.frame, Sample ~ Taxon, value.var = "Measure") ;
  abund.matrix[is.na(abund.matrix)] = 0

  # produce covariate matrix
  cov.abund.matrix = cov(abund.matrix)

  variation.taxa = data.frame("Tax.ID"=names(diag(cov.abund.matrix)),
                              "Variance" = diag(cov.abund.matrix))

  # List correlating microbes above threshold
  melt.frame = melt(cov.abund.matrix) ;
  melt.frame = melt.frame[melt.frame$Var1 != melt.frame$Var2,] ;
  melt.frame = melt.frame[melt.frame$value >= Correlate.threshold,] ;
  colnames(melt.frame) = c("Tax.ID.1", "Tax.ID.2", "Value")

  return(list("Taxon.Variation" = variation.taxa,
              "Taxon.Correlates" = melt.frame,
              "Samples" = samples.save,
              "BatchName" = BatchName,
              "BatchID" = BatchID)) ;
  }
}
