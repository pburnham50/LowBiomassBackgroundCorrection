#' BatchVariationCollect
#' This function looks for correlation between a particular taxon abundance and biomass.
#' @param AbundanceObject Subset object with abundance
#' @param AbundanceVariable Character value on which to produce correlation.
#' @param MetaDataObject Metadata object containing biomass.
#' @param BatchName Character value for Batch column. NULL by DEFAULT.
#' @param BatchIDList Vector of the collection of all bacthes. NULL by DEFAULT.
#' @param Correlate.threshold Float value of lowerbound of correlation between taxa. 10**-4 by DEFAULT
#' @param BiomassInputAdjust Boolean value to declare adjustment by biomass.FALSE by DEFAULT
#' @param BiomassVariable string to describe column name of biomass vector in metadata frame. NULL by DEFAULT
#' @importFrom reshape2 acast
#' @keywords collect-variance
#' @export
#' @examples
#' BatchVariationCollect()
#'


BatchVariationCollect = function(AbundanceObject,AbundanceVariable, MetaDataObject,
                                 BiomassInputAdjust = F, BiomassVariable = NULL,BatchName = NULL,
                                 BatchIDList,Correlate.threshold = 10**-4){

  batch.group.correlates = c() ;
  
  for (i in BatchIDList){
    tmp.cor.df = WithinBatchVariation(AbundanceObject = AbundanceObject,AbundanceVariable = AbundanceVariable,
                                      MetaDataObject = MetaDataObject,
                                      BiomassInputAdjust = BiomassInputAdjust,BiomassVariable = BiomassVariable,
                                      BatchName = BatchName, BatchID = i,
                                      Correlate.threshold = Correlate.threshold) ;
    
    if(!is.null(tmp.cor.df)){cor.var.df = tmp.cor.df$Taxon.Variation ;
    cor.sam.df = tmp.cor.df$Samples ;
    cor.df = merge(cor.var.df, cor.sam.df) ;
    colnames(cor.df)[3] = "Sample" ;
    cor.df$BatchName = BatchName
    cor.df$BatchID = i
    batch.group.correlates = rbind(batch.group.correlates, cor.df)}
  }
  
  batch.group.correlates = batch.group.correlates[!is.na(batch.group.correlates$Variance),] ;
  
  batch.group.correlates$SampleTax = paste(batch.group.correlates$Sample,batch.group.correlates$Tax.ID,sep = "-")
  
  return(batch.group.correlates)
  
}
