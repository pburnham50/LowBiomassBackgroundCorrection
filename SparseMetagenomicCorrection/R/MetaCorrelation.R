#' MetaCorrelation
#' This function looks for correlation between a particular taxon abundance and biomass.
#' @param AbundanceObject Subset object with abundance
#' @param AbundanceVariable Character value on which to produce correlation.
#' @param MetaDataObject Metadata object containing biomass.
#' @param BiomassVariable Character value for biomass column.
#' @param tax.ID NCBI taxification identifier
#' @param tax.ID NCBI taxification identifier.
#' @keywords biomass
#' @export
#' @examples
#' MetaCorrelation()
#'


MetaCorrelation = function(tax.ID, AbundanceObject,AbundanceVariable, MetaDataObject,FeatureVariable){
  # Clean frames
  sub.data.frame = subset.data.frame(AbundanceObject, subset = (Taxon == tax.ID), select = c("Sample",AbundanceVariable )) ;
  colnames(sub.data.frame) = c("Sample", "Measure") ;
  non.empty = nrow(sub.data.frame) ;

  sub.meta.frame = subset.data.frame(MetaDataObject, select =  c("Sample",FeatureVariable)) ;
  colnames(sub.meta.frame) = c("Sample", "Feature") ;

  # merge frames
  merge.df = merge(sub.data.frame, sub.meta.frame, by="Sample",all.y=T)
  merge.df[is.na(merge.df$Measure),]$Measure = min(merge.df$Measure,na.rm = T)

  if (non.empty > 1){
        # determine correlation
        cor.values = cor.test(merge.df$Feature,merge.df$Measure, method = "spearman")

        # summarize in dataframe
        summary.frame = data.frame("Taxon" = tax.ID,
                                   "Feature" = FeatureVariable,
                                   "Measure" = AbundanceVariable,
                                   "Correlation" = cor.values$estimate,
                                   "P.value" = cor.values$p.value,
                                   "Number.non.empty" = non.empty)
  }else{
        # summarize in dataframe
        summary.frame = data.frame("Taxon" = tax.ID,
                                   "Feature" = FeatureVariable,
                                   "Measure" = AbundanceVariable,
                                   "Correlation" = NA,
                                   "P.value" = NA,
                                   "Number.non.empty" = non.empty)
      }
  return(summary.frame) ;
}
