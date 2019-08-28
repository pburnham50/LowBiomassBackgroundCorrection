#' BinaryScore
#' This function determines a score from linear combination of weights with respect to TP, TN, FP, FN, and untested counts.
#' @param AbundanceObject Subset object with abundance
#' @param MetaDataObject Metadata object containing biomass.
#' @param AbundanceTaxName Character value for vector name of clinical taxa identified.
#' @param MetaDataTaxName Character value for vector name of metagenomic taxa identified.
#' @param MetaDataTaxName String to identify negative tests. "Negative" by DEFAULT.
#' @param CIBacteriaVector Character vector of taxa tested. Enterococcus and Escherichia genera by DEFAULT.
#' @param WeightsVector Numeric vector of weights for (in order): true positive, true negative, false postive, false negative, and untested. DEFAULT vector provided.
#' @keywords score
#' @export
#' @examples
#' BinaryScore()
#'

BinaryScore = function(AbundanceObject, MetaDataObject, AbundanceTaxName, MetaDataTaxName, NegativeName = "Negative",
                                CIBacteriaVector = c("Escherichia","Enterococcus"), WeightsVector = c(1,2,-1,-3,-0.1)){
  # Clean frames

  abund.tmp.frame = subset.data.frame(AbundanceObject, select = c("Sample", AbundanceTaxName, "Measurement"))
  meta.tmp.frame = subset.data.frame(MetaDataObject, select = c("Sample", MetaDataTaxName))
  colnames(abund.tmp.frame)[2] = "MS.Tax"
  colnames(meta.tmp.frame)[2] = "CI.Tax"

  merge.frame = merge(abund.tmp.frame, meta.tmp.frame, by = "Sample", all = T)
  merge.frame[is.na(merge.frame$MS.Tax),]$MS.Tax = NegativeName
  merge.frame[is.na(merge.frame$Measurement),]$Measurement = 0

  merge.frame$Match = (merge.frame$MS.Tax == merge.frame$CI.Tax)

  Total.tests.count = nrow(meta.tmp.frame)
  TotalPositive.tests.count = nrow(meta.tmp.frame[meta.tmp.frame$CI.Tax != NegativeName,])
  TotalNegative.tests.count = nrow(meta.tmp.frame[meta.tmp.frame$CI.Tax == NegativeName,])


  TruePositive.count = nrow(subset.data.frame(merge.frame, (Match == T)&(CI.Tax != NegativeName))) ;
  FalsePositive.count = length(unique(subset.data.frame(merge.frame, (MS.Tax%in% CIBacteriaVector)&(CI.Tax == NegativeName))$Sample)) ;

  TrueNegative.count = (TotalNegative.tests.count - FalsePositive.count) ;
  FalseNegative.count = (TotalPositive.tests.count - TruePositive.count) ;

  Nontested.count = nrow(subset.data.frame(merge.frame,  !(MS.Tax%in% CIBacteriaVector))) ;

  count.vector = c(TruePositive.count,TrueNegative.count,FalsePositive.count,FalseNegative.count,Nontested.count)

  Score = sum(WeightsVector * count.vector) / sum(WeightsVector[1:2]*c(TotalPositive.tests.count,TotalNegative.tests.count)) ;

  return(Score) ;
}
