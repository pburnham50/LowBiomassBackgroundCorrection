#' SubsetAbundance
#' This function allows specification of taxon level if applicable and will sum data to that level, depending on specific column being measured.
#' @param AbundanceObject Abundance matrix to be loaded.
#' @param tax.level Column name corresponding to taxification on which data will be aggregated. "species" by default.
#' @param measurement Column name corresponding to abundance measurement.
#' @param sample.column Column number corresponding to sample.
#' @param include.samples Names of samples to be included for subsetting.
#' @param exclude.samples Names of samples to be excluded for subsetting.
#' @param proportional.abundance Determine the abundance with respect to other measurements within the same Sample. TRUE by default.
#' @keywords subset
#' @export
#' @examples
#' SubsetAbundance()

SubsetAbundance <- function(AbundanceObject,
                            tax.level = "species",
                            measurement,
                            sample.column,
                            include.samples=c(),
                            exclude.samples=c(),
                            proportional.abundance=T){
  #subset abundance
  colnames(AbundanceObject)[sample.column] = "Sample"
  AbundanceObject.subset = subset.data.frame(x=AbundanceObject, select = c("Sample",tax.level, measurement)) ;
  colnames(AbundanceObject.subset) = c("Sample", "Taxon", "Measurement")

  #include or exclude samples
  if(length(include.samples)>0){
    AbundanceObject.subset = AbundanceObject.subset[AbundanceObject.subset$Sample %in% include.samples,]
  }
  if(length(exclude.samples)>0){
    AbundanceObject.subset = AbundanceObject.subset[!(AbundanceObject.subset$Sample %in% exclude.samples),]
  }

  #aggregate data
  AbundanceObject.subset.aggregate = aggregate(as.numeric(as.matrix(Measurement)) ~ (Sample + Taxon), AbundanceObject.subset, FUN = sum) ;
  colnames(AbundanceObject.subset.aggregate) = c("Sample", "Taxon", "Measurement")

  #determnine proportion of abundance from in regards to each sample.
  if(proportional.abundance ==T){
    AbundanceObject.subset.aggregate$Proportion = with(AbundanceObject.subset.aggregate, Measurement / ave(Measurement, Sample, FUN = sum))
  }

  return(AbundanceObject.subset.aggregate)
}
