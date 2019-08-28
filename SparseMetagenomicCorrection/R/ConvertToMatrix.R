#' ConvertToMatrix
#' This function allows specification of taxon level if applicable and will sum data to that level, depending on specific column being measured.
#' @param AbundanceObject Abundance matrix to be loaded. Must be filtered to tri-column form.
#' @param CastingValue Character value, reflecting column name, to fill in matrix. "Measurement" is DEFAULT.
#' @param Samples.filter Character vector of samples to remove. NULL by DEFAULT.
#' @param Taxon.filter Character vector of taxa to remove. NULL by DEFAULT.
#' @keywords convert
#' @importFrom reshape2 acast
#' @export
#' @examples
#' ConvertToMatrix()

ConvertToMatrix <- function(AbundanceObject, CastingValue, Samples.filter = NULL, Taxon.filter = NULL){
  #cast into matrix
  abund.matrix = acast(AbundanceObject, Taxon ~ Sample, value.var = CastingValue) ;
  abund.matrix[is.na(abund.matrix)] = 0

  if(!is.null(Samples.filter)){
    abund.matrix = abund.matrix[,!(colnames(abund.matrix) %in% Samples.filter)]
  }
  if(!is.null(Taxon.filter)){
    abund.matrix = abund.matrix[!(row.names(abund.matrix) %in% Taxon.filter),]
  }


  return(abund.matrix)
}
