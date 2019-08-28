#' GenerateStats
#' This function allows specification of taxon level if applicable and will sum data to that level, depending on specific column being measured.
#' @param AbundanceObject Abundance matrix to be loaded. Must be filtered to tri-column form.
#' @keywords stats
#' @importFrom reshape2 acast
#' @export
#' @examples
#' GenerateStats()

GenerateStats <- function(AbundanceObject){
  #cast into matrix
  abund.matrix = acast(AbundanceObject, Taxon ~ Sample, value.var = "Measurement") ;
  abund.matrix[is.na(abund.matrix)] = 0

  total.microbial.abundance = colSums(abund.matrix)
  abund.matrix[abund.matrix>0] = 1
  total.microbes = colSums(abund.matrix)

  nTaxa = dim(abund.matrix)[1]
  nSamples = dim(abund.matrix)[2]

  return(list(data.frame("Tax.number" = nTaxa,
                    "Sample.number" = nSamples),
              data.frame("Sample" = names(total.microbial.abundance),
                         "Total.abundance" = total.microbial.abundance,
                         "Total.microbes" = total.microbes)
              )
         )

}
