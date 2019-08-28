#' LoadAbundance
#'
#' This function reads in an abundance matrix. Matrix should have rows corresponding to samples and microbes, e.g. row 1: Sample #1, Escherichia coli, 0.78
#' @param dir Directory containing microbial abundance matrix. DEFAULT is current directory.
#' @param file Abundance matrix filename.
#' @keywords load
#' @export
#' @examples
#' LoadAbundance()

LoadAbundance <- function(dir = "./", file){
  if(!file.exists(paste0(dir,file))){
    stop("File not found")
  }
  return(fread(paste0(dir,file), sep="\t",header=T))
}
