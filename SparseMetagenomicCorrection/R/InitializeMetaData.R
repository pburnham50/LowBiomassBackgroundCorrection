#' InitializeMetaData
#' This function initializes metadata matrix. It should at least include batch information.
#' @param SampleVector Character vector of sample names.
#' @param BatchVector Character vector of names corresponding to batches in sample extraction.
#' @param BatchName Character string of batch.
#' @keywords meta
#' @importFrom reshape2 acast
#' @export
#' @examples
#' InitializeMetaData()
#'


InitializeMetaData = function(SampleVector, BatchVector, BatchName){
  init.data.frame = data.frame("Sample" = SampleVector, BatchVector) ;
  colnames(init.data.frame)[2] = BatchName ;
  return(init.data.frame) ;
}
