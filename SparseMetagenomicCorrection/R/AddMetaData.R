#' AddMetaData
#' This function adds metadata to the original matrix. It can include up to one parameter and one batch identifiers.
#' @param SampleVector Character vector of sample names.
#' @param BatchFrame Character vector of names corresponding to batches in sample extraction.
#' @param BatchName Character string of batch.
#' @param ParameterFrame Vector corresponding to information about samples.
#' @param ParameterName Character string to name parameter
#' @keywords addmeta
#' @export
#' @examples
#' AddMetaData()
#'


AddMetaData = function(MetaDataObject,BatchFrame=NULL, BatchName=NULL, ParameterFrame=NULL, ParameterName=NULL){
  original.frame = MetaDataObject ;

  if(!is.null(BatchFrame)){
    colnames(BatchFrame) = c("Sample", BatchName) ;
    new.frame = merge(original.frame, BatchFrame, by="Sample", all.x=T) ;
  }
  if(!is.null(ParameterFrame)){
    colnames(ParameterFrame) = c("Sample", ParameterName) ;
    new.frame = merge(original.frame, ParameterFrame, by="Sample", all.x=T) ;
  }
  if(!is.null(BatchFrame) & !is.null(ParameterFrame)){
    colnames(BatchFrame) = c("Sample", BatchName) ;
    colnames(ParameterFrame) = c("Sample", ParameterName) ;
    new.frame = merge(merge(original.frame, BatchFrame, by="Sample", all.x=T), ParameterFrame, by="Sample", all.x =T) ;
  }
  return(new.frame) ;
}
