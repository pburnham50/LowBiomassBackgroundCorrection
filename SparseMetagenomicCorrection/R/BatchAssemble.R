#' BatchAssemble
#' This function is used to assemble a dataframe containing the variance of each taxa in samples.
#' @param BatchList List of batch objects.
#' @param metric "mean" or "variation". "variation" by DEFAULT.
#' @importFrom reshape2 acast
#' @keywords batch
#' @export
#' @examples
#' BatchAssemble()
#'


BatchAssemble = function(BatchList, metric = "variation"){

  n.batchframes = length(BatchList) ;
  batchnames = as.character(sapply(1:n.batchframes,
                                   function(x) paste0(BatchList[[x]]$BatchName,"-",
                                                      BatchList[[x]]$BatchID))) ;

  if( metric == "variation"){
    Batch.Var.Frame = BatchList[[1]]$Taxon.Variation
    colnames(Batch.Var.Frame)[2] = batchnames[1]

    for ( i in 2:n.batchframes){
      Batch.Var.Frame = merge(Batch.Var.Frame, BatchList[[i]]$Taxon.Variation, by ="Tax.ID", all=T)
      colnames(Batch.Var.Frame)[i+1] = batchnames[i]
    }
  }else if( metric == "mean"){
    Batch.Var.Frame = BatchList[[1]]$Mean.Batch.Abundance
    colnames(Batch.Var.Frame)[2] = batchnames[1]

    for ( i in 2:n.batchframes){
      Batch.Var.Frame = merge(Batch.Var.Frame, BatchList[[i]]$Mean.Batch.Abundance, by ="Tax.ID", all=T)
      colnames(Batch.Var.Frame)[i+1] = batchnames[i]
    }
  }


  colnames(Batch.Var.Frame)[which(colnames(Batch.Var.Frame) == "NoBatch-NoBatch")] = "All"

  return(Batch.Var.Frame) ;
}
