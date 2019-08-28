#' AbundanceProportionVariation
#' This function uses the total number of raw reads and BLAST alignment info to determine the proportion of reads from each taxa. This is used for negative control analysis.
#' @param AbundanceObject Matrix including BLAST hits for each sample.
#' @param MetaDataObject Metadata matrix to describe biomass and batches.
#' @param TotalReadsFrame Data frame including sample names and total reads.
#' @param VariationMin Minimum variation to use as a filter. Numeric. 0 by default.
#' @param BatchName Variable name identifying different batches. Character. "all" by default.
#' @param BatchID Batch identifier. Character. NULL by default.
#' @param TaxLevel Taxonomic level to perform analysis. Character. "species" by default.
#' @keywords mass-proportions
#' @import data.table
#' @export
#' @examples
#' AbundanceProportionVariation()


AbundanceProportionVariation <- function(AbundanceObject,
                                         MetaDataObject,
                                         TotalReadsFrame,
                                         VariationMin = 0,
                                         BatchName = "all",
                                         BatchID = NULL,
                                         TaxLevel = "species"){

  # Biomass data frame
  biomass.vec = subset.data.frame(MetaDataObject, select = c("Sample","Biomass")) ;

  # Total reads data frame
  read.frame = TotalReadsFrame ;

  MetaDataObject = merge(MetaDataObject, read.frame, by = "Sample") ;

  # Subset batch
  if(BatchName != "all"){
    meta.sub.frame = subset.data.frame(MetaDataObject,
                                       select = c("Sample", "Biomass", "Total.Reads", BatchName)) ;
    colnames(meta.sub.frame)[4] = "Batch" ;

    merged.matrix = unique(merge(AbundanceObject, meta.sub.frame, by = "Sample")) ;
    merged.matrix$AdjMass = merged.matrix$Biomass * (merged.matrix$AdjustedBlast / merged.matrix$Total.Reads) ;

    selected.matrix = subset.data.frame(x = merged.matrix,
                                        select = c("Sample", "AdjMass", "Biomass", "Batch", TaxLevel)) ;
    colnames(selected.matrix)[5] = "Tax.ID" ;

    selected.matrix = subset.data.frame(selected.matrix, subset = (Batch == BatchID));
    submeta = MetaDataObject[MetaDataObject[[BatchName]] == BatchID,]
  }else{
    meta.sub.frame = subset.data.frame(MetaDataObject,
                                       select = c("Sample", "Biomass", "Total.Reads")) ;

    merged.matrix = unique(merge(AbundanceObject, meta.sub.frame, by = "Sample")) ;
    merged.matrix$AdjMass = merged.matrix$Biomass * (merged.matrix$AdjustedBlast / merged.matrix$Total.Reads) ;

    selected.matrix = subset.data.frame(x = merged.matrix,
                                        select = c("Sample", "AdjMass", "Biomass", TaxLevel)) ;
    colnames(selected.matrix)[4] = "Tax.ID" ;
    submeta = MetaDataObject
  }

  # Aggregate matrix and include zero abundance taxa
  agg.matrix = aggregate(AdjMass ~(Sample + Tax.ID + Biomass) ,selected.matrix, FUN = sum) ;
  full.combo.matrix = data.frame(expand.grid(unique(submeta$Sample), unique(agg.matrix$Tax.ID))) ;
  full.combo.matrix$Var3 = 0
  full.combo.matrix$Var4 = 0
  colnames(full.combo.matrix) = colnames(agg.matrix) ;
  all.matrix = rbind(full.combo.matrix, agg.matrix) ;
  all.matrix = aggregate(AdjMass ~ ( Sample + Tax.ID), all.matrix, FUN = sum) ;

  # Determine covariation matrix
  abund.matrix = acast(all.matrix, Sample ~ Tax.ID, value.var = "AdjMass") ;
  cov.abund.matrix = cov(abund.matrix) ;
  variation.taxa = data.frame("Tax.ID"=names(diag(cov.abund.matrix)),
                              "Variance" = diag(cov.abund.matrix)) ;

  # Specify taxa that need to be dropped and concatenate to sample in Batch
  variation.taxa$Drop = (variation.taxa$Variance < VariationMin)
  taxa.drop = as.numeric(as.matrix(variation.taxa[variation.taxa$Drop,]$Tax.ID))
  sample.drop = unique(all.matrix$Sample)

  #Return list of sample-taxa
  drop.df = data.frame(expand.grid(sample.drop,taxa.drop))
  return(paste0(drop.df[,1],"-",drop.df[,2]));
}

