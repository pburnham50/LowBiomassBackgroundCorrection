#' AlignmentStatisticsforTax
#' This function generates alignment information and statistics for a specific microbe at taxonomic level.
#' @param Sample Abundance matrix to be loaded. Must be filtered to tri-column form.
#' @param tax.id ncbi taxon identifier
#' @param tax.level ncbi taxification level. "species" by DEFAULT.
#' @param Merge.data.frame Data frame including alignment and abundance information
#' @param bin.size Size of bin in bp. Default is 1000 bp.
#' @keywords stats
#' @importFrom ineq Gini
#' @export
#' @examples
#' AlignmentStatisticsforTax()
#'
#
AlignmentStatisticsforTax = function(sample,tax.id, tax.level = "species", bin.size = 10**3, Merge.data.frame){

  # choose by species
  species.data.frame = Merge.data.frame[Merge.data.frame[[tax.level]] == tax.id,]

  # take average genome length and total RGE for the species (in case aggregating across strains)
  compact.strain = unique(data.frame(species.data.frame$Taxid,species.data.frame$Length,species.data.frame$RelCoverage, species.data.frame$AdjustedBlast))
  max.genome.length = max(as.numeric(as.matrix(compact.strain$species.data.frame.Length)))
  min.genome.length = min(as.numeric(as.matrix(compact.strain$species.data.frame.Length)))
  mean.genome.length = round(sum(compact.strain$species.data.frame.Length*(compact.strain$species.data.frame.AdjustedBlast / sum(compact.strain$species.data.frame.AdjustedBlast))))
  RGE = sum(as.numeric(as.matrix(compact.strain$species.data.frame.RelCoverage)))
  AdjBlast = sum(as.numeric(as.matrix(compact.strain$species.data.frame.AdjustedBlast)))

  if(mean.genome.length > 5*bin.size){
    # Find starting positions of reads and bin the genome
    start.positions = as.numeric(as.matrix(species.data.frame$SeqStart)) ;
    initiate.genome = table(factor(round(start.positions,-1*log10(bin.size)),levels = seq(0,mean.genome.length+bin.size,bin.size))) ;
    aggregate.coverage.df = data.frame("Bin"=seq(0,mean.genome.length+bin.size,bin.size),
                                       "Count"=as.numeric(as.matrix(initiate.genome))) ;


    # Calculate coverage stats
    gini.coefficient = ineq::Gini(aggregate.coverage.df$Count) ;
    mean.binned.coverage = mean(aggregate.coverage.df$Count) ;
    coeff.variation = sd(aggregate.coverage.df$Count)/mean.binned.coverage ;

    # Determine simulated stats
    div.genome.size = ceiling(mean.genome.length/bin.size)
    simulated.covering = as.numeric(as.matrix(unlist(sapply(c(0:10000),function(x) rep(x,round(dpois(x,lambda = length(start.positions)/div.genome.size)*(div.genome.size)))))))
    simulated.covar = (sd(simulated.covering,na.rm = T)/mean(simulated.covering,na.rm = T))
    simulated.gini = ineq::Gini(simulated.covering)

    return(data.frame("Sample" = sample,
                      "Tax.lvl" = tax.level,
                      "Tax.ID" = tax.id,
                      "Max.Genome.Size" = max.genome.length,
                      "Mean.Genome.Size" = mean.genome.length,
                      "Min.Genome.Size" = min.genome.length,
                      "Bin.Size" = bin.size,
                      "Gini" = gini.coefficient,
                      "Simul.Gini" = simulated.gini,
                      "CV" = coeff.variation,
                      "Simul.CV" = simulated.covar,
                      "Adjusted.Blast" = AdjBlast,
                      "Rel.Coverage" = RGE,
                      "Total.Hits" = length(start.positions)))
  }
}
