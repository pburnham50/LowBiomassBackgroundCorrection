#' AlignStatsTax
#' This function generates alignment information and statistics for a specific microbe at taxonomic level.
#' @param Sample Sample name to load in tblat.1 file
#' @param Tax ID number for taxon. 
#' @param BinSize Size of bin in bp. Default is 1000 bp. Numeric. 1000 bp by default.
#' @param TblatGIObject Tblat / GI/Tax table path. Character.
#' @keywords stats
#' @importFrom ineq Gini
#' @importFrom data.table fread
#' @export
#' @examples
#' AlignStatsTax()
#'
#

AlignStatsTax = function(Sample, Tax, BinSize = 10**3, TblatGIObject){

  # Subset to specific taxon
  Sub.df = subset.data.frame(TblatGIObject, Tax.ID == Tax) ;
  
  # calculate genome length
  Genome.length = unique(Sub.df$Tax.LENGTH)[1]
  
  # Bin genome and aggregate
  Start.positions = as.numeric(as.matrix(Sub.df$Tax.POS)) ;
  Init.genome = table(factor(round(Start.positions,-1*log10(BinSize)),
                                 levels = seq(0, Genome.length + BinSize, BinSize))) ;
  Agg.df = data.frame("Bin"=seq(0,Genome.length + BinSize, BinSize),
                      "Count"=as.numeric(as.matrix(Init.genome))) ;
  
  # Simulate genome with same number of hits
  Div.Genome.length = ceiling(Genome.length/BinSize) ;
  Sim.df = as.numeric(as.matrix(unlist(sapply(c(0:10000),
                      function(x) rep(x,round(dpois(x,
                      lambda = length(Start.positions)/Div.Genome.length)*(Div.Genome.length))))))) ;
  
  # Calculate statistics
  Mean.cov = mean(Agg.df$Count) ;
  Obs.cv = sd(Agg.df$Count)/Mean.cov ;
  Sim.cv = sd(Sim.df,na.rm = T)/mean(Sim.df,na.rm = T) ;
  
  return(data.frame("Sample" = Sample,
                    "Tax.ID" = Tax,
                    "Genome.Size" = Genome.length,
                    "Bin.Size" = BinSize,
                    "Hits.per.Bin" = Mean.cov,
                    "Obs.CV" = Obs.cv,
                    "Sim.CV" = Sim.cv,
                    "Total.Hits" = length(Start.positions))) ;
  
}
