#' AlignStatsSample
#' This function generates alignment information and statistics for a specific microbe at taxonomic level.
#' @param Sample Sample name to load in tblat.1 file
#' @param Tax ID number for taxon.
#' @param BinSize Size of bin in bp. Default is 1000 bp. Numeric. 1000 bp by default.
#' @param TblatPath Path to tblat file directory. Character. Current directory by default.
#' @param GITable GI/Tax table path. Character.
#' @param MinHits Minimum number of hits to consider. Numeric. 5 by default.
#' @param OutPath Path to write out tables.
#' @param Overwrite Overwrite existing file. Logical. FALSE by default.
#' @keywords stats
#' @importFrom ineq Gini
#' @importFrom data.table fread
#' @importFrom taxize ncbi_get_taxon_summary
#' @export
#' @examples
#' AlignStatsSample()
#'
#

AlignStatsSample = function(Sample, BinSize = 10**3, TblatPath = "./", GITable = "gi_tax_info.tab",
                            MinHits = 5, OutPath = "./", Overwrite = F){

  if(file.exists( paste0(OutPath,"/", Sample, ".cv_aln_stats.tab")) & !(Overwrite) ){

    Final.df = data.frame(read.table(paste0(OutPath,"/", Sample, ".cv_aln_stats.tab") , sep = "\t", header = T)) ;
  }else{

    # Load in tables
    Gi.table = fread(GITable,sep = "\t",header = T) ;
    Tblat.table = fread(paste0(TblatPath, "/", Sample, ".tblat.1"), sep = "\t",header = F)[,c(2,9,13)] ;
    colnames(Tblat.table) = c("GI.NAME", "GI.POS", "Tax.ID") ;
    Tblat.table$GI.ID = as.numeric(matrix(unlist(strsplit(Tblat.table$GI.NAME,split = "\\|")),ncol = 4,byrow = T)[,2]) ;

    # Merge tables into one object, subset to specific taxon
    Merge.df = merge(Gi.table, Tblat.table) ;
    Merge.df$Tax.POS = Merge.df$GI.OFFSET + Merge.df$GI.POS ;

    # Get list of all taxa IDs
    TaxTable = data.frame(table(Merge.df$Tax.ID)) ;
    colnames(TaxTable) = c("Tax.ID", "Hits") ;
    TaxList = as.numeric(as.matrix(TaxTable[ TaxTable$Hits >= MinHits,]$Tax.ID));

    Final.df = data.frame(matrix(unlist(lapply(TaxList, function (x) AlignStatsTax(Sample = Sample, Tax = x, BinSize = 10**3,
                                                             TblatGIObject = Merge.df))), nrow = length(TaxList), byrow = T)) ;

    colnames(Final.df) = colnames(AlignStatsTax(Sample = Sample, Tax = TaxList[1], BinSize = 10**3, TblatGIObject = Merge.df)) ;
    Final.df$Delta.CV = (Final.df$Obs.CV - Final.df$Sim.CV) ;
    Final.df$Sample = Sample ;

    write.table(Final.df, file = paste0(OutPath,"/", Sample, ".cv_aln_stats.tab") , sep = "\t", quote = F, row.names = F, col.names = T) ;
  }

  return(Final.df) ;
}
