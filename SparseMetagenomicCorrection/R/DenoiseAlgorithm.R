#' DenoiseAlgorithm
#' This function uses the total number of raw reads and BLAST alignment info to determine the proportion of reads from each taxa. This is used for negative control analysis.
#' @param AbundanceObject Matrix including BLAST hits for each sample.
#' @param MetaDataObject Metadata matrix to describe biomass and batches.
#' @param NegativeObject Negative control matrix.
#' @param ReadAbundMatrix Matrix listing samples and read abundance.
#' @param CV.Filter Filter matrix based on CV. Logical. TRUE by default.
#' @param MassVar.Filter Filter matrix based on mass variation. Logical. TRUE by default.
#' @param NegCtrl.Filter Filter matrix based on negative control taxa. Logical. TRUE by default.
#' @param deltaCV.Param CV based parameter for filtering difference (max value). Numeric. 10**5 by default.
#' @param MassVar.Param Mass variation based parameter for filtering (min value). Numeric. 0 by default.
#' @param NegCtrl.Param Negative control based parameter for filtering (min value). Numeric. 0 by default.
#' @param TaxLevel Taxonomic level to perform analysis. Character. "species" by default.
#' @param FastqPath Path to raw, gzipped fastqs. Current directory by default.
#' @param TblatPath Path to tblat files. Current directory by default.
#' @param TablePath Path to stored proportion tables. Current directory by default.
#' @param AlnStatsPath Path to CV and other genome stats. Current directory by default.
#' @param Measurement Name of column in abundance object containing measurements. Character. "RelCoverage" by default.
#' @param ReturnContams Return contaminants instead of filtered matrix. Logical. FALSE by default.
#' @keywords denoise
#' @import data.table
#' @import taxize
#' @export
#' @examples
#' DenoiseAlgorithm()


DenoiseAlgorithm <- function(AbundanceObject, MetaDataObject, NegativeObject, ReadAbundMatrix,
                                         CV.Filter = T, MassVar.Filter = T, NegCtrl.Filter = T,
                                         deltaCV.Param = 10**5, MassVar.Param = 0, NegCtrl.Param = 0,
                                         TaxLevel = "genus", FastqPath ="./", TblatPath="./", TablePath="./", AlnStatsPath = "./",
                                         Measurement = "RelCoverage", ReturnContams = F, GITable = "gi_tax_info.tab"){

  sample.list = as.character(as.matrix(MetaDataObject$Sample)) ;

  if(CV.Filter){

    CV.filtered.list = c() ;
    for (j in sample.list){
      CV.filtered.list = rbind(AlignStatsSample(Sample = j,BinSize = 10**3,
                                                TblatPath = TblatPath,
                                                GITable = GITable,
                                                MinHits = 5,
                                                OutPath = AlnStatsPath),
                               CV.filtered.list) ;
    }
    if(nrow(CV.filtered.list[is.na(CV.filtered.list$Delta.CV),]) > 0){
      CV.filtered.list[is.na(CV.filtered.list$Delta.CV),]$Delta.CV = 0 ;
    }


    species_super.df = data.frame(unique(subset.data.frame(AbundanceObject,select = c("Taxid","superkingdom")))) ;
    colnames(species_super.df)[1] = "Tax.ID" ;

    # compute all taxa to include
    CV.filtered.tab = merge(CV.filtered.list,species_super.df, by ="Tax.ID") ;
    CV.filtered.in.tab = subset.data.frame(CV.filtered.tab, (Delta.CV <= (deltaCV.Param))) ;
    CV.filtered.in.tab$sample_tax = paste(CV.filtered.in.tab$Sample, CV.filtered.in.tab$Tax.ID,sep = "-") ;

    AbundanceObject$sample_tax = paste(AbundanceObject$Sample, AbundanceObject$Taxid,sep = "-") ;
    ContamObject.1 = subset.data.frame(AbundanceObject, !(sample_tax %in% as.character(as.matrix(CV.filtered.in.tab$sample_tax)))) ;
    AbundanceObject = subset.data.frame(AbundanceObject, (sample_tax %in% as.character(as.matrix(CV.filtered.in.tab$sample_tax)))) ;

  }

  if(MassVar.Filter){

    ex.batch.list = c()
    ex.batch.names = data.frame(table((MetaDataObject$EX.batch)))
    for (j in as.character(as.matrix(ex.batch.names[ex.batch.names$Freq >= 2,1]))){
      ex.batch.list = c(AbundanceProportionVariation(AbundanceObject = AbundanceObject,
                                                     MetaDataObject = MetaDataObject,
                                                     VariationMin = 10**MassVar.Param,
                                                     BatchName = "EX.batch",
                                                     BatchID = j,
                                                     TotalReadsFrame = ReadAbundMatrix,
                                                     TaxLevel = TaxLevel),
                        ex.batch.list) ;
    }

    lp.batch.list = c()
    lp.batch.names = data.frame(table((MetaDataObject$LP.batch)))
    for (j in as.character(as.matrix(lp.batch.names[lp.batch.names$Freq >= 2,1]))){
      lp.batch.list = c(AbundanceProportionVariation(AbundanceObject = AbundanceObject,
                                                     MetaDataObject = MetaDataObject,
                                                     VariationMin = 10**MassVar.Param,
                                                     BatchName = "LP.batch",
                                                     BatchID = j,
                                                     TotalReadsFrame = ReadAbundMatrix,
                                                     TaxLevel = TaxLevel),
                        lp.batch.list) ;
    }

    all.batch.list = AbundanceProportionVariation(AbundanceObject = AbundanceObject,
                                                   MetaDataObject = MetaDataObject,
                                                   VariationMin = 10**MassVar.Param,
                                                   BatchName = "all",
                                                   TotalReadsFrame = ReadAbundMatrix,
                                                   TaxLevel = TaxLevel)


    filter.batch.list = unique(c(ex.batch.list,lp.batch.list,all.batch.list))

    AbundanceObject$sample_tax = paste(AbundanceObject$Sample, AbundanceObject[[TaxLevel]],sep = "-")
    ContamObject.2 = AbundanceObject[(AbundanceObject$sample_tax %in% filter.batch.list), ]
    AbundanceObject = AbundanceObject[!(AbundanceObject$sample_tax %in% filter.batch.list), ]

  }

  if(NegCtrl.Filter){

    negative.filters = c()
    for (j in sample.list){
      negative.filters = rbind(negative.filters,
                               FilterWithNegativeControl(sample = j,
                                                         contam.set = NegativeObject,out.path = TablePath,
                                                         factor.contam = NegCtrl.Param,
                                                         raw.data.path = FastqPath,
                                                         ReadAbundMatrix = ReadAbundMatrix,
                                                         tblat.path = TblatPath))
    }

    negative.filters.samspec = paste0(negative.filters$Sample,"-",negative.filters$Tax.ID)
    AbundanceObject$samtax = paste0(AbundanceObject$Sample,"-",AbundanceObject$tax_id)
    ContamObject.3 = AbundanceObject[(AbundanceObject$samtax %in% negative.filters.samspec),]
    AbundanceObject = AbundanceObject[!(AbundanceObject$samtax %in% negative.filters.samspec),]
  }

  AbundanceObject$Measurement = AbundanceObject[[Measurement]];
  AbundanceObject$SpecTax = AbundanceObject[[TaxLevel]];

  AggAbundanceObject = aggregate( Measurement ~ (Sample + SpecTax),AbundanceObject, FUN = sum)
  AggAbundanceObject$Name = ncbi_get_taxon_summary(AggAbundanceObject$SpecTax)$name

  if(ReturnContams == F){
    return(AggAbundanceObject);
  }else if(ReturnContams == T & CV.Filter == T & MassVar.Filter == T & NegCtrl.Filter == T) {
    return(list(ContamObject.1,ContamObject.2,ContamObject.3))
  }else{
    return(list(c(),c(),ContamObject.3))
  }


}

