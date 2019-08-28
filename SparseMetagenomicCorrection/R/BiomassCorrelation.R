#' BiomassCorrelation
#' This function looks for correlation between a particular taxon abundance and biomass.
#' For negatively correlated taxa, biomass-proportions are transformed and scaled.
#' @param AbundanceObject Subset object with abundance
#' @param MetaDataObject Metadata object containing biomass.
#' @param BiomassVariable string to describe column name of biomass vector in metadata frame. "Biomass" by DEFAULT.
#' @param Cor.signif.threshold Maximum value for p value of correlated used to initialize Box-cox. 10**-2 by DEFAULT.
#' @param Zscore.threshold Maximum threshold to reject values. 1.65 by DEFAULT (roughly 95%).
#' @importFrom MASS boxcox
#' @keywords biomass
#' @export
#' @examples
#' BiomassCorrelation()
#'


BiomassCorrelation = function(AbundanceObject, MetaDataObject,
                              BiomassVariable = "Biomass", Cor.signif.threshold = 10**-2,
                              Zscore.threshold = 1.65){
  
  sub.abund.df = subset.data.frame(AbundanceObject,select = c("Sample","Taxon","Proportion"))
  sub.meta.df = subset.data.frame(MetaDataObject,select = c("Sample","Biomass"))
  
  # list all taxa
  taxa.list = as.numeric(unique(sub.abund.df$Taxon))
  
  taxa.correlation = c() #initialize vector containing correlation of taxa with biomass
  taxa.reject.list = c() #initialize vector from scaled box-cox transforms
  
  for (j in taxa.list){
    sub.ba.taxon = subset.data.frame(sub.abund.df,Taxon == j, select = c("Sample","Proportion"))
    
    if ( nrow(sub.ba.taxon) > 1){
      merge.ba.taxon = merge(sub.ba.taxon,sub.meta.df,all.y = T) ;
      merge.ba.taxon[is.na(merge.ba.taxon)] = 0 ;
      
      correlation = cor.test(merge.ba.taxon$Proportion,merge.ba.taxon$Biomass, method="spearman")
      correlation.df = data.frame("Taxa" = j, "Samples.with.taxa" = nrow(sub.ba.taxon),
                                  "P.value" = correlation$p.value, "Cor.value"=correlation$estimate) ;
      taxa.correlation = rbind(taxa.correlation, correlation.df) ;
      
      if ( (correlation$p.value <= Cor.signif.threshold) & (correlation$estimate < 0)){
          merge.ba.taxon$ScaleNorm = merge.ba.taxon$Proportion*merge.ba.taxon$Biomass
      
          testing = merge.ba.taxon$ScaleNorm ;
          testing[testing == 0] = min(testing[testing != 0]) ;
          
          #Transform data using Box-Cox power transformation
          Box = boxcox(testing~1, lambda = seq(-6,6,0.01),plotit = F) ;
          Cox = data.frame(Box$x, Box$y) ;
          Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] ;
          lambda = Cox2[1, "Box.x"] ;
          T_box = (testing**(lambda - 1))/lambda ;
          merge.ba.taxon$BCtransformed = T_box ;
          merge.ba.taxon$ScaleNorm = scale(merge.ba.taxon$BCtransformed)[,1] ;

          rejected =  subset.data.frame(merge.ba.taxon, (ScaleNorm <= Zscore.threshold) & (Proportion > 0)) ;
          taxa.reject.list = c(taxa.reject.list, paste(rejected$Sample,j,sep = "-")) ;
      }
    }
  }
  return(list( "Biomass.Correlation" = taxa.correlation,
               "Rejected.Sample-Tax" = taxa.reject.list)) ;
}
