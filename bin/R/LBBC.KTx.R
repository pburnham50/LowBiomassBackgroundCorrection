### Repruding figure 1 and some of the statistics discussed in the text
# Filename: LBBC.KTx.R
# Author: Philip Burnham
# Year: 2019

### Libraries -----------
source("bin/R/load_packages.R") # install and load required packages

# !!! set path to the main directory of the LBBC Github Repo
#setwd("/path/to/GitHub/LowBiomassBackgroundCorrection/")      # <------- Fix this to proper directory, uncomment and run.

# install and load the package to perform filtering
install("./SparseMetagenomicCorrection") 
library(SparseMetagenomicCorrection)

# set relavant path information (relative to LBBC cloned directory)
path.grammy = "grammys/"   # path to grammy files
path.metadata = "metadata/" # path to metadata files
path.reads = "total_reads/" # path to tables with total number of reads for each sample
path.stats = "aln_stats/" # path to directory where CV tables are stored.
path.tblat = "tblats/" # path to where BLAST alignment output for each sample is stored.
path.figs = "figs/" # path to figure

### Load abundance matrix -------------------------------------------------------------------------
KT.abundance = LoadAbundance(dir = path.grammy, file = "KTx.SMA.grammy.tab") ;
KT.abundance = subset.data.frame(KT.abundance, superkingdom == 2) ; # here we select only bacteria
KT.abundance = subset.data.frame(KT.abundance, RelCoverage >= 10**-6)
colnames(KT.abundance)[2] = "Sample"
KT.abundance$Measurement = KT.abundance$RelCoverage

### Load metadata matrices ------------------------------------------------------------------------
clinical.metadata = data.frame(read.table(paste0(path.metadata,"KTx_SMA.metadata.upd.tab"), header = T, sep = "\t")) ;
lab.metadata = data.frame(read.table(paste0(path.metadata,"cfDNA_012519.csv"),header = T, sep = ",",fill = T)) ;
lab.metadata = lab.metadata[lab.metadata$Study_name == "KTx",] ;

colnames(lab.metadata)[1] = "Sample"
lab_meta_data = merge(lab.metadata,clinical.metadata, by="Sample")

KT.meta = InitializeMetaData(SampleVector = as.character(as.matrix(lab_meta_data$Sample)),
                             BatchVector = paste0(as.character(as.matrix(lab_meta_data$Date_extraction)),"_",
                                                  as.character(as.matrix(lab_meta_data$Batch_extraction))),
                             BatchName = "EX.batch") ;

lp.meta = subset.data.frame(lab.metadata,select = c("Sample","Date_libprep","Batch_libprep"))
lp.meta$Batch = paste0(lp.meta$Date_libprep,"_",lp.meta$Batch_libprep)

KT.meta = AddMetaData(MetaDataObject = KT.meta,BatchFrame = lp.meta[,c(1,4)],BatchName = "LP.batch") ;

KT.meta = AddMetaData(MetaDataObject = KT.meta,
                      ParameterFrame = subset.data.frame(lab_meta_data, select = c("Sample", "Biomass_at_libprep_ng")),
                      ParameterName = "Biomass") ;

KT.meta = AddMetaData(MetaDataObject = KT.meta,
                      ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Rec.gender")),
                      ParameterName = "R.gender") ;

KT.meta = AddMetaData(MetaDataObject = KT.meta,
                      ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Status")),
                      ParameterName = "Status") ;

KT.meta = AddMetaData(MetaDataObject = KT.meta,
                      ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Bacteria.genus")),
                      ParameterName = "CI.genus") ;

### Load abundance matrix -------------------------------------------------------------------------
Read.Abund.Matrix = TotalReadsGen(KT.meta,
                                  TotalReadsOutput = paste0(path.reads,"KTmeta.totalreads.tab"),
                                  RawDataPath = "./") ; # raw data path is where fastqs are stored. not necessart for this.

colnames(KT.abundance)[2] = "Sample" ;


### Load taxa from negative controls --------------------------------------------------------------
negatives = SetNegativeControl(sample.vector = paste0("MC",LETTERS[c(1:9,11:14,16:20)]),
                               raw.data.path = "./",table.path = path.reads,
                               tblat.path = path.tblat)


### Optimization of parameters --------------------------------------------------------------------
#Here we impliment the optimized parameters

weights = c(4, 2,-1,-2,-0.2) ; # weights for true pos., true neg., etc...

# function to produce a score which can be optimized
Testing.Score = function(ParamVector = c(2,-6),WeightVector = weights){
  score = BinaryScore(AbundanceObject =  DenoiseAlgorithm(AbundanceObject = KT.abundance,MetaDataObject = KT.meta,
                                                          NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                                                          CV.Filter = T, MassVar.Filter = T,NegCtrl.Filter = T,
                                                          deltaCV.Param = ParamVector[1],
                                                          MassVar.Param = ParamVector[2],
                                                          NegCtrl.Param = 10,
                                                          TaxLevel = "genus", FastqPath = "./",TablePath = path.reads, TblatPath = path.tblat, AlnStatsPath = path.stats,
                                                          GITable = "lookups/gi_tax_info.tab"),
                      MetaDataObject = KT.meta,
                      AbundanceTaxName = "Name", MetaDataTaxName = "CI.genus",NegativeName = "Negative",
                      CIBacteriaVector = c("Escherichia","Enterococcus"), WeightsVector =WeightVector) ;

  return(score) ;
}

f <- function(x) {(1-Testing.Score(ParamVector = c(x[1],x[2])))}
optimized_parameters = optim(c(2,-5),f) # seed at CV = 2, Batch = 10 pg**2


### Denoise the abundance matrix ------------------------------------------------------------------
withfilt.tab = DenoiseAlgorithm(AbundanceObject = KT.abundance,MetaDataObject = KT.meta,
                       NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = T, MassVar.Filter = T,NegCtrl.Filter = T,
                       deltaCV.Param = optimized_parameters$par[1], MassVar.Param = optimized_parameters$par[2], 
                      NegCtrl.Param = 10,TaxLevel = "genus",FastqPath = "./",TablePath = path.reads, TblatPath = path.tblat, AlnStatsPath = path.stats,
                      GITable = "lookups/gi_tax_info.tab")

final.withfilt.tab = merge(withfilt.tab, KT.meta, "Sample")

withoutfilt.tab = DenoiseAlgorithm(AbundanceObject = KT.abundance,MetaDataObject = KT.meta,
                       NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = F, MassVar.Filter = F,NegCtrl.Filter = F,
                       TaxLevel = "genus",FastqPath = "./",TablePath = path.reads, TblatPath = path.tblat, AlnStatsPath = path.stats,
                       GITable = "lookups/gi_tax_info.tab")
final.withoutfilt.tab = merge(withoutfilt.tab, KT.meta, "Sample")



### Plot  -----------------------------------------------------------------------------------------
final.withfilt.tab$UTI = !(final.withfilt.tab$Name %in% c("Escherichia", "Enterococcus"))
final.withoutfilt.tab$UTI = !(final.withoutfilt.tab$Name %in% c("Escherichia", "Enterococcus"))
KT.meta$Name = final.withfilt.tab$Name[1]
KT.meta$UTI = F

high_disease = paste0("GU",c(74,81,83,467,1209,5,14,65,69,1125,1495,2709))

withfilt = ggplot(final.withfilt.tab,aes(Sample,Name))+
                  facet_grid(UTI~paste0(R.gender ,"\n",substr(CI.genus,1,3)), scales = "free", space = "free")+
                  geom_tile(aes(fill=log10(Measurement)),col="black")+
                  geom_blank(data = KT.meta,aes(Sample,Name))+
                  scale_fill_gradient2(low = "blue", mid = "grey",high = "red",midpoint = 0, limits = c(-4,4))+
                  theme_bw()+ xlab("Samples")+
                  theme(axis.text.x=element_blank(), axis.title.y=element_blank(), 
                        axis.title.x=element_text(family="Helvetica",size = 8),
                        axis.text.y=element_text(family="Helvetica",size = 8),
                        strip.background = element_blank(),
                        strip.text.y = element_blank(),
                        strip.text.x = element_text(family="Helvetica",size = 8),
                        legend.title = element_text(family="Helvetica",size = 8),
                        legend.text = element_text(family="Helvetica",size = 8))

withoutfilt = ggplot(final.withoutfilt.tab[final.withoutfilt.tab$Name %in% final.withfilt.tab$Name,],aes(Sample,Name))+
                  facet_grid(UTI~paste0(R.gender,"\n",substr(CI.genus,1,3)), scales = "free", space = "free")+
                  geom_tile(aes(fill=log10(Measurement)),col="black")+
                  geom_blank(data = KT.meta,aes(Sample,Name))+
                  scale_fill_gradient2(low = "blue", mid = "grey",high = "red",midpoint = 0, limits = c(-4,4))+
                  theme_bw()+ xlab("Samples")+
                  theme(axis.text.x=element_blank(), axis.title.y=element_blank(), 
                        axis.title.x=element_text(family="Helvetica",size = 8),
                        axis.text.y=element_text(family="Helvetica",size = 8),
                        strip.background = element_blank(),
                        strip.text.y = element_blank(),
                        strip.text.x = element_text(family="Helvetica",size = 8),
                        legend.title = element_text(family="Helvetica",size = 8),
                        legend.text = element_text(family="Helvetica",size = 8))

arrange.fig = ggarrange(withoutfilt, withfilt, common.legend = T, ncol = 2, nrow = 1) ;
save_eps =F
if (save_eps){
  pdf(file=paste0("/home/psb84/LBBCpaper.Fig1bc.061519.eps"),
      width=6.5,height=3.5, paper="special",bg="white",
      fonts="Helvetica", colormodel="cmyk", pointsize = 1)}

arrange.fig

if (save_eps){ dev.off()}
###


###
# show the total number of genera with and without filtering
length(unique(final.withfilt.tab$Name))
length(unique(final.withoutfilt.tab$Name))

######
# Here we add in some functions and steps for calling the confusion matrix.

# function for the confusion matrix.
call.it = function(sample, abundance.dataframe, included.taxa = c("Escherichia", "Enterococcus")){
  returned.df = subset.data.frame(abundance.dataframe,
                                  (Sample == sample)&(Name %in% included.taxa) ) ;
  if(nrow(returned.df) == 0){
    decided = "Negative" ;
  }else if(nrow(returned.df) == 1){
    decided = as.character(as.matrix(returned.df$Name)) ;
  }else if(nrow(returned.df) == 2){
    decided = "Both" ;
  }
  
  return(decided) ;
} 


KT.abundance = LoadAbundance(dir = path.grammy, file = "KTx.SMA.grammy.tab") ;
KT.abundance = subset.data.frame(KT.abundance, superkingdom == 2) ; # here we select only bacteria
KT.abundance = subset.data.frame(KT.abundance, RelCoverage >= 10**-9)
colnames(KT.abundance)[2] = "Sample"
KT.abundance$Measurement = KT.abundance$RelCoverage

final.withoutfilt.tab = DenoiseAlgorithm(AbundanceObject = KT.abundance,MetaDataObject = KT.meta,
                       NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = F, MassVar.Filter = F,NegCtrl.Filter = F,
                       TaxLevel = "genus",FastqPath = "./",TablePath = path.reads, TblatPath = path.tblat, AlnStatsPath = path.stats,
                       GITable = "lookups/gi_tax_info.tab")

meta.sub.df = subset.data.frame(KT.meta, select = c("Sample","CI.genus"))
colnames(meta.sub.df) = c("Sample","Reference")
levels(meta.sub.df$Reference) = c(levels(meta.sub.df$Reference), "Both")


meta.sub.df$Call = as.character(sapply(meta.sub.df$Sample, 
                                       function(x) call.it(sample = x, 
                                                           abundance.dataframe = final.withoutfilt.tab )))
levels(meta.sub.df$Call) = levels(meta.sub.df$Reference)
final.meta.df = data.frame(meta.sub.df)


#Confusion matrix
ref = factor(as.character(as.matrix(final.meta.df$Reference)),
             levels = levels(final.meta.df$Reference))
measure = factor(as.character(as.matrix(final.meta.df$Call)),
                 levels =  levels(final.meta.df$Call))

olc = confusionMatrix(measure, ref)

###------------------------------------------------------
# looking at the mass of contaminants

arm = DenoiseAlgorithm(AbundanceObject = KT.abundance,MetaDataObject = KT.meta,
                             NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                             CV.Filter = T, MassVar.Filter = T,NegCtrl.Filter = T,
                       deltaCV.Param = optimized_parameters$par[2],
                       MassVar.Param = optimized_parameters$par[3],
                             NegCtrl.Param = 10,TaxLevel = "genus",FastqPath = "./",TablePath = path.reads, TblatPath = path.tblat, AlnStatsPath = path.stats,
                       GITable = "lookups/gi_tax_info.tab")

step1.contam = merge(arm[1], list(Read.Abund.Matrix,KT.meta), by = "Sample")
step1.contam$taxamass = step1.contam$Biomass * step1.contam$AdjustedBlast / step1.contam$Total.Reads
aggsample.step1.contam = aggregate(taxamass ~ Sample,step1.contam, FUN = sum )
colnames(aggsample.step1.contam)[2] = "Stage1.Mass"

step2.contam = merge(arm[2], list(Read.Abund.Matrix,KT.meta), by = "Sample")
step2.contam$taxamass = step2.contam$Biomass * step2.contam$AdjustedBlast / step2.contam$Total.Reads
aggsample.step2.contam = aggregate(taxamass ~ Sample,step2.contam, FUN = sum )
colnames(aggsample.step2.contam)[2] = "Stage2.Mass"

step3.contam = merge(arm[3], list(Read.Abund.Matrix,KT.meta), by = "Sample")
step3.contam$taxamass = step3.contam$Biomass * step3.contam$AdjustedBlast / step3.contam$Total.Reads
aggsample.step3.contam = aggregate(taxamass ~ Sample,step3.contam, FUN = sum )
colnames(aggsample.step3.contam)[2] = "Stage3.Mass"

mdf = merge(aggsample.step1.contam, aggsample.step2.contam, by="Sample", all = T)
mdf = merge(mdf, aggsample.step3.contam, by="Sample", all = T)
mdf[is.na(mdf)] = 0

mdf$Physical.Mass = rowSums(subset.data.frame(mdf, select = c("Stage2.Mass","Stage3.Mass")))
mdf.clin = merge(mdf, KT.meta, by="Sample")



arm = DenoiseAlgorithm(AbundanceObject = KT.abundance,MetaDataObject = KT.meta,
                       NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = T, MassVar.Filter = F,NegCtrl.Filter = T,
                       CV.Param = optimized_parameters$par[1],
                       deltaCV.Param = optimized_parameters$par[2],
                       MassVar.Param = optimized_parameters$par[3],
                       NegCtrl.Param = 10,
                       TaxLevel = "genus", FastqPath = "/workdir/Data/KTx/BKVN/", 
                       TblatPath = "~/Grammy_backup/tblat_tables/",ReturnContams = T)


step3.contam = merge(arm[3], list(Read.Abund.Matrix,KT.meta), by = "Sample")
step3.contam$taxamass = step3.contam$Biomass * step3.contam$AdjustedBlast / step3.contam$Total.Reads
aggsample.step3.contam = aggregate(taxamass ~ Sample,step3.contam, FUN = sum )
colnames(aggsample.step3.contam)[2] = "Stage3.Mass"

mdf = merge(aggsample.step1.contam, aggsample.step2.contam, by="Sample", all = T)
mdf = merge(mdf, aggsample.step3.contam, by="Sample", all = T)
mdf[is.na(mdf)] = 0

mdf$Physical.Mass = rowSums(subset.data.frame(mdf, select = c("Stage2.Mass","Stage3.Mass")))
mdf.clin = merge(mdf, KT.meta, by="Sample")
