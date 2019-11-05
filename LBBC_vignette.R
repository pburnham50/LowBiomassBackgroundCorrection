### Vignette for Low Biomass Background Correction method.
# Filename: LBBC_vignette.R
# Author: Philip Burnham
# Year: 2019

### Libraries ------------------------------------------------------------------------------------------------
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

### Parameters -----------------------------------------------------------------------------------------------

# print plot in window by setting to FALSE
export.eps = F

# set parameters
deltaCV.maximum = 2 # controls filtering based on genome sequencing homogeneity
Batch.var.log.min = -5.5 # controls filtering based on batch effects
Negative.ctrl.coef.max = 10 # controls filtering based on comparison to negative control

# determine taxonomic level for aggregation
tax.level = "genus" # other useful options are "family" and "species"

### Load abundance matrix ------------------------------------------------------------------------------------

KT.abundance = LoadAbundance(dir = path.grammy, file = "KTx.SMA.grammy.tab") ;
KT.abundance = subset.data.frame(KT.abundance, superkingdom == 2) ; # here we select only bacteria

colnames(KT.abundance)[2] = "Sample" ;
KT.abundance$Measurement = KT.abundance$RelCoverage ; # setting the measurement of interest to the relative genomic coverage.

### Load metadata matrices -----------------------------------------------------------------------------------
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
                      ParameterFrame = subset.data.frame(clinical.metadata,
                                                         select = c("Sample", "Rec.gender")),
                      ParameterName = "R.gender") ;

KT.meta = AddMetaData(MetaDataObject = KT.meta,
                      ParameterFrame = subset.data.frame(clinical.metadata,
                                                         select = c("Sample", "Status")),
                      ParameterName = "Status") ;

KT.meta = AddMetaData(MetaDataObject = KT.meta,
                      ParameterFrame = subset.data.frame(clinical.metadata,
                                                         select = c("Sample", "Bacteria.genus")),
                      ParameterName = "CI.genus") ;

### Load abundance matrix ------------------------------------------------------------------------------------
Read.Abund.Matrix = TotalReadsGen(KT.meta,
                                  TotalReadsOutput = paste0(path.reads,"KTmeta.totalreads.tab"),
                                  RawDataPath = "./") ; # raw data path is where fastqs are stored. not necessart for this.

colnames(KT.abundance)[2] = "Sample" ;

### Load taxa from negative controls --------------------------------------------------------------
negatives = SetNegativeControl(sample.vector = paste0("MC",LETTERS[c(1:9,11:14,16:20)]),
                               raw.data.path = "./",table.path = path.reads,
                               tblat.path = path.tblat)



### Denoise the abundance matrix ------------------------------------------------------------------
with_all_filters.tab = DenoiseAlgorithm(AbundanceObject = KT.abundance, MetaDataObject = KT.meta,
                       NegativeObject = negatives, ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = T, MassVar.Filter = T, NegCtrl.Filter = T,TaxLevel = tax.level,
                       deltaCV.Param = deltaCV.maximum, MassVar.Param = Batch.var.log.min,NegCtrl.Param = Negative.ctrl.coef.max,
                        FastqPath = "./",TablePath = path.reads, TblatPath = path.tblat, AlnStatsPath = path.stats,
                       GITable = "lookups/gi_tax_info.tab")

final.withfilt.tab = merge(with_all_filters.tab, KT.meta, "Sample")


# we can also look at what happens if you don't apply any filter
with_no_filters.tab = DenoiseAlgorithm(AbundanceObject = KT.abundance,MetaDataObject = KT.meta,
                       NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = F, MassVar.Filter = F,NegCtrl.Filter = F,
                       TaxLevel = tax.level, GITable = "lookups/gi_tax_info.tab")

final.withoutfilt.tab = merge(with_no_filters.tab, KT.meta, "Sample")

### Plot  -----------------------------------------------------------------------------------------
# here we plot a comparison of the resulting filtered microbiome across samples, with and without filtering.

final.withfilt.tab$UTI = !(final.withfilt.tab$Name %in% c("Escherichia", "Enterococcus"))
final.withoutfilt.tab$UTI = !(final.withoutfilt.tab$Name %in% c("Escherichia", "Enterococcus"))
KT.meta$Name = final.withfilt.tab$Name[1]
KT.meta$UTI = F

high_disease = paste0("GU",c(74,81,83,467,1209,5,14,65,69,1125,1495,2709))

withfilt = ggplot(final.withfilt.tab,aes(Sample,Name))+
                  facet_grid(UTI~paste0(R.gender ,"\n",substr(CI.genus,1,3)),
                             scales = "free", space = "free")+
                  geom_tile(aes(fill=log10(Measurement)),col="black")+
                  geom_blank(data = KT.meta,aes(Sample,Name))+
                  scale_fill_gradient2(low = "blue", mid = "grey",high = "red",
                                       midpoint = 0, limits = c(-4,4))+
                  theme_bw()+ xlab("Samples -- AFTER")+
                  theme(axis.text.x=element_blank(), axis.title.y=element_blank(),
                        axis.title.x=element_text(family="Helvetica",size = 8),
                        axis.text.y=element_text(family="Helvetica",size = 8),
                        strip.background = element_blank(),
                        strip.text.y = element_blank(),
                        strip.text.x = element_text(family="Helvetica",size = 8),
                        legend.title = element_text(family="Helvetica",size = 8),
                        legend.text = element_text(family="Helvetica",size = 8))

withoutfilt = ggplot(final.withoutfilt.tab[final.withoutfilt.tab$Name %in% final.withfilt.tab$Name,],aes(Sample,Name))+
                  facet_grid(UTI~paste0(R.gender,"\n",substr(CI.genus,1,3)),
                             scales = "free", space = "free")+
                  geom_tile(aes(fill=log10(Measurement)),col="black")+
                  geom_blank(data = KT.meta,aes(Sample,Name))+
                  scale_fill_gradient2(low = "blue", mid = "grey",high = "red",
                                       midpoint = 0, limits = c(-4,4))+
                  theme_bw()+ xlab("Samples -- BEFORE")+
                  theme(axis.text.x=element_blank(), axis.title.y=element_blank(),
                        axis.title.x=element_text(family="Helvetica",size = 8),
                        axis.text.y=element_text(family="Helvetica",size = 8),
                        strip.background = element_blank(),
                        strip.text.y = element_blank(),
                        strip.text.x = element_text(family="Helvetica",size = 8),
                        legend.title = element_text(family="Helvetica",size = 8),
                        legend.text = element_text(family="Helvetica",size = 8))

arrange.fig = ggarrange(withoutfilt, withfilt, common.legend = T, ncol = 2, nrow = 1) ;


if (export.eps){
  pdf(file=paste0(path.figs, "LBBC_vignette.eps"),
      width=6.5,height=3.5, paper="special",bg="white",
      fonts="Helvetica", colormodel="cmyk", pointsize = 1)}

arrange.fig
if (export.eps){ dev.off()}
