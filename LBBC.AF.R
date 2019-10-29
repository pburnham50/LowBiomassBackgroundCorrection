### Reproducing amniotic fluid filtered microbiome abundance
# Filename: LBBC.AF.R
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

export.eps = T

# set parameters (change for bacteria and virus)
bact.deltaCV.maximum = 2 
bact.Batch.var.log.min = -6 
bact.Negative.ctrl.coef.max = 10 
virus.deltaCV.maximum = 2 # though we will turn this off.
virus.Batch.var.log.min = -9 
virus.Negative.ctrl.coef.max = 10


# determine taxonomic level for aggregation (different for bacteria and virus)
bact.tax.level = "species" 
virus.tax.level = "family" 

Coverage.min.threshold = 10**-9

### Load metadata matrices ------------------------------------------------------------------------

clinical.metadata = data.frame(read.table(paste0(path.metadata,"AFclinical.csv"),header = T, sep = ","))
lab.metadata = data.frame(read.table(paste0(path.metadata,"cfDNA_020519.csv"),header = T, sep = ",",fill = T))
lab.metadata = lab.metadata[lab.metadata$Study_name == "Romero",]
colnames(lab.metadata)[1] = "Sample"
lab_meta_data = merge(lab.metadata,clinical.metadata, by="Sample")


AF.meta = InitializeMetaData(SampleVector = as.character(as.matrix(lab_meta_data$Sample)),
                             BatchVector = paste0(as.character(as.matrix(lab_meta_data$Date_extraction)),"_",
                                                  as.character(as.matrix(lab_meta_data$Batch_extraction))),
                             BatchName = "EX.batch") ;

lp.meta = subset.data.frame(lab.metadata,select = c("Sample","Date_libprep","Batch_libprep"))
lp.meta$Batch = paste0(lp.meta$Date_libprep,"_",lp.meta$Batch_libprep)

AF.meta = AddMetaData(MetaDataObject = AF.meta,
                      BatchFrame = lp.meta[,c(1,4)],
                      BatchName = "LP.batch") ;

AF.meta = AddMetaData(MetaDataObject = AF.meta,
                      ParameterFrame = subset.data.frame(lab_meta_data, select = c("Sample", "Biomass_at_libprep_ng")),
                      ParameterName = "Biomass") ;

AF.meta = AddMetaData(MetaDataObject = AF.meta,
                      ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Sex")),
                      ParameterName = "Sex") ;

AF.meta = AddMetaData(MetaDataObject = AF.meta,
                      ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Cohort2")),
                      ParameterName = "Status") ;
AF.meta = AF.meta[AF.meta$Status != "PTLNegative",]

### Load abundance matrix -------------------------------------------------------------------------
Read.Abund.Matrix = TotalReadsGen(AF.meta, 
                                  TotalReadsOutput = paste0(path.reads,"AFmeta.totalreads.tab"),
                                  RawDataPath = "./") ;

### Load taxa from negative controls --------------------------------------------------------------
negatives = SetNegativeControl(sample.vector = c(paste0(LETTERS[c(1,3:4)],"M.clip75.AF"),
                                                 "BM.AF","AM6.clip75.AF","AM5.clip75.AF"),
                               raw.data.path = "./",table.path = path.reads,
                               tblat.path = path.tblat)

### Denoise the abundance matrix ------------------------------------------------------------------
AF.abundance = LoadAbundance(dir = path.grammy, file = "AF.grammy.tab")
AF.bacteria.abundance = subset.data.frame(AF.abundance, superkingdom == 2)
AF.bacteria.abundance = subset.data.frame(AF.bacteria.abundance, RelCoverage >= Coverage.min.threshold)
colnames(AF.bacteria.abundance)[2] = "Sample"
AF.bacteria.abundance$Measurement = AF.bacteria.abundance$RelCoverage

bacteria = DenoiseAlgorithm(AbundanceObject = AF.bacteria.abundance, MetaDataObject = AF.meta,
                       NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = T, MassVar.Filter = T,NegCtrl.Filter = T,
                       deltaCV.Param = bact.deltaCV.maximum, MassVar.Param = bact.Batch.var.log.min, 
                       NegCtrl.Param = bact.Negative.ctrl.coef.max, TaxLevel = bact.tax.level,
                       FastqPath = "./",TablePath = path.reads, TblatPath = path.tblat, AlnStatsPath = path.stats,
                       GITable = "lookups/gi_tax_info.tab")

final.withfilt.tab = merge(bacteria, AF.meta, "Sample")
AF.meta$Name = final.withfilt.tab$Name[1]
AF.meta$Status2 = factor(AF.meta$Status, levels=c("NormalNegative","ChorioNegative","ChorioPositive"), 
                         labels=c("NormalNegative","ChorioNegative","ChorioPositive")) 



# now doing the same process for virus. Here we turn off the CV filter due to the size of viral genomes.

AF.virus.abundance = subset.data.frame(AF.abundance, superkingdom == 10239)
AF.virus.abundance = subset.data.frame(AF.virus.abundance, RelCoverage >= Coverage.min.threshold)
colnames(AF.virus.abundance)[2] = "Sample"
AF.virus.abundance$Measurement = AF.virus.abundance$RelCoverage

virus = DenoiseAlgorithm(AbundanceObject = AF.virus.abundance, MetaDataObject = AF.meta,
                       NegativeObject = negatives, ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = F, MassVar.Filter = T,NegCtrl.Filter = T,
                       deltaCV.Param = virus.deltaCV.maximum, MassVar.Param = virus.Batch.var.log.min, 
                       NegCtrl.Param = virus.Negative.ctrl.coef.max, TaxLevel = virus.tax.level,
                       FastqPath = "./",TablePath = path.reads, TblatPath = path.tblat, AlnStatsPath = path.stats,
                       GITable = "lookups/gi_tax_info.tab")

final.withfilt.virus.tab = merge(virus, AF.meta[,1:6], "Sample")


### Plot  -----------------------------------------------------------------------------------------

AF.meta$Name = final.withfilt.tab$Name[1]

clin.tmp.1 = subset.data.frame(clinical.metadata, select = c("Sample","IBIS_Bacteria_1"))
clin.tmp.1 = clin.tmp.1[clin.tmp.1$IBIS_Bacteria_1 != "",]
clin.tmp.2 = subset.data.frame(clinical.metadata, select = c("Sample","IBIS_Bacteria_2"))
clin.tmp.2 = clin.tmp.2[clin.tmp.2$IBIS_Bacteria_2 != "",]
colnames(clin.tmp.1) = colnames(clin.tmp.2)  = c("Sample","Name") 

clinical.bact.calls = rbind(clin.tmp.1, clin.tmp.2)
clinical.bact.calls$Superkingdom = "Bacteria"
clinical.calls = merge(clinical.bact.calls, subset.data.frame(AF.meta,select = c(Sample,Status)))


clinical.calls$Measurement = 10


AF.meta$Status2 = factor(AF.meta$Status, levels=c("NormalNegative","ChorioNegative","ChorioPositive"), 
                         labels=c("NormalNegative","ChorioNegative","ChorioPositive")) 
final.withfilt.tab$Status2 = factor(final.withfilt.tab$Status, levels=c("NormalNegative","ChorioNegative","ChorioPositive"), 
                                    labels=c("NormalNegative","ChorioNegative","ChorioPositive")) 
clinical.calls$Status2 = factor(clinical.calls$Status, levels=c("NormalNegative","ChorioNegative","ChorioPositive"), 
                                labels=c("NormalNegative","ChorioNegative","ChorioPositive")) 
final.withfilt.virus.tab$Status2 = factor(final.withfilt.virus.tab$Status, levels=c("NormalNegative","ChorioNegative","ChorioPositive"), 
                                          labels=c("NormalNegative","ChorioNegative","ChorioPositive")) 

final.withfilt.tab$SK  = "Bacteria"
final.withfilt.virus.tab$SK = "Virus"
AF.meta$SK = "Bacteria"
clinical.calls$SK = "Bacteria"

all.final = rbind(final.withfilt.tab, final.withfilt.virus.tab)

fig = ggplot(all.final,aes(Sample,Name))+
  facet_grid(SK~Status2, scales = "free", space = "free")+
  geom_tile(aes(fill=log10(Measurement)),col="black")+
  geom_blank(data = AF.meta,aes(Sample,Name))+
  geom_point(data= clinical.calls,aes(Sample,Name),shape=13,size=2)+
  scale_fill_gradient2(low = "blue", mid = "grey",high = "red",midpoint = 0, limits = c(-4,4))+
  theme_bw()+ xlab("Samples")+
  theme(legend.position = "top",
        axis.text.x=element_blank(), axis.title.y=element_blank(), 
        axis.title.x=element_text(family="Helvetica",size = 8),
        axis.text.y=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        strip.text = element_text(family="Helvetica",size = 8),
        legend.title = element_text(family="Helvetica",size = 8),
        legend.text = element_text(family="Helvetica",size = 8))


if (export.eps){
  pdf(file=paste0(path.figs,"LBBC.AF.eps"),
      width=6.5,height=2, paper="special",bg="white",
      fonts="Helvetica", colormodel="cmyk", pointsize = 1)}

fig

if (export.eps){ dev.off()}

###
