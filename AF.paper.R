### Script for training

### Libraries -----------
library(devtools); library(ggplot2); library(reshape2);
library(MASS); library(taxize); library(ggpubr) ;
install("/workdir/psb84/SparseMetagenomics/SparseMetagenomicCorrection")
library(SparseMetagenomicCorrection)

### Load metadata matrices ------------------------------------------------------------------------

clinical.metadata = data.frame(read.table("/workdir/psb84/KTx/AFclinical.csv",header = T, sep = ","))
lab.metadata = data.frame(read.table("~/cfDNA_020519.csv",header = T, sep = ",",fill = T))
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
                                  TotalReadsOutput = "/workdir/psb84/SparseMetagenomics/AFmeta.totalreads.tab",
                                  RawDataPath = "/workdir/Data/AF/") ;

### Load taxa from negative controls --------------------------------------------------------------
negatives = SetNegativeControl(sample.vector = c(paste0(LETTERS[c(1,3:4)],"M.clip75.AF"),
                                                 "BM.AF","AM6.clip75.AF","AM5.clip75.AF"),
                               raw.data.path = "/workdir/Data/AF/",
                               tblat.path = "~/Grammy_backup/tblat_tables/")

### Denoise the abundance matrix ------------------------------------------------------------------
Coverage.min.threshold = 10**-9
AF.abundance = LoadAbundance(dir = "/workdir/psb84/KTx/", file = "AF.grammy.tab")
AF.abundance = subset.data.frame(AF.abundance, superkingdom == 2)
AF.abundance = subset.data.frame(AF.abundance, RelCoverage >= Coverage.min.threshold)
colnames(AF.abundance)[2] = "Sample"
AF.abundance$Measurement = AF.abundance$RelCoverage

abc = DenoiseAlgorithm(AbundanceObject = AF.abundance, MetaDataObject = AF.meta,
                       NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = T, MassVar.Filter = T,NegCtrl.Filter = T,
                       deltaCV.Param = 2, MassVar.Param = -6, NegCtrl.Param = 10,
                       TaxLevel = "species", FastqPath = "/workdir/Data/AF/", 
                       TblatPath = "~/Grammy_backup/tblat_tables/")
final.withfilt.tab = merge(abc, AF.meta, "Sample")
AF.meta$Name = final.withfilt.tab$Name[1]


### Plot  -----------------------------------------------------------------------------------------


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

####
Coverage.min.threshold = 10**-9
AF.abundance = LoadAbundance(dir = "/workdir/psb84/KTx/", file = "AF.grammy.tab")
AF.abundance = subset.data.frame(AF.abundance, superkingdom == 10239)
AF.abundance = subset.data.frame(AF.abundance, RelCoverage >= Coverage.min.threshold)
colnames(AF.abundance)[2] = "Sample"
AF.abundance$Measurement = AF.abundance$RelCoverage

abc = DenoiseAlgorithm(AbundanceObject = AF.abundance,MetaDataObject = AF.meta,
                       NegativeObject = negatives,ReadAbundMatrix = Read.Abund.Matrix,
                       CV.Filter = F, MassVar.Filter = T,NegCtrl.Filter = T,
                       deltaCV.Param = 2, MassVar.Param = -9, NegCtrl.Param = 10,
                       TaxLevel = "family", FastqPath = "/workdir/Data/AF/", 
                       TblatPath = "~/Grammy_backup/tblat_tables/")
final.withfilt.virus.tab = merge(abc, AF.meta[,1:6], "Sample")
final.withfilt.virus.tab$Status2 = factor(final.withfilt.virus.tab$Status, levels=c("NormalNegative","ChorioNegative","ChorioPositive"), 
                                          labels=c("NormalNegative","ChorioNegative","ChorioPositive")) 

final.withfilt.tab$SK  = "Bacteria"
final.withfilt.virus.tab$SK = "Virus"
AF.meta$SK = "Bacteria"
clinical.calls$SK = "Bacteria"

all.final = rbind(final.withfilt.tab, final.withfilt.virus.tab)

fig.bacteria = ggplot(all.final,aes(Sample,Name))+
  facet_grid(SK~Status2, scales = "free", space = "free")+
  geom_tile(aes(fill=log10(Measurement)),col="black")+
  geom_blank(data = AF.meta,aes(Sample,Name))+
  geom_point(data= clinical.calls,aes(Sample,Name),shape=13,size=2)+
  scale_fill_gradient2(low = "blue", mid = "grey",high = "red",midpoint = 0, limits = c(-4,4))+
  theme_bw()+ xlab("Samples")+
  theme(legend.position = "none",
        axis.text.x=element_blank(), axis.title.y=element_blank(), 
        axis.title.x=element_text(family="Helvetica",size = 8),
        axis.text.y=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        strip.text = element_text(family="Helvetica",size = 8),
        legend.title = element_text(family="Helvetica",size = 8),
        legend.text = element_text(family="Helvetica",size = 8))


save_eps=T
if (save_eps){
  pdf(file=paste0("/home/psb84/LBBCpaper.Fig2.061519.eps"),
      width=6.5,height=2, paper="special",bg="white",
      fonts="Helvetica", colormodel="cmyk", pointsize = 1)}

fig.bacteria

if (save_eps){ dev.off()}

###
