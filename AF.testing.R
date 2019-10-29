### Script to look at KTx select samples

### Libraries -----------
library(devtools); library(ggplot2); library(reshape2); library(MASS); library(taxize)
#install("/workdir/psb84/SparseMetagenomics/SparseMetagenomicCorrection")
library(SparseMetagenomicCorrection)


### Labels for groups -----
group_names <- list(
  'NormalNegative'="Full-term\nChorio Negative",
  'ChorioNegative'="Full-term\nChorio positive\nCulture negative",
  'PTLNegative'="Full-term\nChorio positive\nCulture positive",
  'PTLNegative'="Pre-term\nChorio positive\nCulture negative"
)

plot_labeller <- function(variable,value){
  if (variable=='facet1') {
    return(facet1_names[value])
  } else if (variable=='Status2') {
    return(group_names[value])
  } else {
    return(as.character(value))
  }
}



clinical.metadata = data.frame(read.table("/workdir/psb84/KTx/AFclinical.csv",header = T, sep = ","))
lab.metadata = data.frame(read.table("~/cfDNA_020519.csv",header = T, sep = ",",fill = T))
lab.metadata = lab.metadata[lab.metadata$Study_name == "Romero",]
colnames(lab.metadata)[1] = "Sample"
lab_meta_data = merge(lab.metadata,clinical.metadata, by="Sample")

#AF.subset.abundance$AdjustedMeasurement = AF.subset.abundance$Measurement * AF.subset.abundance$Biomass ;

# Read in metadata

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

# AF.meta = AddMetaData(MetaDataObject = AF.meta,
#                       ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Bacteria.genus")),
#                       ParameterName = "CI.genus") ;

# Generate statistics on the microbiome data

sample.list = as.character(as.matrix(AF.meta$Sample))
abc = c()
for (j in sample.list){
  abc = rbind(SparseMetagenomicCorrection::AlignmentStatisticsforSample(sample = j,bin.size = 10**3, file.align.dir = "~/Grammy_backup/tblat_tables/",tax.level = "species",file.out.dir = "/workdir/psb84/KTx/",
                                                                        file.abund.dir = "~/Grammy_backup/grammy_tables/"),abc)  
}

### Look at the characteristics


####
#### Pipeline ----------------
####




weights = c(4, 2,-1,-2,-0.25)

####-
####-
####-
####--------------------------


CV.max.threshold = 5;
simulCV.max.threshold = 2.5;
Coverage.min.threshold = 10**-4 ;
batch.var.max.threshold = 10**1;
weights.vector = weights;

#######-----------------

AF.abundance = LoadAbundance(dir = "/workdir/psb84/KTx/", file = "AF.grammy.tab")
AF.abundance = subset.data.frame(AF.abundance, superkingdom == 2)
AF.abundance = subset.data.frame(AF.abundance, RelCoverage >= Coverage.min.threshold)

species_super.df = data.frame(unique(subset.data.frame(AF.abundance,select = c("species","superkingdom"))))
colnames(species_super.df)[1] = "Tax.ID"

CV.filtered.tab = merge(abc,species_super.df, by ="Tax.ID")
CV.filtered.in.tab = subset.data.frame(CV.filtered.tab, (CV <= CV.max.threshold)&((CV-Simul.CV) <= simulCV.max.threshold))
CV.filtered.in.tab$sample_species = paste(CV.filtered.in.tab$Sample, CV.filtered.in.tab$Tax.ID,sep = "-")

### filter out from normal set



AF.abundance$sample_species = paste(AF.abundance$SAMPLE, AF.abundance$species,sep = "-")
AF.abundance = subset.data.frame(AF.abundance, (sample_species %in% as.character(as.matrix(CV.filtered.in.tab$sample_species))))

AF.subset.abundance = SubsetAbundance(AbundanceObject = AF.abundance, sample.column = 2,measurement = "RelCoverage", tax.level = "species")
#AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Bacteria.genus")),ParameterName = "CI.genus")
# AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Bacteria.species")),ParameterName = "CI.species")
AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Cohort2")),ParameterName = "Status")
# AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Rec.gender")),ParameterName = "R.gender")


### Ex-batch-aggregate
ex.batch.id.list = unique(as.character(as.matrix(AF.meta$EX.batch))) ;

ex.batch.group.correlates = BatchVariationCollect(AbundanceObject = AF.subset.abundance,AbundanceVariable = "Measurement",
                                                  MetaDataObject = AF.meta, BatchName = "EX.batch",BiomassInputAdjust = T,BiomassVariable = "Biomass",
                                                  BatchIDList = ex.batch.id.list ,Correlate.threshold = 0)

### Lp-batch-aggregate
lp.batch.id.list = unique(as.character(as.matrix(AF.meta$LP.batch)))

lp.batch.group.correlates = BatchVariationCollect(AbundanceObject = AF.subset.abundance,AbundanceVariable = "Measurement",
                                                  MetaDataObject = AF.meta, BatchName = "LP.batch",BiomassInputAdjust = T,BiomassVariable = "Biomass",
                                                  BatchIDList = lp.batch.id.list ,Correlate.threshold = 0)

ex.batch.exclude = subset.data.frame(ex.batch.group.correlates, Variance < batch.var.max.threshold, select = c("SampleTax"))
lp.batch.exclude = subset.data.frame(lp.batch.group.correlates, Variance < batch.var.max.threshold, select = c("SampleTax"))
bc.exclude = BiomassCorrelation(AbundanceObject = AF.subset.abundance, MetaDataObject = AF.meta, Cor.signif.threshold = 0.05)$`Rejected.Sample-Tax`

AF.subset.abundance$SampleTax = paste(AF.subset.abundance$Sample, AF.subset.abundance$Taxon,sep = "-")
AF.subset.abundance = subset.data.frame(AF.subset.abundance, !(SampleTax %in% as.character(as.matrix(ex.batch.exclude$SampleTax))))
AF.subset.abundance = subset.data.frame(AF.subset.abundance, !(SampleTax %in% as.character(as.matrix(lp.batch.exclude$SampleTax))))
AF.subset.abundance = subset.data.frame(AF.subset.abundance, !(SampleTax %in% as.character(as.matrix(bc.exclude))))

AF.subset.abundance$Name = ncbi_get_taxon_summary(AF.subset.abundance$Taxon)$name
Bacteria.subset.abundance = AF.subset.abundance

















#############################################################################
# Add in nonbacteria

batch.var.max.threshold = 10**-3;

AF.abundance = LoadAbundance(dir = "/workdir/psb84/KTx/", file = "AF.grammy.tab")
AF.abundance = subset.data.frame(AF.abundance, superkingdom == 2759)
AF.abundance = subset.data.frame(AF.abundance, RelCoverage >= Coverage.min.threshold)

species_super.df = data.frame(unique(subset.data.frame(AF.abundance,select = c("species","superkingdom"))))
colnames(species_super.df)[1] = "Tax.ID"

### filter out from normal set


AF.subset.abundance = SubsetAbundance(AbundanceObject = AF.abundance, sample.column = 2,measurement = "RelCoverage", tax.level = "species")
#AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Bacteria.genus")),ParameterName = "CI.genus")
# AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Bacteria.species")),ParameterName = "CI.species")
AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Cohort2")),ParameterName = "Status")
# AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Rec.gender")),ParameterName = "R.gender")


### Ex-batch-aggregate
ex.batch.id.list = unique(as.character(as.matrix(AF.meta$EX.batch))) ;

ex.batch.group.correlates = BatchVariationCollect(AbundanceObject = AF.subset.abundance,AbundanceVariable = "Measurement",
                                                  MetaDataObject = AF.meta, BatchName = "EX.batch",BiomassInputAdjust = T,BiomassVariable = "Biomass",
                                                  BatchIDList = ex.batch.id.list ,Correlate.threshold = 0)

### Lp-batch-aggregate
lp.batch.id.list = unique(as.character(as.matrix(AF.meta$LP.batch)))

lp.batch.group.correlates = BatchVariationCollect(AbundanceObject = AF.subset.abundance,AbundanceVariable = "Measurement",
                                                  MetaDataObject = AF.meta, BatchName = "LP.batch",BiomassInputAdjust = T,BiomassVariable = "Biomass",
                                                  BatchIDList = lp.batch.id.list ,Correlate.threshold = 0)

ex.batch.exclude = subset.data.frame(ex.batch.group.correlates, Variance < batch.var.max.threshold, select = c("SampleTax"))
lp.batch.exclude = subset.data.frame(lp.batch.group.correlates, Variance < batch.var.max.threshold, select = c("SampleTax"))
bc.exclude = BiomassCorrelation(AbundanceObject = AF.subset.abundance, MetaDataObject = AF.meta, Cor.signif.threshold = 0.05)$`Rejected.Sample-Tax`

AF.subset.abundance$SampleTax = paste(AF.subset.abundance$Sample, AF.subset.abundance$Taxon,sep = "-")
AF.subset.abundance = subset.data.frame(AF.subset.abundance, !(SampleTax %in% as.character(as.matrix(ex.batch.exclude$SampleTax))))
AF.subset.abundance = subset.data.frame(AF.subset.abundance, !(SampleTax %in% as.character(as.matrix(lp.batch.exclude$SampleTax))))
AF.subset.abundance = subset.data.frame(AF.subset.abundance, !(SampleTax %in% as.character(as.matrix(bc.exclude))))
AF.subset.abundance$Name = ncbi_get_taxon_summary(AF.subset.abundance$Taxon)$name


Euk.subset.abundance = AF.subset.abundance[AF.subset.abundance$Name == "Trichomonas vaginalis",]


########


batch.var.max.threshold = 10**1;

AF.abundance = LoadAbundance(dir = "/workdir/psb84/KTx/", file = "AF.grammy.tab")
AF.abundance = subset.data.frame(AF.abundance, superkingdom == 10239)
AF.abundance = subset.data.frame(AF.abundance, RelCoverage >= Coverage.min.threshold)

species_super.df = data.frame(unique(subset.data.frame(AF.abundance,select = c("species","superkingdom"))))
colnames(species_super.df)[1] = "Tax.ID"

### filter out from normal set


AF.subset.abundance = SubsetAbundance(AbundanceObject = AF.abundance, sample.column = 2,measurement = "RelCoverage", tax.level = "family")
#AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Bacteria.genus")),ParameterName = "CI.genus")
# AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Bacteria.species")),ParameterName = "CI.species")
AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Cohort2")),ParameterName = "Status")
# AF.subset.abundance = AddMetaData(AF.subset.abundance, ParameterFrame = subset.data.frame(clinical.metadata, select = c("Sample", "Rec.gender")),ParameterName = "R.gender")


### Ex-batch-aggregate
ex.batch.id.list = unique(as.character(as.matrix(AF.meta$EX.batch))) ;

ex.batch.group.correlates = BatchVariationCollect(AbundanceObject = AF.subset.abundance,AbundanceVariable = "Measurement",
                                                  MetaDataObject = AF.meta, BatchName = "EX.batch",BiomassInputAdjust = T,BiomassVariable = "Biomass",
                                                  BatchIDList = ex.batch.id.list ,Correlate.threshold = 0)

### Lp-batch-aggregate
lp.batch.id.list = unique(as.character(as.matrix(AF.meta$LP.batch)))

lp.batch.group.correlates = BatchVariationCollect(AbundanceObject = AF.subset.abundance,AbundanceVariable = "Measurement",
                                                  MetaDataObject = AF.meta, BatchName = "LP.batch",BiomassInputAdjust = T,BiomassVariable = "Biomass",
                                                  BatchIDList = lp.batch.id.list ,Correlate.threshold = 0)

ex.batch.exclude = subset.data.frame(ex.batch.group.correlates, Variance < batch.var.max.threshold, select = c("SampleTax"))
lp.batch.exclude = subset.data.frame(lp.batch.group.correlates, Variance < batch.var.max.threshold, select = c("SampleTax"))
bc.exclude = BiomassCorrelation(AbundanceObject = AF.subset.abundance, MetaDataObject = AF.meta, Cor.signif.threshold = 0.05)$`Rejected.Sample-Tax`

AF.subset.abundance$SampleTax = paste(AF.subset.abundance$Sample, AF.subset.abundance$Taxon,sep = "-")
AF.subset.abundance = subset.data.frame(AF.subset.abundance, !(SampleTax %in% as.character(as.matrix(ex.batch.exclude$SampleTax))))
AF.subset.abundance = subset.data.frame(AF.subset.abundance, !(SampleTax %in% as.character(as.matrix(lp.batch.exclude$SampleTax))))
AF.subset.abundance = subset.data.frame(AF.subset.abundance, !(SampleTax %in% as.character(as.matrix(bc.exclude))))
AF.subset.abundance$Name = ncbi_get_taxon_summary(AF.subset.abundance$Taxon)$name


Virus.subset.abundance = AF.subset.abundance

Bacteria.subset.abundance$Superkingdom = "Bacteria"
Euk.subset.abundance$Superkingdom = "Eukaryote"
Virus.subset.abundance$Superkingdom = "Virus"

AF.subset.abundance = rbind(Bacteria.subset.abundance, Euk.subset.abundance, Virus.subset.abundance)


AF.meta$Taxon = AF.subset.abundance$Taxon[1]
AF.meta$Name = ncbi_get_taxon_summary(AF.meta$Taxon)$name
AF.meta$Superkingdom = "Bacteria"

AF.meta$Status2 = factor(AF.meta$Status, levels=c("NormalNegative","ChorioNegative","ChorioPositive","PTLNegative"), 
                         labels=c("NormalNegative","ChorioNegative","ChorioPositive","PTLNegative")) 
AF.subset.abundance$Status2 = factor(AF.subset.abundance$Status, levels=c("NormalNegative","ChorioNegative","ChorioPositive","PTLNegative"), 
                                     labels=c("NormalNegative","ChorioNegative","ChorioPositive","PTLNegative")) 

clin.tmp.1 = subset.data.frame(clinical.metadata, select = c("Sample","IBIS_Bacteria_1"))
clin.tmp.1 = clin.tmp.1[clin.tmp.1$IBIS_Bacteria_1 != "",]
clin.tmp.2 = subset.data.frame(clinical.metadata, select = c("Sample","IBIS_Bacteria_2"))
clin.tmp.2 = clin.tmp.2[clin.tmp.2$IBIS_Bacteria_2 != "",]
clin.tmp.3 = subset.data.frame(clinical.metadata, select = c("Sample","IBIS_Euk_1"))
clin.tmp.3 = clin.tmp.3[clin.tmp.3$IBIS_Euk_1 != "",]
colnames(clin.tmp.1) = colnames(clin.tmp.2) = colnames(clin.tmp.3) = c("Sample","Name") 

clinical.bact.calls = rbind(clin.tmp.1, clin.tmp.2)
clinical.bact.calls$Superkingdom = "Bacteria"
clin.tmp.3$Superkingdom = "Eukaryote"
clinical.calls = merge(rbind(clinical.bact.calls, clin.tmp.3), subset.data.frame(AF.meta,select = c(Sample,Status2)))


clinical.calls$Measurement = 10


save_eps=T
if (save_eps){  
  pdf(file=paste0("/home/psb84/AF.micriobiome.array.020519.eps"),
      width=16.5/2.54,height=10/2.54, paper="special",bg="white",
      fonts="Helvetica", colormodel="cmyk", pointsize = 1)}

ggplot(AF.subset.abundance,aes(Sample,Name))+
  facet_grid(Superkingdom~Status2, scales = "free", space = "free", labeller = plot_labeller)+
  geom_tile(aes(fill=log10(Measurement)),col="black")+
  geom_point(data= clinical.calls,aes(Sample,Name),shape=13,size=2)+
  geom_blank(data=AF.meta)+xlab("Samples")+
  scale_fill_gradient2(low = "lightblue", mid = "darkgrey",high = "red",midpoint = 0)+
  theme_bw()+ 
  theme(legend.position = "none",axis.text.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.title.x=element_text(family="Helvetica",size = 8), 
        axis.text.y=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        strip.text = element_text(family="Helvetica",size = 8))

if(save_eps){dev.off()}

species.of.interest = AF.subset.abundance$Name

write.table(x = AF.subset.abundance,
            file = "/home/psb84/AF.afterfiltering.tab",
            sep="\t",quote = F, row.names = F,col.names = T)
