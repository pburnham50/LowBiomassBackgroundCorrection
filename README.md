# Low Biomass Background Correction

Author: Philip Burnham (phil.burnham.50@gmail.com)

Link to preprint: https://www.biorxiv.org/content/10.1101/734756v1

## Outline

The following is a guide to reduce the background contamination of microbial reads in a cell-free DNA sequencing set. We will show how to use the pipeline and its functions to denoise a dataset from urine samples taken from patients who had received kidney transplants at the New York Presbyterian Hospital in NYC, USA. The data includes samples with a matched same day diagnosis of either E. coli or Enterococcus urinary tract infection. The dataset also includes patients who had no observable bacteria in the urinary tract, and who had no UTI during the time they were monitored (this is our control group).

The pipeline is implemented in R, and takes advantage of two observations: 1) the variation in coverage across microbial genomes, 2) how the abundance of bacteria changes with the total mass of DNA added to library preparation. In (1), microbial genomes with high coverage variability (measured by the coefficient of variation in coverage) are interpreted as digital noise. In (2), if the representation of a microbe in a batch decreases with increasing amounts of input biomass, it is likely a contaminant, presenting itself at low abundance.

## Procedure

### Step 0: Collecting the needed files

We have implemented a bioinformatics pipeline that: 1) removes low quality DNA reads, 2) aligns reads to the human reference genome (UCSC hg19 build), 3) aligns human-unmapped reads to a large microbial genome reference database (NCBI blast), and (4) estimates the abundance of bacteria based on alignment statistics. This pipeline has been used extensively by our lab in a variety of contexts including plasma, urine, amniotic fluid, and peritoneal dialysis effluent (paper coming).

From this pipeline we collect two important files:

The alignment statistics of each nonhuman read to a microbe using NCBI BLAST (\*.tblat.1).

The phylogeny and abundance table that shows the microbiome estimated from all nonhuman reads (\*.grammy.tab).


For the purposes of this tutorial, we avoid providing raw FASTQ files for privacy and storage concerns. In practice, the LBBC algorithm will calculate the number of sequencing reads to look at batch covariation. This can be bypassed by directly providing a \*.tab file with the number of sequencing reads.

We also take note of the batches that our samples are in, and the estimated biomass (in ng) of the cell-free DNA we are using in our library preparation. This is our metadata file.

*If you are missing the .tblat.1 file, don’t worry! If you didn’t measure your DNA before the library preparation, don’t worry! The pipeline can be altered to ignore these steps (though the corresponding step in filtering will be ignored).*

Start by cloning the LBBC GitHub repository to your local path:

```
git clone https://github.com/pburnham50/LowBiomassBackgroundCorrection
```

### Step 1: Prepare required files

In this tutorial we provide the files used to generate the figures in the linked preprint.

The grammy file (table of all microbes in all samples, and their relative abundances) is given as 'grammy/KTx.SMA.grammy.tab'.
In general, if you don’t have this at this point, copy all \*.grammy.tab files into one folder and run the following from the command line from that directory:

```
$ cat *.grammy.tab | sort | uniq > /path/to/LBBC/Project.grammy.tab ;
```

We also need the tblat files (genome positions of all microbial reads). These files are too large to be hosted on this repository but can be downloaded from the De Vlaminck lab's open access Dropbox storage using the following:

```
$ wget ... tblats/ ;
```

### Step 2: Installing the LBBC package.

Start a new R session and open up the file LBBC_vignette.R in RStudio.

Set working directory to the top-level of the LBBC directory.

```
> setwd("path/to/LBBC/") ;
```

Install the following packages if not already present:
ineq, ggplot2, ggpubr, roxygen2, MASS, devtools, reshape2, taxize.

In your R session this can be achieved by running the following line of code:

```
> source('extra/load_packages.R') ;
```

### Step 3: Running the LBBC package on a urinary cell-free DNA dataset

At this point you are ready to follow along the R script to properly format your data and run the workhorse function of the LBBC package "DenoiseAlgorithm".

Running the full script will produce a plot (similar to Fig 1b in the preprint) of the side-by-side comparison of the filtered and unfiltered datasets.

Explore the effects of varying the three filtering parameters by changing their value at the beginning of the script.

To reproduce the filtered two microbiome arrays (urine from kidney transplant recipients and amniotic fluid from a cohort of women suffering infection/inflammation), run the from the command line (for urine):

```
$ Rscript bin/R/LBBC.KTx.R
```

Or, for amniotic fluid:
```
$ Rscript bin/R/LBBC.AF.R
```
The figure files will appear in the directory "figs/" in .eps format.


### Acknowledgments

Iwijn De Vlaminck (Cornell University) - development of methodology.

Alexandre Pellan Cheng (Cornell University) - software validation.
