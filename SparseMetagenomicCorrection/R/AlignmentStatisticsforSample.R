#' AlignmentStatisticsforSample
#' This function generates alignment information and statistics for all microbes above threshold in a sample.
#' It will also print this information to a file.
#' @param Sample Abundance matrix to be loaded. Must be filtered to tri-column form.
#' @param file.align.dir Directory including file
#' @param file.align.extension extension of alignment file. DEFAULT is *.tblat.1
#' @param file.abund.dir Directory including file
#' @param file.abund.extension extension of abundance file. DEFAULT is *.grammy.tab
#' @param file.out.dir Directory including file
#' @param file.out.extension extension of abundance file. DEFAULT is *.microbes_align_stats.tab
#' @param bin.size Size of bin in bp. Default is 1000 bp.
#' @param tax.level ncbi taxification level. "species" by DEFAULT.
#' @param Sample Abundance matrix to be loaded. Must be filtered to tri-column form.
#' @keywords stats
#' @importFrom data.table fread
#' @export
#' @examples
#' AlignmentStatisticsforSample()
#'
#


AlignmentStatisticsforSample = function(sample, bin.size = 10**3,
                           tax.level = "species",
                           file.align.dir,
                           file.align.extension = ".tblat.1",
                           file.abund.dir,
                           file.abund.extension = ".grammy.tab",
                           file.out.dir = "./",
                           file.out.extension = ".microbes_align_stats.tab"){

  if(!(file.exists(paste0(file.align.dir,sample,file.align.extension)))){
    stop("No alignment file present.")
  }
  if(!(file.exists(paste0(file.abund.dir,sample,file.abund.extension)))){
    stop("No abundance file present.")
  }
  if((file.exists(paste0(file.out.dir,sample,".",bin.size,"bp.",tax.level,".genus",file.out.extension)))){
    species.finalstats.df = data.frame(fread(paste0(file.out.dir,sample,".",bin.size,"bp.",tax.level,".genus",file.out.extension), header = T, sep = "\t"))
  }else{

        # import alignment table
        align.df = data.frame(fread(paste0(file.align.dir,sample,file.align.extension), sep ="\t", header = F))
        align.names = c("ReadPaired","GI","PercIdent","ReadLength","Mismatch","GapOpen","QueryStart","QueryEnd",
                        "SeqStart","SeqEnd","EValue","BitScore","Taxid","TaxName","TaxKingdom")
        colnames(align.df) = align.names
        align.df$ReadID = substr(align.df$ReadPaired,1,nchar(align.df$ReadPaired)-2)

        # import abundance table
        abund.df =  data.frame(fread(paste0(file.abund.dir,sample,file.abund.extension), sep ="\t", header = T))

        # merge grammy and tblat.1
        merged.df = merge(align.df,abund.df, by="Taxid",all.x=T)
        merged.df = merged.df[!is.na(merged.df$species),]

        # remove species for which there are fewer than 10 hits, and list all species.
        species.count.df = data.frame(table(merged.df[[tax.level]]))
        colnames(species.count.df) = c("species","count")
        species.list = as.numeric(as.matrix(species.count.df[species.count.df$count >= 10,]$species))
        species.finalstats.df=c()

        if(length(species.list) == 0){
          stop("No species exceeding minimum number (10) alignments.")
        }else{
          species.finalstats.df=c()
          for (i in species.list){
            species.finalstats.df = rbind(species.finalstats.df,
                                          AlignmentStatisticsforTax(sample=sample,
                                                                    tax.id = i,
                                                                    tax.level = tax.level,
                                                                    bin.size = bin.size,
                                                                    Merge.data.frame = merged.df))}
        }
        write.table(x = species.finalstats.df,
                    file =paste0(file.out.dir,sample,".",bin.size,"bp.",tax.level,".genus",file.out.extension),
                    quote = F, sep="\t",row.names = F, col.names = T) ;
  }
  return(species.finalstats.df) ;
}
