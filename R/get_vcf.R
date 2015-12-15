#' Get data for a genomic region from a remote VCF file.
#' 
#' Returns a list with three dataframes for individuals, SNPs, and genotypes.
#' 
#' Currently, this is hard-coded to access 1000 Genomes phase3 data hosted by
#' Brian Browning (author of BEAGLE):
#' 
#' \url{http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/}
#' 
#' This implementation discards multi-allelic markers that have a "," in the
#' ALT column.
#' 
#' The \code{pop} can be any of: ACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN,
#' FIN, GBR, GIH, GWD, IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU,
#' TSI, YRI. It can also be any super-population: AFR, AMR, EAS, EUR, SAS.
#'
#' Find more details here:
#' \url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
#' 
#' @param chrom a chromosome name (1-22,X) without "chr"
#' @param start a positive integer indicating the start of a genomic region
#' @param end a positive integer indicating the end of a genomic region
#' @param pop the name of a 1000 Genomes population (AMR,AFR,ASN,EUR,...)
#' @return A list with two dataframes:
#' \describe{
#'   \item{meta}{First 8 columns of the VCF file: CHROM, POS, ID, REF, ALT,
#'    QUAL, FILTER, INFO}
#'   \item{geno}{Columns 10 onward of the VCF file. All genotypes are converted
#'    to 0s and 1s representing REF and ALT alleles. This dataframe has two
#'    columns for each individual.}
#' }
#' 
#' @export
get_vcf <- function(chrom, start, end, pop = "EUR") {
  
  # Hard-coded superpopulations for each individual.
  superpops <- rep(
    c("EUR","EAS","AMR","EAS","AMR","EAS","AMR","EAS","AMR","EUR","AMR",
      "EUR","AMR","SAS","EAS","EUR","EAS","AFR","AMR","AFR","AMR","AFR",
      "AMR","AFR","AMR","AFR","EAS","AFR","EAS","AMR","AFR","AMR","AFR",
      "EAS","AFR","AMR","EAS","EUR","EAS","AMR","AFR","AMR","AFR","AMR",
      "AFR","AMR","AFR","AMR","EAS","AFR","AMR","AFR","SAS","AFR","EAS",
      "AFR","SAS","AFR","SAS","AFR","SAS","AFR","SAS","AFR","SAS","AFR",
      "SAS","AFR","SAS","AFR","SAS","AFR","SAS","EUR","AFR","EAS","AFR",
      "EAS","AFR","EAS","AFR","AMR","AFR","AMR","AFR","EUR","SAS"),
    c(185,42,3,35,4,28,9,15,111,1,68,24,9,4,6,73,56,8,2,5,24,2,13,5,7,5,
      17,4,21,2,1,3,3,17,2,3,20,14,1,2,2,11,4,8,3,1,13,2,35,2,1,22,4,7,
      4,27,5,16,10,6,13,6,9,9,12,58,10,84,6,74,3,25,310,99,19,103,29,70,
      18,34,142,19,8,45,52,107,103)
  )
  
  # Hard-coded populations for each individual.
  pops <- rep(
    c("GBR","FIN","GBR","FIN","CHS","PUR","CHS","PUR","CHS","PUR","CDX","PUR",
      "CLM","PUR","CLM","PUR","GBR","CLM","PUR","CLM","IBS","CLM","PEL","PJL",
      "KHV","IBS","GBR","CDX","KHV","ACB","PEL","ACB","PEL","ACB","PEL","ACB",
      "PEL","ACB","KHV","ACB","KHV","PEL","ACB","PEL","ACB","KHV","ACB","PEL",
      "CDX","GBR","IBS","CDX","PEL","ACB","PEL","ACB","PEL","ACB","PEL","ACB",
      "PEL","CDX","ACB","PEL","ACB","GWD","ACB","PJL","ACB","KHV","ACB","GWD",
      "ACB","GWD","PJL","GWD","PJL","GWD","PJL","GWD","PJL","GWD","PJL","GWD",
      "ESN","GWD","BEB","PJL","GWD","MSL","ESN","MSL","PJL","GWD","ESN","MSL",
      "PJL","ESN","GWD","MSL","BEB","PJL","STU","PJL","STU","PJL","STU","ITU",
      "STU","ITU","STU","PJL","ITU","BEB","STU","ITU","STU","BEB","STU","ITU",
      "STU","ITU","STU","ITU","STU","ITU","STU","ITU","STU","ITU","BEB","ITU",
      "STU","ITU","STU","ITU","CEU","YRI","CHB","YRI","JPT","LWK","JPT","YRI",
      "LWK","ASW","MXL","ASW","MXL","ASW","TSI","GIH"),
    c(55,17,31,82,42,3,35,4,28,9,15,41,18,26,16,10,1,27,11,30,24,3,6,4,6,70,
      3,22,34,8,2,5,24,2,13,5,7,5,17,4,21,2,1,3,3,17,2,3,20,1,13,1,2,2,11,4,
      8,3,1,13,2,35,2,1,9,4,9,4,7,4,10,7,2,8,5,16,10,6,13,6,9,9,12,37,19,2,4,
      6,10,27,43,4,6,6,29,39,3,8,2,15,13,8,5,7,19,6,1,11,5,1,12,3,20,20,13,15,
      13,21,9,12,8,2,2,10,7,8,1,7,4,1,28,5,1,7,2,3,99,19,103,29,70,18,34,60,
      81,1,19,8,45,52,107,103)
  )
  
  if (!is.na(pop) && !pop %in% pops && !pop %in% superpops) {
    stop(
      "Invalid pop=", pop,
      "\nMust be a pop: ", paste(unique(pops), collapse = " "),
      "\nOr a superpop: ", paste(unique(superpops), collapse = " ")
    )
  }
  
  if (!is.numeric(start) || start < 1) {
    stop("Invalid start=", start, "\nMust be a positive integer.")
  }
  
  if (!is.numeric(end) || end < start) {
    stop("Invalid end=", end, "\nMust be an integer >= start.")
  }
  
  if (!chrom %in% c(1:22, "X")) {
    stop(
      "Invalid chrom=", chrom,
      "\nMust be one of: ", paste(c(1:22, "X"), collapse = " ")
    )
  }
  
  # These are variants filtered by Brian Browning, the developer of BEAGLE.
  data_url = paste(sep = "",
    "http://tabix.iobio.io/?cmd=-h%20%27",
    "http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/",
    "individual_chromosomes/chr", chrom, ".1kg.phase3.v5a.vcf.gz",
    "%27%20", chrom, ":", start, "-", end
  )
  
  # Download the data from the server.
  txt <- RCurl::getURL(data_url)
  
  # Extract the sample identifiers from the VCF header.
  sample_ids <- strsplit(
    # FIXME: Should detect number of comment lines.
    # Get the 5th line and split it.
    strsplit(txt, "\n", fixed = TRUE)[[1]][5],
    "\t", fixed = TRUE
  )[[1]]
  
  # Discard the "CHROM,POS,...,INFO,V9" columns.
  sample_ids <- sample_ids[10:length(sample_ids)]
  
  # Read the body of the data into a dataframe.
  vcf <- read.delim(
    text = txt, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  
  # Assign the standard column names and sample identifiers.
  colnames(vcf) <- c(
    "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "V9",
    sample_ids
  )
  
  # Discard multi-allelic markers.
  vcf <- vcf[grep(",", vcf$ALT, invert = TRUE), ]
  
  # Select the genotype columns that belong to a particular population.
  if (pop %in% pops) {
    vcf <- vcf[,c(rep(TRUE, 9), pops == pop)]
  } else if (pop %in% superpops) {
    vcf <- vcf[,c(rep(TRUE, 9), superpops == pop)]
  }
  
  retval <- list()
  
  #load("data/sysdata.rda")
  retval$ind <- proxysnps::ind[colnames(vcf)[10:ncol(vcf)],]
  
  # Separate the metadata from the genotypes.
  retval$meta <- vcf[,1:8]
  retval$geno <- vcf[,10:ncol(vcf)]
  
  # Convert the genotypes to a numeric matrix.
  retval$geno <- t(apply(retval$geno, 1, function(row) {
    as.numeric(do.call(cbind, strsplit(row, "|", fixed = TRUE)))
  }))
  
  rownames(retval$geno) <- retval$meta$ID
  colnames(retval$geno) <- rep(retval$ind$Individual.ID, each = 2)
  
  return(retval)  
}