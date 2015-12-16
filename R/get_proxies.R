#' Get proxy SNPs for a SNP at a given genomic position.
#' 
#' Returns a dataframe with proxy SNPs.
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
#' @param pos a positive integer indicating the position of a SNP
#' @param window_size a positive integer indicating the size of the window
#' @param pop the name of a 1000 Genomes population (AMR,AFR,ASN,EUR,...). Set
#'        this to NA to use all populations.
#' @return A dataframe with the following columns:
#'   \describe{
#'     \item{CHROM}{Chromosome name, e.g. "1"}
#'     \item{POS}{Position, e.g. 583090}
#'     \item{ID}{Identifier, e.g. "rs11063140"}
#'     \item{REF}{Reference allele, e.g. "A"}
#'     \item{ALT}{Alternative allele, e.g. "G"}
#'     \item{MAF}{Minor allele frequency, e.g. 0.1}
#'     \item{R.squared}{Squared Pearson correlation coefficient, e.g. 1.0}
#'     \item{D.prime}{D prime value, e.g. 1.0}
#'     \item{CHOSEN}{Binary indicator set to TRUE for the SNP of interest}
#'   }
#' 
#' @examples
#' d <- get_proxies(chrom = "12", pos = 583090, window_size = 1e5, pop = "AFR")
#' head(d)
#'  
#' @export
get_proxies <- function(
  chrom = NA, pos = NA, query = NA, window_size = 1e5, pop = NA) {
  
#   chrom = "12"
#   pos = 583090
#   window_size = 1e5
#   pop = "AFR"
  
  # Try to get the position from myvariant.info if there's a query.
  if (is.na(chrom) && is.na(pos) && !is.na(query)) {
    v <- myvariant::queryVariant(query)
    if (v$total > 0) {
      chrom <- v$hits$dbsnp$chrom[1]
      pos <- v$hits$dbsnp$hg19$start[1]
    } else {
      stop(sprintf("[myvariant.info] No hits for query '%s'", query))
    }
  }
  
  # Get VCF data for this genomic region.
  vcf <- get_vcf(
    chrom = chrom, 
    start = floor(pos - window_size / 2),
    end =  floor(pos + window_size / 2),
    pop = pop)
  
  # Separate the metadata from the genotypes.
  meta <- vcf$meta
  geno <- vcf$geno
  
  # Drop unused columns.
  meta$QUAL <- NULL
  meta$FILTER <- NULL
  meta$INFO <- NULL
  
  # Compute the minor allele frequencies.
  maf <- sapply(rowSums(geno) / ncol(geno), function(x) min(x, 1 - x))
  
  # Select the SNP of interest.
  chosen <- meta$CHROM == chrom & meta$POS == pos
  
  # If the position does not match any SNPs, choose the middle one.
  if (!length(which(chosen))) {
    chosen <- rep(FALSE, nrow(meta))
    chosen[floor(nrow(meta) / 2)] <- TRUE
  }
  
  # Compute R.squared and D.prime with binary (0/1) markers.
  ld <- apply(geno, 1, function(y) {
    compute_ld(geno[which(chosen),], y)
  })
  
  # Create a dataframe.
  retval <- cbind(
    meta,
    MAF = maf,
    R.squared = round(sapply(ld, "[[", "R.squared"), 6),
    D.prime = round(sapply(ld, "[[", "D.prime"), 6),
    CHOSEN = chosen
  )
  
  return(retval[order(retval$R.squared, decreasing = TRUE),])
}

#' Compute two commonly used linkage disequilibrium statistics.
#' 
#' Compute R.squared and D.prime for two binary numeric vectors.
#' 
#' Find more details here:
#' \url{https://en.wikipedia.org/wiki/Linkage_disequilibrium}
#' 
#' @param x a numeric vector of ones and zeros
#' @param y a numeric vector of ones and zeros
#' @return A list with two items:
#' \describe{
#'  \item{R.squared}{Squared Pearson correlation coefficient.}
#'  \item{D.prime}{Coefficient of linkage disequilibrium D divided by the
#'   theoretical maximum.}
#' }
#' 
#' @examples
#' compute_ld(c(0,0,0,1,1,1), c(1,1,1,1,0,0))
#' 
#' @export
compute_ld <- function(x, y) {
  stopifnot(all(names(table(x)) %in% c("0", "1")))
  stopifnot(all(names(table(y)) %in% c("0", "1")))
  stopifnot(length(x) == length(y))
  
  # Returned list.
  retval <- list("D.prime" = 0, "R.squared" = 0)
  
  # Allele frequencies.
  px <- max(sum(x) / length(x), .Machine$double.eps)
  py <- max(sum(y) / length(x), .Machine$double.eps)
  pxy <- sum(x & y) / length(x)
  
  # D prime.
  d <- pxy - px * py
  dmax <- 0
  if (d < 0) {
    dmax <- min(px * py, (1 - px) * (1 - py))
  } else {
    dmax <- min(px * (1 - py), (1 - px) * py)
  }
  if (dmax != 0) {
    retval[["D.prime"]] <- d / dmax
  }
  
  # Squared Pearson correlation coefficient.
  retval[["R.squared"]] <- (d / sqrt(px * (1 - px) * py * (1 - py))) ^ 2
  
  return(retval)
}
