% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_proxies.R
\name{get_proxies}
\alias{get_proxies}
\title{Get proxy SNPs for a SNP at a given genomic position.}
\usage{
get_proxies(chrom = NA, pos = NA, query = NA, window_size = 1e+05,
  pop = NA)
}
\arguments{
\item{chrom}{a chromosome name (1-22,X) without "chr"}

\item{pos}{a positive integer indicating the position of a SNP}

\item{window_size}{a positive integer indicating the size of the window}

\item{pop}{the name of a 1000 Genomes population (AMR,AFR,ASN,EUR,...). Set
this to NA to use all populations.}
}
\value{
A dataframe with the following columns:
  \describe{
    \item{CHROM}{Chromosome name, e.g. "1"}
    \item{POS}{Position, e.g. 583090}
    \item{ID}{Identifier, e.g. "rs11063140"}
    \item{REF}{Reference allele, e.g. "A"}
    \item{ALT}{Alternative allele, e.g. "G"}
    \item{MAF}{Minor allele frequency, e.g. 0.1}
    \item{R.squared}{Squared Pearson correlation coefficient, e.g. 1.0}
    \item{D.prime}{D prime value, e.g. 1.0}
    \item{CHOSEN}{Binary indicator set to TRUE for the SNP of interest}
  }
}
\description{
Returns a dataframe with proxy SNPs.
}
\details{
Currently, this is hard-coded to access 1000 Genomes phase3 data hosted by
Brian Browning (author of BEAGLE):

\url{http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/}

This implementation discards multi-allelic markers that have a "," in the
ALT column.

The \code{pop} can be any of: ACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN,
FIN, GBR, GIH, GWD, IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU,
TSI, YRI. It can also be any super-population: AFR, AMR, EAS, EUR, SAS.

Find more details here:
\url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
}
\examples{
d <- get_proxies(chrom = "12", pos = 583090, window_size = 1e5, pop = "AFR")
head(d)
 
}

