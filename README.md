# proxysnps

`proxysnps` is an R package that implements functions to get proxy SNPs
in [linkage disequilibrium][LD] (LD) with a [SNP] in the [1000 Genomes
Project][1000genomes].

```r
library(proxysnps)
d <- get_proxies(query = "rs42")
plot(d$POS, d$R.squared, main="rs42", xlab="Position", ylab=bquote("R"^2))
```

![rs42](https://github.com/slowkow/proxysnps/blob/master/rs42.png)

## Usage

See the [vignette] for more usage examples.

## Installation

```r
install.packages("devtools")
devtools::install_github("slowkow/proxysnps")
```

## Data

This package provides easy access to 1000 Genomes Project VCF files that have
been filtered by [Brian Browning][browning], available
[here][beagle_1000genomes].

[browning]: https://faculty.washington.edu/browning/
[beagle_1000genomes]: http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/

## Contributing

Please [submit an issue][issues] to report bugs or ask questions.

Please contribute bug fixes or new features with a [pull request][pull] to this
repository.

[issues]: https://github.com/slowkow/proxysnps/issues
[pull]: https://help.github.com/articles/using-pull-requests/

## Related work

[HaploReg][haploreg]

> HaploReg is a tool for exploring annotations of the noncoding genome at
> variants on haplotype blocks, such as candidate regulatory SNPs at
> disease-associated loci. Using [LD] information from the [1000 Genomes
> Project][1000genomes], linked SNPs and small indels can be visualized along
> with chromatin state and protein binding annotation from the [Roadmap
> Epigenomics][roadmap] and [ENCODE] projects, sequence conservation across
> mammals, the effect of SNPs on regulatory motifs, and the effect of SNPs on
> expression from [eQTL] studies.

[LDheatmap]

> Produces a graphical display, as a heat map, of measures of pairwise linkage
> disequilibria between SNPs. Users may optionally include the physical
> locations or genetic map distances of each SNP on the plot.

Also see the [Statistical Genetics CRAN Task View][statgen] for additional
R packages.

[LocusZoom][locuszoom]

> LocusZoom is a tool to plot regional association results from genome-wide
> association scans or candidate gene studies.

[SNAP][snap]

> SNAP is a computer program and web-based service for the rapid retrieval of
> linkage disequilibrium proxy SNP results given input of one or more query
> SNPs and based on empirical observations from the [International HapMap
> Project][hapmap] and the [1000 Genomes Project][1000genomes].

[SNPSNAP][snpsnap]

> The SNPsnap Web server enables SNP-based enrichment analysis by providing
> matched sets of SNPs that can be used to calibrate background expectations.
> Specifically, SNPsnap efficiently identifies sets of randomly drawn SNPs
> that are matched to a set of query SNPs based on allele frequency, number
> of SNPs in LD, distance to nearest gene and gene density.

[Tagger][tagger]

> Tagger is a tool for the selection and evaluation of tag SNPs from genotype
> data such as that from the [International HapMap Project][hapmap]. It
> combines the simplicity of pairwise tagging methods with the efficiency
> benefits of multimarker haplotype approaches.

[LD]: https://en.wikipedia.org/wiki/Linkage_disequilibrium
[SNP]: https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism
[eQTL]: https://en.wikipedia.org/wiki/Expression_quantitative_trait_loci

[vignette]: https://github.com/slowkow/proxysnps/blob/master/vignettes/proxysnps.md

[hapmap]: http://www.hapmap.org/
[1000genomes]: http://www.1000genomes.org/
[GRCh37]: http://www.1000genomes.org/faq/which-reference-assembly-do-you-use
[roadmap]: http://www.roadmapepigenomics.org/
[ENCODE]: https://www.encodeproject.org/

[haploreg]: http://www.broadinstitute.org/mammals/haploreg/
[LDheatmap]: https://cran.r-project.org/web/packages/LDheatmap/index.html
[statgen]: https://cran.r-project.org/web/views/Genetics.html
[locuszoom]: http://locuszoom.sph.umich.edu/locuszoom/
[snap]: http://www.broadinstitute.org/mpg/snap/
[snpsnap]: https://data.broadinstitute.org/mpg/snpsnap/about.html
[tagger]: https://www.broadinstitute.org/mpg/tagger/
