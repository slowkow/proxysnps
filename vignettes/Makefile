proxysnps.md : proxysnps.Rmd
	R --slave -e 'knitr::knit("proxysnps.Rmd")'
	perl -i -pe 's{\(figures/}{(https://github.com/slowkow/proxysnps/blob/master/vignettes/figures/}' proxysnps.md
