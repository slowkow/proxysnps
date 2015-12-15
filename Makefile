README.md : README.Rmd
	R --slave -e 'knitr::knit("README.Rmd")'
	perl -i -pe 's{\(figures/}{(https://github.com/slowkow/proxysnps/blob/master/figures/}' README.md
