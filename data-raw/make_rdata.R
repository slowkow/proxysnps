ind <- read.delim(
  file = "data-raw/20130606_g1k.ped.gz",
  stringsAsFactors = FALSE
)

# Delete unused columns.
for (col in colnames(ind)) {
  if (length(table(ind[[col]])) == 1) {
    ind[[col]] <- NULL
  }
}

rownames(ind) <- ind$Individual.ID

ind$Gender <- ifelse(ind$Gender == 2, "female", "male")

ind$SuperPopulation[ind$Population == "CHB"] <- "EAS"
ind$SuperPopulation[ind$Population == "JPT"] <- "EAS"
ind$SuperPopulation[ind$Population == "CHS"] <- "EAS"
ind$SuperPopulation[ind$Population == "CDX"] <- "EAS"
ind$SuperPopulation[ind$Population == "CHB"] <- "EAS"
ind$SuperPopulation[ind$Population == "KHV"] <- "EAS"

ind$SuperPopulation[ind$Population == "CEU"] <- "EUR"
ind$SuperPopulation[ind$Population == "TSI"] <- "EUR"
ind$SuperPopulation[ind$Population == "FIN"] <- "EUR"
ind$SuperPopulation[ind$Population == "GBR"] <- "EUR"
ind$SuperPopulation[ind$Population == "IBS"] <- "EUR"

ind$SuperPopulation[ind$Population == "YRI"] <- "AFR"
ind$SuperPopulation[ind$Population == "LWK"] <- "AFR"
ind$SuperPopulation[ind$Population == "GWD"] <- "AFR"
ind$SuperPopulation[ind$Population == "MSL"] <- "AFR"
ind$SuperPopulation[ind$Population == "ESN"] <- "AFR"
ind$SuperPopulation[ind$Population == "ASW"] <- "AFR"
ind$SuperPopulation[ind$Population == "ACB"] <- "AFR"

ind$SuperPopulation[ind$Population == "MXL"] <- "AMR"
ind$SuperPopulation[ind$Population == "PUR"] <- "AMR"
ind$SuperPopulation[ind$Population == "CLM"] <- "AMR"
ind$SuperPopulation[ind$Population == "PEL"] <- "AMR"

ind$SuperPopulation[ind$Population == "GIH"] <- "SAS"
ind$SuperPopulation[ind$Population == "PJL"] <- "SAS"
ind$SuperPopulation[ind$Population == "BEB"] <- "SAS"
ind$SuperPopulation[ind$Population == "STU"] <- "SAS"
ind$SuperPopulation[ind$Population == "ITU"] <- "SAS"

save(list = c("ind"), file = "data/ind.rda")
