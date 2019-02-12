#### install biomaRt
conda install -c r r-essentials
conda install r-xml
biocLite("biomaRt")
BiocManager::install("biomaRt", version = "3.8")

#### get the data

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = "jan2019.archive.ensembl.org")
  # host = "ensembl.org")
ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
  "ensembl_gene_id", "external_gene_name", "description",
  "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))
head(ttg)

#### serialize the data object
saveRDS(ttg, file = 'ttg.rds')

### load the serialization:
ttg <- readRDS("ttg.rds")

# Write the data to a CSV in R
write.csv(ttg, file = "ttg.csv")

###biomaRt::listMarts()
