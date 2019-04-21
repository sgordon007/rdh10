# General prep of R

#### activate R
R

#### install
source("http://bioconductor.org/biocLite.R")
biocLite("devtools")    # only if devtools not yet installed
biocLite("pachterlab/sleuth")
biocLite("cowplot")

#### activate the library within R
library('sleuth')
library('cowplot')

# Experiment 1, two variables, 'genotype' and sex.  Run the two sexes independently

## males

#### load the serialization:
ttg <- readRDS("ttg.rds")

metadata <- read.table('experiment1_males_matadata_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

# make the sleuth object
so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

### fit the data with intercept model vs model including genotype effect:
so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')

 models(so)
[  genotype  ]
formula:  ~genotype 
data modeled:  obs_counts 
transform sync'ed:  TRUE 
coefficients:
	(Intercept)
 	genotypeWT
[  intercept  ]
formula:  ~1 
data modeled:  obs_counts 
transform sync'ed:  TRUE 
coefficients:
	(Intercept)

#### LRT test for DEG
so <- sleuth_lrt(so, 'intercept', 'genotype')
intercept_genotype <- sleuth_results(so, 'intercept:genotype', 'lrt')
intercept_genotype_table_gene <- dplyr::filter(intercept_genotype, qval <= 0.05)
head(intercept_genotype_table_gene, 20)
write.csv(intercept_genotype_table_gene, file = "experiment1_males_LRT_intercept_genotype_table_gene.csv")

#### Wald test for DEG
##### several genes with Gtpase activity
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_males_Wald_wt_genotypeWT_results_table_significant05.csv")

sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, pval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_males_Wald_wt_genotypeWT_results_table_Pval_significant05.csv")

## females

metadata <- read.table('experiment1_females_matadata_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

# make the sleuth object
so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

### fit the data with intercept model vs model including genotype effect:
so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')


#### LRT test for DEG females
so <- sleuth_lrt(so, 'intercept', 'genotype')
intercept_genotype <- sleuth_results(so, 'intercept:genotype', 'lrt')
intercept_genotype_table_gene <- dplyr::filter(intercept_genotype, qval <= 0.05)
head(intercept_genotype_table_gene, 20)
write.csv(intercept_genotype_table_gene, file = "experiment1_females_LRT_intercept_genotype_table_gene.csv")

#### Wald test for DEG females
##### type 2 diabetes: Madd:  https://uswest.ensembl.org/Mus_musculus/Gene/Phenotype?db=core;g=ENSMUSG00000040687;r=2:91137360-91183837
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_females_Wald_wt_genotypeWT_results_table_significant05.csv")

sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, pval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_females_Wald_wt_genotypeWT_results_table_Pval_significant05.csv")



# Experiment 2, one variable, 'genotype'.  The reset of the variable are the same.

metadata <- read.table('experiment2_matadata_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

# make the sleuth object
so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

### fit the data with intercept model vs model including genotype effect:
so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')

#### LRT test for DEG
so <- sleuth_lrt(so, 'intercept', 'genotype')
intercept_genotype <- sleuth_results(so, 'intercept:genotype', 'lrt')
intercept_genotype_table_gene <- dplyr::filter(intercept_genotype, qval <= 0.05)
head(intercept_genotype_table_gene, 20)
write.csv(intercept_genotype_table_gene, file = "experiment2_LRT_intercept_genotype_table_gene.csv")

#### Wald test for DEG
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment2_Wald_wt_genotypeWT_results_table_significant05.csv")

sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, pval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment2_Wald_wt_genotypeWT_results_table_Pval_significant05.csv")




# Experiment 3, two variables, 'genotype' and sex.  Run the two sexes independently

## males

#### load the serialization:
ttg <- readRDS("ttg.rds")

metadata <- read.table('experiment3_males_matadata_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 10)

# make the sleuth object
so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

### fit the data with intercept model vs model including genotype effect:
so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')


#### LRT test for DEG
so <- sleuth_lrt(so, 'intercept', 'genotype')
intercept_genotype <- sleuth_results(so, 'intercept:genotype', 'lrt')
intercept_genotype_table_gene <- dplyr::filter(intercept_genotype, qval <= 0.05)
head(intercept_genotype_table_gene, 20)
write.csv(intercept_genotype_table_gene, file = "experiment3_males_LRT_intercept_genotype_table_gene.csv")

#### Wald test for DEG
#####  tumor suppressor essential for controlling cell proliferation: https://www.ncbi.nlm.nih.gov/gene/2195
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment3_males_Wald_wt_genotypeWT_results_table_significant05.csv")

sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, pval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment3_males_Wald_wt_genotypeWT_results_table_Pval_significant05.csv")

## females

metadata <- read.table('experiment3_females_matadata_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 10)

# make the sleuth object
so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

### fit the data with intercept model vs model including genotype effect:
so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')


#### LRT test for DEG females
so <- sleuth_lrt(so, 'intercept', 'genotype')
intercept_genotype <- sleuth_results(so, 'intercept:genotype', 'lrt')
intercept_genotype_table_gene <- dplyr::filter(intercept_genotype, qval <= 0.05)
head(intercept_genotype_table_gene, 20)
write.csv(intercept_genotype_table_gene, file = "experiment3_females_LRT_intercept_genotype_table_gene.csv")

#### Wald test for DEG females
#####
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment3_females_Wald_wt_genotypeWT_results_table_significant05.csv")

sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, pval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment3_females_Wald_wt_genotypeWT_results_table_Pval_significant05.csv")

### no overlap between DEG sets derived from qval thresholding
##### Slc7a2: lipid metabolism: http://www.informatics.jax.org/marker/MGI:99828
##### Slc26a10: lipid metabolism: http://www.informatics.jax.org/marker/MGI:2143920

##### Car3 ENSMUSG00000027559: obesity:  http://uswest.ensembl.org/Mus_musculus/Gene/Phenotype?db=core;g=ENSMUSG00000027559;r=3:14863512-14872523
http://uswest.ensembl.org/Mus_musculus/Gene/Phenotype?db=core;g=ENSMUSG00000027559;r=3:14863512-14872523

##### nervous system: Madd:  http://uswest.ensembl.org/Mus_musculus/Gene/Phenotype?db=core;g=ENSMUSG00000040687;r=2:91137360-91183837

#### in one of the sets it seemed like there might entrichment for ubiquitin processes.



