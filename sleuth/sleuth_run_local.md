Background

This html has info on making heatmap with sleuth live:  
file:///Users/sgordon/Documents/project_management/marta/sleuth.html



# Experiment 2, single 'genotype' variable

#### activate the library within R
library('sleuth')
library('cowplot')

metadata <- read.table('experiment2_matadata_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

### load the serialization:
ttg <- readRDS("ttg.rds")

so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

#### plot PCA of the sample and label samples by their genotype
plot_pca(so, color_by = 'genotype')

new_position_theme <- theme(legend.position = c(0.80, 0.90))
plot_pca(so, color_by = 'genotype', text_labels = TRUE) +
  new_position_theme
  
plot_loadings(so, pc_input = 1)

##### get error when using plot_bootstrap:
plot_bootstrap(so, 'ENSMUSG00000096887', color_by = 'genotype') +
  new_position_theme


### fit the data with intercept model vs model including genotype effect:
so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')
models(so)

so <- sleuth_lrt(so, 'intercept', 'genotype')
intercept_genotype <- sleuth_results(so, 'intercept:genotype', 'lrt')
intercept_genotype_table_gene <- dplyr::filter(intercept_genotype, qval <= 0.05)
head(intercept_genotype_table_gene, 20)
write.csv(intercept_genotype_table_gene, file = "experiment2_intercept_genotype_table_gene.csv")

##### again getting error with following get and plot bootstrap:
get_bootstrap_summary(so, "ENSMUSG00000068457", 'tpm')

plot_bootstrap(so, "ENSMUSG00000068457", color_by = 'genotype') +
  new_position_theme

plot_bootstrap(so, 'ENSMUST00000082420.1', color_by = 'genotype') +
  new_position_theme

#### Try Wald test
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment2_sleuth_wt_genotypeWT_results_table_significant05.csv")

sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.1)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment2_sleuth_wt_genotypeWT_results_table_significant1.csv")



# Move onto Experiment 1

### Experiment 1, two variables, 'genotype' and sex

#### activate the library within R
library('sleuth')
library('cowplot')

metadata <- read.table('experiment1_matadata_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

### load the serialization:
ttg <- readRDS("ttg.rds")

so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

#### does this work
plot_pca(so, color_by = 'genotype')
plot_pca(so, color_by = 'sex')

new_position_theme <- theme(legend.position = c(0.80, 0.90))
plot_pca(so, color_by = 'genotype', text_labels = TRUE) +
  new_position_theme
  
plot_loadings(so, pc_input = 1)

plot_bootstrap(so, 'ENSMUSG00000096887', color_by = 'genotype') +
  new_position_theme

vignette('intro', package = 'sleuth')
?plot_volcano
plot_volcano(so, test, test_type = "lrt", which_model = "sex",
       sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
       highlight = NULL)
       
plot_volcano(so, color_by = 'genotype')

so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~sex, 'sex')
so <- sleuth_fit(so, ~1, 'intercept')
models(so)
tests(so)

so <- sleuth_lrt(so, 'intercept', 'genotype')
intercept_genotype <- sleuth_results(so, 'intercept:genotype', test_type = 'lrt')
intercept_genotype_table_gene <- dplyr::filter(intercept_genotype, qval <= 0.05)
head(intercept_genotype_table_gene, 20)
their were no significant genes at default cutoff with LRT test

#### Try Wald test
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
head(sleuth_wt_genotypeWT_results_table)
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
head(sleuth_wt_genotypeWT_results_table_significant, 20)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_sleuth_wt_genotypeWT_results_table_significant.csv")


so <- sleuth_lrt(so, 'intercept', 'sex')
intercept_sex <- sleuth_results(so, 'intercept:sex', test_type = 'lrt')
intercept_sex_table_gene <- dplyr::filter(intercept_sex, qval <= 0.05)
head(intercept_sex_table_gene, 20)
write.csv(intercept_sex_table_gene, file = "experiment1_intercept_sex_table_gene.csv")

#### what about Experiment 1 comparing single factor model to multifactor model:

so <- sleuth_fit(so, ~genotype + sex, 'full')

models(so)

so <- sleuth_lrt(so, 'sex', 'full')
sex_full <- sleuth_results(so, 'sex:full', test_type = 'lrt')
sex_full <- sleuth_results(so, 'sex:full', test_type = 'wald')

sex_full_table_gene <- dplyr::filter(sex_full, qval <= 0.05)
head(intercept_genotype_table_gene, 20)
their were no significant genes at default cutoff

#### Try Wald test
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
head(sleuth_wt_genotypeWT_results_table)
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
head(sleuth_wt_genotypeWT_results_table_significant, 20)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_sleuth_wt_genotypeWT_results_table_significant.csv")

so <- sleuth_wt(so,'genotypeWT','full')
so <- sleuth_wt(so,'sexmale','full')

full_sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'full')
head(full_sleuth_wt_genotypeWT_results_table)
full_sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(full_sleuth_wt_genotypeWT_results_table, qval <= 0.05)


tests(so)

so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~sex, 'sex')

dim(sleuth_table) and dim(sleuth_significant), you will see that 2237/28095 transcripts (7.96%) are significantly DE



# Run Experiment 1, but do male and female separately and leave out samples L and A

# Experiment 1, two variables, 'genotype' and sex

#### make a clean slate, a new so object for male subjects:

#### activate the library within R
library('sleuth')
library('cowplot')

metadata <- read.table('experiment1_matadata_h5_male_NO-L-A.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

### load the serialization:
ttg <- readRDS("ttg.rds")

so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

#### does this work
plot_pca(so, color_by = 'genotype')
plot_pca(so, color_by = 'sex')

new_position_theme <- theme(legend.position = c(0.80, 0.90))
plot_pca(so, color_by = 'genotype', text_labels = TRUE) +
  new_position_theme
  
plot_loadings(so, pc_input = 1)


so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')
models(so)
tests(so)

so <- sleuth_lrt(so, 'intercept', 'genotype')
sleuth_lrt_intercept_genotype <- sleuth_results(so, 'intercept:genotype', test_type = 'lrt')
sleuth_lrt_intercept_genotype_table_gene <- dplyr::filter(sleuth_lrt_intercept_genotype, qval <= 0.1)
head(sleuth_lrt_intercept_genotype_table_gene, 20)


#### Try Wald test
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
head(sleuth_wt_genotypeWT_results_table_significant, 20)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_matadata_h5_male_NO-L-A_wt_genotypeWT_results_table_significant05.csv")
##### test at 0.1 q-value:
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.1)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_matadata_h5_male_NO-L-A_wt_genotypeWT_results_table_significant1.csv")

Gene: Rdh10 ENSMUSG00000025921


#### -----

#### make a clean slate, a new so object for male subjects:

#### activate the library within R
library('sleuth')
library('cowplot')

metadata <- read.table('experiment1_matadata_h5_female_NO-L-A.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

### load the serialization:
ttg <- readRDS("ttg.rds")

so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

#### does this work
plot_pca(so, color_by = 'genotype')
  
plot_loadings(so, pc_input = 1)

so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')
models(so)
tests(so)

#### Try Wald test
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
head(sleuth_wt_genotypeWT_results_table_significant, 20)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_matadata_h5_female_NO-L-A_wt_genotypeWT_results_table_significant05.csv")
##### test at 0.1 q-value:
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.1)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_matadata_h5_female_NO-L-A_wt_genotypeWT_results_table_significant1.csv")

Gene: Rdh10 ENSMUSG00000025921


#### -----

#### make a clean slate, a new so object for male subjects:

#### activate the library within R
library('sleuth')
library('cowplot')

metadata <- read.table('experiment1_matadata_h5_NO-L-A.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

### load the serialization:
ttg <- readRDS("ttg.rds")

so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

#### does this work
plot_pca(so, color_by = 'genotype')
  
plot_loadings(so, pc_input = 1)

so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')
models(so)
tests(so)

#### Try Wald test
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
dim(sleuth_wt_genotypeWT_results_table_significant)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_matadata_h5_NO-L-A_wt_genotypeWT_results_table_significant05.csv")
##### test at 0.1 q-value:
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.1)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment1_matadata_h5_NO-L-A_wt_genotypeWT_results_table_significant1.csv")

Gene: Rdh10 ENSMUSG00000025921



# Experiment 3: sex and genotype within the liver

#### make a clean slate, a new so object for all subjects:

#### activate the library within R
library('sleuth')
library('cowplot')

metadata <- read.table('experiment3_matadata_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

### load the serialization:
ttg <- readRDS("ttg.rds")

so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

#### does this work
plot_pca(so, color_by = 'genotype')
  
plot_loadings(so, pc_input = 1)

so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')
models(so)
tests(so)

#### Try Wald test
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
dim(sleuth_wt_genotypeWT_results_table_significant)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment3_matadata_h5_wt_genotypeWT_results_table_significant05.csv")
##### test at 0.1 q-value:
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.1)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "experiment3_matadata_h5_wt_genotypeWT_results_table_significant1.csv")

Gene: Rdh10 ENSMUSG00000025921

## try lrt reduced vs full model

so <- sleuth_fit(so, ~sex, 'reduced')

so <- sleuth_fit(so, ~sex + genotype, 'full')

so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_gene05 <- dplyr::filter(sleuth_table_gene, qval <= 0.05)
write.csv(sleuth_table_gene05, file = "experiment3_matadata_h5_LRT_ConditionedOnSex_genotypeWT_results_table_significant05.csv")

#### transcript level
sleuth_table_transcript <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
sleuth_table_transcript05 <- dplyr::filter(sleuth_table_gene, qval <= 0.05)
write.csv(sleuth_table_transcript05, file = "experiment3_transcriptLevel_matadata_h5_LRT_ConditionedOnSex_genotypeWT_results_table_significant05.csv")



########################################## Full experiment

# Do some sanity checks to make sure that sleuth is working



#### activate the library within R
library('sleuth')
library('cowplot')

metadata <- read.table('Rdh10KO_liver_BAT_RNAseq_SampleID_MV_2019-02-10_more_data_abundance_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

head(metadata, 5)

### load the serialization:
ttg <- readRDS("ttg.rds")

so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

head(so$target_mapping)

#### 
plot_pca(so, color_by = 'sex')
dev.copy(pdf,'all_samples_labeled_by_sex_pca.pdf')
dev.off()

#### 
plot_pca(so, color_by = 'tissue')
dev.copy(pdf,'all_samples_labeled_by_tissue_pca.pdf')
dev.off()

#### 
plot_pca(so, color_by = 'genotype')
dev.copy(pdf,'all_samples_labeled_by_genotype_pca.pdf')
dev.off()

  
so <- sleuth_fit(so, ~sex, 'sex')
so <- sleuth_fit(so, ~tissue, 'tissue')

so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')
models(so)
tests(so)

#### Try Wald tests for genotype
tests(so)
so <- sleuth_wt(so,'genotypeWT','genotype')
sleuth_wt_genotypeWT_results_table <- sleuth_results(so, 'genotypeWT', test_type = 'wt', 'genotype')
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "all_samples_metadata_h5_wt_genotypeWT_results_table_significant05.csv")
##### test at 0.1 q-value:
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.1)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "all_samples_metadata_h5_wt_genotypeWT_results_table_significant1.csv")

#### Try Wald tests for tissue
tests(so)
so <- sleuth_wt(so,'tissueliver','tissue')
results_table <- sleuth_results(so, 'tissueliver', test_type = 'wt', 'tissue')
significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.05)
write.csv(significant, file = "all_samples_metadata_h5_wt_tissueliver_results_table_significant05.csv")
##### test at 0.1 q-value:
significant <- dplyr::filter(results_table, qval <= 0.1)
write.csv(significant, file = "all_samples_metadata_h5_wt_tissueliver_results_table_significant1.csv")

#### Try Wald tests for sex
tests(so)
so <- sleuth_wt(so,'sexmale','sex')
results_table <- sleuth_results(so, 'sexmale', test_type = 'wt', 'sex')
significant <- dplyr::filter(results_table, qval <= 0.05)
write.csv(significant, file = "all_samples_metadata_h5_wt_sexmale_results_table_significant05.csv")
##### test at 0.1 q-value:
significant <- dplyr::filter(results_table, qval <= 0.1)
write.csv(significant, file = "all_samples_metadata_h5_wt_sexmale_results_table_significant1.csv")



Gene: Rdh10 ENSMUSG00000025921

## try lrt between naive and single parameter model
so <- sleuth_lrt(so, 'intercept', 'sex')
sleuth_table_gene <- sleuth_results(so, 'intercept:sex', 'lrt', show_all = FALSE)
sleuth_table_gene05 <- dplyr::filter(sleuth_table_gene, qval <= 0.1)
write.csv(sleuth_table_gene05, file = "experiment3_matadata_h5_LRT_ConditionedOnSex_genotypeWT_results_table_significant05.csv")


##### test at 0.1 q-value:
sleuth_wt_genotypeWT_results_table_significant <- dplyr::filter(sleuth_wt_genotypeWT_results_table, qval <= 0.1)
write.csv(sleuth_wt_genotypeWT_results_table_significant, file = "all_samples_metadata_h5_wt_genotypeWT_results_table_significant1.csv")



## try lrt reduced vs full model

so <- sleuth_fit(so, ~sex, 'reduced')

so <- sleuth_fit(so, ~sex + genotype, 'full')

so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_gene05 <- dplyr::filter(sleuth_table_gene, qval <= 0.05)
write.csv(sleuth_table_gene05, file = "experiment3_matadata_h5_LRT_ConditionedOnSex_genotypeWT_results_table_significant05.csv")


