## activate conda

source $BSCRATCH/mus_2019_01_29/conda/anaconda2/bin/activate

# activate conda env
source activate sleuth

# activate R

# load libs
library('sleuth')
library('cowplot')

metadata <- read.table('Rdh10KO_liver_BAT_RNAseq_SampleID_MV_2019-02-10_more_data_abundance_h5.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

### load the serialization:
ttg <- readRDS("ttg.rds")

so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)

                 ens_gene         ext_gene             target_id
1      ENSMUSG00000064372            mt-Tp    ENSMUST00000082423

33325  ENSMUSG00000008489           Elavl2  ENSMUST00000147611.1
33326  ENSMUSG00000008489           Elavl2  ENSMUST00000176362.1
33327  ENSMUSG00000008489           Elavl2  ENSMUST00000144769.1
33328  ENSMUSG00000008489           Elavl2  ENSMUST00000176151.1
33329  ENSMUSG00000008489           Elavl2  ENSMUST00000107124.9
33330  ENSMUSG00000111090          Gm30698    ENSMUST00000214789
33331  ENSMUSG00000010045          Tmem115  ENSMUST00000010189.2
33332  ENSMUSG00000093913       Obox4-ps28  ENSMUST00000175654.1

#### does this work
plot_pca(so, color_by = 'genotype')
plot_pca(so, color_by = 'tissue')


plot_pca(so, text_labels = TRUE, color_by = 'center')


plot_scatter(so, 'DRR002306', 'DRR002319')

so <- sleuth_fit(so, ~sex + tissue + genotype, 'full')

so <- sleuth_fit(so, ~sex + tissue, 'reduced')


so <- sleuth_fit(so, ~tissue + sex + genotype, 'full')

so <- sleuth_fit(so, ~tissue + sex, 'reduced')


so <- sleuth_fit(so, ~tissue + genotype, 'full')

so <- sleuth_fit(so, ~tissue, 'reduced')


so <- sleuth_fit(so, ~tissue, 'tissue')
so <- sleuth_fit(so, ~1, 'intercept')
so <- sleuth_lrt(so, 'intercept', 'tissue')
intercept_tissue <- sleuth_results(so, 'intercept:tissue', 'lrt')
intercept_tissue <- dplyr::filter(intercept_tissue, qval <= 0.05)


so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~1, 'intercept')
so <- sleuth_lrt(so, 'intercept', 'genotype')
intercept_tissue <- sleuth_results(so, 'intercept:genotype', 'lrt')
sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)





#### download the TPM data to plot locally
nohup tar -cvzf results_download.tar.gz results_download > nohup.results_download.2019-02-16.compressed &



#### The likelihood ratio test (lrt) is performed with
so <- sleuth_lrt(so, 'reduced', 'full')

#### When running the command ‘sleuth_results,’ sleuth uses the p-values from comparing transcripts to make a gene-level determination and perform gene differential expression.

sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)

#### The most significantly differential genes are
head(sleuth_table_gene, 20)


#### Because gene results are built on transcript results, the gene and transcript results are entirely consistent and compatible with each other. To visualize the transcript results that led to the gene results above, one merely runs sleuth_results again but this time setting the flag ‘pval_aggregate’ to FALSE.
sleuth_table_tx <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
sleuth_table_tx <- dplyr::filter(sleuth_table_tx, qval <= 0.05)
head(sleuth_table_tx, 20)

sleuth_live(so)
