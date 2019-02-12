library('sleuth')
library('cowplot')

metadata <- read.table('Rdh10KO_liver_BAT_RNAseq_SampleID_MV_2019-02-10_more_data.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

### load the serialization:
ttg <- readRDS("ttg.rds")

so <- sleuth_prep(metadata, target_mapping = ttg,
  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)

33325  ENSMUSG00000008489           Elavl2  ENSMUST00000147611.1
33326  ENSMUSG00000008489           Elavl2  ENSMUST00000176362.1
33327  ENSMUSG00000008489           Elavl2  ENSMUST00000144769.1
33328  ENSMUSG00000008489           Elavl2  ENSMUST00000176151.1
33329  ENSMUSG00000008489           Elavl2  ENSMUST00000107124.9
33330  ENSMUSG00000111090          Gm30698    ENSMUST00000214789
33331  ENSMUSG00000010045          Tmem115  ENSMUST00000010189.2
33332  ENSMUSG00000093913       Obox4-ps28  ENSMUST00000175654.1


so <- sleuth_fit(so, ~sex + tissue + genotype, 'full')

so <- sleuth_fit(so, ~sex + tissue, 'reduced')

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
