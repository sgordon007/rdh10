http://ccb.jhu.edu/software/hisat2/manual.shtml#the-hisat2-build-indexer

cd mus_2019_01_29/reference_sequences

gunzip GCF_000001635.26_GRCm38.p6_genomic.fna.gz

grep ">" GCF_000001635.26_GRCm38.p6_genomic.fna

stats.sh GCF_000001635.26_GRCm38.p6_genomic.fna


module load hisat/2.1.0

hisat2-build [options]* <reference_in> <ht2_base>

hisat2-build GCF_000001635.26_GRCm38.p6_genomic.fna GCF_000001635_26

### rebuild the index using known splice sites
hisat2-build --ss splicesites.txt --exon exons.txt GCF_000001635.26_GRCm38.p6_genomic.fna GCF_000001635_26_wsplices




--known-splicesite-infile <path>

--known-splicesite-infile <path>
With this mode, you can provide a list of known splice sites, which HISAT2 makes use of to align reads with small anchors.
You can create such a list using python hisat2_extract_splice_sites.py genes.gtf > splicesites.txt, where hisat2_extract_splice_sites.py is included in the HISAT2 package, 
genes.gtf is a gene annotation file, and splicesites.txt is a list of splice sites with which you provide HISAT2 in this mode. 
Note that it is better to use indexes built using annotated transcripts (such as genome_tran or genome_snp_tran), 
which works better than using this option. It has no effect to provide splice sites that are already included in the indexes.

conda create -n gffread -c bioconda gffread

conda create -n hisat2 -c bioconda hisat2

source activate hisat2

python ~/conda_envs/.conda/envs/hisat2/bin/hisat2_extract_splice_sites.py GCF_000001635.26_GRCm38.p6_genomic.gtf > splicesites.txt

splicesites.txt

python /global/projectb/scratch/sgordon/conda_envs/.conda/envs/hisat2/bin/hisat2_extract_exons.py GCF_000001635.26_GRCm38.p6_genomic.gtf > exons.txt










usage: hisat2_extract_splice_sites.py [-h] [-v] [gtf_file]

Extract splice junctions from a GTF file

positional arguments:
  gtf_file       input GTF file (use "-" for stdin)

optional arguments:
  -h, --help     show this help message and exit
  -v, --verbose  also print some statistics to stderr

## create conda env for gffread to reformat to GTF
conda create -n gffread -c bioconda gffread

GCF_000001635.26_GRCm38.p6_genomic.gff

source activate gffread

gffread -E GCF_000001635.26_GRCm38.p6_genomic.gff -T -o GCF_000001635.26_GRCm38.p6_genomic.gtf

NC_000067.6     Gnomon  exon    3199733 3207317 .       -       .       transcript_id "rna0"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  exon    3213439 3216968 .       -       .       transcript_id "rna0"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  exon    3421702 3421901 .       -       .       transcript_id "rna0"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  exon    3670552 3672278 .       -       .       transcript_id "rna0"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  CDS     3216022 3216968 .       -       2       transcript_id "rna0"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  CDS     3421702 3421901 .       -       1       transcript_id "rna0"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  CDS     3670552 3671348 .       -       0       transcript_id "rna0"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     BestRefSeq      exon    3214482 3216968 .       -       .       transcript_id "rna1"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     BestRefSeq      exon    3421702 3421901 .       -       .       transcript_id "rna1"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     BestRefSeq      exon    3670552 3671498 .       -       .       transcript_id "rna1"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     BestRefSeq      CDS     3216022 3216968 .       -       2       transcript_id "rna1"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     BestRefSeq      CDS     3421702 3421901 .       -       1       transcript_id "rna1"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     BestRefSeq      CDS     3670552 3671348 .       -       0       transcript_id "rna1"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  exon    3322816 3323760 .       -       .       transcript_id "rna2"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  exon    3421702 3421901 .       -       .       transcript_id "rna2"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  exon    3670552 3672278 .       -       .       transcript_id "rna2"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  CDS     3323747 3323760 .       -       2       transcript_id "rna2"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  CDS     3421702 3421901 .       -       1       transcript_id "rna2"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  CDS     3670552 3671348 .       -       0       transcript_id "rna2"; gene_id "gene0"; gene_name "Xkr4";
NC_000067.6     Gnomon  exon    3361552 3366505 .       +       .       transcript_id "rna3"; gene_id "gene2"; gene_name "Gm39585";
NC_000067.6     Gnomon  exon    3368342 3377812 .       +       .       transcript_id "rna3"; gene_id "gene2"; gene_name "Gm39585";
NC_000067.6     Gnomon  exon    3671524 3672264 .       -       .       transcript_id "gene4"; gene_id "gene4"; gene_name "LOC108167595";
NC_000067.6     Gnomon  exon    4115673 4115948 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4119668 4119712 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4120015 4120073 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4142612 4142766 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4147812 4147963 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4148612 4148744 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4163855 4163941 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4170205 4170404 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4197534 4197641 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4206660 4206837 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4226611 4226823 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4228443 4228619 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4231053 4231144 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4233184 4233287 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  exon    4233436 4233728 .       -       .       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  CDS     4115935 4115948 .       -       2       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  CDS     4119668 4119712 .       -       2       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  CDS     4120015 4120073 .       -       1       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  CDS     4142612 4142766 .       -       0       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  CDS     4147812 4147963 .       -       2       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  CDS     4148612 4148744 .       -       0       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";
NC_000067.6     Gnomon  CDS     4163855 4163941 .       -       0       transcript_id "rna5"; gene_id "gene5"; gene_name "Gm38717";



hisat2 -f -x $HISAT2_HOME/example/index/22_20-21M_snp -1 $HISAT2_HOME/example/reads/reads_1.fa -2 $HISAT2_HOME/example/reads/reads_2.fa -S eg2.sam

  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]

hisat2 --known-splicesite-infile splicesites.txt -x ../reference_sequences/GCF_000001635_26 -1 ../data/JNMV1A_S18/JNMV1A_S18_L3_L4_R1_001.fastq.gz -2 ../data/JNMV1A_S18/JNMV1A_S18_L3_L4_R2_001.fastq.gz -S results/JNMV1A_S18/JNMV1A_S18_L3_L4.sam

--downstream-transcriptome-assembly

nano run_hisat.sh
while IFS= read -r var
do
    echo "$var"
    mkdir -p results/"$var"
    hisat2 --threads 3 --downstream-transcriptome-assembly -x ../reference_sequences/GCF_000001635_26_wsplices -1 ../data/"$var"/"$var"_L3_L4_R1_001.fastq.gz -2 ../data/"$var"/"$var"_L3_L4_R2_001.fastq.gz -S results/"$var"/"$var"_L3_L4_R1_001.sam
done < "sample_prefixes.txt"


# convert the SAMs to sorted BAMs
module load samtools/1.3.1

nano run_samtools.sh
while IFS= read -r var
do
    echo "$var"
    samtools view -Su results/"$var"/"$var"_L3_L4_R1_001.sam | samtools sort - > results/"$var"/"$var"_L3_L4_R1_001.sorted
done < "sample_prefixes.txt"

nohup ./run_samtools.sh 



# Running StringTie

module unload python/2.7.4
module load python/3.6-anaconda_4.3.0
conda create --prefix /global/projectb/scratch/sgordon/conda_envs/stringtie_cf201901 -c bioconda/label/cf201901 stringtie 


### could not get conda env to install, try manual install
git clone https://github.com/gpertea/stringtie
cd stringtie
make release


### ran stringtie on server with output formatted for ballgown
nano run_stringtie.sh
while IFS= read -r var
do
    echo "$var"
    ./stringtie-1.3.5.Linux_x86_64/stringtie results/"$var"/"$var"_L3_L4_R1_001.sorted.bam -G ../reference_sequences/GCF_000001635.26_GRCm38.p6_genomic.gtf -B -e -p 3 -o "$var"_stringtie_results/"$var"_stringtie.gtf
done < "sample_prefixes.txt"

nohup ./run_stringtie.sh &


### ran stringtie on server with output formatted for cufflinks..in progress, can't find online code
nano run_stringtie_for_cufflinks.sh
while IFS= read -r var
do
    echo "$var"
    ./stringtie-1.3.5.Linux_x86_64/stringtie results/"$var"/"$var"_L3_L4_R1_001.sorted.bam -G ../reference_sequences/GCF_000001635.26_GRCm38.p6_genomic.gtf -B -e -p 3 -o "$var"_stringtie_results/"$var"_stringtie.gtf
done < "sample_prefixes.txt"

nohup ./run_stringtie.sh &


# Running Ballgown

module load R/3.3.1

source("http://bioconductor.org/biocLite.R")
biocLite("ballgown")
library(ballgown)

module unload python/2.7.4
module load python/3.6-anaconda_4.3.0
conda create --prefix /global/projectb/scratch/sgordon/conda_envs/bioconductor-ballgown -c bioconda bioconductor-ballgown

conda install -c bioconda 

# Could not get the remote install of ballgown working, do local install

source("http://bioconductor.org/biocLite.R")
biocLite("ballgown")
library(ballgown)

## Example data

library(ballgown)
data_directory = system.file('extdata', package='ballgown') # automatically finds ballgown's installation directory
# examine data_directory:
data_directory

# make the ballgown object:
bg = ballgown(dataDir=data_directory, samplePattern='sample', meas='all')
bg

library(ballgown)
data_directory = system.file('extdata', package='ballgown') 
bg = ballgown(dataDir=data_directory, samplePattern='sample', meas='all')
save(bg, file='bg.rda')

structure(bg)$exon

## Real data, just load simplest experiment

# Loading all the data and saving object
bg = ballgown(dataDir=data_directory, samplePattern='J', meas='all')
save(bg, file='bg_all.rda')
sampleNames(bg)

load(bg_all.rda)

# loading just exp2

library(ballgown)
library(dplyr)
data_directory = "stringtie_results"
bg = ballgown(dataDir=data_directory, samplePattern='J', meas='all')
save(bg, file='bg_exp2.rda')
sampleNames(bg)

[1] "JNMV2A_S30_stringtie_results" "JNMV2B_S31_stringtie_results"
[3] "JNMV2C_S32_stringtie_results" "JNMV2D_S33_stringtie_results"
[5] "JNMV2E_S34_stringtie_results" "JNMV2F_S35_stringtie_results"
[7] "JNMV2G_S36_stringtie_results" "JNMV2H_S37_stringtie_results"

pData(bg) = data.frame(id=sampleNames(bg), group=c(0,1,1,1,1,0,0,0))


head(pData(bg), 3)

320	NC_000067.6	+	16105882	16132550	rna287	6	2295	gene136	Rdh10	12.937691	3.255514


plotMeans('XLOC_000454', bg, groupvar='group', meas='FPKM', colorby='transcript')

stat_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group')
head(stat_results)
?stattest

source("http://bioconductor.org/biocLite.R")
biocLite("dplyr")
library(dplyr)
stat_results_sig <- dplyr::filter(stat_results, qval <= 0.05)

gene10042

plotTranscripts(gene='gene10042', gown=bg, samples='sample12', 
    meas='FPKM', colorby='transcript', 
    main='transcripts from gene XLOC_000454: sample 12, FPKM')
    
plotMeans('gene10042', bg, groupvar='group', meas='FPKM', colorby='transcript')

stat_results_sig <- dplyr::filter(stat_results, pval <= 0.01)
# Write CSV in R
write.csv(stat_results_sig, file = "exp2_pval_stat_results_pval0_01.csv")


load(bg_all.rda)


########## re-load experiment 2 data

library(ballgown)
library(dplyr)
data_directory = "stringtie_results"

load('bg_exp2.rda')
sampleNames(bg)

[1] "JNMV2A_S30_stringtie_results" "JNMV2B_S31_stringtie_results"
[3] "JNMV2C_S32_stringtie_results" "JNMV2D_S33_stringtie_results"
[5] "JNMV2E_S34_stringtie_results" "JNMV2F_S35_stringtie_results"
[7] "JNMV2G_S36_stringtie_results" "JNMV2H_S37_stringtie_results"

pData(bg) = data.frame(id=sampleNames(bg), group=c(0,1,1,1,1,0,0,0))

head(pData(bg), 3)


320	NC_000067.6	+	16105882	16132550	rna287	6	2295	gene136	Rdh10	12.937691	3.255514

plotMeans('gene136', bg, groupvar='group', meas='FPKM', colorby='transcript')
plotMeans('rna287', bg, groupvar='group', meas='FPKM', colorby='transcript')

plotTranscripts('gene136', bg, 
    samples=c('JNMV2A_S30_stringtie_results', 'JNMV2F_S35_stringtie_results', 'JNMV2G_S36_stringtie_results', 'JNMV2H_S37_stringtie_results', 'JNMV2B_S31_stringtie_results', 'JNMV2C_S32_stringtie_results', 'JNMV2D_S33_stringtie_results', 'JNMV2E_S34_stringtie_results'), 
    meas='FPKM', colorby='transcript')

results/JNMV2A_S30
results/JNMV2B_S31
results/JNMV2C_S32
results/JNMV2D_S33
results/JNMV2E_S34
results/JNMV2F_S35
results/JNMV2G_S36
results/JNMV2H_S37

 
stat_results = stattest(bg, feature='gene', getFC=TRUE, meas='FPKM', covariate='group')

library(dplyr)
stat_results_sig <- dplyr::filter(stat_results, pval <= 0.01)
# Write CSV in R
write.csv(stat_results_sig, file = "exp2_gene_pval_stat_results_pval0_01.csv")

## do the same for transcripts
stat_results = stattest(bg, feature='transcript', getFC=TRUE, meas='FPKM', covariate='group')
# Write CSV in R
write.csv(stat_results_sig, file = "exp2_transcript_stat_results_pval0_01.csv")





# loading just exp1, both sexes

library(ballgown)
library(dplyr)
data_directory = "stringtie_results"


#### could not get this to work
# load the sample information
list.files(path = ".",)
pheno_data <- read.csv("stringtie_results/exp1_ballgown_pheno_data.csv")
head(pheno_data)


# create a ballgown object
bg <- ballgown(dataDir = "data_directory",
                    samplePattern = "J",
                    meas='all',
                    pData = pheno_data)
                

# use this
bg = ballgown(dataDir=data_directory, samplePattern='J', meas='all')
save(bg, file='bg_exp1.rda')
sampleNames(bg)

 [1] "JNMV1A_S18_stringtie_results" "JNMV1B_S19_stringtie_results"
 [3] "JNMV1C_S20_stringtie_results" "JNMV1D_S21_stringtie_results"
 [5] "JNMV1E_S22_stringtie_results" "JNMV1F_S23_stringtie_results"
 [7] "JNMV1G_S24_stringtie_results" "JNMV1H_S25_stringtie_results"
 [9] "JNMV1I_S26_stringtie_results" "JNMV1J_S27_stringtie_results"
[11] "JNMV1K_S28_stringtie_results" "JNMV1L_S29_stringtie_results"


# load rda
library(ballgown)
library(dplyr)
data_directory = "stringtie_results"

load('bg_exp1.rda')
sampleNames(bg)

pData(bg) = data.frame(id=sampleNames(bg), genotype=c(0,1,0,0,1,0,1,1,1,0,0,1),sex=c(0,0,1,1,1,1,1,1,0,0,0,0))

head(pData(bg), 3)

320	NC_000067.6	+	16105882	16132550	rna287	6	2295	gene136	Rdh10	12.937691	3.255514


plotMeans('gene136', bg, groupvar='genotype', meas='FPKM', colorby='transcript')

# filter on variance before calc DE
library(genefilter)
bg_filt <- subset(bg, "rowVars(texpr(bg)) >1", genomesubset=TRUE)
#### ballgown instance with 13903 transcripts and 12 samples


# test on transcripts
results_transcripts <- stattest(bg_filt,
                                feature="transcript",
                                covariate="genotype",
                                adjustvars = c("sex"),
                                getFC=TRUE, meas="FPKM")

table(results_transcripts$pval<0.01)

# the order is the same so we can simply combine the information
results_transcripts <- data.frame(geneNames = geneNames(bg_filt),
                                  geneIDs = geneIDs(bg_filt),
                                  transcriptIDs = transcriptNames(bg_filt),
                                  results_transcripts)
stat_results_sig <- dplyr::filter(results_transcripts, pval <= 0.01)

# Write CSV in R for genotype with geneID and transcriptID, filtered on pval
write.csv(stat_results_sig, file = "genotype_ballgown_exp1_transcripts_stat_results_pval0_01.csv")

stat_results_sig_fc <- stat_results_sig %>% filter(fc < 0.5 | fc > 2)
# Write CSV in R for genotype with geneID and transcriptID, filtered on pval and FC
write.csv(stat_results_sig_fc, file = "genotype_ballgown_exp1_transcripts_stat_results_pval0_01_FC2.csv")

plotMeans('gene29382', bg, groupvar='genotype', meas='FPKM', colorby='transcript')
http://uswest.ensembl.org/Mus_musculus/Gene/Family?db=core;g=ENSMUSG00000025144;r=11:120710942-120713738

plotMeans('gene39548', bg, groupvar='genotype', meas='FPKM', colorby='transcript')

plotMeans('gene42935', bg, groupvar='genotype', meas='FPKM', colorby='transcript')





library(ggplot2)
library(cowplot)

results_transcripts$mean <- rowMeans(texpr(bg_filt))
 
ggplot(results_transcripts, aes(log2(mean), log2(fc), colour = pval<0.01)) +
  scale_color_manual(values=c("#999999", "#FF0000")) +
  geom_point() +
  geom_hline(yintercept=0)
  
  

           
# test on genes for genotype
results_genes <- stattest(bg,
                                feature="gene",
                                covariate="genotype",
                                adjustvars = c("sex"),
                                getFC=TRUE, meas="FPKM")

table(results_genes$pval<0.05)
library(dplyr)
stat_results_sig <- dplyr::filter(results_genes, pval <= 0.01)

rdh10 <- dplyr::filter(results_genes, id == "gene136")

  feature      id        fc      pval     qval
1    gene gene136 0.8243779 0.2584735 0.830791


## Plot all the samples individually, grouped by genotype
plotTranscripts('gene136', bg, 
    samples=c('JNMV1A_S18_stringtie_results', 'JNMV1C_S20_stringtie_results', 'JNMV1D_S21_stringtie_results', 'JNMV1F_S23_stringtie_results', 'JNMV1J_S27_stringtie_results', 'JNMV1K_S28_stringtie_results', 'JNMV1B_S19_stringtie_results', 'JNMV1E_S22_stringtie_results', 'JNMV1G_S24_stringtie_results', 'JNMV1H_S25_stringtie_results', 'JNMV1I_S26_stringtie_results', 'JNMV1L_S29_stringtie_results'), 
    meas='FPKM', colorby='transcript')

stat_results_sig <- dplyr::filter(results_genes, pval <= 0.01)

# Write CSV in R
write.csv(stat_results_sig, file = "ballgown_exp1_genes_stat_results_pval0_01.csv")

stat_results_sig <- dplyr::filter(results_transcripts, pval <= 0.01)

# Write CSV in R
write.csv(stat_results_sig, file = "ballgown_exp1_transcripts_stat_results_pval0_01.csv")


# test on genes for sex
results_genes <- stattest(bg,
                                feature="gene",
                                covariate="sex",
                                adjustvars = c("genotype"),
                                getFC=TRUE, meas="FPKM")
stat_results_sig <- dplyr::filter(results_genes, qval <= 0.05)
dim(stat_results_sig)
[1] 332   5

# Write CSV in R for sex
write.csv(stat_results_sig, file = "sex_ballgown_exp1_gene_stat_results_pval0_01.csv")

                  
# test on transcripts for sex
results_transcripts <- stattest(bg,
                                feature="transcript",
                                covariate="sex",
                                adjustvars = c("genotype"),
                                getFC=TRUE, meas="FPKM")
                                

# the order is the same so we can simply combine the information
results_transcripts <- data.frame(geneNames = geneNames(bg),
                                  geneIDs = geneIDs(bg),
                                  results_transcripts)
stat_results_sig <- dplyr::filter(results_transcripts, qval <= 0.05)

# Write CSV in R for sex
write.csv(stat_results_sig, file = "sex_ballgown_exp1_transcripts_stat_results_pval0_01.csv")


library(ggplot2)
library(cowplot)


bg_filt <- subset(bg, "rowVars(texpr(bg)) >1", genomesubset=TRUE)
#### ballgown instance with 13903 transcripts and 12 samples

results_transcripts$mean <- rowMeans(texpr(bg_filt))
 
ggplot(results_transcripts, aes(log2(mean), log2(fc), colour = qval<0.05)) +
  scale_color_manual(values=c("#999999", "#FF0000")) +
  geom_point() +
  geom_hline(yintercept=0)




# loading just exp1, male data only

library(ballgown)
library(dplyr)
data_directory = "stringtie_results"
                
# load the data and save
bg = ballgown(dataDir=data_directory, samplePattern='J', meas='all')
save(bg, file='bg_exp1_male_only.rda')
sampleNames(bg)

# re-load rda later
library(ballgown)
library(dplyr)
data_directory = "stringtie_results"

load('bg_exp1.rda')
sampleNames(bg)

[1] "JNMV1A_S18_stringtie_results" "JNMV1B_S19_stringtie_results"
[3] "JNMV1I_S26_stringtie_results" "JNMV1J_S27_stringtie_results"
[5] "JNMV1K_S28_stringtie_results" "JNMV1L_S29_stringtie_results"

# genotype of male only
pData(bg) = data.frame(id=sampleNames(bg), genotype=c(0,1,1,0,0,1))

head(pData(bg), 3)

320	NC_000067.6	+	16105882	16132550	rna287	6	2295	gene136	Rdh10	12.937691	3.255514


plotMeans('gene136', bg, groupvar='genotype', meas='FPKM', colorby='transcript')

library(genefilter)
bg_filt <- subset(bg, "rowVars(texpr(bg)) >1", genomesubset=TRUE)

# test on transcripts
results_transcripts <- stattest(bg_filt,
                                feature="transcript",
                                covariate="genotype",
                                getFC=TRUE, meas="FPKM")

results_gene <- stattest(bg_filt,
                                feature="gene",
                                covariate="genotype",
                                getFC=TRUE, meas="FPKM")

table(results_transcripts$qval<0.05)
table(results_gene$qval<0.05)

table(results_transcripts$pval<0.01)

# the order is the same so we can simply combine the information
results_transcripts <- data.frame(geneNames = geneNames(bg_filt),
                                  geneIDs = geneIDs(bg_filt),
                                  transcriptIDs = transcriptNames(bg_filt),
                                  results_transcripts)
stat_results_sig <- dplyr::filter(results_transcripts, pval <= 0.01)

# Write CSV in R for genotype with geneID and transcriptID, filtered on pval
write.csv(stat_results_sig, file = "genotype_ballgown_exp1_male_only_transcripts_stat_results_pval0_01.csv")

stat_results_sig_fc <- stat_results_sig %>% filter(fc < 0.5 | fc > 2)
# Write CSV in R for genotype with geneID and transcriptID, filtered on pval and FC
write.csv(stat_results_sig_fc, file = "genotype_ballgown_exp1_male_only_transcripts_stat_results_pval0_01_FC2.csv")



# loading just exp1, female data only

library(ballgown)
library(dplyr)
data_directory = "stringtie_results"
                
# load the data and save
bg = ballgown(dataDir=data_directory, samplePattern='J', meas='all')
save(bg, file='bg_exp1_female_only.rda')
sampleNames(bg)
[1] "JNMV1C_S20_stringtie_results" "JNMV1D_S21_stringtie_results"
[3] "JNMV1E_S22_stringtie_results" "JNMV1F_S23_stringtie_results"
[5] "JNMV1G_S24_stringtie_results" "JNMV1H_S25_stringtie_results"

# re-load rda later
library(ballgown)
library(dplyr)
data_directory = "stringtie_results"
load('bg_exp1.rda')
sampleNames(bg)

# genotype of female only
pData(bg) = data.frame(id=sampleNames(bg), genotype=c(0,0,1,0,1,1))

head(pData(bg), 3)

320	NC_000067.6	+	16105882	16132550	rna287	6	2295	gene136	Rdh10	12.937691	3.255514


plotMeans('gene136', bg, groupvar='genotype', meas='FPKM', colorby='transcript')

library(genefilter)
bg_filt <- subset(bg, "rowVars(texpr(bg)) >1", genomesubset=TRUE)

# test on transcripts
results_transcripts <- stattest(bg,
                                feature="transcript",
                                covariate="genotype",
                                getFC=TRUE, meas="FPKM")

results_gene <- stattest(bg,
                                feature="gene",
                                covariate="genotype",
                                getFC=TRUE, meas="FPKM")

table(results_transcripts$qval<0.05)
table(results_gene$qval<0.05)

table(results_transcripts$pval<0.01)

# the order is the same so we can simply combine the information
results_transcripts <- data.frame(geneNames = geneNames(bg),
                                  geneIDs = geneIDs(bg),
                                  transcriptIDs = transcriptNames(bg),
                                  results_transcripts)
stat_results_sig <- dplyr::filter(results_transcripts, pval <= 0.01)

# Write CSV in R for genotype with geneID and transcriptID, filtered on pval
write.csv(stat_results_sig, file = "genotype_ballgown_exp1_female_only_transcripts_stat_results_pval0_01.csv")

stat_results_sig_fc <- stat_results_sig %>% filter(fc < 0.5 | fc > 2)
# Write CSV in R for genotype with geneID and transcriptID, filtered on pval and FC
write.csv(stat_results_sig_fc, file = "genotype_ballgown_exp1_female_only_transcripts_stat_results_pval0_01_FC2.csv")


