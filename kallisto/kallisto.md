## History on working with Kallisto for RNA-Seq

## Install

## activate environment

source mus_2019_01_29/conda/anaconda2/bin/activate

source activate kallisto

Run Analysis:

##### probably not what we want: kallisto index -i GCF_000001635_transcripts.idx GCF_000001635.26_GRCm38.p6_rna.fna.gz


#### this one: 
kallisto index -i Mus_musculus_cdna.idx Mus_musculus.GRCm38.cdna.all.fa.gz

kallisto quant -i Mus_musculus_cdna.idx -o results/JNMV1A_S18 -b 100 -t 8 ../data/JNMV1A_S18_L003_R1_001.fastq.gz ../data/JNMV1A_S18_L003_R2_001.fastq.gz

### Each sample is run separately through kallisto quant.  Do this in a shell script:

cd $BSCRATCH/mus_2019_01_29/kallisto_run

cp ../data/sample_prefixes.txt ./

nano run_kallisto.sh

while IFS= read -r var
do
    echo "$var"
    mkdir -p results/"$var"
    kallisto quant -i ../reference_sequences/Mus_musculus_cdna.idx -o results/"$var" -b 100 -t 8 ../data/"$var"/"$var"_L3_L4_R1_001.fastq.gz ../data/"$var"/"$var"_L3_L4_R2_001.fastq.gz
done < "sample_prefixes.txt"

nohup ./run_kallisto.sh


The results of the main quantification, i.e. the abundance estimate using kallisto on the data is in the abundance.tsv file. 
Abundances are reported in “estimated counts” (est_counts) and in Transcripts Per Million (TPM). 


The file is tab delimited so that it can easily parsed. The output can also be analyzed with the sleuth tool.

The run_info.json file contains a summary of the run, including data on the number targets used for quantification, 
the number of bootstraps performed, the version of the program used and how it was called. 

GCF_000001635.26_GRCm38.p6_rna.fna.gz

