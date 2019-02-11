## History on working with Kallisto for RNA-Seq

## Install

## first install conda

wget -c https://repo.continuum.io/archive/Anaconda2-5.0.1-Linux-x86_64.sh

bash Anaconda2-5.0.1-Linux-x86_64.sh

source $BSCRATCH/mus_2019_01_29/conda/anaconda2/bin/activate

conda create --name kallisto

source activate kallisto

conda install kallisto

## Tutorial
https://pachterlab.github.io/kallisto/starting.html

kallisto

### Uploaded test data from source download to test install

Run the test:

kallisto index -i transcripts.idx transcripts.fasta.gz

kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq.gz reads_2.fastq.gz

The results of the main quantification, i.e. the abundance estimate using kallisto on the data is in the abundance.tsv file. 
Abundances are reported in “estimated counts” (est_counts) and in Transcripts Per Million (TPM). 


The file is tab delimited so that it can easily parsed. The output can also be analyzed with the sleuth tool.

The run_info.json file contains a summary of the run, including data on the number targets used for quantification, 
the number of bootstraps performed, the version of the program used and how it was called. 

## More Information

https://pachterlab.github.io/kallisto-sleuth-workshop-2016/