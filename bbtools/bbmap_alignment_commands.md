maxindel=16000          Don't look for indels longer than this. Lower is faster.
                        Set to >=100k for RNAseq with long introns like mammals.
                        
maxindel=160000

slow=t

samtools/1.4




### Each sample is run separately through kallisto quant.  Do this in a shell script:

cd $BSCRATCH/mus_2019_01_29/kallisto_run

cp ../data/sample_prefixes.txt ./

nano run_bbmap.sh


./run_bbmap.sh