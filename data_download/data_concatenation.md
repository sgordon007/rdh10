# Cat fastq files spread across the two lanes:

#### get a list of sample prefixes
for f in *L003_R1_001.fastq.gz; do
    echo "${f}"
    STEM=$(basename "${f}" _L003_R1_001.fastq.gz)
    echo "${STEM}" >> sample_prefixes.txt
done

32 sample_prefixes.txt

JNMV1A_S18_L003_R1_001.fastq.gz

### run the job for R1
while IFS= read -r var
do
    echo "$var"
    mkdir "$var"
    cat "$var"_L003_R1_001.fastq.gz "$var"_L004_R1_001.fastq.gz > ./"$var"/"$var"_L3_L4_R1_001.fastq.gz
done < "sample_prefixes.txt"

JNMV1A_S18_L003_R2_001.fastq.gz

### run the job for R2
while IFS= read -r var
do
    echo "$var"
    mkdir "$var"
    cat "$var"_L003_R2_001.fastq.gz "$var"_L004_R2_001.fastq.gz > ./"$var"/"$var"_L3_L4_R2_001.fastq.gz
done < "sample_prefixes.txt"

