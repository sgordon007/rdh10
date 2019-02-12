
Here is the sample ID table that you need. (attached excel file)
1. For experiment 1 (1A to 1L), we need comparison WT vs KO for all mice, and then separately for males and females. 
This means the comparison is first: 1A+1J+1K + 1C+1D+1F    vs   1B +1I+1L +1E +1G +1H

Then: 1A+1J+1K  vs  1B +1I+1L

Then:  1C+1D+1F  vs  1E +1G +1H

2. For Experiment 2 (2A to 2H) we need WT vs KO comparison: 2A +2F+2F+2H  vs  2B+2C+2D+2E

3. For experiment 3 (3A to 3L) we need comparison WT vs KO for all mice, and then for males and females separately:

3D+3E+3F+ 3G+3H+3I   vs  3A+3B+3C+3J+3K+3L

then: 3D+3E+3F  vs.  3A+3B+3C

and then: 3G+3H+3I  vs.  3J+3K+3L



metadata <- read.table('Rdh10KO_liver_BAT_RNAseq_SampleID_MV_2019-02-10_more_data.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

