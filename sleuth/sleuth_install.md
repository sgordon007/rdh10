## activate environment

source $BSCRATCH/mus_2019_01_29/conda/anaconda2/bin/activate

conda create --name sleuth

source activate sleuth

conda install --channel bioconda r-sleuth

source("http://bioconductor.org/biocLite.R")
biocLite("devtools")    # only if devtools not yet installed
biocLite("pachterlab/sleuth")

#### activate the library
library('sleuth')
library('cowplot')

#### install biomaRt
https://bioconductor.org/packages/release/bioc/html/biomaRt.html


We recommend starting with the vignette:

vignette('intro', package = 'sleuth')


## Tutorial
https://github.com/pachterlab/sleuth

https://pachterlab.github.io/sleuth/walkthroughs


https://pachterlab.github.io/sleuth_walkthroughs/pval_agg/analysis.html



