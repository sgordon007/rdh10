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

We recommend starting with the vignette:

vignette('intro', package = 'sleuth')


## Tutorial
https://github.com/pachterlab/sleuth

https://pachterlab.github.io/sleuth/walkthroughs




