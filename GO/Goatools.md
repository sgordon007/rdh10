curl http://snapshot.geneontology.org/ontology/go-basic.obo > go-basic.obo

# also downloaded goslim_mouse.obo

conda create -n goatools goatools

conda install goatools

source activate goatools


# Get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
from goatools.base import download_ncbi_associations
gene2go = download_ncbi_associations()

