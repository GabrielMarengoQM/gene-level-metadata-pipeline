# Import web-files (accessed manually from web)
library(utils)

# DepMap
models <- utils::read.csv("./data/raw/web-files/Model.csv")
gene_effect <- utils::read.csv("./data/raw/web-files/CRISPRGeneEffect.csv")

# GTEx
gtex_expression <- utils::read.delim("./data/raw/web-files/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz", skip=2)

# OGEE/STRING
ppi_connectivity <- read.delim("./data/raw/web-files/connectivity.txt.gz")

# dbNSFP/GWAS/HI
dbNSFP <- read_delim('./data/raw/web-files/dbNSFP5.1_gene')

# gnomAD
gnomad <- read.table("./data/raw/web-files/gnomad.v4.1.constraint_metrics.tsv", head = TRUE)

# Protein atlas bulk expression/subcellular location/protein class
hpa_bulk_expression <- read_tsv("./data/raw/web-files/normal_tissue.tsv.zip")
hpa_protein_classes <- read_tsv("./data/raw/web-files/proteinatlas.tsv.zip")
