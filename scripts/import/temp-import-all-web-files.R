# Import web-files (accessed manually from web)
library(utils)
library(readr)

# DepMap
models <- utils::read.csv("./data/raw/web-files/Model.csv")
gene_effect <- utils::read.csv("./data/raw/web-files/CRISPRGeneEffect.csv")

# GTEx
gtex_expression <- utils::read.delim("./data/raw/web-files/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz", skip=2)

# OGEE/STRING
ppi_connectivity <- utils::read.delim("./data/raw/web-files/connectivity.txt.gz")

# dbNSFP/GWAS/HI
dbNSFP <- readr::read_delim('./data/raw/web-files/dbNSFP5.1_gene')

# gnomAD
gnomad <- read.table("./data/raw/web-files/gnomad.v4.1.constraint_metrics.tsv", head = TRUE)

# Protein atlas; bulk expression / subcellular location / protein class / 
hpa_bulk_expression <- readr::read_tsv("./data/raw/web-files/normal_tissue.tsv.zip")
hpa_atlas <- readr::read_tsv("./data/raw/web-files/proteinatlas.tsv.zip")

# lymphoblastoid_time_series_expression across 7 organs
lymphoblastoid_time_series_expression <- readr::read_tsv("./data/raw/web-files/E-MTAB-6814-query-results.tpms.tsv", skip = 4)

# mouse lethal genes
mp_lethal_terms <- utils::read.csv("./data/raw/web-files/lethal_terms.csv")

# impc window of lethality
wol <- read.delim("./data/raw/web-files/sup_file_1.txt")

# human phenotype ontology
hpo_genes_to_phenotype <- readr::read_tsv("./data/raw/web-files/genes_to_phenotype.txt") %>% 
  dplyr::select(gene_symbol, hpo_id, hpo_name)
hpo_genes_to_dissease <- readr::read_tsv("./data/raw/web-files/genes_to_disease.txt") %>% 
  dplyr::select(gene_symbol, association_type, disease_id)

