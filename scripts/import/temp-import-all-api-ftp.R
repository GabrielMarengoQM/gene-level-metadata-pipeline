# Import API/FTP data 
library(biomaRt)
library(dplyr)
library(PANTHER.db)
library(STRINGdb)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)
library(RJSONIO)
library(tidyverse)
library(arrow)

# Keys / tokens / version numbers
impc_release_version <- "22.1"
omim_api_key <- NA # hidden env / gitignore
lethal_genes_token <- "" # hidden env / gitignore
g2p_folder_file <- "2025_02_28/DDG2P_2025-02-28.csv.gz"
panelapp_max <- 338

# Logs
log_errors <- list()

# Genomic & Sequence Features ----

# HGNC genes file; Gene IDs / Prev/alias names / gene name / gene group
tryCatch({
  protein.coding.genes <- utils::read.delim('https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/locus_types/gene_with_protein_product.txt')
  arrow::write_parquet(protein.coding.genes, './data/raw/api-ftp/protein.coding.genes.parquet')
}, error = function(e) {
  log_errors[["HGNC_protein_coding_genes"]] <<- conditionMessage(e)
})

# chromosome, positions, length, GC content
tryCatch({
  my_genes <- protein.coding.genes$symbol
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_position_gc_content <- biomaRt::getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "percentage_gene_gc_content"),
    filters = "hgnc_symbol",
    values = my_genes,
    mart = ensembl
  )
  arrow::write_parquet(gene_position_gc_content, './data/raw/api-ftp/gene_position_gc_content.parquet')
}, error = function(e) {
  log_errors[["biomaRt_gene_position_gc"]] <<- conditionMessage(e)
})

# Proteins/Proteomics ----

# PANTHER protein classes
tryCatch({
  PANTHER.db::pthOrganisms(PANTHER.db) <- "HUMAN"
  cols_short <- c("UNIPROT", "CLASS_ID", "CLASS_TERM", "FAMILY_ID", "FAMILY_TERM", "SUBFAMILY_TERM")
  uniprot_ids <- PANTHER.db::keys(PANTHER.db, keytype="UNIPROT")
  panther_data_shortcols <- PANTHER.db::select(PANTHER.db, 
                                               keys=uniprot_ids, 
                                               columns=cols_short, 
                                               keytype="UNIPROT")
  arrow::write_parquet(panther_data_shortcols, './data/raw/api-ftp/panther_data_shortcols.parquet')
}, error = function(e) {
  log_errors[["PANTHER_protein_classes"]] <<- conditionMessage(e)
})

# string ppi interactions
tryCatch({
  hgnc_ensembl <- protein.coding.genes %>% dplyr::select(hgnc_id, ensembl_gene_id)
  string_db <- STRINGdb::STRINGdb$new(version="12.0", species=9606,
                                      score_threshold=700, network_type="full", 
                                      input_directory="")
  input_genes_mapped <- string_db$map(data.frame(hgnc_ensembl), "ensembl_gene_id", removeUnmappedRows = TRUE )
  input_genes_mapped_vector <- input_genes_mapped %>% pull(STRING_id)
  interactions <- string_db$get_interactions(input_genes_mapped_vector)
  arrow::write_parquet(interactions, './data/raw/api-ftp/ppi_interactions.parquet')
}, error = function(e) {
  log_errors[["STRINGdb_interactions"]] <<- conditionMessage(e)
})

# Mouse Perturbation Assays ----

# impc; viability / phenotypes 
tryCatch({
  mouse_viability_impc <- data.table::fread(paste0(
    "http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-",
    impc_release_version,
    "/results/viability.csv.gz"
  ))
  arrow::write_parquet(mouse_viability_impc, './data/raw/api-ftp/mouse_viability_impc.parquet')
}, error = function(e) {
  log_errors[["IMPC_viability"]] <<- conditionMessage(e)
})

tryCatch({
  mouse_phenotypes_impc <- data.table::fread(paste0(
    "http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-",
    impc_release_version,
    "/results/genotype-phenotype-assertions-ALL.csv.gz"
  ))
  arrow::write_parquet(mouse_phenotypes_impc, './data/raw/api-ftp/mouse_phenotypes_impc.parquet')
}, error = function(e) {
  log_errors[["IMPC_phenotypes"]] <<- conditionMessage(e)
})

# mgi; viability / phenotypes 
tryCatch({
  mouse_viability_mgi <- data.table::fread(
    "https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt"
  )
  arrow::write_parquet(mouse_viability_mgi, './data/raw/api-ftp/mouse_viability_mgi.parquet')
}, error = function(e) {
}, error = function(e) {
  log_errors[["MGI_viability"]] <<- conditionMessage(e)
})

tryCatch({
  mouse_phenotypes_mgi <- data.table::fread(
    header = FALSE,
    "https://www.informatics.jax.org/downloads/reports/VOC_MammalianPhenotype.rpt"
  )
  arrow::write_parquet(mouse_phenotypes_mgi, './data/raw/api-ftp/mouse_phenotypes_mgi.parquet')
}, error = function(e) {
  log_errors[["MGI_phenotypes"]] <<- conditionMessage(e)
})

tryCatch({
  mouse_protein_coding_genes_mgi <- data.table::fread(
    "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"
  )
  arrow::write_parquet(mouse_protein_coding_genes_mgi, './data/raw/api-ftp/mouse_protein_coding_genes_mgi.parquet')
}, error = function(e) {
  log_errors[["MGI_protein_coding_genes"]] <<- conditionMessage(e)
})

# Functional Annotation ----

# go
tryCatch({
  go_anot <- toTable(org.Hs.egGO) %>%
    dplyr::select(gene_id, go_id, Ontology) %>%
    distinct()
  keys <- unique(go_anot$go_id)
  go_terms <- AnnotationDbi::select(GO.db, keys = keys, keytype="GOID", columns=c("TERM")) %>%
    dplyr::rename(go_id = GOID, go_term = TERM)
  go_anot_term_raw <- go_anot %>% left_join(go_terms)
  arrow::write_parquet(go_anot_term_raw, './data/raw/api-ftp/go_anot_term_raw.parquet')
}, error = function(e) {
  log_errors[["GO_annotations"]] <<- conditionMessage(e)
})

# reactome
tryCatch({
  genes_universe <- protein.coding.genes %>% dplyr::select(symbol, entrez_id)
  universe <- as.character(genes_universe$entrez_id)
  entrez_pathid <- as.data.frame(reactomeEXTID2PATHID)
  pathid_name <- as.data.frame(reactomePATHID2NAME)
  reactome_annotations <- entrez_pathid %>% full_join(pathid_name)
  arrow::write_parquet(reactome_annotations, './data/raw/api-ftp/reactome_annotations.parquet')
}, error = function(e) {
  log_errors[["Reactome_annotations"]] <<- conditionMessage(e)
})

# Human Disease & Phenotypes ----

# omim
tryCatch({
  mim.Titles <- read.delim(paste0("https://data.omim.org/downloads/", omim_api_key, "/mimTitles.txt"), skip = 2)
  arrow::write_parquet(mim.Titles, './data/raw/api-ftp/mim.Titles.parquet')
}, error = function(e) {
  log_errors[["OMIM_mimTitles"]] <<- conditionMessage(e)
})

tryCatch({
  morbidmap <- read.delim(paste0("https://data.omim.org/downloads/", omim_api_key, "/morbidmap.txt"), skip = 3)
  arrow::write_parquet(morbidmap, './data/raw/api-ftp/morbidmap.parquet')
}, error = function(e) {
  log_errors[["OMIM_morbidmap"]] <<- conditionMessage(e)
})

tryCatch({
  genemap <- read.delim(paste0("https://data.omim.org/downloads/", omim_api_key, "/genemap2.txt"), skip = 3)
  arrow::write_parquet(genemap, './data/raw/api-ftp/genemap.parquet')
}, error = function(e) {
  log_errors[["OMIM_genemap"]] <<- conditionMessage(e)
})

# omim lethal phenotypes
tryCatch({
  lethal.phenotypes <- fread(paste0(
    "https://raw.github.qmul.ac.uk/whri-phenogenomics/catalogue_lethal_genes/refs/heads/master/catalogue_lethal_genes_app/data/omim_curation.tsv?token=",
    lethal_genes_token
  ))
  arrow::write_parquet(lethal.phenotypes, './data/raw/api-ftp/lethal.phenotypes.parquet')
}, error = function(e) {
  log_errors[["OMIM_lethal_phenotypes"]] <<- conditionMessage(e)
})

# g2p
tryCatch({
  dd.g2p <- fread(paste0(
    "http://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/",
    g2p_folder_file
  ))
  arrow::write_parquet(dd.g2p, './data/raw/api-ftp/dd.g2p.parquet')
}, error = function(e) {
  log_errors[["G2P_ddg2p"]] <<- conditionMessage(e)
})

# panelapp
tryCatch({
  jfiles_unlist <- NULL
  for (i in 1:panelapp_max){
    page_number <- i 
    page <- paste0("https://panelapp-aus.org/api/v1/genes/?page=", page_number, "")
    jfiles <- fromJSON(page)
    jfiles_unlist [[i]] <- enframe(unlist(jfiles)) 
    message( i, ' of ', panelapp_max)
  }
  write_rds(jfiles_unlist, './data/raw/api-ftp/jfiles_unlist.rds') # list so rds not parquet which request table format
}, error = function(e) {
  log_errors[["PanelApp_download"]] <<- conditionMessage(e)
})

# Comparative/Evolutionary Genomics ----

# human-mouse orthologs
tryCatch({
  mouse_orthologs <- data.table::fread("https://www.gentar.org/orthology-api/api/ortholog/write_to_tsv_file")
  arrow::write_parquet(mouse_orthologs, './data/raw/api-ftp/mouse_orthologs.parquet')
}, error = function(e) {
  log_errors[["Mouse_orthologs"]] <<- conditionMessage(e)
})

# human paralogs
tryCatch({
  human<- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  atts <- listAttributes(human) %>% 
    filter(stringr::str_detect(name, "paralog"))
  human_paralogs <- getBM(c('ensembl_gene_id', atts$name),
                          filters = "ensembl_gene_id",
                          values = protein.coding.genes$ensembl_gene_id, mart = human, uniqueRows=TRUE)
  arrow::write_parquet(human_paralogs, './data/raw/api-ftp/human_paralogs.parquet')
}, error = function(e) {
  log_errors[["human_paralogs"]] <<- conditionMessage(e)
})

# Final summary ----
if (length(log_errors) > 0) {
  cat("\n\u274C The following steps failed:\n")
  print(log_errors)
} else {
  cat("\n\u2705 All steps completed successfully.\n")
}