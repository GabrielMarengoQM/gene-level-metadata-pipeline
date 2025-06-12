library(dplyr)
library(tidyr)
library(STRINGdb)
library(arrow)
library(RJSONIO)

# Genomic & Sequence Features ----

# HGNC genes file; Gene IDs / Prev/alias names / gene name / gene group / position / length / GC content
protein.coding.genes <- arrow::read_parquet('./data/raw/api-ftp/protein.coding.genes.parquet')

# Gene IDs
gene_ids <- protein.coding.genes %>% 
  dplyr::select(symbol, hgnc_id, entrez_id, ensembl_gene_id, vega_id, ucsc_id, omim_id) 

# Prev/alias names
prev_names <- protein.coding.genes %>% 
  dplyr::select(symbol, prev_symbol) %>% 
  separate_rows(prev_symbol, sep = "\\|")

alias_names <- protein.coding.genes %>% 
  dplyr::select(symbol, alias_symbol) %>% 
  separate_rows(alias_symbol, sep = "\\|")

# human mouse mappings
mouse_ids <- protein.coding.genes %>% 
  dplyr::select(symbol, mgd_id) %>% 
  tidyr::separate_rows(mgd_id, sep = '\\|')

# gene uniprot protein mappings
uniprot_ids <- protein.coding.genes %>% 
  dplyr::select(symbol, uniprot_ids) %>%
  tidyr::separate_rows(uniprot_ids, sep = '\\|')

# gene name
gene_names <- protein.coding.genes %>% 
  dplyr::select(symbol, name)

# gene group
gene_groups <- protein.coding.genes %>% 
  dplyr::select(symbol, gene_group) %>% 
  tidyr::separate_rows(gene_group, sep = "\\|") %>% 
  filter(!is.na(gene_group)) %>% 
  filter(gene_group != "")

# chromosome, positions, length
gene_position_length<- arrow::read_parquet('./data/raw/api-ftp/gene_position_gc_content.parquet') %>% 
  dplyr::mutate(gene_length = end_position - start_position + 1) %>% 
  dplyr::select(-percentage_gene_gc_content)

# GC content
gene_gc_content <- arrow::read_parquet('./data/raw/api-ftp/gene_position_gc_content.parquet') %>% 
  dplyr::select(hgnc_symbol, percentage_gene_gc_content)

# Proteins/Proteomics ----

# PANTHER protein classes
panther_data_shortcols <- arrow::read_parquet('./data/raw/api-ftp/panther_data_shortcols.parquet') %>% 
  left_join(uniprot_ids, by = c("UNIPROT" = "uniprot_ids"))

# string ppi interactions
interactions <- arrow::read_parquet('./data/raw/api-ftp/ppi_interactions.parquet')
hgnc_ensembl <- protein.coding.genes %>%
  dplyr::select(hgnc_id, ensembl_gene_id)
string_db <- STRINGdb$new(version="12.0", species=9606,
                          score_threshold=700, network_type="full", 
                          input_directory="")
input_genes_mapped <- string_db$map(data.frame(hgnc_ensembl), "ensembl_gene_id", removeUnmappedRows = TRUE )
interactions_cleaned <- input_genes_mapped %>% 
  left_join(interactions, join_by(STRING_id == from), 
            relationship = "many-to-many") %>%
  dplyr::select(hgnc_id, STRING_id, to, combined_score) %>% 
  dplyr::rename(protein1_hgnc_id = hgnc_id) %>% 
  dplyr::rename(protein1_string_id = STRING_id) %>%
  dplyr::rename(protein2_string_id = to)
head(interactions_cleaned)

interactions_cleaned_hgnc <- input_genes_mapped %>% 
  left_join(interactions_cleaned, join_by(STRING_id == protein2_string_id), 
            relationship = "many-to-many") %>% 
  dplyr::select(protein1_hgnc_id, protein1_string_id, hgnc_id, STRING_id, combined_score) %>%
  dplyr::rename(protein2_hgnc_id = hgnc_id) %>%
  dplyr::rename(protein2_string_id = STRING_id) %>% 
  dplyr::arrange(protein1_hgnc_id)

gene_symbol_hgnc_mappings <- protein.coding.genes %>%
  dplyr::select(symbol, hgnc_id) %>%
  dplyr::rename(gene_symbol = symbol) %>%
  dplyr::rename(protein1_hgnc_id = hgnc_id)

symbol_hgnc_mappings <- protein.coding.genes %>% 
  dplyr::select(symbol, hgnc_id) %>% 
  dplyr::rename(protein2_hgnc_id = hgnc_id, protein2_gene_symbol = symbol)

string_ppi_df <- interactions_cleaned_hgnc %>%
  mutate(
    protein1_string_id = str_replace(protein1_string_id, "^9606\\.", ""),
    protein2_string_id = str_replace(protein2_string_id, "^9606\\.", "")
  ) %>%
  left_join(gene_symbol_hgnc_mappings) %>% 
  left_join(symbol_hgnc_mappings) %>% 
  filter(!is.na(combined_score)) %>% 
  mutate(combined_score = combined_score/1000) %>% 
  dplyr::select(gene_symbol, protein1_string_id, protein2_string_id, combined_score)

# Mouse Perturbation Assays ----

# impc viability 
mouse_viability_impc <- arrow::read_parquet('./data/raw/api-ftp/mouse_viability_impc.parquet')
human_mouse_genes <- protein.coding.genes %>%
  separate_rows(mgd_id, sep = "\\|") %>%
  dplyr::select(symbol, mgd_id) %>%
  dplyr::rename(gene_symbol = symbol, mgi_id = mgd_id)

mouse.viability.impc2 <- mouse_viability_impc %>%
  dplyr::filter(Comment == "") %>% # Remove conflicting evidence
  dplyr::select('Gene Accession Id', 'Viability Phenotype HOMs/HEMIs') %>%
  dplyr::rename(mgi_id = 'Gene Accession Id',
                impc_viability = 'Viability Phenotype HOMs/HEMIs') %>%
  distinct() %>%
  left_join(human_mouse_genes) %>%  # Join gene symbols to mouse viability
  dplyr::filter(!is.na(gene_symbol)) %>%
  dplyr::select(-mgi_id) %>%
  distinct()

conflicts <- mouse.viability.impc2 %>% # remove conflicts that arise from one2many mappings
  distinct() %>%
  count(gene_symbol) %>%
  filter(n == 1)

mouse.viability.impc3 <- mouse.viability.impc2 %>%
  filter(gene_symbol %in% conflicts$gene_symbol) %>%
  distinct()

# impc phenotypes 
mouse_phenotypes_impc <- arrow::read_parquet('./data/raw/api-ftp/mouse_phenotypes_impc.parquet')

mouse.phenotypes.impc2 <- mouse_phenotypes_impc %>%
  dplyr::select(marker_accession_id, zygosity, mp_term_name, mp_term_id) %>%
  dplyr::rename(mgi_id = marker_accession_id,
                impc_zygosity = zygosity,
                impc_phenotypes = mp_term_name,
                mp_id = mp_term_id) %>%
  dplyr::distinct() %>%
  left_join(human_mouse_genes) %>%
  dplyr::filter(!is.na(gene_symbol))

# mgi
lethal_terms <- mp_lethal_terms$x
mouse_protein_coding_genes_mgi <- arrow::read_parquet('./data/raw/api-ftp/mouse_protein_coding_genes_mgi.parquet')
mouse_proteincoding_genes <- mouse_protein_coding_genes_mgi$`MGI Accession ID`
mouse_viability_mgi <- arrow::read_parquet('./data/raw/api-ftp/mouse_viability_mgi.parquet')
mouse.viability.mgi2 <- mouse_viability_mgi %>%
  dplyr::select(7,5) %>%
  distinct() %>%
  dplyr::rename(mgi_id = V7, mp_term = V5) %>%
  filter(mgi_id %in% mouse_proteincoding_genes) %>%
  mutate(mp_term_lethal = ifelse(mp_term %in% lethal_terms, "y","n")) %>%
  dplyr::select(mgi_id, mp_term_lethal) %>%
  distinct() %>%
  arrange(mgi_id, mp_term_lethal) %>%
  group_by(mgi_id) %>%
  summarise(mgi_lethal_term = paste0(unique(mp_term_lethal), collapse = "|")) %>%
  mutate(viability_mgi = ifelse(mgi_lethal_term == "n","viable","lethal")) %>%
  dplyr::select(mgi_id, viability_mgi) %>%
  distinct() 
mouse.viability.mgi3 <- mouse.viability.mgi2 %>% 
  left_join(human_mouse_genes) %>%
  filter(!is.na(gene_symbol))

# Functional Annotation ----

# go
gene_symbol_entrez_id_mapping <- protein.coding.genes %>%
  dplyr::select(symbol, entrez_id) %>%
  dplyr::rename(gene_symbol = symbol)
gene_symbol_entrez_id_mapping$entrez_id <- as.character(gene_symbol_entrez_id_mapping$entrez_id)

go_anot_term_raw <- arrow::read_parquet('./data/raw/api-ftp/go_anot_term_raw.parquet') 

go_bp <- go_anot_term_raw %>% 
  dplyr::rename(entrez_id = gene_id) %>% 
  left_join(gene_symbol_entrez_id_mapping) %>% 
  filter(Ontology == "BP") %>% 
  pivot_wider(names_from = Ontology,
              values_from = c("go_id","go_term")) %>% 
  unnest(c(go_id_BP, go_term_BP)) %>% 
  dplyr::select(-entrez_id)

go_mf <- go_anot_term_raw %>% 
  dplyr::rename(entrez_id = gene_id) %>% 
  left_join(gene_symbol_entrez_id_mapping) %>% 
  filter(Ontology == "MF") %>% 
  pivot_wider(names_from = Ontology,
              values_from = c("go_id","go_term")) %>% 
  unnest(c(go_id_MF, go_term_MF)) %>% 
  dplyr::select(-entrez_id)

go_cc <- go_anot_term_raw %>% 
  dplyr::rename(entrez_id = gene_id) %>% 
  left_join(gene_symbol_entrez_id_mapping) %>% 
  filter(Ontology == "CC") %>% 
  pivot_wider(names_from = Ontology,
              values_from = c("go_id","go_term")) %>% 
  unnest(c(go_id_CC, go_term_CC)) %>% 
  dplyr::select(-entrez_id)

# reactome
reactome_annotations <- arrow::read_parquet('./data/raw/api-ftp/reactome_annotations.parquet')
reactome_annotations2 <- reactome_annotations %>%
  filter(grepl("Homo sapiens:", path_name)) %>%
  mutate(path_name = sub("Homo sapiens: ", "", path_name)) %>%
  dplyr::rename(entrez_id = 'gene_id')  %>%
  left_join(gene_symbol_entrez_id_mapping) %>%
  dplyr::select(-entrez_id) %>%
  filter(!is.na(gene_symbol))

# Human Disease & Phenotypes ----

# lethal phenotypes
lethal.phenotypes <- arrow::read_parquet('./data/raw/api-ftp/lethal.phenotypes.parquet')
lethal.phenotypes2 <- lethal.phenotypes %>%
  dplyr::select(
    omim_id,
    # omim_phenotype,
    disease_gene_lethal,
    earliest_lethality_category
  ) %>%
  mutate(
    earliest_lethality_category = case_when(
      is.na(earliest_lethality_category) ~ NA_character_,
      earliest_lethality_category == "L1" ~ "Prenatal death; L1",
      earliest_lethality_category == "L2" ~ "Neonatal death; L2",
      earliest_lethality_category == "L3" ~ "Death in infancy; L3",
      earliest_lethality_category == "L4" ~ "Death in childhood; L4",
      earliest_lethality_category == "L5" ~ "Death in adolescence; L5",
      earliest_lethality_category == "L6" ~ "Death in adulthood; L6",
      earliest_lethality_category == "LU" ~ "Not determined; LU",
      earliest_lethality_category == "NL" ~ "Non lethal; NL",
      TRUE ~ earliest_lethality_category
    )) %>%
  mutate_all(., ~ ifelse(. == "-", NA, .)) %>%
  mutate(phenotype_id = as.character(omim_id)) %>%
  dplyr::select(-omim_id)

# panel app
jfiles_unlist <- read_rds('./data/raw/api-ftp/jfiles_unlist.rds')
all_jfiles <- do.call(rbind, jfiles_unlist)

all_jfiles_df <- all_jfiles  %>%
  filter(name %in% c("results.gene_data.gene_symbol", "results.entity_type", 
                     "results.confidence_level", "results.mode_of_inheritance",
                     "results.panel.id", "results.panel.name", 
                     "results.panel.disease_group",
                     "results.panel.disease_sub_group")) %>%
  replace(is.na(.),"-") %>%
  mutate(value = as.character(value))


gene_symbol = all_jfiles_df %>%
  filter(name == "results.gene_data.gene_symbol")

entity_type = all_jfiles_df %>%
  filter(name == "results.entity_type")

confidence_level = all_jfiles_df %>%
  filter(name == "results.confidence_level")

mode_of_inheritance = all_jfiles_df %>%
  filter(name == "results.mode_of_inheritance")

panel_id = all_jfiles_df %>%
  filter(name == "results.panel.id")

panel_name = all_jfiles_df %>%
  filter(name == "results.panel.name")

panel_disease_group = all_jfiles_df %>%
  filter(name == "results.panel.disease_group")

panel_disease_subgroup = all_jfiles_df %>%
  filter(name == "results.panel.disease_sub_group")

panelapp_df <- data.frame(gene_symbol = gene_symbol$value,
                          confidence_level = confidence_level$value,
                          mode_of_inheritance = mode_of_inheritance$value,
                          panel_id = panel_id$value,
                          panel_name = panel_name$value,
                          panel_disease_group = panel_disease_group$value,
                          panel_disease_subgroup = panel_disease_subgroup$value)
panelapp_df2 <- panelapp_df %>% 
  mutate(confidence_level = case_when(
    confidence_level == 3 ~ "green",
    confidence_level == 2 ~ "amber",
    confidence_level == 1 ~ "red",
    TRUE ~ NA_character_  # optional: handles missing or unexpected values
  )) %>% 
  dplyr::select(gene_symbol, panel_disease_group, confidence_level, mode_of_inheritance) %>% 
  filter(panel_disease_group != "")

# Comparative/Evolutionary Genomics ----

# human-mouse orthologs
mouse_orthologs <- arrow::read_parquet('./data/raw/api-ftp/mouse_orthologs.parquet')
orthologs_names <- names(mouse_orthologs)[-1]
orthologs_names <- c(orthologs_names, 'blank')
names(mouse_orthologs) <- orthologs_names
mouse_orthologs <- mouse_orthologs %>% 
  dplyr::select(-blank)
orthologs2 <- mouse_orthologs %>%
  dplyr::select(`Human Gene Symbol`, `Human Category For Threshold`, `Mgi Gene Acc Id`, `Mouse Gene Symbol`) %>%
  dplyr::rename(gene_symbol = `Human Gene Symbol`, ortholog_mapping = `Human Category For Threshold`, mgi_id = `Mgi Gene Acc Id`, mouse_symbol = `Mouse Gene Symbol`)

# human paralogs
human_paralogs <- arrow::read_parquet('./data/raw/api-ftp/human_paralogs.parquet')
gene_symbol_ensg_id_mappings <- protein.coding.genes %>%
  dplyr::select(symbol, ensembl_gene_id) 
human_paralogs2 <- human_paralogs %>% 
  left_join(gene_symbol_ensg_id_mappings) %>% 
  dplyr::select(symbol, hsapiens_paralog_associated_gene_name, hsapiens_paralog_subtype) %>% 
  filter(hsapiens_paralog_associated_gene_name != "") %>% 
  distinct() %>% 
  as.data.frame()
