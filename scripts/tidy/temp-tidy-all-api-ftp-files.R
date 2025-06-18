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

# Prev names
prev_names1 <- protein.coding.genes %>% 
  dplyr::select(symbol, prev_symbol) %>% 
  separate_rows(prev_symbol, sep = "\\|") %>% 
  dplyr::rename(hgnc_gene_symbol = symbol) %>% 
  filter(prev_symbol != "")

prev_names <- hgnc_symbol_template_func() %>% 
  left_join(prev_names1) %>% 
  distinct()

# Alias names
alias_names1 <- protein.coding.genes %>% 
  dplyr::select(symbol, alias_symbol) %>% 
  separate_rows(alias_symbol, sep = "\\|") %>% 
  dplyr::rename(hgnc_gene_symbol = symbol) %>% 
  filter(alias_symbol != "")

alias_names <- hgnc_symbol_template_func() %>% 
  left_join(alias_names1) %>% 
  distinct()

# human mouse mappings
mouse_ids <- protein.coding.genes %>% 
  dplyr::select(symbol, mgd_id) %>% 
  tidyr::separate_rows(mgd_id, sep = '\\|') %>% 
  dplyr::rename(hgnc_gene_symbol = symbol) %>% 
  filter(mgd_id != "")

hgnc_symbol_mgi_id_mappings <- hgnc_symbol_template_func() %>% 
  left_join(mouse_ids) %>% 
  distinct()


# gene uniprot protein mappings
uniprot_ids <- protein.coding.genes %>% 
  dplyr::select(symbol, uniprot_ids) %>%
  tidyr::separate_rows(uniprot_ids, sep = '\\|') %>% 
  dplyr::rename(hgnc_gene_symbol = symbol) %>% 
  filter(uniprot_ids != "")

hgnc_symbol_uniprot_id_mappings <- hgnc_symbol_template_func() %>% 
  left_join(uniprot_ids) %>% 
  distinct()

# gene name
gene_names1 <- protein.coding.genes %>% 
  dplyr::select(symbol, name) %>% 
  dplyr::rename(hgnc_gene_symbol = symbol) %>% 
  filter(name != "")

gene_names <- hgnc_symbol_template_func() %>% 
  left_join(gene_names1) %>% 
  distinct()

# gene group
gene_groups1 <- protein.coding.genes %>% 
  dplyr::select(symbol, gene_group) %>% 
  tidyr::separate_rows(gene_group, sep = "\\|") %>% 
  filter(!is.na(gene_group)) %>% 
  filter(gene_group != "") %>% 
  dplyr::rename(hgnc_gene_symbol = symbol) %>% 
  filter(gene_group != "")

gene_groups <- hgnc_symbol_template_func() %>% 
  left_join(gene_groups1) %>% 
  distinct()

# chromosome, positions, length
chromosomes <- c(as.character(1:22), "X", "Y")
latest_ensg_id <- protein.coding.genes %>% pull(ensembl_gene_id)
gene_position_length1 <- arrow::read_parquet('./data/raw/api-ftp/gene_position_gc_content.parquet') %>% 
  dplyr::mutate(gene_length = end_position - start_position + 1) %>% 
  dplyr::select(-percentage_gene_gc_content) %>% 
  dplyr::rename(hgnc_gene_symbol = hgnc_symbol) %>% 
  filter(ensembl_gene_id %in% latest_ensg_id) %>% 
  filter(chromosome_name != "") %>% 
  filter(chromosome_name %in% chromosomes) %>% 
  filter(start_position != "") %>% 
  filter(end_position != "") %>% 
  filter(gene_length != "") %>% 
  distinct()

gene_position_length <- hgnc_symbol_template_func() %>% 
  left_join(gene_position_length1) %>% 
  distinct()

# GC content
gene_gc_content1 <- arrow::read_parquet('./data/raw/api-ftp/gene_position_gc_content.parquet') %>% 
  dplyr::select(hgnc_symbol, percentage_gene_gc_content, chromosome_name, ensembl_gene_id) %>% 
  dplyr::rename(hgnc_gene_symbol = hgnc_symbol) %>% 
  filter(ensembl_gene_id %in% latest_ensg_id) %>% 
  filter(percentage_gene_gc_content != "") %>% 
  filter(chromosome_name %in% chromosomes) %>% 
  dplyr::select(-chromosome_name, -ensembl_gene_id)

gene_gc_content <- hgnc_symbol_template_func() %>% 
  left_join(gene_gc_content1) %>% 
  distinct()

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

# OMIM genemap 
genemap <- arrow::read_parquet('./data/raw/api-ftp/genemap.parquet')

# Split phenotypes across new rows for one2many gene-pheno mappings
genemap_2 <- genemap %>%
  separate_rows(Phenotypes, sep = "; ")

# Grep mois
moi_keywords <- c("Autosomal recessive", "Autosomal dominant", "Digenic recessive", "Digenic dominant",
                  "Somatic mutation", "Multifactorial", "Isolated cases", "Mitochondrial",
                  "Pseudoautosomal dominant", "Pseudoautosomal recessive", "X-linked recessive",
                  "X-linked dominant", "X-linked", "Y-linked")

genemap_3 <- genemap_2 %>%
  mutate(moi = sapply(str_extract_all(Phenotypes, paste(moi_keywords, collapse = "|")), 
                      paste, collapse = "; ")) %>%
  separate_rows(moi, sep = "; ")

# Grep symbols
genemap_4 <- genemap_3 %>%
  mutate(symbol_key = case_when(
    grepl("^\\{", Phenotypes) ~ "Mutations that contribute to susceptibility to multifactorial disorders",
    grepl("^\\[", Phenotypes) ~ "Nondiseases",
    grepl("^\\?", Phenotypes) ~ "The relationship between the phenotype and gene is provisional",
    TRUE ~ NA_character_  # Handle other cases if needed
  ))

# Grep numbers
genemap_5 <- genemap_4 %>%
  mutate(number_key = case_when(
    grepl("\\(1)", Phenotypes) ~ "The disorder is placed on the map based on its association with a gene, but the underlying defect is not known",
    grepl("\\(2)", Phenotypes) ~ "The disorder has been placed on the map by linkage or other statistical method; no mutation has been found",
    grepl("\\(3)", Phenotypes) ~ "The molecular basis for the disorder is known; a mutation has been found in the gene",
    grepl("\\(4)", Phenotypes) ~ "A contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated causing the phenotype",
    TRUE ~ NA_character_  # Handle other cases if needed
  ))

# Get phenotype id
genemap_6 <- genemap_5 %>%
  mutate(phenotype_id = str_extract(Phenotypes, "\\d{6}"))

# Get Phenotype
genemap_7 <- genemap_6 %>%
  mutate(phenotypes = str_extract(Phenotypes, "^(.*?)(?=\\d{6})")) %>%
  mutate(phenotypes = str_replace_all(phenotypes, "[\\{\\[\\?\\]\\}]", "")) %>%
  mutate(phenotypes = str_replace(phenotypes, ", $", "")) 

genemap8 <- genemap_7 %>%
  dplyr::select(`Approved.Gene.Symbol`, moi, symbol_key, number_key, phenotype_id, phenotypes) %>%
  filter(!is.na(`Approved.Gene.Symbol`))

genemap9 <- genemap8 %>% 
  rename(hgnc_gene_symbol = `Approved.Gene.Symbol`) %>% 
  filter(hgnc_gene_symbol != "") %>% 
  filter(!is.na(phenotypes))



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
