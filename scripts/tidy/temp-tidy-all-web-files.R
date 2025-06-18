library(dplyr)
library(stringr)
library(tidyr)

source('./utils/hgnc_symbol_template_func.R')

# Cellular assays ----

# DepMap
# cell0
# Transform and turn into matrix of 1 and 0 where 1 is essential and 0 is not (-0.5 gene effect threshold)
cell0 <- gene_effect
colnames(cell0) <-  stringr::str_split(names(cell0), "\\ ", simplify=T)[,1]
colnames(cell0)[colnames(cell0) == "...1"] <- "model_id"
cell1  = cell0[,-1]
cell2 <- t(cell1)
cell3 <- cell2
colnames(cell3) = cell0$V1
cell4 = cell3
cell4[cell4 >= -0.5] <- 0
cell4[cell4 != 0] <- 1
cell4 <- data.frame(cell4)
# cell 5 = with n_essential, percentage_essential columns
cell5 <- cell4 %>%
  mutate(n_essential = rowSums(.)) %>%
  mutate(percentage_essential = (n_essential/ncol(cell4))*100)
# cell 7 = mean scores
cell7 <- data.frame(cell3) %>%
  mutate(mean_score_all = rowMeans(.)) %>%
  mutate(gene_symbol = row.names(cell3)) %>%
  dplyr::select(gene_symbol, mean_score_all) %>%
  `rownames<-`( NULL )

depmap.data1 <- cell5 %>%
  mutate(gene_symbol = row.names(.)) %>%
  dplyr::select(gene_symbol, percentage_essential) %>%
  mutate(gene_symbol = gsub("\\.\\..*?\\.", "", gene_symbol))

depmap.data2 <- cell7 %>%
  dplyr::select(gene_symbol, mean_score_all) %>%
  mutate(gene_symbol = gsub("\\.\\..*?\\.", "", gene_symbol))

depmap.data3 <- depmap.data1 %>%
  full_join(depmap.data2, by = "gene_symbol") %>%
  mutate(percentage_essential = round(percentage_essential, 3),
         mean_score_all = round(mean_score_all, 3)) %>% 
  dplyr::rename(hgnc_gene_symbol = gene_symbol)

DepMap <- hgnc_symbol_template_func() %>% 
  left_join(depmap.data3) %>% 
  distinct()

# Transcriptomics ----

# gtex
gene_symbol_ensg_id_mappings <- protein.coding.genes %>%
  dplyr::select(symbol, ensembl_gene_id) 
gtex_expression2 <- gtex_expression %>%
  dplyr::filter(!str_detect(Name, "_PAR_Y")) %>% # Remove PAR Y genes
  dplyr::mutate(Name = str_replace(Name, "\\..*", "")) %>%  # Remove full stop and number
  dplyr::filter(Name %in% protein.coding.genes$ensembl_gene_id) %>% # Only keep current gene symbols and ensembl gene ids in all pcg
  dplyr::rename(ensembl_gene_id = Name) %>% 
  dplyr::select(-Description) %>% 
  left_join(gene_symbol_ensg_id_mappings) %>% 
  dplyr::rename(hgnc_gene_symbol = symbol) %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))

GTEx_expression <- hgnc_symbol_template_func() %>% 
  left_join(gtex_expression2) %>% 
  distinct()

# Protein atlas; bulk expression / subcellular location / protein class / 
hpa_bulk_expression2 <- hpa_bulk_expression %>% 
  filter(Reliability == "Approved") %>% 
  dplyr::select(-Gene, -Reliability) %>% 
  dplyr::rename(hgnc_gene_symbol = `Gene name`)

HPA_bulk_expression <- hgnc_symbol_template_func() %>% 
  left_join(hpa_bulk_expression2) %>% 
  distinct()

# lymphoblastoid_time_series_expression across 7 organs (TPM)
lymphoblastoid_time_series_expression2 <- lymphoblastoid_time_series_expression %>% 
  dplyr::select(-`Gene ID`) %>% 
  dplyr::rename(hgnc_gene_symbol = `Gene Name`)

Lymphoblastoid_time_series_expression <- hgnc_symbol_template_func() %>% 
  left_join(lymphoblastoid_time_series_expression2) %>% 
  distinct()

# Genomic features ----

# subcellular location
HPA_subcellular_location_0 <- hpa_atlas %>% 
  dplyr::select(Gene, `Subcellular location`) %>% 
  tidyr::separate_rows(`Subcellular location`, sep = "\\," ) %>% 
  filter(!is.na(`Subcellular location`)) %>% 
  dplyr::rename(hgnc_gene_symbol = Gene)

HPA_subcellular_location <- hgnc_symbol_template_func() %>% 
  left_join(HPA_subcellular_location_0) %>% 
  distinct()

# Proteins ----

# protein class
HPA_protein_classes_0 <- hpa_atlas %>% 
  dplyr::select(Gene, `Protein class`) %>% 
  tidyr::separate_rows(`Protein class`, sep = "\\," ) %>% 
  mutate(`Protein class` = trimws(`Protein class`)) %>% 
  dplyr::rename(hgnc_gene_symbol = Gene)

HPA_protein_classes <- hgnc_symbol_template_func() %>% 
  left_join(HPA_protein_classes_0) %>% 
  distinct()

# ogee/string connectivity scores
ppi_connectivity2 <- ppi_connectivity %>% 
  filter(gene %in% protein.coding.genes$symbol) %>% 
  dplyr::select(gene, score, connectivity, percentile) %>% 
  dplyr::rename(hgnc_gene_symbol = gene)

OGEE_ppi_connectivity <- hgnc_symbol_template_func() %>% 
  left_join(ppi_connectivity2) %>% 
  distinct()

# Disease/human phenotypes ----

# dbNSFP Clingen HI scores
dbNSFP_hi <- dbNSFP %>% 
  dplyr::select(Gene_name, ClinGen_Haploinsufficiency_Score) %>% 
  filter(ClinGen_Haploinsufficiency_Score != ".") %>% 
  dplyr::rename(hgnc_gene_symbol = Gene_name)

dbNSFP_haplo_insufficiency <- hgnc_symbol_template_func() %>% 
  left_join(dbNSFP_hi) %>% 
  distinct()

# Constraint metrics ----

# gnomAD
gnomad2 <- gnomad[, c("gene", "transcript", "mane_select", "lof.oe_ci.upper", "mis.oe_ci.upper", "constraint_flags")]
names(gnomad2) <- c("gene", "gnomad_transcript", "gnomad_mane_select", 
                      "gnomad_lof_upper_90_ci", "gnomad_mis_upper_90_ci", 'gnomad_constraint_flags')
mane_transcripts <- mane_select_data %>% filter(transcript_mane_select != "") %>% pull(ensembl_transcript_id)
canonical_transcripts <- mane_select_data %>% filter(transcript_is_canonical == 1) %>% pull(ensembl_transcript_id)
gnomad3 <- gnomad2 %>% 
  dplyr::select(gene, gnomad_transcript, gnomad_mane_select) %>% 
  distinct() %>% 
  filter(gnomad_transcript %in% mane_transcripts | gnomad_transcript %in% canonical_transcripts) 
gnomad4 <- gnomad2 %>% 
  filter(gnomad_transcript %in% gnomad3$gnomad_transcript) %>% 
  dplyr::select(gene, gnomad_lof_upper_90_ci, gnomad_mane_select) %>% 
  dplyr::rename(hgnc_gene_symbol = gene) %>% 
  filter(!is.na(hgnc_gene_symbol)) %>% 
  distinct()
dups <- gnomad4 %>% 
  group_by(hgnc_gene_symbol) %>%                        
  filter(n() > 1 & gnomad_mane_select == "true") %>% 
  ungroup()
gnomad5 <- gnomad4 %>% 
  group_by(hgnc_gene_symbol) %>%                        
  filter(n() == 1) %>% 
  ungroup() %>% 
  bind_rows(dups) %>% 
  dplyr::select(-gnomad_mane_select) %>% 
  dplyr::rename(LOEUF = gnomad_lof_upper_90_ci)

gnomAD <- hgnc_symbol_template_func() %>% 
  left_join(gnomad5) %>% 
  distinct()

# Mouse assays ----
# impc window of lethality
wol2 <- wol %>% 
  dplyr::select(hs_gene_symbol, wol) %>% 
  dplyr::rename(hgnc_gene_symbol = hs_gene_symbol) %>% 
  distinct()

WoL <- hgnc_symbol_template_func() %>% 
  left_join(wol2) %>% 
  distinct()

# GWAS
dbNSFP_GWAS_catalog1 <- dbNSFP %>%
  dplyr::select(Gene_name, `Trait_association(GWAS)`) %>% 
  # 1) Remove [ ... ] including anything inside
  mutate(no_brackets = str_remove_all(`Trait_association(GWAS)`, "\\[[^]]*\\]")) %>%
  # 2) Separate rows on semicolons
  separate_rows(no_brackets, sep = ";") %>%
  # 3) Trim whitespace and remove any empty strings
  mutate(no_brackets = str_trim(no_brackets)) %>%
  filter(no_brackets != "") %>% 
  dplyr::mutate(across(everything(), ~na_if(.x, "."))) %>% 
  dplyr::rename(GWAS_trait = no_brackets) %>% 
  dplyr::select(-`Trait_association(GWAS)`) %>% 
  distinct() %>% 
  filter(!is.na(GWAS_trait)) %>% 
  dplyr::rename(hgnc_gene_symbol = Gene_name)

dbNSFP_GWAS_catalog <- hgnc_symbol_template_func() %>% 
  left_join(dbNSFP_GWAS_catalog1) %>% 
  distinct()

