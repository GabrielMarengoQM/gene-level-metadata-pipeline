library(biomaRt)
library(dplyr)

# Constraint metrics ----

# Alphamissense ~13k genes
alpham <- alphamissense_pathogenicity %>%
  dplyr::mutate(transcript_id = sub("\\.[0-9]+$", "", transcript_id)) %>% 
  dplyr::rename(ensembl_transcript_id = transcript_id) 
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl")
mane_select_data <- getBM(
  attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_mane_select", "transcript_is_canonical"),
  filters = "hgnc_symbol",
  values = protein.coding.genes$symbol,
  mart = ensembl
)
alpham2 <- alpham %>%
  left_join(mane_select_data) %>%
  mutate(
    flag = case_when(
      transcript_mane_select != "" & is.na(transcript_is_canonical) ~ "mane only",
      transcript_is_canonical == 1 & transcript_mane_select == "" ~ "canonical only",
      transcript_mane_select != "" & transcript_is_canonical == 1 ~ "mane and canonical",
      transcript_mane_select == "" & is.na(transcript_is_canonical) ~ "neither",
      TRUE ~ NA_character_  
    )
  ) %>% 
  filter(!is.na(hgnc_symbol))
# ~27 genes canoncical only, ~4k neither
# 53 dups
# Keep mane -> if dups then select mane or canoncical if possible -> average any remaining dups

# SCoNeS
scones2 <- scones[, c(1, 17, 19)] %>%
  dplyr::rename(gene_symbol = Gene)
scones3 <- subset(scones2, gene_symbol %in% protein.coding.genes$symbol)
scones4 <- scones3 %>%
  dplyr::rename(scones = SCoNeS) %>% 
  dplyr::select(-DOMINO)

# DOMINO - three dups, removed
domino2 <- domino %>% 
  dplyr::select(`#HGNC ID`, Score) %>%
  dplyr::rename(hgnc_symbol = `#HGNC ID`, domino = Score)
domino_dups <- domino2 %>% count(hgnc_symbol) %>% filter(n > 1)
domino3 <- domino2 %>% filter(!hgnc_symbol %in% domino_dups$hgnc_symbol)

# GISMO & GISMO-mis
gene_symbol_ensg_id_mappings <- protein.coding.genes %>%
  dplyr::select(symbol, ensembl_gene_id) %>%
  dplyr::rename(gene_symbol = symbol) %>%
  dplyr::rename(gene = ensembl_gene_id)
gismo2 <- gismo %>%
  left_join(gene_symbol_ensg_id_mappings) %>%
  dplyr::select(median, gene_symbol, decile) %>%
  dplyr::select(gene_symbol, everything()) %>%
  filter(!is.na(gene_symbol)) %>%
  dplyr::rename(gismo_median = median, gismo_decile = decile)
gismo_mis2 <- gismo_mis %>%
  dplyr::select(mean.comb, genename, dec, gene) %>%
  dplyr::rename(gene_symbol = genename) %>%
  dplyr::select(gene_symbol, everything()) %>%
  dplyr::rename(gismo_mis_mean_comb = mean.comb, gismo_mis_decile = dec, transcript = gene)
gismo_gismo_mis2 <- gismo2 %>%
  full_join(gismo_mis2) %>%
  filter(gene_symbol %in% protein.coding.genes$symbol) %>%
  dplyr::select(-transcript) %>%
  mutate(
    gismo_median = round(gismo_median, 3),
    gismo_mis_mean_comb = round(gismo_mis_mean_comb, 3)
  ) %>%
  distinct()
gismo_gismo_mis2_dups <- gismo_gismo_mis2 %>% count(gene_symbol) %>% filter(n > 1)
gismo_gismo_mis3 <- gismo_gismo_mis2 %>% 
  filter(!gene_symbol %in% gismo_gismo_mis2_dups$gene_symbol) %>% 
  distinct()

# S het
shet_post <- shet[, c(1, 2, 7, 8, 9)]
names(shet_post) <- c("ens_gene_id", "hgnc_id", "shet_post_mean", "shet_post_lower", "shet_post_upper")
shet_post <- shet_post %>%
  mutate(
    shet_post_mean = round(shet_post_mean, 3),
    shet_post_lower = round(shet_post_lower, 3),
    shet_post_upper = round(shet_post_upper, 3)
  ) %>%
  dplyr::select(hgnc_id, shet_post_mean) %>% 
  distinct() 

# Cellular assays ----

# Cell fitness in different culture conditions - Laminin & MEF (Mair 2019)
mef <- Mair_2019_mef[,c(1, 2, 7)]
names(mef) <- c("gene_symbol", "bf_mef", "fdr_mef")
mef <- mef %>%
  mutate(fdr_mef = round(fdr_mef, 3)) 

laminin <- Mair_2019_laminin[,c(1, 2, 7)]
names(laminin) <- c("gene_symbol", "bf_lam", "fdr_lam")
laminin <- laminin %>%
  mutate(fdr_lam = round(fdr_lam, 3)) 

# rosen 2024
rosen24_NE_pluripotency_score2 <- rosen24_NE_pluripotency_score %>%
  dplyr::select(X1, 11) %>% 
  setNames(as.character(.[1,])) %>%  
  dplyr::slice(-1) %>% 
  dplyr::mutate(across(2, ~ round(as.numeric(.), 2))) %>% 
  dplyr::rename()
rosen24_DE_pluripotency_score2 <- rosen24_DE_pluripotency_score %>%
  dplyr::select(X1, 11) %>% 
  setNames(as.character(.[1,])) %>%  
  slice(-1) %>% 
  mutate(across(2, ~ round(as.numeric(.), 2)))
rosen24_E8_self_renewel_score2 <- rosen24_E8_self_renewel_score %>%
  dplyr::select(X1, 11) %>% 
  setNames(as.character(.[1,])) %>%  
  slice(-1) %>% 
  mutate(across(2, ~ round(as.numeric(.), 2))) %>% 
  dplyr::rename(E8_self_renewal_score=2)
rosen24_E6_self_renewel_score2 <- rosen24_E6_self_renewel_score %>%
  dplyr::select(X1, 11) %>% 
  setNames(as.character(.[1,])) %>%  
  slice(-1) %>% 
  mutate(across(2, ~ round(as.numeric(.), 2))) %>% 
  dplyr::rename(E6_self_renewal_score=2)


