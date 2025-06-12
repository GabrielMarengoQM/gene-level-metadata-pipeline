# Import manuscript-files
library(readr)
library(openxlsx)

# AlphaMissense pathogencity scores
alphamissense_pathogenicity <- readr::read_tsv("./data/raw/manuscript-files/AlphaMissense_gene_hg38.tsv.gz", skip = 3)

# Cell fitness in different culture conditions - Laminin & MEF (Mair 2019)
Mair_2019_mef <- openxlsx::read.xlsx("./data/raw/manuscript-files/1-s2.0-S2211124719302128-mmc2.xlsx")
Mair_2019_laminin <- openxlsx::read.xlsx("./data/raw/manuscript-files/1-s2.0-S2211124719302128-mmc4.xlsx")

# SCoNeS (Rapaport 2021)
scones <- openxlsx::read.xlsx("./data/raw/manuscript-files/pnas.2001248118.sd01.xlsx", startRow = 2)

# DOMINO
domino <- readr::read_tsv("./data/raw/manuscript-files/score_all_final_19.02.19.txt")

# Pluripotency and fitness scores (Rosen 2024)
rosen24_NE_pluripotency_score <- openxlsx::read.xlsx("./data/raw/manuscript-files/41467_2024_53284_MOESM4_ESM.xlsx", sheet=1) 
rosen24_DE_pluripotency_score <- openxlsx::read.xlsx("./data/raw/manuscript-files/41467_2024_53284_MOESM4_ESM.xlsx", sheet=2) 
rosen24_E8_self_renewel_score <- openxlsx::read.xlsx("./data/raw/manuscript-files/41467_2024_53284_MOESM4_ESM.xlsx", sheet=3) 
rosen24_E6_self_renewel_score <- openxlsx::read.xlsx("./data/raw/manuscript-files/41467_2024_53284_MOESM4_ESM.xlsx", sheet=4) 
 

# GISMO, GISMO missense
gismo <- openxlsx::read.xlsx("./data/raw/manuscript-files/media-2.xlsx", sheet = 'Supplementary Table 2')
gismo_mis <- openxlsx::read.xlsx("./data/raw/manuscript-files/media-2.xlsx", sheet = 'Supplementary Table 3')

# S het
shet <- readr::read_tsv("./data/raw/manuscript-files/media-1.tsv")

