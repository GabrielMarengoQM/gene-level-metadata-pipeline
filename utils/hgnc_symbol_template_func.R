# Gene symbols for mapping to data tables

hgnc_symbol_template_func <- function() {
  hgnc_symbol_template <- protein.coding.genes %>% 
    dplyr::select(symbol) %>% 
    dplyr::rename(hgnc_gene_symbol = symbol) %>% 
    filter(!is.na(hgnc_gene_symbol)) %>% 
    distinct()
  return(hgnc_symbol_template)
}

