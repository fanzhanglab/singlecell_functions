
# GLM regression: identify genes that are associated with pseudotime from single-cell trajectory analysis
# logit(position) ~ 1 + 
  
markers <- c("THY1", "NOTCH3", "PDGFRA", "PDPN", "CD34", "CD55", "PRG4", "PTPRC")
with(tissue, {
    idx_use <- which(!is.na(meta_data$time))
    glm_res <- Reduce(rbind, lapply(markers, function(marker) {
        message(marker)
        ldata <- data.frame(value = exprs_norm[marker, idx_use], meta_data[idx_use, ])
        suppressWarnings({
            fit_full <- glm(time ~ 1 + value + donor + log(nGene), family = "binomial", data = ldata)
        })
        data.table(symbol = marker, subset(broom::tidy(fit_full), term == "value"))
    })) %>%
        dplyr::mutate(fdr = p.adjust(p.value, "BH")) %>%
        dplyr::arrange(fdr, p.value)
})



# prioritize trajectory associated genes
idx_tissue <- which(meta_data$source == "primary")
tissue <- list()
tissue$meta_data <- meta_data[idx_tissue, ]
tissue$exprs_norm <- exprs_norm[, idx_tissue]

with(tissue, {
  idx_use <- which(meta_data$cell_type %in% c('lining', 'sublining'))
  design <- model.matrix(~time + donor, meta_data[idx_use, ])
  limma_res <<- exprs_norm[, idx_use] %>% 
    limma::lmFit(design) %>% 
    limma::eBayes() %>% 
    limma::topTable(coef = 2, number = 1e4)     
})

limma_res %>% head
limma_res %>% subset(adj.P.Val < 1e-5) %>% nrow
nrow(limma_res)