library(data.table)
library(magrittr)
library(ggplot2)
library(pheatmap)

## now by specific protein level

dt <- data.table::fread("../Data/p53_scores.tsv") %>% .[effect == "Missense"]

dt <- dt[!protein_change %like% "delins"]

dt[, locus := protein_change %>% stringr::str_extract("(?<=p\\.)[A-Z]+[0-9]+")]

wl <- dt[, .N, by = "locus"] %>% .[N > 3, locus %>% unique]

dt <- dt[locus %chin% wl]

m <- dt[, .N, by = c("protein_change", "p53_intensity", "ID")]

m %<>% dcast(p53_intensity ~ protein_change + ID, value.var = "N")

rn <- c("0", "1+", "2+", "3+")

m <- m[, .SD, .SDcols = 2:ncol(m)] %>% as.matrix()

rownames(m) <- rn

m[is.na(m)] <- 0

reorder <- order(colnames(m) %>% stringr::str_extract("[0-9]+[A-z]+"))

m <- m[, reorder]

v2 <- colnames(m)

v <- v2 %>% stringr::str_replace("^p\\.","") %>% stringr::str_replace("_[0-9]+$", "")

anno <- as.data.frame(v)

rownames(anno) <- colnames(m)

colnames(anno) <- "Variant"

colnames(m) %<>% stringr::str_replace("^p\\.", "") %>%
  stringr::str_replace("_.*$", "")

breaks <- seq(from = 1, to = ncol(m), by = 31)

breaks %<>% c(ncol(m) + 1)

pl <- lapply(breaks[1:(length(breaks) - 1)] %>% seq_along, function(n){
  
  m2 <- m[, breaks[n]:(breaks[n + 1] - 1)]
  
  hm <- pheatmap::pheatmap(m2,
                           color = colorRampPalette(c("white", "dodgerblue"))(4),
                           cellwidth = 16, cellheight = 14,
                           cluster_cols = FALSE, cluster_rows = FALSE,
                           display_numbers = TRUE,
                           number_color = "black", number_face = "bold",
                           number_format = "%.0f",
                           main = "Staining intensity for recurrent missense loci",
                           legend = FALSE,
                           angle_col = 45,
                           fontsize_number = 8,
                           fontsize_col = 12,
                           fontsize_row = 12
  )
  
  return(hm$gtable)
  
})

plot <- cowplot::plot_grid(plotlist = pl, ncol = 1)

ggsave(plot = plot, filename = "../Figure 3/scores_heatmap_labeled.pdf", width = 7.5, height = 4)

pl <- lapply(breaks[1:(length(breaks) - 1)] %>% seq_along, function(n){
  
  m2 <- m[, breaks[n]:(breaks[n + 1] - 1)]
  
  colnames(m2) <- v2[breaks[n]:(breaks[n + 1] - 1)] %>%
    stringr::str_replace("^.*_", "Case #")
  
  hm <- pheatmap::pheatmap(m2,
                           color = colorRampPalette(c("white", "dodgerblue"))(4),
                           cellwidth = 16, cellheight = 14,
                           cluster_cols = FALSE, cluster_rows = FALSE,
                           display_numbers = TRUE,
                           number_color = "black", number_face = "bold",
                           number_format = "%.0f",
                           main = "Staining intensity for recurrent missense loci",
                           legend = FALSE,
                           angle_col = 45,
                           fontsize_number = 8,
                           fontsize_col = 12,
                           fontsize_row = 12
  )
  
  return(hm$gtable)
  
})

plot <- cowplot::plot_grid(plotlist = pl, ncol = 1)

ggsave(plot = plot, filename = "Figure_3B_caseID.pdf", width = 7.5, height = 5)


pl <- lapply(breaks[1:(length(breaks) - 1)] %>% seq_along, function(n){
  
  m2 <- m[, breaks[n]:(breaks[n + 1] - 1)]
  
  colnames(m2) %<>% stringr::str_replace("^[A-Z][0-9]+", "")
  
  hm <- pheatmap::pheatmap(m2,
                           color = colorRampPalette(c("white", "dodgerblue"))(4),
                           cellwidth = 16, cellheight = 14,
                           cluster_cols = FALSE, cluster_rows = FALSE,
                           display_numbers = TRUE,
                           number_color = "black", number_face = "bold",
                           number_format = "%.0f",
                           main = "Staining intensity for recurrent missense loci",
                           legend = FALSE,
                           angle_col = 0,
                           fontsize_number = 8,
                           fontsize_col = 12,
                           fontsize_row = 12
  )
  
  return(hm$gtable)
  
})

plot <- cowplot::plot_grid(plotlist = pl, ncol = 1)

ggsave(plot = plot, filename = "Figure_3B.pdf", width = 7.5, height = 3)

