library(data.table)
library(magrittr)
library(ggplot2)
library(pheatmap)

dt <- data.table::fread("../Data/p53_scores.tsv")

consensus <- function(x){
  
  t <- table(x)
  
  score <- ifelse(any(t > 1), names(t)[which.max(t)], as.character(NA))
  
  return(as.numeric(score))
  
}

dt[, con_score := p53_intensity %>% consensus, by = "ID"]

d_all <- dt %>% data.table::copy() %>%
  .[protein_change == "WT", protein_change := "Wild-type"] %>%
  .[effect == "Null", protein_change := "Splice/nonsense/frameshift"] %>%
  .[effect == "Missense", protein_change := "Missense"] %>%
  .[, ID %>% unique %>% length, by = c("protein_change", "con_score")]

m <- d_all

m %<>% dcast(protein_change ~ con_score, value.var = "V1")

rn <- m[, protein_change]

cn <- c("No consensus", "0", "1+", "2+", "3+")

m <- m[, .SD, .SDcols = 2:ncol(m)] %>% as.matrix()

rownames(m) <- rn

m[is.na(m)] <- 0

raw_m <- m

m <- apply(m, 1, function(v) prop.table(v) %>%
             sprintf(fmt = "%.2f") %>% as.numeric) %>% t

colnames(m) <- cn

display_mod <- function(m){
  
  m <-  (m * 100) %>% sprintf(fmt = "%.0f%%") %>%
    matrix(ncol = ncol(m))
  
  m[m == "0%"] <- ""
  
  return(m)
  
}

anno <- paste0(m %>% display_mod, " (", raw_m, ")") %>%
  matrix(nrow = nrow(m))

hm <- pheatmap::pheatmap(m,
                         color = colorRampPalette(c("white", "red"))(100),
                         cellwidth = 60, cellheight = 40,
                         cluster_cols = FALSE, cluster_rows = FALSE,
                         display_numbers = anno,
                         number_color = "black", number_face = "bold",
                         main = "Consensus p53 staining intensity",
                         legend = FALSE,
                         angle_col = 45,
                         fontsize_number = 12,
                         fontsize = 16,
                         fontsize_col = 16,
                         fontsize_row = 16
)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(hm, "Figure_1B.pdf", width = 8.5, height = 4)


