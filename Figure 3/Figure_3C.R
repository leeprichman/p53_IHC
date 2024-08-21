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

dt <- dt[, .SD, .SDcols = c("protein_change", "p53_intensity", "ID", "Training_level")]

m <- dt %>% dcast(protein_change + ID ~ Training_level, value.var = "p53_intensity")

rn <- m[, paste0(protein_change, "_", ID)]

m <- m[, 3:5] %>% as.matrix

rownames(m) <- rn

dm <- dist(m) %>% as.matrix

dm[upper.tri(dm, diag = TRUE)] <- NA

cjt <- data.table::CJ(Var1 = rn, Var2 = rn) %>%
  .[stringr::str_extract(Var1, "[A-Z][0-9]+[A-Z]") ==
      stringr::str_extract(Var2, "[A-Z][0-9]+[A-Z]")] %>%
  .[, edge := 1]

dm %<>% reshape2::melt() %>% na.omit %>% as.data.table()

dm <- merge(dm, cjt, all.x = TRUE, by = c("Var1", "Var2"))

dm[is.na(edge), edge := 0]
# 
# ss <- lapply(dt[, locus %>% unique], function(l){
#   
#   m <- dt[locus == l]
#   
#   m %<>% dcast(protein_change + ID ~ Training_level, value.var = "p53_intensity")
#   
#   rn <- m[, paste0(protein_change, "_", ID)]
#   
#   m <- m[, 3:5] %>% as.matrix
#   
#   rownames(m) <- rn
#   
#   # doubkle counting here with the diagonal
#   dm <- dist(m) %>% as.matrix
#   
#   dm[upper.tri(dm, diag = TRUE)] <- NA
#   
#  dm %<>% reshape2::melt() %>% as.data.table() %>% na.omit %>%
#     .[ Var1 != Var2]
#   
#   dm[, win := ifelse(
#     stringr::str_extract(Var1, "[0-9]+[A-Z]") ==
#       stringr::str_extract(Var2, "[0-9]+[A-Z]"),
#     1, 0)]
#   
#   # average outside edge - average inside edge / max edge distance
# 
#   dm[, locus := l]
#   
#   return(dm)
#   
# }) %>% data.table::rbindlist()

# score <- (mean(ss[win == 0, value]) - mean(ss[win == 1, value])) / ss[, value %>% max]

score <- (mean(dm[edge == 0, value]) - mean(dm[edge == 1, value])) / dm[, value %>% max]

lab <- paste0("Silhouette score = ", score %>% sprintf(fmt = "%.2f"))

labp <- wilcox.test(value ~ edge, data = dm)$p.value %>%
  sprintf(fmt = "%.2f") %>% paste("p =", .)
  

bpp <- ggplot(dm, aes(x  = factor(edge), y = value)) +
  geom_boxplot(fill = "dodgerblue", width = 0.3, size = 0.75) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = unit(6, "points")),
        axis.title = element_text(size = unit(8, "points"), face = "bold"),
        ) +
  annotate("label", label = lab, x = 1.8, y = 6.75, fontface = "bold", size = 2) +
  annotate("segment", x = 1, xend = 2, y = 5.5, yend = 5.5, size = 0.75) +
  annotate("text", y = 5.8, x = 1.5, label = labp, size = 2) +
  scale_x_discrete(labels = c("Different\nmissense", "Identical\nmissense")) +
  ylim(c(0, 7)) +
  xlab(NULL) +
  ylab("Distance between\nIHC scores") 

ggsave(plot = bpp, filename = "Figure_3C.pdf", height = 2, width = 2)

