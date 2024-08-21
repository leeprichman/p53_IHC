library(magrittr)
library(data.table)
library(ggplot2)

scores <- "../Data/p53_scores.tsv" %>% data.table::fread()

scores[, ID %>% unique %>% length, by = "effect"]
#      effect V1
# 1: Missense 80
# 2:       WT 34
# 3:     Null 20

scorers <- scores[, Training_level %>% unique]

sl <- scores[, .(ID, p53_intensity, Training_level)] %>%
  dcast(ID~Training_level, value.var = "p53_intensity")

irr::kappam.fleiss(sl[, 2:4])

m <- matrix(replicate(length(scorers), scorers), nrow = length(scorers), ncol = length(scorers))

m <- matrix(paste(m, t(m), sep = ";"), nrow = length(scorers))

comps <- m[lower.tri(m)]

compt <- lapply(comps, function(n){
  
  v <- n %>% strsplit(split = ";") %>% unlist
  
  t <- data.table::data.table(v[1],v[2])
  
  return(t)
  
}) %>% data.table::rbindlist()

# plot correlation and get spearmans
pl <- lapply(1:nrow(compt), function(n){
  
  x <- scores[Training_level == compt[n, V1], p53_pct %>% as.numeric] / 100
  
  y <- scores[Training_level == compt[n, V2], p53_pct %>% as.numeric] / 100
  
  dt <- data.table::data.table(x, y)
  
  setnames(dt, names(dt), compt[n] %>% unlist)
  
  dt <- na.omit(dt)
  
  ct <- cor.test(dt[, 1] %>% unlist, dt[, 2] %>% unlist, method = "pearson")
  
  #rho <- paste("Spearman's \u2374" , "=" , ct$estimate %>% round(digits = 2))
  rsq <- paste("r\u00B2 =" , ct$estimate^2 %>% round(digits = 2))
  
  pval <- paste("p =", signif(ct$p.value, digits = 2))
  
  p1 <- ggplot(dt, aes(x = .data[[names(dt)[1]]], y = .data[[names(dt)[2]]])) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", level = 0) +
    annotate("label", x = 0.8, y = 0.1, label = paste0(rsq, "\n", pval),
             fontface = "bold", size = 4) +
    theme_bw() +
    scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black", size = 14),
          axis.title = element_text(size = 18, face = "bold"))
  
  # now intensity scores
  x <- scores[Training_level == compt[n, V1], p53_intensity %>% as.numeric]
  
  y <- scores[Training_level == compt[n, V2], p53_intensity %>% as.numeric]
  
  dt <- data.table::data.table(x, y)
  
  setnames(dt, names(dt), compt[n] %>% unlist)
  
  dt <- na.omit(dt)
  
  ct <- psych::cohen.kappa(dt, levels = 0:3)
  
  kappa <- ct$confid["weighted kappa", "estimate"] %>% signif(digits = 2) %>%
    paste("Cohen's weighted \u03BA =", .)
  
  ci <- ct$confid["weighted kappa", c("lower", "upper")] %>% signif(digits = 2) %>%
    paste(collapse = " - ") %>% paste("(95% CI: ", ., ")", sep = "")
  
  st <- scores[, .SD, .SDcols = c("ID", "p53_intensity", "Training_level")] %>%
    dcast(ID ~ Training_level, value.var = "p53_intensity")
      
    c <- st[, .SD, .SDcols = names(dt)]
      
    c %<>% dcast(formula(paste(names(dt)[1], "~", names(dt)[2])))
      
    rn <- c[, 1] %>% unlist
      
    c <- c[, 2:ncol(c)] %>% as.matrix
      
    rownames(c) <- rn
      
    # placeholder for missing columns
    if (!all(0:3 %in% as.numeric(colnames(c)))){
        
      mn <- 0:3 %>% .[!0:3 %in% as.numeric(colnames(c))]
        
      rv <- replicate(nrow(c), 0)
        
      for (m in mn){
          
        c <- cbind(rv, c)
          
        colnames(c)[1] <- m
          
      }
        
    }
      
    if (!all(0:3 %in% as.numeric(rownames(c)))){
        
      mn <- 0:3 %>% .[!0:3 %in% as.numeric(rownames(c))]
        
      rv <- replicate(ncol(c), 0)
        
      for (m in mn){
          
        c <- rbind(rv, c)
          
        rownames(c)[1] <- m
          
      }
        
    }
      
    c <- c[c("0", "1", "2", "3"), c("0", "1", "2", "3")]
    
    c <- c[rev(1:nrow(c)), ]
    
    pheatmap::pheatmap(c,
                       filename = paste(names(dt)[1], "_", names(dt)[2], "_hm.pdf", sep = ""),
                       color = colorRampPalette(c("#FFFFFF", "#99000D"))(100),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       border_color = "black",
                       cellheight = 30,
                       cellwidth = 30,
                       legend = FALSE,
                       main = paste0(kappa, "\n", ci),
                       display_numbers = TRUE,
                       number_color = "black",
                       number_format = "%.0f",
                       angle_col = "0",
                       width = 3,
                       height = 4
                       )

    return(p1)
    
})

pl <- cowplot::plot_grid(plotlist = pl, nrow = 1, scale = 0.9)

ggsave(plot = pl, "Figure_2B.pdf", height = 3.75, width = 12)

