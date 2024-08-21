library(magrittr)
library(data.table)
library(ggplot2)

cv <- c("#c51162", "#aa00ff", "#0091ea", "#64dd17", "#ffab00", "#00b8d4", 
  "#d50000", "#6200ea", "#2962ff", "#a7ffeb", "#00c853", "#ff6d00", 
  "#aeea00", "#dd2c00")

dt <- data.table::fread("../Data/p53_scores.tsv")

dt <- dt[, .SD, .SDcols = c("ID", "effect", "Training_level", "p53_intensity")]

effects <- dt[, effect %>% sort %>% unique]

consensus <- function(x){
  
  t <- table(x)
  
  score <- ifelse(any(t > 1), names(t)[which.max(t)], as.character(NA))
  
  return(as.numeric(score))
  
}

condt <- dt[, p53_intensity %>% consensus, by = c("ID", "effect")]

condt %>% setnames("V1", "p53_intensity")

condt[, Training_level := "Consensus"]

# condt %>% data.table::fwrite("Consensus_ratings.tsv", sep = "\t")

dt <- data.table::rbindlist(list(dt, condt), use.names = TRUE)

posl <- list(3, 0, 1)

pt <- lapply(1:1000, function(bs){
  
  print(bs)
  
  set.seed(bs)
  
  wl <- sample(dt[, ID %>% unique],
               size = floor(dt[, ID %>% unique %>% length] * 0.8))
  
  bt <- dt[ID %chin% wl]
  
  tab <- lapply(dt[, Training_level %>% unique], function(t){
    
    it <- bt[Training_level == t]
    
    ot <- lapply(effects %>% seq_along, function(n){
      
      posv <- posl[[n]]
      
      d <- it %>% data.table::copy()
      
      d[, case := ifelse(effect == effects[n], "case", "control")]
      
      d[, test := ifelse(p53_intensity %in% posv, 1, 0)]
      
      d <- d[, .N, by = c("case", "test")] %>% dcast(case ~ test, value.var = "N")
      
      m <- d[, .SD, .SDcols = 2:ncol(d)] %>% as.matrix
      
      m[is.na(m)] <- 0 
      
      rownames(m) <- d[, case]
      
      # if the bootstrap caught no positives for a small set like NULLs
      if (ncol(m) < 2){
        
        mn <- !c("0", "1") %chin% colnames(m)
        
        mn <- c("0", "1")[mn]
        
        mn <- matrix(c(0,0), ncol = 1, dimnames = list(rownames(m), mn))
        
        m <- cbind(m, mn)
        
      }
      
      tp <- m["case", "1"]
      
      fp <- m["control", "1"]
      
      tn <- m["control", "0"]
      
      fn <- m["case", "0"]
      
      ppv <- tp / (tp + fp)
      
      npv <- tn / (tn + fn)
      
      sens <- tp / (tp + fn)
      
      spec <- tn / (tn + fp)
      
      return(data.table(effect = effects[n], ppv, npv,
                        Sensitivity = sens, Specificity = spec))
      
    }) %>% data.table::rbindlist()
    
    ot[, Training_level := t]
    
    return(ot)  
    
  }) %>% data.table::rbindlist()
  
  tab[, seed := bs]
  
  return(tab)
  
}) %>% data.table::rbindlist()

plott <- pt %>% melt(id.vars = c("Training_level", "seed", "effect"))

gt <- plott[, list(mean(value), sd(value)), by = c("Training_level", "effect", "variable")]

gt[effect == "WT", effect := "Wild-type"]

gt %>% data.table::fwrite("../Figure 1/sensspectable.tsv", sep = "\t")

p_all <- ggplot(gt[variable %in% c("Sensitivity", "Specificity")], aes(x = effect, y = V1)) +
  geom_col(aes(fill = variable), position = "dodge") +
  geom_errorbar(aes(color = effect, ymin = V1 - (2 * V2),
                    ymax = V1 + (2*V2)), position = "dodge",
                    width = 0.4, linewidth = 0.7) +
  facet_wrap(variable~Training_level, ncol = 3, scales = "free_x") +
  scale_color_manual(values = replicate(4, "black")) +
  scale_fill_manual(values = cv) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold", color = "black"),
        axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1.5))
        )

p <- ggplot(gt[Training_level == "Consensus" & 
                 variable %in% c("Sensitivity", "Specificity")],
            aes(x = effect, y = V1)) +
  geom_col(aes(fill = variable), position = "dodge", width = 0.5) +
  geom_errorbar(aes(color = effect, ymin = V1 - (2 * V2),
                    ymax = V1 + (2*V2)), position = "dodge",
                width = 0.3, linewidth = 0.7) +
  facet_wrap(~variable, nrow = 1, scales = "free_x") +
  scale_color_manual(values = replicate(4, "black")) +
  scale_fill_manual(values = cv) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = rel(1)),
        axis.text = element_text(face = "bold", color = "black"),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1))
        ) +
  ggtitle("Consensus p53 staining intensity\nas a diagnostic test")

ggsave(p, filename = "Figure_1C.pdf", width = 6, height = 3)
