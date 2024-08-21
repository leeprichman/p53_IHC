library(magrittr)
library(data.table)
library(ggplot2)

cv <- c("#c51162", "#aa00ff", "#0091ea", "#64dd17", "#ffab00", "#00b8d4", 
        "#d50000", "#6200ea", "#2962ff", "#a7ffeb", "#00c853", "#ff6d00", 
        "#aeea00", "#dd2c00")

dt <- "../Data/p53_scores.tsv" %>% data.table::fread()

dt <- dt[effect == "Missense", .SD, .SDcols = c("ID", "p53_intensity",
                                                    "p53_pct", "maxvaf")]

dt[, p53_pct := p53_pct %>% as.numeric]

dt[, p53s := ifelse(sum(p53_intensity == 3) >= 2, "3+", "0 - 2+"), by = "ID"]

dt2 <- dt[, .SD, .SDcols = c("ID", "maxvaf", "p53s")] %>% unique(by = "ID")

pvaf <- wilcox.test(maxvaf ~ p53s, dt2)$p.value

dt2[, p53s := paste0(p53s, ", n = ", ID %>% unique %>% length), by = "p53s"]

gvaf <- ggplot(dt2, aes(x = p53s, y = maxvaf)) +
  geom_boxplot(fill = "dodgerblue") +
  geom_jitter(color = "black", width = 0.2, height = 0) +
  annotate("text", size = 6, x = 1.5, y = 99, fontface = "bold",
           label = "***") +
  theme_bw() +
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", face = "bold", size = 24),
        panel.grid = element_blank()
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  xlab("\np53 IHC intensity") +
  ylab("VAF\n") +
  ggtitle("IHC intensity by variant\nallele fraction (VAF)\n")

ggsave(plot = gvaf, "Figure_4C.pdf", height = 6, width = 6)

