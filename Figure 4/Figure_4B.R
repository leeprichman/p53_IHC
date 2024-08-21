library(magrittr)
library(data.table)
library(ggplot2)

cv <- c("#c51162", "#aa00ff", "#0091ea", "#64dd17", "#ffab00", "#00b8d4", 
        "#d50000", "#6200ea", "#2962ff", "#a7ffeb", "#00c853", "#ff6d00", 
        "#aeea00", "#dd2c00")

dt <- "../Data/p53_scores.tsv" %>% data.table::fread()

dt <- dt[effect != "WT"]

dt[effect == "Missense", p53s :=
     ifelse(sum(p53_intensity == 3) >= 2, "3+", "0 - 2+"), by = "ID"]

dt[effect == "Null", p53s :=
     ifelse(sum(p53_intensity == 0) >= 2, "0", "1+ - 3+"), by = "ID"]

dt[, Allele_status := LOH]

dt[is.na(Allele_status) | Allele_status == "Noisy", Allele_status := "Unknown"]

dt[Allele_status %chin% c("No CNV", "") & !(is.na(Second) | Second == ""),
   Allele_status := "2nd Missense"]

dt[Allele_status == "Gain", Allele_status := "No CNV"]

dt[Allele_status == "CN LOH", Allele_status := "LOH"]

dt[maxvaf > 50 & Allele_status == "No CNV", Allele_status := "Inferred"]

d <- dt[, .SD %>% unique, .SDcols = c("ID", "effect", "p53s", "Allele_status")]

d <- d[Allele_status != "Unknown"]

d[Allele_status == "No CNV", `2nd allele` := "Intact"]

d[Allele_status %chin% c("2nd Missense", "LOH", "Inferred"), `2nd allele` := "Disrupted"]

d <- d[, ID %>% unique %>% length, by = c("effect", "p53s", "2nd allele")]

d[, V2 := V1/sum(V1), by = "effect"]

d[, effect := paste0(effect, ", n = ", sum(V1)), by = "effect"]

p1 <- ggplot(d, aes(x = p53s, y = V2)) +
  geom_col(aes(fill = `2nd allele`), color = "white") +
  facet_wrap(~effect, scales = "free")  +
  scale_fill_manual(values = cv[c(12,3)]) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1),
        panel.border = element_blank(),
        strip.background = element_rect(linewidth = 1),
        strip.text = element_text(face = "bold", size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", face = "bold", size = 24),
        panel.grid = element_blank(),
        legend.position = c(0.9,0.82)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    xlab("\np53 IHC intensity") +
    ylab("% of cases\n") +
    ggtitle("IHC intensity by 2nd allele status\n")

ggsave(plot = p1, filename ="Figure_4B.pdf", height = 6, width = 8)

#fisher tests
d[effect %like% "Missense", V1, by = c("2nd allele", "p53s")] %>%
  dcast(`2nd allele` ~ p53s, value.var = "V1") %>% .[, .SD, .SDcols = 2:3] %>%
  as.matrix %>%
  fisher.test

d[effect %like% "Null", V1, by = c("2nd allele", "p53s")] %>%
  dcast(`2nd allele` ~ p53s, value.var = "V1") %>% .[, .SD, .SDcols = 2:3] %>%
  as.matrix %>%
  fisher.test

