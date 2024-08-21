library(magrittr)
library(data.table)
library(ggplot2)

scores <- "../Data/p53_scores.tsv" %>% data.table::fread()

# p.H179delinsLN  and p.RH213QY excluded
scores <- scores[!is.na(Transcription)]

scores[, sum := sum(as.numeric(p53_intensity) == 3), by = "ID"]

dt <- scores[, .SD %>% unique, .SDcols = c("ID","sum", "Transcription")]

dt[, pos := ifelse(sum >= 2, "pos", "neg")]
#dt[, pos := ifelse(sum >= 9, "pos", "neg")]

pred <- ROCR::prediction(predictions = dt$Transcription,
                         labels = dt$pos, label.ordering = c("pos", "neg"))
perf <- ROCR::performance(prediction.obj = pred, measure = "auc")

auc <- perf@y.values[[1]]

roct <- data.table::data.table(thresh = pred@cutoffs[[1]], 
          tp = pred@tp[[1]], fp = pred@fp[[1]], tn = pred@tn[[1]], 
          fn = pred@fn[[1]], tpr = pred@tp[[1]]/pred@n.pos[[1]], 
          fpr = pred@fp[[1]]/pred@n.neg[[1]], ppv = pred@tp[[1]]/(pred@tp[[1]] + 
          pred@fp[[1]]), npv = pred@tn[[1]]/(pred@tn[[1]] + 
          pred@fn[[1]]), auc = auc * 100)

dt[, fact := pos %>% factor(levels = c("neg", "pos"), ordered = TRUE)]

sglm <- glm(fact~Transcription, data = dt, family = "binomial") %>% summary

pval <- sglm$coefficients["Transcription", "Pr(>|z|)"] %>% signif(digits = 2)

roct[, tpr2 := -fpr + 1] %>% .[, fpr2 := -tpr + 1]

roc <- ggplot(roct, aes(x = fpr2, y = tpr2)) +
  #ggrepel::geom_text_repel(aes(label = thresh), max.overlaps = 100) +
  #geom_point(color = "dodgerblue", shape = "cross") +
  geom_path(color = "dodgerblue", size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("label", x = 0.8, y = 0.1, size = 6, fontface = "bold",
           label = paste("AUC = ", auc %>% signif(digits = 2), 
                         "\np = ", pval, sep = "")
           ) +
  theme_bw() +
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14)) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  ggtitle("Receiver operating characteristic", subtitle = "p53 residual transcriptional activity as a\npredictor of increased p53 IHC staining (3+)")

ggsave(plot = roc, "Figure_4A.pdf", height = 7, width = 7)
 
