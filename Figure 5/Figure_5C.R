library(magrittr)
library(data.table)
library(ggplot2)

cv <- c("#c51162", "#aa00ff", "#0091ea", "#64dd17", "#ffab00", "#00b8d4", 
  "#d50000", "#6200ea", "#2962ff", "#a7ffeb", "#00c853", "#ff6d00", 
  "#aeea00", "#dd2c00")

dt <- "../Data/Survival.tsv" %>% data.table::fread()

dt <- dt %>% unique(by = "surv_ID")

# null vs all missense, then null vs high TS

co <- quantile(dt[effect == "Missense", Transcription], na.rm = TRUE, 0.75)

# n.b. two NA values for delins and a 2 base sub are excluded cause no PHANTM
hist <- ggplot(dt[effect == "Missense"], aes(x = Transcription)) +
  geom_histogram(bins = 55, color = "black", fill = "#aa00ff",
                 aes(y=after_stat(count/sum(count)))) +
  geom_vline(xintercept = co, linetype = "dashed", linewidth = 1) +
  annotate("text", x = co + 15, y = 0.12, label = "Upper 25%ile\nHigh RTA", size = 6) +
  theme_bw() +
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(linewidth = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", face = "bold", size = 24),
        panel.grid = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  ylab("% of missense mutants\n") +
  xlab("\nTranscriptional activity") +
  ggtitle("Distribution of p53 missense\nresidual transcriptional activity (RTA)")

ggsave(hist, filename = "Figure_5C.pdf", width = 8, height = 6)
