library(survminer)
library(survival)
library(magrittr)
library(data.table)
library(ggplot2)
library(lubridate)

cv <- c("#c51162", "#aa00ff", "#0091ea", "#64dd17", "#ffab00", "#00b8d4", 
  "#d50000", "#6200ea", "#2962ff", "#a7ffeb", "#00c853", "#ff6d00", 
  "#aeea00", "#dd2c00")

cvv <- cv[c(1:3, 4:7, 8, 5)]

names(cvv) <- c("WT", "Missense", "Null", "1+", "2+", "3+", "0",
                "High RTA", "Low RTA")

dt <- "../Data/Survival.tsv" %>% data.table::fread()

dt <- dt %>% unique(by = "surv_ID")

co <- quantile(dt[effect == "Missense", Transcription], na.rm = TRUE, 0.75)

# null vs all missense, then null vs high TS

dt[, stain := paste0(p53_intensity, "+")]

dt[p53_intensity == 0, stain := "0"]

dt[, p53_status := factor(effect, levels = c("WT", "Null", "Missense"))]

dt[, p53_IHC := factor(stain, levels = c("1+", "0", "2+", "3+"))]

dt2 <- dt

dt2[effect == "WT", trans := "WT"]

dt2[effect == "Null", trans := "Null"]

dt2[effect == "Missense" & Transcription >= co, trans := "High RTA"]

dt2[effect == "Missense" & Transcription < co, trans := "Low RTA"]

dt2 <- dt2[!is.na(trans)]

dt2[, p53_activity := factor(trans, levels = c("WT", "Low RTA", "High RTA", "Null"))]

fit_mut <- survfit(Surv(surv, event) ~ p53_status, data = dt)

fit_stain <- survfit(Surv(surv, event) ~ p53_IHC, data = dt)

gstain <- ggsurvplot(fit_stain)

fit_trans <- survfit(Surv(surv, event) ~ p53_activity, data = dt2)

gtrans <- ggsurvplot(fit_trans)

fit_v <-list(fit_stain, fit_mut, fit_trans)

title_v <- c("p53 IHC intensity", "p53 status", "p53 RTA")

var_v <- c("p53_IHC", "p53_status", "p53_activity")

dl <- list(dt,dt,dt2)

gl <- lapply(fit_v %>% seq_along, function(n){
  
  print(n)
  
  g1 <- ggsurvplot(fit_v[[n]], palette = cvv, surv.median.line = "hv",
             linetype = "dashed", conf.int = FALSE,
             legend.labs = names(fit_v[[n]]$strata) %>%
               stringr::str_replace("^.*=", ""),
             ggtheme = theme_survminer() + theme(
               legend.title = element_blank(),
               legend.text = element_text(size = 12),
               axis.title.x = element_text(face = "bold", size = 12),
               axis.title.y = element_text(face = "bold", size = 12),
               axis.text.x = element_text(size = 12),
               axis.text.y = element_text(size = 12),
               plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
             )) +
    ylab("Survival probability\n") +
    xlab("\nOverall survival (days)") +
    ggtitle(title_v[n])
  
  if (title_v[n] == "p53 RTA") g1$plot <-
    g1$plot + guides(color=guide_legend(nrow=2)) 
  
  ff <- formula(paste("Surv(surv, event) ~", var_v[n]))
  
  g2 <- coxph(ff, data = dl[[n]]) %>% ggforest(fontsize = 1, data = dl[[n]])
  
  return(g1$plot)
  
})

fig <- cowplot::plot_grid(plotlist = gl, scale = 1, ncol = 1)

ggsave("Survival plots.pdf", plot = fig, height = 10, width = 4)
