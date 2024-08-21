library(ggplot2)
library(magrittr)
library(data.table)
library(pheatmap)
library(ROCR)
library(psych)
library(irr)

lapply(paste("Figure", 4:5), function(f){
  
  setwd(f)
  
  scripts <- list.files(pattern = "\\.R$")
  
  lapply(scripts, function(s){
    
    print(s)
    
    source(s)
    
  })
  
  setwd("..")
  
})

