#### Code title: run_Survival
#### Description: Feature-wise survival analysis
#### Last log: 23.07.14

## Settings
rm(list = ls())

## Needed libraries
library(tidyverse)
library(survminer)
library(survival)

## Load sample information and datasets
# You can downloaded from supplementary file

s = readxl::read_xlsx('TableS1', sheet = 2)

in.data = list()
in.data$Protein = readxl::read_xlsx('TableS7', sheet = 3)
in.data$Phospho = readxl::read_xlsx('TableS7', sheet = 4)
in.data$Acetyl = readxl::read_xlsx('TableS7', sheet = 5)

in.data$Protein = in.data$Protein %>%
  separate_rows(Gene_name, sep="; ") %>%  filter(!is.na(Gene_name)) %>% dplyr::select(-ID) %>% dplyr::select(ID=Gene_name,everything())

in.data$Phospho = in.data$Phospho %>% 
  separate_rows(Gene_name, sep="; ") %>% filter(!is.na(Gene_name)) %>% mutate(site= do.call(rbind.data.frame,strsplit(ID, '_'))[[2]], ID3=paste(Gene_name,site,sep="_")) %>%
  dplyr::select(-ID,-Gene_name, -site) %>% dplyr::select(ID=ID3,everything()) 

in.data$Acetyl <- in.data$Acetyl %>%
  separate_rows(Gene_name, sep="; ") %>% filter(!is.na(Gene_name)) %>% mutate(site= do.call(rbind.data.frame,strsplit(ID, '_'))[[2]], ID3=paste(Gene_name,site,sep="_")) %>%
  dplyr::select(-ID,-Gene_name, -site) %>% dplyr::select(ID=ID3,everything())

## Run feature-wise survival analysis
out = list()
dts = names(in.data)
for (dt in dts){
  res = data.frame()
  in.data1 = in.data[[dt]]
  
  for (i in seq(1, nrow(in.data1))){
    print (c(dt, i))
    id = in.data1[i,]$ID
    
    mm = in.data[[dt]] %>% filter(ID==id) %>% dplyr::select(contains('-T')) %>% gather(sid, value)
    
    s.high = mm %>% filter(value >= quantile(mm$value, 0.5, na.rm = T)) %>% pull(sid)
    s.low = mm %>% filter(value < quantile(mm$value, 0.5, na.rm = T)) %>% pull(sid)
    a1 = s %>% dplyr::select(sid, NSCLC_death, OS) %>% mutate(NSCLC_death2 = ifelse(NSCLC_death=='Y', 1, 0)) %>%  
      mutate(Group = ifelse(sid %in% s.high, 'high', NA), Group = ifelse(sid %in% s.low, 'low', Group)) %>% 
      filter(!is.na(Group)) %>% mutate(Group=factor(Group, levels=c("high","low")))
    fit = survfit(Surv(OS, NSCLC_death2) ~ Group, data = a1) 
    cox = coxph(Surv(OS, NSCLC_death2) ~ Group, data = a1)
    direction = cox$coefficients
    
    ## Add to result data frame
    res = rbind.data.frame(res, 
                           surv_pvalue(fit) %>%
                             mutate(data_type=dt, ID=id,
                                    cause_high_survival = ifelse(direction<0, 'N','Y')) %>% 
                             dplyr::select(ID, data_type, pval, cause_high_survival))
  }
  
  out[[dt]] <- res
}
