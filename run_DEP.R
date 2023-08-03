#### Code title: run_DEP
#### Description: Differential expression analysis for all variables using this study
#### Last log: 23.07.14

## Settings
rm(list = ls())

## Needed libraries
library(tidyverse)
library(broom)

## Load data
in.data = list()
in.data$Protein = readxl::read_xlsx('files/TableS7_Expression_data.xlsx', sheet = 3)
in.data$Phospho = readxl::read_xlsx('files/TableS7_Expression_data.xlsx', sheet = 4)
in.data$Acetyl = readxl::read_xlsx('files/TableS7_Expression_data.xlsx', sheet = 5)

s = readxl::read_xlsx('files/TableS1_Data_overview.xlsx', sheet = 2)
exp.proteom = readxl::read_excel('files/TMTchannel_info.xlsx', sheet=1) %>% as.data.frame() %>% dplyr::rename(sid=3)

## Imputation
perc_imputation = 0.3 # same as NMF

imp.data = list()
for (dt in names(in.data)){
  d1 = in.data[[dt]] %>% column_to_rownames('ID') %>% select(contains('RE-')) %>% as.data.frame() %>% t()
  d1 = scale(d1)
  d1 = t(d1)
  d1 = as.data.frame(d1)
  d1 = d1 %>% filter(rowSums(is.na(.))<=ncol(.)*perc_imputation)
  mm1 = impute::impute.knn(as.matrix(d1), k=5)
  imp.data[[dt]] = mm1$data
}

## Add variables to test
s1 = s %>% mutate(is_Subtype1 = ifelse(Subtype=='1', 'Y', 'N'),
                  is_Subtype2 = ifelse(Subtype=='2', 'Y', 'N'),
                  is_Subtype3 = ifelse(Subtype=='3', 'Y', 'N'),
                  is_Subtype4 = ifelse(Subtype=='4', 'Y', 'N'),
                  is_Subtype5 = ifelse(Subtype=='5', 'Y', 'N'),
                  is_wgd = ifelse(WGD=='TRUE', 'Y', NA),
                  is_wgd = ifelse(WGD=='FALSE', 'N', is_wgd))
s1 = merge(s1, exp.proteom, by.x='Sample.ID', by.y='sid', all.x=T)
s1 = s1 %>% mutate(Set = as.factor(Set))
s1$sid = s1$Sample.ID

### Function to run
runReg <- function(imp.data1, samples_dt1, dt1, in.samples, in.data1, test_col, to_test, group_subset=F){
  samples_to_use = intersect(samples_dt, in.samples$sid)
  dat = imp.data1[[dt]][,samples_to_use]
  
  if (group_subset){
    res.deg = dat %>% as.data.frame() %>% rownames_to_column('ID') %>% gather(sid, ratio, -ID) %>% 
      merge(., in.samples %>% select(Sex, !!as.name(test_col), Set, sid), by='sid') %>%
      group_by(ID) %>%
      do(fitTest = tidy(lm(ratio ~ !!as.name(test_col) + Set, data = ., na.action = "na.omit"))) %>% 
      unnest(fitTest) %>%
      filter(term==to_test) %>%
      mutate(padj = p.adjust(p.value, method = 'BH', n = nrow(dat)))
    res.deg = res.deg %>% merge(., in.data1 %>% select(ID, gene=`Gene Symbol`), by='ID')
  } else {
    res.deg = dat %>% as.data.frame() %>% rownames_to_column('ID') %>% gather(sid, ratio, -ID) %>% 
      merge(., in.samples %>% select(DX, Sex, !!as.name(test_col), Set, sid, Group), by='sid') %>%
      group_by(ID) %>%
      do(fitTest = tidy(lm(ratio ~ !!as.name(test_col) + Group + Set, 
                           data = ., na.action = "na.omit"))) %>% 
      unnest(fitTest) %>%
      filter(term==to_test) %>%
      mutate(padj = p.adjust(p.value, method = 'BH', n = nrow(dat)))
    res.deg = res.deg %>% merge(., in.data1 %>% select(ID, gene=`Gene Symbol`), by='ID')
  }
  
  print (paste('FDR <5%', sum(res.deg$padj < 0.05)))
  print (paste('FDR <10%', sum(res.deg$padj < 0.1)))
  return(res.deg)
}

## Run regression
dts = names(imp.data)
output.run = list()
for (dt in dts){
  list.run = list()
  samples_dt = colnames(imp.data[[dt]])
  
  ## NAT vs. Tumor
  run_name = 'NAT_Tumor'
  s2 = data.frame(sid = colnames(imp.data[[dt]])) %>% 
    mutate(isTumor = ifelse(grepl('-T', sid), 'Y', 'N'), isTumor = factor(isTumor, levels=c('N', 'Y'))) %>%
    merge(., s1%>%select(sid, Sex, Age), by='sid', all.x=T)
  s2 = merge(s2, exp.proteom, by='sid', all.x=T) %>% mutate(Set = as.factor(Set))
  
  list.run[[run_name]][['sample']] = s2
  list.run[[run_name]][['DEG']] = runReg(imp.data1 = imp.data, samples_dt1 = samples_dt, dt1 = dt, 
                                         in.samples = list.run[[run_name]][['sample']], in.data1=in.data[[dt]], 
                                         test_col = 'isTumor', to_test = 'isTumorY', group_subset=T)
  
  ## NMF Subtype
  ## Subtype1
  run_name = 'Subtype1'
  name1 = 'Y'; name2 = 'N'
  list.run[[run_name]][['sample']] = s1 %>% filter(!is.na(Subtype)) %>% 
    mutate(isSubtype = ifelse(Subtype==1, 'Y', 'N'), isSubtype = factor(isSubtype, levels=c(name2, name1))) %>% 
    filter(sid %in% samples_dt )
  list.run[[run_name]][['DEG']] = runReg(imp.data1 = imp.data, samples_dt1 = samples_dt, dt1 = dt, 
                                         in.samples = list.run[[run_name]][['sample']], 
                                         in.data1=in.data[[dt]], 
                                         test_col = 'isSubtype', to_test = 'isSubtypeY')
  head(list.run[[run_name]][['DEG']])
  
  
  ## Subtype 2
  run_name = 'Subtype2'
  name1 = 'Y'; name2 = 'N'
  list.run[[run_name]][['sample']] = s1 %>% filter(!is.na(Subtype)) %>% 
    mutate(isSubtype = ifelse(Subtype==2, 'Y', 'N'), isSubtype = factor(isSubtype, levels=c(name2, name1))) %>% 
    filter(sid %in% samples_dt )
  list.run[[run_name]][['DEG']] = runReg(imp.data1 = imp.data, samples_dt1 = samples_dt, dt1 = dt, 
                                         in.samples = list.run[[run_name]][['sample']], 
                                         in.data1=in.data[[dt]], 
                                         test_col = 'isSubtype', to_test = 'isSubtypeY')
  head(list.run[[run_name]][['DEG']])
  
  
  ## Subtype 3
  run_name = 'Subtype3'
  name1 = 'Y'; name2 = 'N'
  list.run[[run_name]][['sample']] = s1 %>% filter(!is.na(Subtype)) %>% 
    mutate(isSubtype = ifelse(Subtype==3, 'Y', 'N'), isSubtype = factor(isSubtype, levels=c(name2, name1))) %>% 
    filter(sid %in% samples_dt )
  list.run[[run_name]][['DEG']] = runReg(imp.data1 = imp.data, samples_dt1 = samples_dt, dt1 = dt, 
                                         in.samples = list.run[[run_name]][['sample']], 
                                         in.data1=in.data[[dt]], 
                                         test_col = 'isSubtype', to_test = 'isSubtypeY')
  head(list.run[[run_name]][['DEG']])
  
  
  ## Subtype 4
  run_name = 'Subtype4'
  name1 = 'Y'; name2 = 'N'
  list.run[[run_name]][['sample']] = s1 %>% filter(!is.na(Subtype)) %>% 
    mutate(isSubtype = ifelse(Subtype==4, 'Y', 'N'), isSubtype = factor(isSubtype, levels=c(name2, name1))) %>% 
    filter(sid %in% samples_dt )
  list.run[[run_name]][['DEG']] = runReg(imp.data1 = imp.data, samples_dt1 = samples_dt, dt1 = dt, 
                                         in.samples = list.run[[run_name]][['sample']], 
                                         in.data1=in.data[[dt]], 
                                         test_col = 'isSubtype', to_test = 'isSubtypeY')
  head(list.run[[run_name]][['DEG']])
  
  
  ## Subtype 5
  run_name = 'Subtype5'
  name1 = 'Y'; name2 = 'N'
  list.run[[run_name]][['sample']] = s1 %>% filter(!is.na(Subtype)) %>% 
    mutate(isSubtype = ifelse(Subtype==5, 'Y', 'N'), isSubtype = factor(isSubtype, levels=c(name2, name1))) %>% 
    filter(sid %in% samples_dt )
  list.run[[run_name]][['DEG']] = runReg(imp.data1 = imp.data, samples_dt1 = samples_dt, dt1 = dt, 
                                         in.samples = list.run[[run_name]][['sample']], 
                                         in.data1=in.data[[dt]], 
                                         test_col = 'isSubtype', to_test = 'isSubtypeY')
  head(list.run[[run_name]][['DEG']])
  
  
  ## Whole Genome Doubling
  run_name = 'WGD_Y_vs_N'
  name1 = 'Y'; name2 = 'N'
  list.run[[run_name]][['sample']] = s1 %>% filter(is_wgd %in% c(name1, name2)) %>% 
    mutate(is_wgd = factor(is_wgd, levels=c(name2, name1)), Group=Subtype) %>% filter(sid %in% samples_dt )
  list.run[[run_name]][['DEG']] = runReg(imp.data1 = imp.data, samples_dt1 = samples_dt, dt1 = dt, 
                                         in.samples = list.run[[run_name]][['sample']], 
                                         in.data1=in.data[[dt]], 
                                         test_col = 'is_wgd', to_test = 'is_wgdY')
  
  ## Save
  output.run[[dt]][['DEP']] <- list.run
}

##