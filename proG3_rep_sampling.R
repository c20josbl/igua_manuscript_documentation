##Script for random sampling of progenomes3 representative genomes

#Used to read proGenomes3_specI_lineageGTDB.tab file, downloaded from https://progenomes.embl.de/download.cgi
data=read.table(file.choose(), header = T,sep = "\t")
str(data)
data$family=as.factor(data$family)
data$order=as.factor(data$order)

specI_v4=rownames(data)
data=data.frame(data,specI_v4)

set.seed(1001)
library(dplyr)

stratified=data %>%
  group_by(genus) %>%
  sample_n(1)

set.seed(1001)
stratified1=sample(stratified$specI, 5000)
writeLines(stratified1, "stratified1.txt")

set.seed(1002)
stratified=data %>%
  group_by(genus) %>%
  sample_n(1)

set.seed(1002)
stratified2=sample(stratified$specI, 5000)
writeLines(stratified2, "stratified2.txt")

set.seed(1003)
stratified=data %>%
  group_by(genus) %>%
  sample_n(1)

set.seed(1003)
stratified3=sample(stratified$specI, 5000)
writeLines(stratified3, "stratified3.txt")

set.seed(1004)
stratified=data %>%
  group_by(genus) %>%
  sample_n(1)

set.seed(1004)
stratified4=sample(stratified$specI, 5000)
writeLines(stratified4, "stratified4.txt")


set.seed(1005)
stratified=data %>%
  group_by(genus) %>%
  sample_n(1)

set.seed(1005)
stratified5=sample(stratified$specI, 5000)
writeLines(stratified5, "stratified5.txt")

data2=table(data$genus)
data3=data.frame(data2)
only1=subset(data3,Freq==1)
#about half of genera only have 1 bacterium in them

repeats=intersect(c(stratified1,stratified2,stratified3,stratified4),stratified5)
repeat1=intersect(stratified4,stratified5)
(6481/11302)*100
(3511/5000)*100

#overlap between any two= about 30%
#in total dataset, 57.3 genera only have one sample
