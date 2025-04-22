##HTGCF benchmark parameter script

#Load relevant packages
library(readxl)
library(stringr)
library(tidyr)
library(tidyverse)
library(funtimes)
library(aricode)
library(ggplot2)

#import curated dataset, fix formatting and save what is needed
curated=read_excel("41589_2019_400_MOESM31_ESM.xlsx")
sep_curated=curated %>%
  separate_rows(`MIBiG accession`, sep = ",")
sep_curated$`MIBiG accession`=str_trim(sep_curated$`MIBiG accession`, side=c("both"))
sep_curated=data.frame(sep_curated$`MIBiG accession`,sep_curated$Group)
colnames(sep_curated)=c("mibig","group")

x=seq(0.01,1.0,0.01)
distance=c(0,x)
x[10]=paste0(x[10],"0")
x[20]=paste0(x[20],"0")
x[30]=paste0(x[30],"0")
x[40]=paste0(x[40],"0")
x[50]=paste0(x[50],"0")
x[60]=paste0(x[60],"0")
x[70]=paste0(x[70],"0")
x[80]=paste0(x[80],"0")
x[90]=paste0(x[90],"0")
x[100]=paste0(x[100],".00")
#average_purity=0
average_GCFs=0
#average_ARI=0
n=length(x)
for (i in 1:n){
  benchmark=read.delim(paste0("average_results/benchmark_",x[i],".tsv"))
  benchmark$cluster_id=strtrim(benchmark$cluster_id,10)
  benchmark=benchmark %>%
    group_by(gcf_id) %>%
    filter(n() >=2)
  htgcf=data.frame(benchmark$cluster_id,benchmark$gcf_id)
  colnames(htgcf)=c("mibig","GCF")
  BGCs=c("BGC0000015", "BGC0000095","BGC0000322",
         "BGC0000462","BGC0000629","BGC0000875",
         "BGC0000945", "BGC0001081", "BGC0001147",
         "BGC0001208", "BGC0001232", "BGC0001347")
  b=length(BGCs)
  temp=0
  for (y in 1:b){
    temp=c(temp,c(which(htgcf$mibig==BGCs[y])))
  }
  htgcf_f=htgcf[-temp,]
  merged_data=merge(sep_curated,htgcf_f,by="mibig")
 # ARI_htgcf=AMI(merged_data$group,merged_data$GCF)
 # average_ARI=c(average_ARI,ARI_htgcf)
  merged_data$GCF=as.factor(merged_data$GCF)
  level=nlevels(merged_data$GCF)
  average_GCFs=c(average_GCFs,level)
  #purity_htgcf=purity(merged_data$group,merged_data$GCF)
  #average_purity=c(average_purity,purity_htgcf[["pur"]])
}

#complete_purity=0
complete_GCFs=0
#complete_ARI=0
n=length(x)
for (i in 1:n){
  benchmark=read.delim(paste0("complete_results/benchmark_", x[i],".tsv"))
  benchmark$cluster_id=strtrim(benchmark$cluster_id,10)
  benchmark=benchmark %>%
    group_by(gcf_id) %>%
    filter(n() >=2)
  htgcf=data.frame(benchmark$cluster_id,benchmark$gcf_id)
  colnames(htgcf)=c("mibig","GCF")
  BGCs=c("BGC0000015", "BGC0000095","BGC0000322",
         "BGC0000462","BGC0000629","BGC0000875",
         "BGC0000945", "BGC0001081", "BGC0001147",
         "BGC0001208", "BGC0001232", "BGC0001347")
  b=length(BGCs)
  temp=0
  for (y in 1:b){
    temp=c(temp,c(which(htgcf$mibig==BGCs[y])))
  }
  htgcf_f=htgcf[-temp,]
  merged_data=merge(sep_curated,htgcf_f,by="mibig")
 # ARI_htgcf=AMI(merged_data$group,merged_data$GCF)
  #complete_ARI=c(complete_ARI,ARI_htgcf)
  merged_data$GCF=as.factor(merged_data$GCF)
  level=nlevels(merged_data$GCF)
  complete_GCFs=c(complete_GCFs,level)
  #purity_htgcf=purity(merged_data$group,merged_data$GCF)
  #complete_purity=c(complete_purity,purity_htgcf[["pur"]])
}

#single_purity=0
single_GCFs=0
#single_ARI=0
n=length(x)
for (i in 1:n){
  benchmark=read.delim(paste0("single_results/benchmark_", x[i],".tsv"))
  benchmark$cluster_id=strtrim(benchmark$cluster_id,10)
  benchmark=benchmark %>%
    group_by(gcf_id) %>%
    filter(n() >=2)
  htgcf=data.frame(benchmark$cluster_id,benchmark$gcf_id)
  colnames(htgcf)=c("mibig","GCF")
  BGCs=c("BGC0000015", "BGC0000095","BGC0000322",
         "BGC0000462","BGC0000629","BGC0000875",
         "BGC0000945", "BGC0001081", "BGC0001147",
         "BGC0001208", "BGC0001232", "BGC0001347")
  b=length(BGCs)
  temp=0
  for (y in 1:b){
    temp=c(temp,c(which(htgcf$mibig==BGCs[y])))
  }
  htgcf_f=htgcf[-temp,]
  merged_data=merge(sep_curated,htgcf_f,by="mibig")
 # ARI_htgcf=AMI(merged_data$group,merged_data$GCF)
  #single_ARI=c(single_ARI,ARI_htgcf)
  merged_data$GCF=as.factor(merged_data$GCF)
  level=nlevels(merged_data$GCF)
  single_GCFs=c(single_GCFs,level)
  #purity_htgcf=purity(merged_data$group,merged_data$GCF)
 # single_purity=c(single_purity,purity_htgcf[["pur"]])
}

#weighted_purity=0
weighted_GCFs=0
#weighted_ARI=0
n=length(x)
for (i in 1:n){
  benchmark=read.delim(paste0("weighted_results/benchmark_", x[i],".tsv"))
  benchmark$cluster_id=strtrim(benchmark$cluster_id,10)
  benchmark=benchmark %>%
    group_by(gcf_id) %>%
    filter(n() >=2)
  htgcf=data.frame(benchmark$cluster_id,benchmark$gcf_id)
  colnames(htgcf)=c("mibig","GCF")
  BGCs=c("BGC0000015", "BGC0000095","BGC0000322",
         "BGC0000462","BGC0000629","BGC0000875",
         "BGC0000945", "BGC0001081", "BGC0001147",
         "BGC0001208", "BGC0001232", "BGC0001347")
  b=length(BGCs)
  temp=0
  for (y in 1:b){
    temp=c(temp,c(which(htgcf$mibig==BGCs[y])))
  }
  htgcf_f=htgcf[-temp,]
  merged_data=merge(sep_curated,htgcf_f,by="mibig")
 # ARI_htgcf=AMI(merged_data$group,merged_data$GCF)
  #weighted_ARI=c(weighted_ARI,ARI_htgcf)
  merged_data$GCF=as.factor(merged_data$GCF)
 level=nlevels(merged_data$GCF)
  weighted_GCFs=c(weighted_GCFs,level)
  #purity_htgcf=purity(merged_data$group,merged_data$GCF)
  #weighted_purity=c(weighted_purity,purity_htgcf[["pur"]])
}

#centroid_purity=0
centroid_GCFs=0
#centroid_ARI=0
n=length(x)
for (i in 1:n){
  benchmark=read.delim(paste0("centroid_results/benchmark_", x[i],".tsv"))
  benchmark$cluster_id=strtrim(benchmark$cluster_id,10)
  benchmark=benchmark %>%
    group_by(gcf_id) %>%
    filter(n() >=2)
  htgcf=data.frame(benchmark$cluster_id,benchmark$gcf_id)
  colnames(htgcf)=c("mibig","GCF")
  BGCs=c("BGC0000015", "BGC0000095","BGC0000322",
         "BGC0000462","BGC0000629","BGC0000875",
         "BGC0000945", "BGC0001081", "BGC0001147",
         "BGC0001208", "BGC0001232", "BGC0001347")
  b=length(BGCs)
  temp=0
  for (y in 1:b){
    temp=c(temp,c(which(htgcf$mibig==BGCs[y])))
  }
  htgcf_f=htgcf[-temp,]
  merged_data=merge(sep_curated,htgcf_f,by="mibig")
 # ARI_htgcf=AMI(merged_data$group,merged_data$GCF)
  #centroid_ARI=c(centroid_ARI,ARI_htgcf)
 merged_data$GCF=as.factor(merged_data$GCF)
  level=nlevels(merged_data$GCF)
  centroid_GCFs=c(centroid_GCFs,level)
 # purity_htgcf=purity(merged_data$group,merged_data$GCF)
 # centroid_purity=c(centroid_purity,purity_htgcf[["pur"]])
}

#median_purity=0
median_GCFs=0
#median_ARI=0
n=length(x)
for (i in 1:n){
  benchmark=read.delim(paste0("median_results/benchmark_", x[i],".tsv"))
  benchmark$cluster_id=strtrim(benchmark$cluster_id,10)
  benchmark=benchmark %>%
    group_by(gcf_id) %>%
    filter(n() >=2)
  htgcf=data.frame(benchmark$cluster_id,benchmark$gcf_id)
  colnames(htgcf)=c("mibig","GCF")
  BGCs=c("BGC0000015", "BGC0000095","BGC0000322",
         "BGC0000462","BGC0000629","BGC0000875",
         "BGC0000945", "BGC0001081", "BGC0001147",
         "BGC0001208", "BGC0001232", "BGC0001347")
  b=length(BGCs)
  temp=0
  for (y in 1:b){
    temp=c(temp,c(which(htgcf$mibig==BGCs[y])))
  }
  htgcf_f=htgcf[-temp,]
  merged_data=merge(sep_curated,htgcf_f,by="mibig")
 # ARI_htgcf=AMI(merged_data$group,merged_data$GCF)
  #median_ARI=c(median_ARI,ARI_htgcf)
  merged_data$GCF=as.factor(merged_data$GCF)
  level=nlevels(merged_data$GCF)
  median_GCFs=c(median_GCFs,level)
 # purity_htgcf=purity(merged_data$group,merged_data$GCF)
 # median_purity=c(median_purity,purity_htgcf[["pur"]])
}

#ward_purity=0
#ward_ARI=0
ward_GCFs=0
n=length(x)
for (i in 1:n){
  benchmark=read.delim(paste0("ward_results/benchmark_", x[i],".tsv"))
  benchmark$cluster_id=strtrim(benchmark$cluster_id,10)
  benchmark=benchmark %>%
    group_by(gcf_id) %>%
    filter(n() >=2)
  htgcf=data.frame(benchmark$cluster_id,benchmark$gcf_id)
  colnames(htgcf)=c("mibig","GCF")
  BGCs=c("BGC0000015", "BGC0000095","BGC0000322",
         "BGC0000462","BGC0000629","BGC0000875",
         "BGC0000945", "BGC0001081", "BGC0001147",
         "BGC0001208", "BGC0001232", "BGC0001347")
  b=length(BGCs)
  temp=0
  for (y in 1:b){
    temp=c(temp,c(which(htgcf$mibig==BGCs[y])))
  }
  htgcf_f=htgcf[-temp,]
  merged_data=merge(sep_curated,htgcf_f,by="mibig")
 # ARI_htgcf=AMI(merged_data$group,merged_data$GCF)
  #ward_ARI=c(ward_ARI,ARI_htgcf)
  merged_data$GCF=as.factor(merged_data$GCF)
  level=nlevels(merged_data$GCF)
  ward_GCFs=c(ward_GCFs,level)
  #purity_htgcf=purity(merged_data$group,merged_data$GCF)
  #ward_purity=c(ward_purity,purity_htgcf[["pur"]])
}
#Remove 0s used to start the vectors

distance=distance[-1]
complete_GCFs=complete_GCFs[-1]
single_GCFs=single_GCFs[-1]
average_GCFs=average_GCFs[-1]
weighted_GCFs=weighted_GCFs[-1]
centroid_GCFs=centroid_GCFs[-1]
median_GCFs=median_GCFs[-1]
ward_GCFs=ward_GCFs[-1]


#Make final dataframe with purity values

distance=rep(distance,7)
group=c(rep("complete",100), rep("single",100),
        rep("average",100), rep("weighted",100), rep("centroid",100),
        rep("median",100), rep("ward",100))
#purity=c(complete_purity, single_purity, average_purity,
 #        weighted_purity, centroid_purity, median_purity, ward_purity)

GCFs=c(complete_GCFs, single_GCFs, average_GCFs,
         weighted_GCFs, centroid_GCFs, median_GCFs, ward_GCFs)
#data=data.frame(distance,group,GCFs)
#ARIs=c(complete_ARI, single_ARI, average_ARI,
 #        weighted_ARI, centroid_ARI, median_ARI, ward_ARI)
#data=data.frame(distance,group,GCFs)

#Making the plot

data=data.frame(distance,group,GCFs)
write_delim(data, "./GCFs_benchmark_data.tsv")

