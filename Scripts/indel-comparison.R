rm(list = ls())
library(stringr)
library(dplyr)
setwd("/Users/adm-loc/Documents/test-benoit/F75-F78/")

my_file1 = "indels_F75.csv"
my_file2 = "indels_F78.csv"
reference = "reference_data.txt"

#preparation data
my_data = readr::read_csv(my_file1)
my_data1 = as.data.frame(my_data)
my_data = readr::read_csv(my_file2)
my_data2 = as.data.frame(my_data)
reference = readr::read_tsv(reference, col_names = F)
reference = as.data.frame(reference)

#combine data
combine_data = cbind(my_data1, my_data2$`Reads Count`, my_data2$`Plus Count`,
                     my_data2$`Minus Count`, my_data2$`Dominant Symbol`, my_data2$`Dominant Proportion (%)`, 
                     str_to_upper(reference$X3))

#>=10X
combine_data = combine_data[combine_data$`Reads Count`>=10,]
combine_data = combine_data[combine_data$`my_data2$\`Reads Count\``>=10,]



#keep only different sites
combine_data <- combine_data[combine_data[, 6] != combine_data[, 11], ]


#retain only prevalent indels
combine_data$to_keep = 0
for (i in 1:nrow(combine_data)){
  if(combine_data$`Dominant Proportion (%)`[i]>= 50){
    combine_data$to_keep[i] = 1
  }
  if(combine_data$`my_data2$\`Dominant Proportion (%)\``[i]>= 50){
    combine_data$to_keep[i] = 1
  }
}

combine_data = combine_data[combine_data$to_keep == 1,]

#focusing on the core genome
combine_data$to_keep = 0
for(i in 1:nrow(combine_data)){
  if(combine_data$Chromosome[i] == "Pf3D7_01_v3"){
    if(combine_data$Position[i] >= 91183 && combine_data$Position[i] <= 591736){
      combine_data$to_keep[i] = 1
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_02_v3"){
    if(combine_data$Position[i] >= 67438 && combine_data$Position[i] <= 880726){
      combine_data$to_keep[i] = 1
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_03_v3"){
    if(combine_data$Position[i] >= 63099 && combine_data$Position[i] <= 1016045){
      combine_data$to_keep[i] = 1
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_04_v3"){
    if(combine_data$Position[i] >= 86540 && combine_data$Position[i] <= 1145399){
      combine_data$to_keep[i] = 1
      if(combine_data$Position[i] >= 546226 && combine_data$Position[i] <= 601217){
        combine_data$to_keep[i] = 0
      }
      if(combine_data$Position[i] >= 936358 && combine_data$Position[i] <= 979881){
        combine_data$to_keep[i] = 0
      }
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_05_v3"){
    if(combine_data$Position[i] >= 37149 && combine_data$Position[i] <= 1333478){
      combine_data$to_keep[i] = 1
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_06_v3"){
    if(combine_data$Position[i] >= 39483 && combine_data$Position[i] <= 1317743){
      combine_data$to_keep[i] = 1
      if(combine_data$Position[i] >= 724169 && combine_data$Position[i] <= 742098){
        combine_data$to_keep[i] = 0
      }
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_07_v3"){
    if(combine_data$Position[i] >= 76561 && combine_data$Position[i] <= 1386108){
      combine_data$to_keep[i] = 1
      if(combine_data$Position[i] >= 512373 && combine_data$Position[i] <= 591896){
        combine_data$to_keep[i] = 0
      }
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_08_v3"){
    if(combine_data$Position[i] >= 57343 && combine_data$Position[i] <= 1391320){
      combine_data$to_keep[i] = 1
      if(combine_data$Position[i] >= 430803 && combine_data$Position[i] <= 466095){
        combine_data$to_keep[i] = 0
      }
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_09_v3"){
    if(combine_data$Position[i] >= 79027 && combine_data$Position[i] <= 1448081){
      combine_data$to_keep[i] = 1
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_10_v3"){
    if(combine_data$Position[i] >= 49338 && combine_data$Position[i] <= 1600156){
      combine_data$to_keep[i] = 1
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_11_v3"){
    if(combine_data$Position[i] >= 77654 && combine_data$Position[i] <= 2003815){
      combine_data$to_keep[i] = 1
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_12_v3"){
    if(combine_data$Position[i] >= 49282 && combine_data$Position[i] <= 2175030){
      combine_data$to_keep[i] = 1
      if(combine_data$Position[i] >= 768014 && combine_data$Position[i] <= 775380){
        combine_data$to_keep[i] = 0
      }
      if(combine_data$Position[i] >= 1695106 && combine_data$Position[i] <= 1741581){
        combine_data$to_keep[i] = 0
      }
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_13_v3"){
    if(combine_data$Position[i] >= 72312 && combine_data$Position[i] <= 2857354){
      combine_data$to_keep[i] = 1
    }
  }else if(combine_data$Chromosome[i] == "Pf3D7_14_v3"){
    if(combine_data$Position[i] >= 32983 && combine_data$Position[i] <= 3259777){
      combine_data$to_keep[i] = 1
    }
  }
}

combine_data = combine_data[combine_data$to_keep == 1,]




combine_data$totalS1 = 0
combine_data$totalS2 = 0
combine_data$total_freq1 = 0
combine_data$total_freq2 = 0


#consider only very different variants between the two samples
#at least 70% of difference in one allele
combine_data$to_keep = 0
for(i in 1:nrow(combine_data)){
  two_step = 0
  sample1_allele = combine_data$`Dominant Symbol`[i]
  sample2_allele = combine_data$`my_data2$\`Dominant Symbol\``[i]
  
  sample1_allele_freq = combine_data$`Dominant Proportion (%)`[i]
  if(sample1_allele == "+"){
    sample2_allele_freq = combine_data$`my_data2$\`Plus Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }else if(sample1_allele == "-"){
    sample2_allele_freq = combine_data$`my_data2$\`Minus Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }else{
    sample2_allele_freq = 0
  }
  
  if(abs(sample1_allele_freq - sample2_allele_freq) >= 50){
    two_step = two_step + 1
  }
  
  if(two_step >= 1){
    combine_data$to_keep[i] = 1
  }
} 

combine_data = combine_data[combine_data$to_keep == 1,]


combine_data$to_keep = 0
for(i in 1:nrow(combine_data)){
  two_step = 0
  sample1_allele = combine_data$`Dominant Symbol`[i]
  sample2_allele = combine_data$`my_data2$\`Dominant Symbol\``[i]
  sample1_allele_freq = 0
  sample2_allele_freq = 0
  
  
  sample2_allele_freq = combine_data$`my_data2$\`Dominant Proportion (%)\``[i]
  if(sample2_allele == "+"){
    sample1_allele_freq = combine_data$`Plus Count`[i] / combine_data$`Reads Count`[i] * 100
  }else if(sample2_allele == "-"){
    sample1_allele_freq = combine_data$`Minus Count`[i] / combine_data$`Reads Count`[i] * 100
  }else{
    sample1_allele_freq = 0
  }
  
  if(abs(sample1_allele_freq - sample2_allele_freq) >= 50){
    two_step = two_step + 1
  }
  
  
  if(two_step >= 1){
    combine_data$to_keep[i] = 1
  }
}


combine_data = combine_data[combine_data$to_keep == 1,]

#consider only variants in second round not present in first round and vice versa
#very stringent
combine_data$to_keep = 0
for(i in 1:nrow(combine_data)){
  two_step = 0
  sample1_allele = combine_data$`Dominant Symbol`[i]
  sample2_allele = combine_data$`my_data2$\`Dominant Symbol\``[i]
  
  sample1_allele_freq = combine_data$`Dominant Proportion (%)`[i]
  if(sample1_allele == "+"){
    sample2_allele_freq = combine_data$`my_data2$\`Plus Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }else if(sample1_allele == "-"){
    sample2_allele_freq = combine_data$`my_data2$\`Minus Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }
  
  if(sample2_allele_freq <= 10){
    two_step = two_step + 1
  }
  
  if(combine_data$to_keep[i] == 0){
    
    sample2_allele_freq = combine_data$`my_data2$\`Dominant Proportion (%)\``[i]
    if(sample2_allele == "+"){
      sample1_allele_freq = combine_data$`Plus Count`[i] / combine_data$`Reads Count`[i] * 100
    }else if(sample2_allele == "-"){
      sample1_allele_freq = combine_data$`Minus Count`[i] / combine_data$`Reads Count`[i] * 100
    }
    
    if(sample1_allele_freq <= 10){
      two_step = two_step + 1
    }
  }
  
  if(two_step >= 1){
    combine_data$to_keep[i] = 1
  }
  
}

combine_data = combine_data[combine_data$to_keep == 1,]



