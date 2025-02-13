rm(list = ls())
library(stringr)
library(dplyr)


#my_file1 = "snp-S178_D17.csv"
#my_file2 = "snp-S178J24.csv"
my_file1 = "pileup_all-75.csv"
my_file2 = "pileup_all-78.csv"
reference = "reference_data.txt"

#Q341, 150 SNPs


#preparation data
my_data = readr::read_csv(my_file1)
my_data1 = as.data.frame(my_data)
my_data = readr::read_csv(my_file2)
my_data2 = as.data.frame(my_data)
reference = readr::read_tsv(reference, col_names = F)
reference = as.data.frame(reference)


#combine data
combine_data = cbind(my_data1, my_data2$`Reads Count`, my_data2$`A Count`,
                     my_data2$`T Count`, my_data2$`C Count`, my_data2$`G Count`,
                     my_data2$`Dominant Base`, my_data2$`Dominant Proportion (%)`, 
                     str_to_upper(reference$X3))

#>=10X
combine_data = combine_data[combine_data$`Reads Count`>=10,]
combine_data = combine_data[combine_data$`my_data2$\`Reads Count\``>=10,]

#retain only sites different to references
combine_data = combine_data[combine_data$`Dominant Base` != combine_data$`str_to_upper(reference$X3)` | combine_data$`my_data2$\`Dominant Base\`` != combine_data$`str_to_upper(reference$X3)`,]


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


#keep only different sites
combine_data <- combine_data[combine_data[, 8] != combine_data[, 15], ]

combine_data$totalS1 = 0
combine_data$totalS2 = 0
combine_data$total_freq1 = 0
combine_data$total_freq2 = 0

#remove indels, associated with more than 150% in freq
for(i in 1:nrow(combine_data)){
  combine_data$totalS1[i] = sum(combine_data$`A Count`[i], combine_data$`T Count`[i],
                                combine_data$`C Count`[i], combine_data$`G Count`[i])
  combine_data$totalS2[i] = sum(combine_data$`my_data2$\`A Count\``[i], combine_data$`my_data2$\`T Count\``[i],
                                combine_data$`my_data2$\`C Count\``[i], combine_data$`my_data2$\`G Count\``[i])
  combine_data$total_freq1[i] = combine_data$totalS1[i] / combine_data$`Reads Count`[i] * 100
  combine_data$total_freq2[i] = combine_data$totalS2[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
}

combine_data = combine_data[combine_data$total_freq1 <= 150,]
combine_data = combine_data[combine_data$total_freq2 <= 150,]


#consider only very different variants between the two samples
#at least 50% of difference in one allele
combine_data$to_keep = 0
for(i in 1:nrow(combine_data)){
  two_step = 0
  sample1_allele = combine_data$`Dominant Base`[i]
  sample2_allele = combine_data$`my_data2$\`Dominant Base\``[i]
  sample1_allele_freq = 0
  sample2_allele_freq = 0
  
  sample1_allele_freq = combine_data$`Dominant Proportion (%)`[i]
  if(sample1_allele == "A"){
    sample2_allele_freq = combine_data$`my_data2$\`A Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }else if(sample1_allele == "T"){
    sample2_allele_freq = combine_data$`my_data2$\`T Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }else if(sample1_allele == "C"){
    sample2_allele_freq = combine_data$`my_data2$\`C Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }else{
    sample2_allele_freq = combine_data$`my_data2$\`G Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
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
  sample1_allele = combine_data$`Dominant Base`[i]
  sample2_allele = combine_data$`my_data2$\`Dominant Base\``[i]
  sample1_allele_freq = 0
  sample2_allele_freq = 0
  sample2_allele_freq = combine_data$`my_data2$\`Dominant Proportion (%)\``[i]
  if(sample2_allele == "A"){
    sample1_allele_freq = combine_data$`A Count`[i] / combine_data$`Reads Count`[i] * 100
  }else if(sample2_allele == "T"){
    sample1_allele_freq = combine_data$`T Count`[i] / combine_data$`Reads Count`[i] * 100
  }else if(sample2_allele == "C"){
    sample1_allele_freq = combine_data$`C Count`[i] / combine_data$`Reads Count`[i] * 100
  }else{
    sample1_allele_freq = combine_data$`G Count`[i] / combine_data$`Reads Count`[i] * 100
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
  sample1_allele = combine_data$`Dominant Base`[i]
  sample2_allele = combine_data$`my_data2$\`Dominant Base\``[i]
  sample1_allele_freq = 0
  sample2_allele_freq = 0
  
  sample1_allele_freq = combine_data$`Dominant Proportion (%)`[i]
  if(sample1_allele == "A"){
    sample2_allele_freq = combine_data$`my_data2$\`A Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }else if(sample1_allele == "T"){
    sample2_allele_freq = combine_data$`my_data2$\`T Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }else if(sample1_allele == "C"){
    sample2_allele_freq = combine_data$`my_data2$\`C Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }else{
    sample2_allele_freq = combine_data$`my_data2$\`G Count\``[i] / combine_data$`my_data2$\`Reads Count\``[i] * 100
  }
  
  if(sample2_allele_freq <= 10){
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
  sample1_allele = combine_data$`Dominant Base`[i]
  sample2_allele = combine_data$`my_data2$\`Dominant Base\``[i]
  sample1_allele_freq = 0
  sample2_allele_freq = 0
  
  sample2_allele_freq = combine_data$`my_data2$\`Dominant Proportion (%)\``[i]
  if(sample2_allele == "A"){
    sample1_allele_freq = combine_data$`A Count`[i] / combine_data$`Reads Count`[i] * 100
  }else if(sample2_allele == "T"){
    sample1_allele_freq = combine_data$`T Count`[i] / combine_data$`Reads Count`[i] * 100
  }else if(sample2_allele == "C"){
    sample1_allele_freq = combine_data$`C Count`[i] / combine_data$`Reads Count`[i] * 100
  }else{
    sample1_allele_freq = combine_data$`G Count`[i] / combine_data$`Reads Count`[i] * 100
  }
  
  if(sample1_allele_freq <= 10){
    two_step = two_step + 1
  }
  
  
  if(two_step >= 1){
    combine_data$to_keep[i] = 1
  }
  
}

combine_data = combine_data[combine_data$to_keep == 1,]

write.table(combine_data, "result-Q341.txt", quote = F, sep = "\t", row.names = F,
            col.names = T)




