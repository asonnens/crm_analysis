require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 
require(kernlab)
require(ROCR)
library(lattice)
library(plyr)
require(ROSE)
library(ggplot2)
library(reshape2)
require(caTools)

#maybe two functions
#is this an enhancer, if yes, what's its expression pattern

setwd("shadow_enhancers/machine_learning_input/October")

known_data <- read.table("Known_enhancers.tsv", header = TRUE)
indel <- read.table("known_enhancers_indel.tsv", header = TRUE)
PC <- known_data$Average_PC
plot(PC ~ known_data$Dist_TSS)
barplot(PC, col = c('blue', 'blue', 'darkgreen', 'darkgreen', 'darkgoldenrod', 'darkgoldenrod', 'coral2', 'pink', 'darkorange', 'darkorange', 'chocolate4', 'red', 'red', 'lightgreen', 'lightgreen', 'cornsilk1', 'cornsilk2', 'cornsilk2', 'cornsilk3', 'cornsilk4', 'black'), main = "Conservation of Known Enhancers", ylab = "Phastcons score")
indel_nums <- indel$indel
barplot(indel_nums, col = c('blue', 'blue', 'darkgreen', 'darkgreen', 'darkgoldenrod', 'darkgoldenrod', 'coral2', 'pink', 'darkorange', 'darkorange', 'chocolate4', 'red', 'red', 'lightgreen', 'lightgreen', 'cornsilk1', 'cornsilk2', 'cornsilk2', 'cornsilk3', 'cornsilk4', 'black'), main = "Indel Presence in Known Enhancers", ylab = "Percentage genomes with indels within 100bp window")


Mac_Dl_train_cons <- read.table("cons_train_MacArthur_Dorsal.tsv", header = TRUE)
Mac_Dl_train_lax <- read.table("lax_train_MacArthur_Dorsal.tsv", header = TRUE)
Mac_Dl_train_rand <- read.table("rand_train_MacArthur_Dorsal.tsv", header = TRUE)
Mac_Dl_train_all <- read.table("all_train_MacArthur_Dorsal.tsv", header = TRUE)
Mac_Dl_put <- read.table("Dorsal_putative_Oct1.tsv", header = TRUE)

RDl_train_cons <- read.table("cons_train_Rushlow_Dorsal.tsv", header = TRUE)
RDl_train_lax <- read.table("lax_train_Rushlow_Dorsal.tsv", header = TRUE)
RDl_train_rand <- read.table("rand_train_Rushlow_Dorsal.tsv", header = TRUE)
RDl_train_all <- read.table("all_train_Rushlow_Dorsal.tsv", header = TRUE)
RDl_put <- read.table("RDl_putative_Oct1.tsv", header = TRUE)

FSna_train_cons <- read.table("cons_train_Furlong_Snail.tsv", header = TRUE)
FSna_train_lax <- read.table("lax_train_Furlong_Snail.tsv", header = TRUE)
FSna_train_rand <- read.table("rand_train_Furlong_Snail.tsv", header = TRUE)
FSna_train_all <- read.table("all_train_Furlong_Snail.tsv", header = TRUE)
FSna_put <- read.table("Furlong_Sna_putative_Oct1.tsv", header = TRUE)
FSna_put500 <- read.table("Furlong_Sna_putative500_Oct1.tsv", header = TRUE)
FSna_prosp <- read.table("FSna_put.tsv", header = TRUE)
FSna_put500 <- rename(FSna_put500, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  
FSna_prosp <- rename(FSna_prosp, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  


FTwi_train_cons <- read.table("cons_train_Furlong_Twist.tsv", header = TRUE)
FTwi_train_lax <- read.table("lax_train_Furlong_Twist.tsv", header = TRUE)
FTwi_train_rand <- read.table("rand_train_Furlong_Twist.tsv", header = TRUE)
FTwi_train_all <- read.table("all_train_Furlong_Twist.tsv", header = TRUE)
#FTwi_put <- read.table("Furlong_Twi_putative_Oct1.tsv", header = TRUE)
#FTwi_put500 <- read.table("Furlong_Twi_putative500_Oct1.tsv", header = TRUE)

Zld_train_cons <- read.table("cons_train_Zelda180.tsv", header = TRUE)
Zld_train_lax <- read.table("lax_train_Zelda180.tsv", header = TRUE)
Zld_train_rand <- read.table("rand_train_Zelda180.tsv", header = TRUE)
Zld_train_all <- read.table("all_train_Zelda180.tsv", header = TRUE)
#Zld_put <- read.table("Furlong_Twi_putative_Oct1.tsv", header = TRUE)
#Zld_put500 <- read.table("Furlong_Twi_putative500_Oct1.tsv", header = TRUE)

express <- read.table("expression_training_Oct.tsv", header = TRUE)

express <- rename(express, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  


FSna_train_cons <- FSna_train_cons
FSna_all <- FSna_train_all
current_put <- FSna_put500



FSna_train_all <- FSna_train_all[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]
FSna_train_cons <- FSna_train_cons[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]

FSna_train_all <- rename(FSna_train_all, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  
FSna_train_cons <- rename(FSna_train_cons, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  

tab_FSna_cons <- table(FSna_train_cons$status)
tab_FSna_all <- table(FSna_train_all$status)

FTwi_train_all <- FTwi_train_all[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]
FTwi_train_cons <- FTwi_train_cons[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]

FTwi_train_all <- rename(FTwi_train_all, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  
FTwi_train_cons <- rename(FTwi_train_cons, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  

tab_FTwi_cons <- table(FTwi_train_cons$status)
tab_FTwi_all <- table(FTwi_train_all$status)

Zld_train_all <- Zld_train_all[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]
Zld_train_cons <- Zld_train_cons[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]

Zld_train_all <- rename(Zld_train_all, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  
Zld_train_cons <- rename(Zld_train_cons, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  

tab_Zld_cons <- table(Zld_train_cons$status)
tab_Zld_all <- table(Zld_train_all$status)

RDl_train_all <- RDl_train_all[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]
RDl_train_cons <- RDl_train_cons[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]

RDl_train_all <- rename(RDl_train_all, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  
RDl_train_cons <- rename(RDl_train_cons, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  

tab_RDl_cons <- table(RDl_train_cons$status)
tab_RDl_all <- table(RDl_train_all$status)

Mac_Dl_train_all <- Mac_Dl_train_all[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]
Mac_Dl_train_cons <- Mac_Dl_train_cons[c("name","status","Zelda180", "Rushlow_Dorsal", "MacArthur_Dorsal", "Furlong_Snail", "MacArthur_Snail", 
                             "Zeitlinger1_Twist", "Furlong_Twist", "MacArthur_Twist", "Eisen2013_Bicoid", "MacArthur_Bicoid", 
							 "Eisen2010_Caudal", "MacArthur_Caudal", "Eisen2013_Hunchback", "MacArthur_Hunchback", "Eisen2013_Giant", "MacArthur_Giant", 
							 "MacArthur_Hairy", "Eisen2010_Knirps", "Eisen2013_Kruppel", "KK_H3K27ac", "White_H3K27ac", "KK_H3K4me1", "White_p300")]

Mac_Dl_train_all <- rename(Mac_Dl_train_all, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  
Mac_Dl_train_cons <- rename(Mac_Dl_train_cons, c("Zelda180" = "Zelda_seq", "Rushlow_Dorsal" = "Dorsal_seq", "MacArthur_Dorsal"= "Dorsal_chip", "Furlong_Snail" = "Snail_2014_chip",
                      "MacArthur_Snail" = "Snail_2009_chip", "Zeitlinger1_Twist" = "Twist_seq", "Furlong_Twist" = "Twist_2014_chip", "MacArthur_Twist" = "Twist_2009_chip",
                      "Eisen2013_Bicoid" = "Bicoid_seq", "MacArthur_Bicoid" = "Bicoid_chip", "Eisen2010_Caudal"= "Caudal_seq", "MacArthur_Caudal" = "Caudal_chip",
                      "Eisen2013_Hunchback" = "Hunchback_seq", "MacArthur_Hunchback" = "Hunchback_chip", "Eisen2013_Giant" = "Giant_seq", "MacArthur_Giant" = "Giant_chip",
                      "MacArthur_Hairy" = "Hairy_chip", "Eisen2010_Knirps" = "Knirps_seq", "Eisen2013_Kruppel" = "Kruppel_seq", "KK_H3K27ac" = "H3K27ac_2015_seq",
                      "White_H3K27ac" = "H3K27ac_2010_seq", "KK_H3K4me1" = "H3K4me1_2015_seq", "White_p300" = "p300_2010_seq"))  

tab_Mac_Dl_cons <- table(Mac_Dl_train_cons$status)
tab_Mac_Dl_all <- table(Mac_Dl_train_all$status)

all_features <- c(3:25)
my_features <- all_features

data1 <- current_all[,3:25]
data1_scale <- scale(data1)
status <- current_all
status <- as.numeric(status$status)
data3 <- cbind(status, data1)
data3$status[data3$status == 1] <- -1
data3$status[data3$status == 2] <- 1
data_scale <- scale(data3)
q3 <- qplot(x=Var1, y=Var2, data=melt(cor(data_scale, use="p")), fill=value, geom="tile") +
   scale_fill_gradient2(limits=c(-1, 1))
q3 + theme(axis.text.x = element_text(angle = 90, hjust = 1))

q3 <- qplot(x=Var1, y=Var2, data=melt(cor(Twist_vals, use="p")), fill=value, geom="tile") +
   scale_fill_gradient2(limits=c(-1, 1))
q3 + theme(axis.text.x = element_text(angle = 90, hjust = 1))

data_data =cor(data_scale, use="p")
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}

upper_tri <- get_upper_tri(data_data)

reordered_data <- reorder_cormat(data_data)
upper_reorderd <- get_upper_tri(reordered_data)

melted_data <- melt(upper_tri)
#melted_data <- melt(upper_reorderd)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "seashell", high = "steelblue4", mid = "cadetblue3", 
   midpoint = 0.375, limit = c(0,0.75), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()

strata_2var <- function(dataset, variable_name, variable_table, variable_pos, dataframe_length = 14, second_val = 7, samples = 300){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
      size = c(samples,samples))
  my_data <- dataset[, c(1,variable_pos,start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

strata_Percent_var <- function(dataset, variable_name, variable_table, variable_pos, dataframe_length = 14, second_val = 7, samplesize){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
      #size = c(100,40))
      size = c(round(variable_table[2]* (2/3)), round(variable_table[1]*(2/3))))
  my_data <- dataset[, c(1,variable_pos,start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}



resample_default <- function(function_name, dataset, datatable, variable_name, outfilename, description="README!", reps, samplesize, stratafun){
    temp_dir <- paste("../September/", outfilename, sep = "")
    dir.create(temp_dir)
	info_file_name <- paste(temp_dir, "input_data_summary.txt", sep = "/")
	input_data_name <- paste(temp_dir, "input_data.csv", sep = "/")
	readme_file_name <- paste(temp_dir, "README.txt", sep = "/")
	write(description, file = readme_file_name)
	write.csv(dataset, file = input_data_name)
	#summary_info <- str(dataset)
	#write.table(datatable, informational_file)
	#writeLines(summary_info, informational_file)
	#close(informational_file)
	datathing <- data.frame()
	count <- 0
    for(i in 1:reps){
	    temp_file_name <- paste(outfilename, i, sep = "_")
		temp_file_name <- paste(temp_file_name, ".csv", sep = "")
		temp_file_name <- paste(temp_dir, temp_file_name, sep = "/")
        myresult <- suppressWarnings(function_name(dataset, datatable, temp_file_name, info_file_name, samplesize, stratafun, count))
			datathing[i,"rf500_acc"] <- myresult[[1]][1]
			datathing[i,"rf500_tpos"] <- myresult[[2]][1]
			datathing[i,"rf500_tneg"] <- myresult[[3]][1]
			datathing[i,"rf500_fpos"] <- myresult[[4]][1]
			datathing[i,"rf500_fneg"] <- myresult[[5]][1]
			datathing[i,"rf1000_acc"] <- myresult[[1]][2]
			datathing[i,"rf1000_tpos"] <- myresult[[2]][2]
			datathing[i,"rf1000_tneg"] <- myresult[[3]][2]
			datathing[i,"rf1000_fpos"] <- myresult[[4]][2]
			datathing[i,"rf1000_fneg"] <- myresult[[5]][2]
			datathing[i,"rf1500_acc"] <- myresult[[1]][3]
			datathing[i,"rf1500_tpos"] <- myresult[[2]][3]
			datathing[i,"rf1500_tneg"] <- myresult[[3]][3]
			datathing[i,"rf1500_fpos"] <- myresult[[4]][3]
			datathing[i,"rf1500_fneg"] <- myresult[[5]][3]
			
     }      
	return(datathing)            
}   

resample_rf500 <- function(function_name, dataset, datatable, variable_name, outfilename, description="README!", reps, samplesize, stratafun, graph){
    temp_dir <- paste("../September/", outfilename, sep = "")
    dir.create(temp_dir)
	info_file_name <- paste(temp_dir, "input_data_summary.txt", sep = "/")
	input_data_name <- paste(temp_dir, "input_data.csv", sep = "/")
	readme_file_name <- paste(temp_dir, "README.txt", sep = "/")
	write(description, file = readme_file_name)
	write.csv(dataset, file = input_data_name)
	#summary_info <- str(dataset)
	#write.table(datatable, informational_file)
	#writeLines(summary_info, informational_file)
	#close(informational_file)
	datathing <- data.frame()
	par(mfrow=c(1,1))
	count <- 0
    for(i in 1:reps){
		count <- count + 1
	    temp_file_name <- paste(outfilename, i, sep = "_")
		temp_file_name <- paste(temp_file_name, ".csv", sep = "")
		temp_file_name <- paste(temp_dir, temp_file_name, sep = "/")
        myresult <- suppressWarnings(function_name(dataset, datatable, temp_file_name, info_file_name, samplesize, stratafun = strata_Percent_var, count, graph))
			datathing[i,"rf500_acc"] <- myresult[1]
			datathing[i,"rf500_tpos"] <- myresult[2]
			datathing[i,"rf500_tneg"] <- myresult[3]
			datathing[i,"rf500_fpos"] <- myresult[4]
			datathing[i,"rf500_fneg"] <- myresult[5]
			datathing[i, "rf500_AUC"] <- myresult[6]
			datathing[i, "rf500_PRAUC"] <- myresult[7]
     }      
	return(datathing)            
}   



enhancer_class_status_scale <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize, stratafun, count, graph){
    enhancers_set <- stratafun(dataset, "status", datatable,2, 25, 22, samplesize)   
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_training_vars <- enhancers_training2[,c(my_features)]
	enhancers_data_trainingscale <- scale(enhancers_training_vars)
    enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_trainingscale)
	enhancers_training <- enhancers_trainingscale
	enhancers_test2 <- enhancers_set[[2]]
	enhancers_test_vars <- enhancers_test2[,c(my_features)]
	enhancers_data_testscale <- scale(enhancers_test_vars)
    enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_test <- enhancers_testscale 
	actual <- enhancers_test$status
	#random forest 500 trees                  
	enhancer_rf_500 <- randomForest(status ~., data = enhancers_training[,2:length(enhancers_training)],  ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$status, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	rf_500_false_pos <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[1,])
	predictions=as.vector(enhancer_rf_500$votes[,2])
	pred1=prediction(predictions,enhancers_training$status)
	perf1 <- performance(pred1, "prec", "rec")
	perf_ROC=performance(pred1,"tpr","fpr")
	perf_AUC=performance(pred1,"auc")
	AUC <- perf_AUC@y.values[[1]]
	precision <- performance(pred1, "prec")@y.values[[1]]
	precision[1] <- 0
	recall <- performance(pred1, "rec")@y.values[[1]]
	prec_reca <- trapz(recall, precision)
	if (count == 1 && graph == "recall") {plot(perf1, col = "darkgreen", main = "Precision Recall")} else if (count !=1 && graph == "recall") {plot(perf1, add = TRUE, col = "darkgreen")}
	if (count == 1 && graph == "ROC") {plot(perf_ROC, col = "darkgoldenrod", main = "ROC plot")} else if (count !=1 && graph == "ROC") {plot(perf_ROC, add = TRUE, col = "darkgoldenrod")}
	write.csv(enhancers_test, file = summaryfile)
	#close(filename)
	return_list <- list(rf_500_prediction_success, rf_500_actual_on, rf_500_true_off, rf_500_false_on, rf_500_false_off, AUC, prec_reca)
	return(return_list)             
}       




Mac_Dl_all_plot_ROC <- resample_rf500(enhancer_class_status_scale, Mac_Dl_train_all, tab_Mac_Dl_all, "MacDlAll", description = "MacDlAll", "MacDlAll", 100, 300, strata_Percent_var, "ROC")
Mac_Dl_all_plot_Recall <- resample_rf500(enhancer_class_status_scale, Mac_Dl_train_all, tab_Mac_Dl_all, "MacDlAll", description = "MacDlAll", "MacDlAll", 100, 300, strata_Percent_var, "recall")
Mac_Dl_cons_plot_ROC <- resample_rf500(enhancer_class_status_scale, Mac_Dl_train_cons, tab_Mac_Dl_cons, "MacDlcons", description = "MacDlcons", "MacDlcons", 100, 300, strata_2var, "ROC")
Mac_Dl_cons_plot_Recall <- resample_rf500(enhancer_class_status_scale, Mac_Dl_train_cons, tab_Mac_Dl_cons, "MacDlcons", description = "MacDlcons", "MacDlcons", 100, 300, strata_2var, "recall")

RDl_all_plot_ROC <- resample_rf500(enhancer_class_status_scale, RDl_train_all, tab_RDl_all, "RDlAll", description = "RDlAll", "RDlAll", 100, 100, strata_Percent_var, "ROC")
RDl_all_plot_Recall <- resample_rf500(enhancer_class_status_scale, RDl_train_all, tab_RDl_all, "RDlAll", description = "RDlAll", "RDlAll", 100, 100, strata_Percent_var, "recall")
RDl_cons_plot_ROC <- resample_rf500(enhancer_class_status_scale, RDl_train_cons, tab_RDl_cons, "RDlcons", description = "RDlcons", "RDlcons", 100, 100, strata_2var, "ROC")
RDl_cons_plot_Recall <- resample_rf500(enhancer_class_status_scale, RDl_train_cons, tab_RDl_cons, "RDlcons", description = "RDlcons", "RDlcons", 100, 100, strata_2var, "recall")

Zld_all_plot_ROC <- resample_rf500(enhancer_class_status_scale, Zld_train_all, tab_Zld_all, "ZldAll", description = "ZldAll", "ZldAll", 100, 70, strata_Percent_var, "ROC")
Zld_all_plot_Recall <- resample_rf500(enhancer_class_status_scale, Zld_train_all, tab_Zld_all, "ZldAll", description = "ZldAll", "ZldAll", 100, 70, strata_Percent_var, "recall")
Zld_cons_plot_ROC <- resample_rf500(enhancer_class_status_scale, Zld_train_cons, tab_Zld_cons, "Zldcons", description = "Zldcons", "RDlcons", 100, 70, strata_2var, "ROC")
Zld_cons_plot_Recall <- resample_rf500(enhancer_class_status_scale, Zld_train_cons, tab_Zld_cons, "Zldcons", description = "Zldcons", "Zldcons", 100, 70, strata_2var, "recall")

FTwi_all_plot_ROC <- resample_rf500(enhancer_class_status_scale, FTwi_train_all, tab_FTwi_all, "FTwiAll", description = "FTwiAll", "FTwiAll", 100, 250, strata_Percent_var, "ROC")
FTwi_all_plot_Recall <- resample_rf500(enhancer_class_status_scale, FTwi_train_all, tab_FTwi_all, "FTwiAll", description = "FTwiAll", "FTwiAll", 100, 250, strata_Percent_var, "recall")
FTwi_cons_plot_ROC <- resample_rf500(enhancer_class_status_scale, FTwi_train_cons, tab_FTwi_cons, "FTwicons", description = "FTwicons", "FTwicons", 100, 250, strata_2var, "ROC")
FTwi_cons_plot_Recall <- resample_rf500(enhancer_class_status_scale, FTwi_train_cons, tab_FTwi_cons, "FTwicons", description = "FTwicons", "FTwicons", 100, 250, strata_2var, "recall")

FSna_all_plot_ROC <- resample_rf500(enhancer_class_status_scale, FSna_train_all, tab_FSna_all, "FSnaAll", description = "FSnaAll", "FSnaAll", 100, 250, strata_Percent_var, "ROC")
FSna_all_plot_Recall <- resample_rf500(enhancer_class_status_scale, FSna_train_all, tab_FSna_all, "FSnaAll", description = "FSnaAll", "FSnaAll", 100, 250, strata_Percent_var, "recall")
FSna_cons_plot_ROC <- resample_rf500(enhancer_class_status_scale, FSna_train_cons, tab_FSna_cons, "FSnacons", description = "FSnacons", "FSnacons", 100, 250, strata_2var, "ROC")
FSna_cons_plot_Recall <- resample_rf500(enhancer_class_status_scale, FSna_train_cons, tab_FSna_cons, "FSnacons", description = "FSnacons", "FSnacons", 100, 250, strata_2var, "recall")

Mac_Dl_all_noPlot <- resample_rf500(enhancer_class_status_scale, Mac_Dl_train_all, tab_Mac_Dl_all, "MacDlAll", description = "MacDlAll", "MacDlAll", 100, 300, strata_Percent_var, "none")
Mac_Dl_all_means <- apply(Mac_Dl_all_noPlot, MARGIN = 2, FUN = mean)
Mac_Dl_all_sd <- apply(Mac_Dl_all_noPlot, MARGIN=2, FUN=sd)
Mac_Dl_cons_noPlot <- resample_rf500(enhancer_class_status_scale, Mac_Dl_train_cons, tab_Mac_Dl_cons, "MacDlcons", description = "MacDlcons", "MacDlcons", 100, 300, strata_2var, "none")
Mac_Dl_cons_means <- apply(Mac_Dl_cons_noPlot, MARGIN = 2, FUN = mean)
Mac_Dl_cons_sd <- apply(Mac_Dl_cons_noPlot, MARGIN=2, FUN=sd)
FSna_all_noPlot <- resample_rf500(enhancer_class_status_scale, FSna_train_all, tab_FSna_all, "FSnaAll", description = "FSnaAll", "FSnaAll", 100, 250, strata_Percent_var, "none")
FSna_all_means <- apply(FSna_all_noPlot, MARGIN = 2, FUN = mean)
FSna_all_sd <- apply(FSna_all_noPlot, MARGIN=2, FUN=sd)
FSna_cons_noPlot <- resample_rf500(enhancer_class_status_scale, FSna_train_cons, tab_FSna_cons, "FSnacons", description = "FSnacons", "FSnacons", 100, 250, strata_2var, "none")
FSna_cons_means <- apply(FSna_cons_noPlot, MARGIN = 2, FUN = mean)
FSna_cons_sd <- apply(FSna_cons_noPlot, MARGIN=2, FUN=sd)
FTwi_all_noPlot <- resample_rf500(enhancer_class_status_scale, FTwi_train_all, tab_FTwi_all, "FTwiAll", description = "FTwiAll", "FTwiAll", 100, 250, strata_Percent_var, "none")
FTwi_all_means <- apply(FTwi_all_noPlot, MARGIN = 2, FUN = mean)
FTwi_all_sd <- apply(FTwi_all_noPlot, MARGIN=2, FUN=sd)
FTwi_cons_noPlot <- resample_rf500(enhancer_class_status_scale, FTwi_train_cons, tab_FTwi_cons, "FTwicons", description = "FTwicons", "FTwicons", 100, 250, strata_2var, "none")
FTwi_cons_means <- apply(FTwi_cons_noPlot, MARGIN = 2, FUN = mean)
FTwi_cons_sd <- apply(FTwi_cons_noPlot, MARGIN=2, FUN=sd)
Zld_all_noPlot <- resample_rf500(enhancer_class_status_scale, Zld_train_all, tab_Zld_all, "ZldAll", description = "ZldAll", "ZldAll", 100, 70, strata_Percent_var, "none")
Zld_all_means <- apply(Zld_all_noPlot, MARGIN = 2, FUN = mean)
Zld_all_sd <- apply(Zld_all_noPlot, MARGIN=2, FUN=sd)
Zld_cons_noPlot <- resample_rf500(enhancer_class_status_scale, Zld_train_cons, tab_Zld_cons, "Zldcons", description = "Zldcons", "RDlcons", 100, 70, strata_2var, "none")
Zld_cons_means <- apply(Zld_cons_noPlot, MARGIN = 2, FUN = mean)
Zld_cons_sd <- apply(Zld_cons_noPlot, MARGIN=2, FUN=sd)
RDl_all_noPlot <- resample_rf500(enhancer_class_status_scale, RDl_train_all, tab_RDl_all, "RDlAll", description = "RDlAll", "RDlAll", 100, 100, strata_Percent_var, "none")
RDl_all_means <- apply(RDl_all_noPlot, MARGIN = 2, FUN = mean)
RDl_all_sd <- apply(RDl_all_noPlot, MARGIN=2, FUN=sd)
RDl_cons_noPlot <- resample_rf500(enhancer_class_status_scale, RDl_train_cons, tab_RDl_cons, "RDlcons", description = "RDlcons", "RDlcons", 100, 100, strata_2var, "none")
RDl_cons_means <- apply(RDl_cons_noPlot, MARGIN = 2, FUN = mean)
RDl_cons_sd <- apply(RDl_cons_noPlot, MARGIN=2, FUN=sd)

AUC_means <- as.vector(cbind(Mac_Dl_all_means[6], Mac_Dl_cons_means[6], FSna_all_means[6], FSna_cons_means[6], FTwi_all_means[6], FTwi_cons_means[6], Zld_all_means[6], Zld_cons_means[6], RDl_all_means[6], RDl_cons_means[6]))
AUC_sd <- as.vector(cbind(Mac_Dl_all_sd[6], Mac_Dl_cons_sd[6], FSna_all_sd[6], FSna_cons_sd[6], FTwi_all_sd[6], FTwi_cons_sd[6], Zld_all_sd[6], Zld_cons_sd[6], RDl_all_sd[6], RDl_cons_sd[6]))

PRAUC_means <- as.vector(cbind(Mac_Dl_all_means[7], Mac_Dl_cons_means[7], FSna_all_means[7], FSna_cons_means[7], FTwi_all_means[7], FTwi_cons_means[6], Zld_all_means[6], Zld_cons_means[6], RDl_all_means[6], RDl_cons_means[6]))
PRAUC_sd <- as.vector(cbind(Mac_Dl_all_sd[7], Mac_Dl_cons_sd[7], FSna_all_sd[6], FSna_cons_sd[6], FTwi_all_sd[6], FTwi_cons_sd[6], Zld_all_sd[6], Zld_cons_sd[6], RDl_all_sd[6], RDl_cons_sd[6]))

barplot(AUC_means, ylim= c(-0.1, 1), col = c("darkred", "coral3", "darkgoldenrod4", "goldenrod", "darkgreen","darkolivegreen4", "antiquewhite4", "azure4", "darkblue", "cornflowerblue" ), xaxt= "n",main = "AUC various models, on vs off")
arrows(seq(0.7,18,1.2), (AUC_means - AUC_sd), seq(.7, 18, 1.2), (AUC_means + AUC_sd), length=.05, angle=90, code=3)
text(y=-.01, x = 0.5, srt=45, "uDlchip")
text(y=-.01, x = 1.5,srt= 45, "bDlchip")
text(y=-.01, x = 2.5, srt=45, "uSnachip")
text(y=-.01, x = 3.6,srt= 45, "bSnachip")
text(y=-.01, x = 5, srt=45, "uTwichip")
text(y=-.01, x = 6.3,srt= 45, "bTwichip")
text(y=-.01, x = 7.5, srt=45, "uZldseq")
text(y=-.01, x = 8.6,srt= 45, "bZldseq")
text(y=-.01, x = 9.8, srt=45, "uDlseq")
text(y=-.01, x = 11,srt= 45, "bDlseq")

plot(PRAUC_means, ylim= c(-0.1, 1), col = c("darkred", "coral3", "darkgoldenrod4", "goldenrod", "darkgreen","darkolivegreen4", "antiquewhite4", "azure4", "darkblue", "cornflowerblue"),type = "l", xaxt= "n",main = "Precision/Recall AUC various models, on vs off")

barplot(PRAUC_means, ylim= c(-0.1, 1), col = c("darkred", "coral3", "darkgoldenrod4", "goldenrod", "darkgreen","darkolivegreen4", "antiquewhite4", "azure4", "darkblue", "cornflowerblue"),  xaxt= "n",main = "Precision/Recall AUC various models, on vs off")
arrows(seq(0.7,18,1.2), (PRAUC_means - PRAUC_sd), seq(.7, 18, 1.2), (PRAUC_means + PRAUC_sd), length=.05, angle=90, code=3)
text(y=-.01, x = 0.5, srt=45, "uDlchip")
text(y=-.01, x = 1.5,srt= 45, "bDlchip")
text(y=-.01, x = 2.5, srt=45, "uSnachip")
text(y=-.01, x = 3.6,srt= 45, "bSnachip")
text(y=-.01, x = 5, srt=45, "uTwichip")
text(y=-.01, x = 6.3,srt= 45, "bTwichip")
text(y=-.01, x = 7.5, srt=45, "uZldseq")
text(y=-.01, x = 8.6,srt= 45, "bZldseq")
text(y=-.01, x = 9.8, srt=45, "uDlseq")
text(y=-.01, x = 11,srt= 45, "bDlseq")


current_all <- FSna_train_all
current_tab <- tab_FSna_all
    enhancers_set <- strata_Percent_var(current_all, "status", tab_all,2, 25, 23, 100)   
	enhancers_training2 <- enhancers_set[[1]]

current_all_rf_500 <- randomForest(status ~., data = current_all[,2:length(current_all)], ntree = 500)
varImpPlot(current_all_rf_500, main = "current all, status", lcolor = "gray")
target_labels=as.vector(current_all$status)
target_labels[target_labels=="on"]="T"
target_labels[target_labels=="off"]="F"
predictions=as.vector(current_all_rf_500$votes[,2])
pred=prediction(predictions,current_all$status)
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]

perf_ROC=performance(pred,"tpr","fpr")
plot(perf_ROC, main="ROC plot")

perf1 <- performance(pred, "prec", "rec")
plot(perf1, add = TRUE)
perf_what <- perf1@y.values[[1]]
precis <- performance(pred, "prec")@y.values[[1]]
precisnew <- precis
precisnew[1] <- 0
reca <- performance(pred, "rec")@y.values[[1]]
recanew <- reca
plot(recanew, precisnew, type = "l")
auctest <- trapz(recanew, precisnew)
plot(perf1, main = "Precision Recall")



resample_all <- resample_rf500(enhancer_class_status_scale, current_all, tab_all, "tempwhat", description = "tempwhat", "sts", 100, 250, strata_Percent_var)

current_all_rf_1500 <- randomForest(status ~., data = current_all[,2:length(current_all)], ntree = 1500)
temp2 <- importance(current_all_rf_1500) 
varImpPlot(current_all_rf_1500, main = "current all, status", lcolor = "gray")

rename(temp_frame, c("trial1", "trial2"))

Mac_dorsal_rf_500_unbalanced <- randomForest(status ~., data = Mac_Dl_train_all[,2:length(Mac_Dl_train_all)], ntree = 500)
varImpPlot(Mac_dorsal_rf_500_unbalanced, main = "Dl 2009 unbalanced, status", lcolor = "gray")



Rushlow_dorsal_rf_500 <- randomForest(status ~., data = FSna_train_cons[,2:length(FSna_train_cons)], ntree = 500)
varImpPlot(Rushlow_dorsal_rf_500, main = "Dl 2015 balanced, status", lcolor = "gray")

pred_FSna_FSna_rf_500 <- predict(Rushlow_dorsal_rf_500, newdata = FSna_put, type = "class")
FSna_FSna_on_frame <- cbind(FSna_put[,1:2], pred_FSna_FSna_rf_500)
write.table(FSna_FSna_on_frame, "FSna_FSna_on_predict.tsv", sep="\t") 

Rushlow_dorsal_rf_500_unbalanced <- randomForest(status ~., data = FSna_train_all[,2:length(FSna_train_all)], ntree = 500)
varImpPlot(Rushlow_dorsal_rf_500_unbalanced, main = "Dl 2015 unbalanced, status", lcolor = "gray")

Mac_dorsal_rf_500_unbalanced <- randomForest(status ~., data = Mac_Dl_train_all[,2:length(Mac_Dl_train_all)], ntree = 500)
varImpPlot(Mac_dorsal_rf_500_unbalanced, main = "Dl 2009 unbalanced, status", lcolor = "gray")

prediction_MDl <- predict(Mac_dorsal_rf_500, newdata = Mac_Dl_put, type = "class")
Macdl_on_frame <- cbind(Mac_Dl_put[,1:2], prediction_MDl)

prediction_MDl_unbalanced <- predict(Mac_dorsal_rf_500_unbalanced, newdata = Mac_Dl_put, type = "class")
Macdl_on_frame_unbalanced <- cbind(Mac_Dl_put[,1:2], prediction_MDl_unbalanced)

FSna_balanced <- resample_graph(enhancer_class_status_scale, FSna_train_cons, tab_FSna_cons, "FSnaCons", "FSnaCons", "FSnaCons", 20, 60, strata_2var)
FSna_unbalanced <- resample_graph(enhancer_class_status_scale, FSna_train_all, tab_FSna_all, "FSna_all", "FSna_all", "FSna_all", 20, 60, strata_Percent_var)


FSna_on_un <- randomForest(status ~., data = current_all[,2:length(current_all)], ntree = 500)
FSna_on_bal <- randomForest(status ~., data = FSna_train_cons[,2:length(FSna_train_cons)], ntree = 500)
prediction_FSna_unbalanced <- predict(FSna_on_un, newdata = FSna_put, type = "class")
prediction_FSna_balanced <- predict(FSna_on_bal, newdata = FSna_put, type = "class")
FSna_on_un_frame <- cbind(FSna_put[,1:2], prediction_FSna_unbalanced)
FSna_on_bal_frame <- cbind(FSna_put[,1:2], prediction_FSna_balanced)
write.table(FSna_on_un_frame, "Furlong_Snail_on_unbalanced_pred.tsv", sep = "\t")
write.table(FSna_on_bal_frame, "Furlong_Snail_on_balanced_pred.tsv", sep = "\t")

prediction_FSna_unbalanced500 <- predict(FSna_on_un, newdata = FSna_put500, type = "class")
FSna_on_un_frame500 <- cbind(FSna_put500[,1:2], prediction_FSna_unbalanced500)
write.table(FSna_on_un_frame500, "Furlong_Snail_on_unbalanced_pred500.tsv", sep = "\t")

prediction_FSna_prosp <- predict(FSna_on_un, newdata = FSna_prosp, type = "class")
FSna_on_prosp <- cbind(FSna_prosp[,1:2], prediction_FSna_prosp)
write.table(FSna_on_prosp, "Furlong_Snail_on_prosp.tsv", sep = "\t")

write.table(Macdl_on_frame, "Mac_Dl_on_predict.tsv", sep="\t") 
write.table(Macdl_on_frame_unbalanced, "Mac_Dl_on_predict_unbalanced.tsv", sep = "\t")

Rdorsal_rf_500 <- randomForest(status ~., data = FSna_train_cons[,2:length(FSna_train_cons)], ntree = 500)

prediction_FSna <- predict(Rdorsal_rf_500, newdata = FSna_put, type = "class")
FSna_on_frame <- cbind(FSna_put[,1:2], prediction_FSna)

prediction_FSna_unbalanced <- predict(FSna_on_un, newdata = FSna_put, type = "class")
FSna_on_frame_unbalanced <- cbind(FSna_put[,1:2], prediction_FSna_unbalanced)


write.table(FSna_on_frame, "FSna_on_predict.tsv", sep="\t")
write.table(FSna_on_frame_unbalanced, "FSna_on_predict_unbalanced.tsv", sep="\t")


enhancer_class_Ectoderm <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize, stratafun, count, graph){
    enhancers_set <- stratafun(dataset, "Ectoderm", datatable, 8, 33, 22, samplesize) 
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_training_vars <- enhancers_training2[,c(my_features)]
	enhancers_data_trainingscale <- scale(enhancers_training_vars)
    enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_trainingscale)
	enhancers_training <- enhancers_trainingscale
	enhancers_test2 <- enhancers_set[[2]]
	enhancers_test_vars <- enhancers_test2[,c(my_features)]
	enhancers_data_testscale <- scale(enhancers_test_vars)
    enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_test <- enhancers_testscale
	enhancers_weights <- table(enhancers_test$Ectoderm)
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(Ectoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$Ectoderm ,predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	rf_500_false_pos <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[1,])
	predictions=as.vector(enhancer_rf_500$votes[,2])
	pred1=prediction(predictions,enhancers_training$Ectoderm)
	perf1 <- performance(pred1, "prec", "rec")
	perf_ROC=performance(pred1,"tpr","fpr")
	perf_AUC=performance(pred1,"auc")
	AUC <- perf_AUC@y.values[[1]]
	precision <- performance(pred1, "prec")@y.values[[1]]
	precision[1] <- 0
	recall <- performance(pred1, "rec")@y.values[[1]]
	prec_reca <- trapz(recall, precision)
	if (count == 1 && graph == "recall") {plot(perf1, col = "darkgreen", main = "Precision Recall")} else if (count !=1 && graph == "recall") {plot(perf1, add = TRUE, col = "darkgreen")}
	if (count == 1 && graph == "ROC") {plot(perf_ROC, col = "darkgoldenrod", main = "ROC plot")} else if (count !=1 && graph == "ROC") {plot(perf_ROC, add = TRUE, col = "darkgoldenrod")}
	write.csv(enhancers_test, file = summaryfile)
	#close(filename)
	return_list <- list(rf_500_prediction_success, rf_500_actual_on, rf_500_true_off, rf_500_false_on, rf_500_false_off, AUC, prec_reca)
	return(return_list)                        
}        


enhancer_class_Mesoderm <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize, stratafun,count, graph){
    enhancers_set <- stratafun(dataset, "Mesoderm", datatable, 10, 33, 22, samplesize) 
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_training_vars <- enhancers_training2[,c(my_features)]
	enhancers_data_trainingscale <- scale(enhancers_training_vars)
    enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_trainingscale)
	enhancers_training <- enhancers_trainingscale
	enhancers_test2 <- enhancers_set[[2]]
	enhancers_test_vars <- enhancers_test2[,c(my_features)]
	enhancers_data_testscale <- scale(enhancers_test_vars)
    enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_test <- enhancers_testscale
	enhancers_weights <- table(enhancers_test$Mesoderm)
	actual <- enhancers_test$Mesoderm
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(Mesoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$Mesoderm, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	predictions=as.vector(enhancer_rf_500$votes[,2])
	pred1=prediction(predictions,enhancers_training$Mesoderm)
	perf1 <- performance(pred1, "prec", "rec")
	perf_ROC=performance(pred1,"tpr","fpr")
	perf_AUC=performance(pred1,"auc")
	AUC <- perf_AUC@y.values[[1]]
	precision <- performance(pred1, "prec")@y.values[[1]]
	precision[1] <- 0
	recall <- performance(pred1, "rec")@y.values[[1]]
	prec_reca <- trapz(recall, precision)
	if (count == 1 && graph == "recall") {plot(perf1, col = "darkgreen", main = "Precision Recall")} else if (count !=1 && graph == "recall") {plot(perf1, add = TRUE, col = "darkgreen")}
	if (count == 1 && graph == "ROC") {plot(perf_ROC, col = "darkgoldenrod", main = "ROC plot")} else if (count !=1 && graph == "ROC") {plot(perf_ROC, add = TRUE, col = "darkgoldenrod")}
	write.csv(enhancers_test, file = summaryfile)
	#close(filename)
	return_list <- list(rf_500_prediction_success, rf_500_actual_on, rf_500_true_off, rf_500_false_on, rf_500_false_off, AUC, prec_reca)
	return(return_list)             
}     

enhancer_class_Endoderm <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize,stratafun, count, graph){
    enhancers_set <- stratafun(dataset, "Endoderm", datatable, 9, 33, 22, samplesize) 
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_training_vars <- enhancers_training2[,c(my_features)]
	enhancers_data_trainingscale <- scale(enhancers_training_vars)
    enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_trainingscale)
	enhancers_training <- enhancers_trainingscale
	enhancers_test2 <- enhancers_set[[2]]
	enhancers_test_vars <- enhancers_test2[,c(my_features)]
	enhancers_data_testscale <- scale(enhancers_test_vars)
    enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_test <- enhancers_testscale
	enhancers_test <- enhancers_testscale
	enhancers_weights <- table(enhancers_test$Endoderm)
	actual <- enhancers_test$Endoderm
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(Endoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$Endoderm, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	predictions=as.vector(enhancer_rf_500$votes[,2])
	pred1=prediction(predictions,enhancers_training$Endoderm)
	perf1 <- performance(pred1, "prec", "rec")
	perf_ROC=performance(pred1,"tpr","fpr")
	perf_AUC=performance(pred1,"auc")
	AUC <- perf_AUC@y.values[[1]]
	precision <- performance(pred1, "prec")@y.values[[1]]
	precision[1] <- 0
	recall <- performance(pred1, "rec")@y.values[[1]]
	prec_reca <- trapz(recall, precision)
	if (count == 1 && graph == "recall") {plot(perf1, col = "darkgreen", main = "Precision Recall")} else if (count !=1 && graph == "recall") {plot(perf1, add = TRUE, col = "darkgreen")}
	if (count == 1 && graph == "ROC") {plot(perf_ROC, col = "darkgoldenrod", main = "ROC plot")} else if (count !=1 && graph == "ROC") {plot(perf_ROC, add = TRUE, col = "darkgoldenrod")}
	write.csv(enhancers_test, file = summaryfile)
	#close(filename)
	return_list <- list(rf_500_prediction_success, rf_500_actual_on, rf_500_true_off, rf_500_false_on, rf_500_false_off, AUC, prec_reca)
	return(return_list)             
}    
enhancer_class_A <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize, stratafun, count, graph){
    enhancers_set <- stratafun(dataset, "A", datatable, 5, 33, 22, samplesize) 
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_training_vars <- enhancers_training2[,c(my_features)]
	enhancers_data_trainingscale <- scale(enhancers_training_vars)
    enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_trainingscale)
	enhancers_training <- enhancers_trainingscale
	enhancers_test2 <- enhancers_set[[2]]
	enhancers_test_vars <- enhancers_test2[,c(my_features)]
	enhancers_data_testscale <- scale(enhancers_test_vars)
    enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_test <- enhancers_testscale
	enhancers_weights <- table(enhancers_test$A)
	actual <- enhancers_test$A
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(A ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$A, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	predictions=as.vector(enhancer_rf_500$votes[,2])
	pred1=prediction(predictions,enhancers_training$A)
	perf1 <- performance(pred1, "prec", "rec")
	perf_ROC=performance(pred1,"tpr","fpr")
	perf_AUC=performance(pred1,"auc")
	AUC <- perf_AUC@y.values[[1]]
	precision <- performance(pred1, "prec")@y.values[[1]]
	precision[1] <- 0
	recall <- performance(pred1, "rec")@y.values[[1]]
	prec_reca <- trapz(recall, precision)
	if (count == 1 && graph == "recall") {plot(perf1, col = "darkgreen", main = "Precision Recall")} else if (count !=1 && graph == "recall") {plot(perf1, add = TRUE, col = "darkgreen")}
	if (count == 1 && graph == "ROC") {plot(perf_ROC, col = "darkgoldenrod", main = "ROC plot")} else if (count !=1 && graph == "ROC") {plot(perf_ROC, add = TRUE, col = "darkgoldenrod")}
	write.csv(enhancers_test, file = summaryfile)
	#close(filename)
	return_list <- list(rf_500_prediction_success, rf_500_actual_on, rf_500_true_off, rf_500_false_on, rf_500_false_off, AUC, prec_reca)
	return(return_list)             
}    
enhancer_class_C <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize, stratafun, count, graph){
    enhancers_set <- stratafun(dataset, "C", datatable, 6, 33, 22, samplesize) 
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_training_vars <- enhancers_training2[,c(my_features)]
	enhancers_data_trainingscale <- scale(enhancers_training_vars)
    enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_trainingscale)
	enhancers_training <- enhancers_trainingscale
	enhancers_test2 <- enhancers_set[[2]]
	enhancers_test_vars <- enhancers_test2[,c(my_features)]
	enhancers_data_testscale <- scale(enhancers_test_vars)
    enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_test <- enhancers_testscale
	enhancers_weights <- table(enhancers_test$C)
	actual <- enhancers_test$C
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(C ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$C, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	predictions=as.vector(enhancer_rf_500$votes[,2])
	pred1=prediction(predictions,enhancers_training$C)
	perf1 <- performance(pred1, "prec", "rec")
	perf_ROC=performance(pred1,"tpr","fpr")
	perf_AUC=performance(pred1,"auc")
	AUC <- perf_AUC@y.values[[1]]
	precision <- performance(pred1, "prec")@y.values[[1]]
	precision[1] <- 0
	recall <- performance(pred1, "rec")@y.values[[1]]
	prec_reca <- trapz(recall, precision)
	if (count == 1 && graph == "recall") {plot(perf1, col = "darkgreen", main = "Precision Recall")} else if (count !=1 && graph == "recall") {plot(perf1, add = TRUE, col = "darkgreen")}
	if (count == 1 && graph == "ROC") {plot(perf_ROC, col = "darkgoldenrod", main = "ROC plot")} else if (count !=1 && graph == "ROC") {plot(perf_ROC, add = TRUE, col = "darkgoldenrod")}
	write.csv(enhancers_test, file = summaryfile)
	#close(filename)
	return_list <- list(rf_500_prediction_success, rf_500_actual_on, rf_500_true_off, rf_500_false_on, rf_500_false_off, AUC, prec_reca)
	return(return_list)             
}    

enhancer_class_P <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize, stratafun,count, graph){
    enhancers_set <- stratafun(dataset, "P", datatable, 7, 33, 22, samplesize) 
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_training_vars <- enhancers_training2[,c(my_features)]
	enhancers_data_trainingscale <- scale(enhancers_training_vars)
    enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_trainingscale)
	enhancers_training <- enhancers_trainingscale
	enhancers_test2 <- enhancers_set[[2]]
	enhancers_test_vars <- enhancers_test2[,c(my_features)]
	enhancers_data_testscale <- scale(enhancers_test_vars)
    enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_test <- enhancers_testscale
	enhancers_weights <- table(enhancers_test$P)
	actual <- enhancers_test$P
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(P ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$P, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	predictions=as.vector(enhancer_rf_500$votes[,2])
	pred1=prediction(predictions,enhancers_training$P)
	perf1 <- performance(pred1, "prec", "rec")
	perf_ROC=performance(pred1,"tpr","fpr")
	perf_AUC=performance(pred1,"auc")
	AUC <- perf_AUC@y.values[[1]]
	precision <- performance(pred1, "prec")@y.values[[1]]
	precision[1] <- 0
	recall <- performance(pred1, "rec")@y.values[[1]]
	prec_reca <- trapz(recall, precision)
	if (count == 1 && graph == "recall") {plot(perf1, col = "darkgreen", main = "Precision Recall")} else if (count !=1 && graph == "recall") {plot(perf1, add = TRUE, col = "darkgreen")}
	if (count == 1 && graph == "ROC") {plot(perf_ROC, col = "darkgoldenrod", main = "ROC plot")} else if (count !=1 && graph == "ROC") {plot(perf_ROC, add = TRUE, col = "darkgoldenrod")}
	write.csv(enhancers_test, file = summaryfile)
	#close(filename)
	return_list <- list(rf_500_prediction_success, rf_500_actual_on, rf_500_true_off, rf_500_false_on, rf_500_false_off, AUC, prec_reca)
	return(return_list)             
}    


tab_ect <- table(express$Ectoderm)
tab_ect
tab_end <- table(express$Endoderm)
tab_end
tab_mes <- table(express$Mesoderm)
tab_mes
tab_A <- table(express$A)
tab_A
tab_C <- table(express$C)
tab_C 
tab_P <- table(express$P)
tab_P              

balanced_ect <- ovun.sample(Ectoderm~., data = express, p = 0.5, seed = 3, method = "both")$data
balanced_ect_table <- table(balanced_ect$Ectoderm)
balanced_ect_table
balanced_mes <- ovun.sample(Mesoderm~., data = express, p = 0.5, seed =3, method = "both")$data
balanced_mes_table <- table(balanced_mes$Mesoderm)
balanced_mes_table
balanced_end <- ovun.sample(Endoderm~., data = express, p = 0.5, seed =3, method = "both")$data
balanced_end_table <- table(balanced_end$Endoderm)
balanced_end_table
balanced_A <- ovun.sample(A~., data = express, p = 0.5, seed = 3, method = "both")$data
balanced_A_table <- table(balanced_A$A)
balanced_A_table
balanced_C <- ovun.sample(C~., data = express, p = 0.5, seed =3, method = "both")$data
balanced_C_table <- table(balanced_C$C)
balanced_C_table
balanced_P <- ovun.sample(P~., data = express, N = nrow(express), p = 0.5, seed = 3, method = "both")$data
balanced_P_table <- table(balanced_P$P)
balanced_P_table

Ectoderm_unbalanced_noplot <- resample_rf500(enhancer_class_Ectoderm, express, tab_ect, "tempwhat", description = "tempwhat", "sts2", 100, 80, strata_Percent_var, "none")
Ectoderm_balanced_noplot <- resample_rf500(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "tempwhat", description = "tempwhat", "sts2", 100, 100, strata_2var, "none")
Endoderm_unbalanced_noplot <- resample_rf500(enhancer_class_Endoderm, express, tab_end, "tempwhat", description = "tempwhat", "sts2", 100, 80, strata_Percent_var, "none")
Endoderm_balanced_noplot <- resample_rf500(enhancer_class_Endoderm, balanced_end, balanced_end_table, "tempwhat", description = "tempwhat", "sts2", 100, 100, strata_2var, "none")
Mesoderm_unbalanced_noplot <- resample_rf500(enhancer_class_Mesoderm, express, tab_mes, "tempwhat", description = "tempwhat", "sts2", 100, 80, strata_Percent_var, "none")
Mesoderm_balanced_noplot <- resample_rf500(enhancer_class_Mesoderm, balanced_mes, balanced_mes_table, "tempwhat", description = "tempwhat", "sts2", 100, 100, strata_2var, "none")
A_unbalanced_noplot <- resample_rf500(enhancer_class_A, express, tab_A, "tempwhat", description = "tempwhat", "sts2", 100, 80, strata_Percent_var, "none")
A_balanced_noplot <- resample_rf500(enhancer_class_A, balanced_A, balanced_A_table, "tempwhat", description = "tempwhat", "sts2", 100, 100, strata_2var, "none")
C_unbalanced_noplot <- resample_rf500(enhancer_class_C, express, tab_C, "tempwhat", description = "tempwhat", "sts2", 100, 80, strata_Percent_var, "none")
C_balanced_noplot <- resample_rf500(enhancer_class_C, balanced_C, balanced_C_table, "tempwhat", description = "tempwhat", "sts2", 100, 100, strata_2var, "none")
P_unbalanced_noplot <- resample_rf500(enhancer_class_P, express, tab_P, "tempwhat", description = "tempwhat", "sts2", 100, 80, strata_Percent_var, "none")
P_balanced_noplot <- resample_rf500(enhancer_class_P, balanced_P, balanced_P_table, "tempwhat", description = "tempwhat", "sts2", 100, 100, strata_2var, "none")

Ectoderm_unbalanced_means <- apply(Ectoderm_unbalanced_noplot, MARGIN = 2, FUN = mean)
Ect_un_sd <- apply(Ectoderm_unbalanced_noplot, MARGIN=2, FUN=sd)
Ectoderm_balanced_means <- apply(Ectoderm_balanced_noplot, MARGIN = 2, FUN = mean)
Ect_bal_sd <- apply(Ectoderm_balanced_noplot, MARGIN=2, FUN=sd)
Endoderm_unbalanced_means <- apply(Endoderm_unbalanced_noplot, MARGIN = 2, FUN = mean)
End_un_sd <- apply(Endoderm_unbalanced_noplot, MARGIN=2, FUN=sd)
Endoderm_balanced_means <- apply(Endoderm_balanced_noplot, MARGIN = 2, FUN = mean)
End_bal_sd <- apply(Endoderm_balanced_noplot, MARGIN=2, FUN=sd)
Mesoderm_unbalanced_means <- apply(Mesoderm_unbalanced_noplot, MARGIN = 2, FUN = mean)
Mes_un_sd <- apply(Mesoderm_unbalanced_noplot, MARGIN=2, FUN=sd)
Mesoderm_balanced_means <- apply(Mesoderm_balanced_noplot, MARGIN = 2, FUN = mean)
Mes_bal_sd <- apply(Mesoderm_balanced_noplot, MARGIN=2, FUN=sd)

DV_AUC_means <- as.vector(cbind(Ectoderm_unbalanced_means[6], Ectoderm_balanced_means[6], Endoderm_unbalanced_means[6], Endoderm_balanced_means[6], Mesoderm_unbalanced_means[6], Mesoderm_balanced_means[6]))
DV_AUC_sd <- as.vector(cbind(Ect_un_sd[6], Ect_bal_sd[6], End_un_sd[6], Ect_bal_sd[6], Mes_un_sd[6], Mes_bal_sd[6]))
DV_PRAUC_means <- as.vector(cbind(Ectoderm_unbalanced_means[7], Ectoderm_balanced_means[7], Endoderm_unbalanced_means[7], Endoderm_balanced_means[7], Mesoderm_unbalanced_means[7], Mesoderm_balanced_means[7]))
DV_PRAUC_sd <-  as.vector(cbind(Ect_un_sd[7], Ect_bal_sd[7], End_un_sd[7], End_bal_sd[7], Mes_un_sd[7], Mes_bal_sd[7]))

A_unbalanced_means <- apply(A_unbalanced_noplot, MARGIN = 2, FUN = mean)
A_un_sd <- apply(A_unbalanced_noplot, MARGIN=2, FUN=sd)
A_balanced_means <- apply(A_balanced_noplot, MARGIN = 2, FUN = mean)
A_bal_sd <- apply(A_balanced_noplot, MARGIN=2, FUN=sd)
C_unbalanced_means <- apply(C_unbalanced_noplot, MARGIN = 2, FUN = mean)
C_un_sd <- apply(C_unbalanced_noplot, MARGIN=2, FUN=sd)
C_balanced_means <- apply(C_balanced_noplot, MARGIN = 2, FUN = mean)
C_bal_sd <- apply(C_balanced_noplot, MARGIN=2, FUN=sd)
P_unbalanced_means <- apply(P_unbalanced_noplot, MARGIN = 2, FUN = mean)
P_un_sd <- apply(P_unbalanced_noplot, MARGIN=2, FUN=sd)
P_balanced_means <- apply(P_balanced_noplot, MARGIN = 2, FUN = mean)
P_bal_sd <- apply(P_balanced_noplot, MARGIN=2, FUN=sd)

AP_AUC_means <- as.vector(cbind(A_unbalanced_means[6], A_balanced_means[6], C_unbalanced_means[6], C_balanced_means[6], P_unbalanced_means[6], P_balanced_means[6]))
AP_AUC_sd <- as.vector(cbind(A_un_sd[6], A_bal_sd[6], C_un_sd[6], C_bal_sd[6], P_un_sd[6], P_bal_sd[6]))
AP_PRAUC_means <- as.vector(cbind(A_unbalanced_means[7], A_balanced_means[7], C_unbalanced_means[7], C_balanced_means[7], P_unbalanced_means[7], P_balanced_means[7]))
AP_PRAUC_sd <-  as.vector(cbind(A_un_sd[7], A_bal_sd[7], C_un_sd[7], C_bal_sd[7], P_un_sd[7], P_bal_sd[7]))


barplot(DV_AUC_means, ylim= c(-0.4, 1), col = c("darkred", "coral3", "darkgoldenrod4", "goldenrod", "darkgreen","darkolivegreen4", "antiquewhite4", "azure4", "darkblue", "cornflowerblue" ), xaxt= "n",main = "AUC Dorsal-Ventral")
arrows(seq(0.7,18,1.2), (DV_AUC_means - DV_AUC_sd), seq(.7, 18, 1.2), (DV_AUC_means + DV_AUC_sd), length=.05, angle=90, code=3)
text(y=-.16, x = 0.5, srt=45, "unbalanced")
text(y=-.16, x = 1.6,srt= 45, "balanced")
text(y = -.32, x = 1, "Ectoderm")
text(y=-.16, x = 2.8, srt=45, "unbalanced")
text(y=-.16, x = 4,srt= 45, "balanced")
text(y = -.32, x = 3.4, "Endoderm")
text(y=-.16, x = 5.2, srt=45, "unbalanced")
text(y=-.16, x = 6.4,srt= 45, "balanced")
text(y = -.32, x = 5.8, "Mesoderm")

barplot(DV_PRAUC_means, ylim= c(-0.4, 1), col = c("darkred", "coral3", "darkgoldenrod4", "goldenrod", "darkgreen","darkolivegreen4", "antiquewhite4", "azure4", "darkblue", "cornflowerblue" ), xaxt= "n",main = "PR-AUC Dorsal-Ventral")
arrows(seq(0.7,18,1.2), (DV_PRAUC_means - DV_PRAUC_sd), seq(.7, 18, 1.2), (DV_PRAUC_means + DV_PRAUC_sd), length=.05, angle=90, code=3)
text(y=-.16, x = 0.5, srt=45, "unbalanced")
text(y=-.16, x = 1.6,srt= 45, "balanced")
text(y = -.32, x = 1, "Ectoderm")
text(y=-.16, x = 2.8, srt=45, "unbalanced")
text(y=-.16, x = 4,srt= 45, "balanced")
text(y = -.32, x = 3.4, "Endoderm")
text(y=-.16, x = 5.2, srt=45, "unbalanced")
text(y=-.16, x = 6.4,srt= 45, "balanced")
text(y = -.32, x = 5.8, "Mesoderm")

barplot(AP_AUC_means, ylim= c(-0.4, 1), col = c("darkred", "coral3", "darkgoldenrod4", "goldenrod", "darkgreen","darkolivegreen4", "antiquewhite4", "azure4", "darkblue", "cornflowerblue" ), xaxt= "n",main = "AUC Anterior-Posterior")
arrows(seq(0.7,18,1.2), (AP_AUC_means - AP_AUC_sd), seq(.7, 18, 1.2), (AP_AUC_means + AP_AUC_sd), length=.05, angle=90, code=3)
text(y=-.16, x = 0.5, srt=45, "unbalanced")
text(y=-.16, x = 1.6,srt= 45, "balanced")
text(y = -.32, x = 1, "Anterior")
text(y=-.16, x = 2.8, srt=45, "unbalanced")
text(y=-.16, x = 4,srt= 45, "balanced")
text(y = -.32, x = 3.4, "Central")
text(y=-.16, x = 5.2, srt=45, "unbalanced")
text(y=-.16, x = 6.4,srt= 45, "balanced")
text(y = -.32, x = 5.8, "Posterior")

barplot(AP_PRAUC_means, ylim= c(-0.4, 1), col = c("darkred", "coral3", "darkgoldenrod4", "goldenrod", "darkgreen","darkolivegreen4", "antiquewhite4", "azure4", "darkblue", "cornflowerblue" ), xaxt= "n",main = "PR-AUC Anterior-Posterior")
arrows(seq(0.7,18,1.2), (AP_PRAUC_means - AP_PRAUC_sd), seq(.7, 18, 1.2), (AP_PRAUC_means + AP_PRAUC_sd), length=.05, angle=90, code=3)
text(y=-.16, x = 0.5, srt=45, "unbalanced")
text(y=-.16, x = 1.6,srt= 45, "balanced")
text(y = -.32, x = 1, "Anterior")
text(y=-.16, x = 2.8, srt=45, "unbalanced")
text(y=-.16, x = 4,srt= 45, "balanced")
text(y = -.32, x = 3.4, "Central")
text(y=-.16, x = 5.2, srt=45, "unbalanced")
text(y=-.16, x = 6.4,srt= 45, "balanced")
text(y = -.32, x = 5.8, "Posterior")

colnames(FSna_put)[2] <- "Ectoderm"
balanced_ect_comp <- express[,c(8,11:length(express))]
FSna_Ectoderm_rf_500 <- randomForest(Ectoderm ~., data = balanced_ect_comp, ntree = 500)
prediction_FSna_Ectoderm <- predict(FSna_Ectoderm_rf_500, newdata = FSna_put, type = "class")
FSna_Ectoderm_frame <- cbind(FSna_put[,1:2], prediction_FSna_Ectoderm)


#varImpPlot(Rdorsal_Ectoderm_rf_500, main = "Balanced Ectoderm", lcolor = "gray")



colnames(FSna_put)[2] <- "Endoderm"
balanced_end_comp <- express[,c(9,11:length(express))]
Rdorsal_Endoderm_rf_500 <- randomForest(Endoderm ~., data = balanced_end_comp, ntree = 500)
prediction_FSna_Endoderm <- predict(Rdorsal_Endoderm_rf_500, newdata = FSna_put, type = "class")
FSna_Endoderm_frame <- cbind(FSna_put[,1:2], prediction_FSna_Endoderm)

#varImpPlot(Rdorsal_Endoderm_rf_500, main = "Balanced Endoderm", lcolor = "gray")

colnames(FSna_put)[2] <- "Mesoderm"
balanced_mes_comp <- express[,c(10:length(express))]
Rdorsal_Mesoderm_rf_500 <- randomForest(Mesoderm ~., data = balanced_mes_comp, ntree = 500)
prediction_FSna_Mesoderm <- predict(Rdorsal_Mesoderm_rf_500, newdata = FSna_put, type = "class")
FSna_Mesoderm_frame <- cbind(FSna_put[,1:2], prediction_FSna_Mesoderm)

ect_unbalanced <- express[,c(8,11:length(express))]
end_unbalanced <- express[,c(9,11:length(express))]
mes_unbalanced <- express[,c(10:length(express))]
A_unbalanced <- express[,c(5,11:length(express))]
C_unbalanced <- express[,c(6,11:length(express))]
P_unbalanced <- express[,c(7,11:length(express))]
#varImpPlot(Rdorsal_Mesoderm_rf_500, main = "Balanced Mesoderm", lcolor = "gray")
FSna_Ectoderm_rf_500 <- randomForest(Ectoderm ~., data = ect_unbalanced, ntree = 500)
prediction_FSnaEct_prosp <- predict(FSna_Ectoderm_rf_500, newdata = FSna_prosp, type = "class")
FSna_Endoderm_rf_500 <- randomForest(Endoderm ~., data = end_unbalanced, ntree = 500)
prediction_FSnaEnd_prosp <- predict(FSna_Endoderm_rf_500, newdata = FSna_prosp, type = "class")
FSna_Mesoderm_rf_500 <- randomForest(Mesoderm ~., data = mes_unbalanced, ntree = 500)
prediction_FSnaMes_prosp <- predict(FSna_Mesoderm_rf_500, newdata = FSna_prosp, type = "class")
FSna_A_rf_500 <- randomForest(A ~., data = A_unbalanced, ntree = 500)
prediction_FSnaA_prosp <- predict(FSna_A_rf_500, newdata = FSna_prosp, type = "class")
FSna_C_rf_500 <- randomForest(C ~., data = C_unbalanced, ntree = 500)
prediction_FSnaC_prosp <- predict(FSna_C_rf_500, newdata = FSna_prosp, type = "class")
FSna_P_rf_500 <- randomForest(P ~., data = P_unbalanced, ntree = 500)
prediction_FSnaP_prosp <- predict(FSna_P_rf_500, newdata = FSna_prosp, type = "class")

FSna_exp_prosp <- cbind(FSna_prosp[,1:2], prediction_FSnaEct_prosp, prediction_FSnaEnd_prosp, prediction_FSnaMes_prosp, prediction_FSnaA_prosp, prediction_FSnaC_prosp, prediction_FSnaP_prosp)

write.table(FSna_exp_prosp, "Furlong_Snail_500_prospective.tsv", sep = "\t")

write.table(FSna_Ect_prosp, "Furlong_Snail_Ect_prosp.tsv", sep = "\t")

write.table(FSna_Ectoderm_frame, "FSna_Ectoderm_predict.tsv", sep="\t")
write.table(FSna_Endoderm_frame, "FSna_Endoderm_predict.tsv", sep="\t")
write.table(FSna_Mesoderm_frame, "FSna_Mesoderm_predict.tsv", sep="\t")

colnames(FSna_put)[2] <- "A"
balanced_A_comp <- express[,c(5,11:length(express))]
Rdorsal_A_rf_500 <- randomForest(A ~., data = balanced_A_comp, ntree = 500)
prediction_FSna_A <- predict(Rdorsal_A_rf_500, newdata = FSna_put, type = "class")
FSna_A_frame <- cbind(FSna_put[,1:2], prediction_FSna_A)

varImpPlot(Rdorsal_A_rf_500, main = "Balanced Anterior", lcolor = "gray")

colnames(FSna_put)[2] <- "C"
balanced_C_comp <- express[,c(6,11:length(express))]
Rdorsal_C_rf_500 <- randomForest(C ~., data = balanced_C_comp, ntree = 500)
prediction_FSna_C <- predict(Rdorsal_C_rf_500, newdata = FSna_put, type = "class")
FSna_C_frame <- cbind(FSna_put[,1:2], prediction_FSna_C)

varImpPlot(Rdorsal_C_rf_500, main = "Balanced Center", lcolor = "gray")

colnames(FSna_put)[2] <- "P"
balanced_P_comp <- express[c(7,11:length(express))]
Rdorsal_P_rf_500 <- randomForest(P ~., data = balanced_P_comp, ntree = 500)
prediction_FSna_P <- predict(Rdorsal_P_rf_500, newdata = FSna_put, type = "class")
FSna_P_frame <- cbind(FSna_put[,1:2], prediction_FSna_P)

varImpPlot(Rdorsal_P_rf_500, main = "Balanced Posterior", lcolor = "gray")

write.table(FSna_A_frame, "FSna_A_predict.tsv", sep="\t")
write.table(FSna_C_frame, "FSna_C_predict.tsv", sep="\t")
write.table(FSna_P_frame, "FSna_P_predict.tsv", sep="\t")




colnames(Mac_Dl_put)[2] <- "Ectoderm"
balanced_ect_comp <- balanced_ect[,c(8,11:length(balanced_ect))]
Mdorsal_Ectoderm_rf_500 <- randomForest(Ectoderm ~., data = balanced_ect_comp, ntree = 500)
prediction_MDl_Ectoderm <- predict(Mdorsal_Ectoderm_rf_500, newdata = Mac_Dl_put, type = "class")
Mac_Dl_Ectoderm_frame <- cbind(Mac_Dl_put[,1:2], prediction_MDl_Ectoderm)

colnames(Mac_Dl_put)[2] <- "Mesoderm"
balanced_mes_comp <- balanced_mes[,c(10,11:length(balanced_ect))]
Mdorsal_Mesoderm_rf_500 <- randomForest(Mesoderm ~., data = balanced_mes_comp, ntree = 500)
prediction_MDl_Mesoderm <- predict(Mdorsal_Mesoderm_rf_500, newdata = Mac_Dl_put, type = "class")
Mac_Dl_Mesoderm_frame <- cbind(Mac_Dl_put[,1:2], prediction_MDl_Mesoderm)

colnames(Mac_Dl_put)[2] <- "Endoderm"
balanced_end_comp <- balanced_end[,c(9,11:length(balanced_end))]
Rdorsal_Endoderm_rf_500 <- randomForest(Endoderm ~., data = balanced_end_comp, ntree = 500)
prediction_Mac_Dl_Endoderm <- predict(Rdorsal_Endoderm_rf_500, newdata = Mac_Dl_put, type = "class")
Mac_Dl_Endoderm_frame <- cbind(Mac_Dl_put[,1:2], prediction_Mac_Dl_Endoderm)

colnames(Mac_Dl_put)[2] <- "A"
balanced_A_comp <- balanced_A[,c(5,11:length(balanced_A))]
Rdorsal_A_rf_500 <- randomForest(A ~., data = balanced_A_comp, ntree = 500)
prediction_Mac_Dl_A <- predict(Rdorsal_A_rf_500, newdata = Mac_Dl_put, type = "class")
Mac_Dl_A_frame <- cbind(Mac_Dl_put[,1:2], prediction_Mac_Dl_A)


colnames(Mac_Dl_put)[2] <- "C"
balanced_C_comp <- balanced_C[,c(6,11:length(balanced_C))]
Rdorsal_C_rf_500 <- randomForest(C ~., data = balanced_C_comp, ntree = 500)
prediction_Mac_Dl_C <- predict(Rdorsal_C_rf_500, newdata = Mac_Dl_put, type = "class")
Mac_Dl_C_frame <- cbind(Mac_Dl_put[,1:2], prediction_Mac_Dl_C)

colnames(Mac_Dl_put)[2] <- "P"
balanced_P_comp <- balanced_P[,c(7,11:length(balanced_P))]
Rdorsal_P_rf_500 <- randomForest(P ~., data = balanced_P_comp, ntree = 500)
prediction_Mac_Dl_P <- predict(Rdorsal_P_rf_500, newdata = Mac_Dl_put, type = "class")
Mac_Dl_P_frame <- cbind(Mac_Dl_put[,1:2], prediction_Mac_Dl_P)


write.table(Mac_Dl_Ectoderm_frame, "Mac_Dl_Ectoderm_predict.tsv", sep="\t")
write.table(Mac_Dl_Endoderm_frame, "Mac_Dl_Endoderm_predict.tsv", sep="\t")
write.table(Mac_Dl_Mesoderm_frame, "Mac_Dl_Mesoderm_predict.tsv", sep="\t")

write.table(Mac_Dl_A_frame, "Mac_Dl_A_predict.tsv", sep="\t")
write.table(Mac_Dl_C_frame, "Mac_Dl_C_predict.tsv", sep="\t")
write.table(Mac_Dl_P_frame, "Mac_Dl_P_predict.tsv", sep="\t")

colnames(Mac_Dl_put)[2] <- "Ectoderm"
ect_unbalanced <- express[,c(8,11:length(express))]
Mdorsal_Ectoderm_rf_500_unbalanced <- randomForest(Ectoderm ~., data = ect_unbalanced, ntree = 500)
prediction_MDl_Ectoderm_unbalanced <- predict(Mdorsal_Ectoderm_rf_500_unbalanced, newdata = Mac_Dl_put, type = "class")
Mdl_Ectoderm_frame_unbalanced <- cbind(Mac_Dl_put[,1:2], prediction_MDl_Ectoderm_unbalanced)


colnames(Mac_Dl_put)[2] <- "A"
balanced_A_comp <- balanced_A[,c(5,11:length(balanced_A))]
Rdorsal_A_rf_500 <- randomForest(A ~., data = balanced_A_comp, ntree = 500)
prediction_Mac_Dl_A <- predict(Rdorsal_A_rf_500, newdata = Mac_Dl_put, type = "class")
Mac_Dl_A_frame <- cbind(Mac_Dl_put[,1:2], prediction_Mac_Dl_A)

varImpPlot(Rdorsal_A_rf_500, main = "Balanced Anterior", lcolor = "gray")

colnames(Mac_Dl_put)[2] <- "C"
balanced_C_comp <- balanced_C[,c(6,11:length(balanced_C))]
Rdorsal_C_rf_500 <- randomForest(C ~., data = balanced_C_comp, ntree = 500)
prediction_Mac_Dl_C <- predict(Rdorsal_C_rf_500, newdata = Mac_Dl_put, type = "class")
Mac_Dl_C_frame <- cbind(Mac_Dl_put[,1:2], prediction_Mac_Dl_C)

varImpPlot(Rdorsal_C_rf_500, main = "Balanced Center", lcolor = "gray")

colnames(Mac_Dl_put)[2] <- "P"
balanced_P_comp <- balanced_P[,c(7,11:length(balanced_P))]
Rdorsal_P_rf_500 <- randomForest(P ~., data = balanced_P_comp, ntree = 500)
prediction_Mac_Dl_P <- predict(Rdorsal_P_rf_500, newdata = Mac_Dl_put, type = "class")
Mac_Dl_P_frame <- cbind(Mac_Dl_put[,1:2], prediction_Mac_Dl_P)

varImpPlot(Rdorsal_P_rf_500, main = "Balanced Posterior", lcolor = "gray")

write.table(Mac_Dl_A_frame, "Mac_Dl_A_predict.tsv", sep="\t")
write.table(Mac_Dl_C_frame, "Mac_Dl_C_predict.tsv", sep="\t")
write.table(Mac_Dl_P_frame, "Mac_Dl_P_predict.tsv", sep="\t")

write.table(Mdl_Ectoderm_frame_unbalanced, "Mac_Dl_Ect_unbalanced.tsv", sep = "\t")


resample_graph <- function(function_choice, dataset, tab_choice, run_name, readme_info, graph_title, reps, samplenum, stratafun){
    mymodel <- resample_rf500(function_choice, dataset, tab_choice, run_name, readme_info, description = "README!", reps, samplenum, stratafun, "none")
	my_means <- apply(mymodel, MARGIN = 2, FUN = mean)
	my_sd <- apply(mymodel, MARGIN=2, FUN=sd)
	par(mgp = c(0,1,0))
	barplot(my_means, ylim=c(-6,100), col=c(33,"darkgreen","lightgreen","orange", "darkred"), , xaxt="n", main = graph_title)
	arrows(seq(.7, 18, 1.2), (my_means - my_sd), seq(.7, 18, 1.2), (my_means + my_sd), length=.05, angle=90, code=3)
	text(y=-4, x=3, "rf_500")
	text(y=-4, x=9, "rf_1000")
	text(y=-4, x=15, "rf_1500")
	return_list <- list(my_means[12:15], my_sd[12:15])
	return_vals <- c(return_list[[1]], return_list[[2]])
	return(return_vals)
}

FSna_balanced <- resample_graph(enhancer_class_status_scale, FSna_train_cons, tab_FSna_cons, "FSnaCons", "FSnaCons", "FSnaCons", 20, 60, strata_2var)
FSna_unbalanced <- resample_graph(enhancer_class_status_scale, current_all, tab_all, "FSna_all", "FSna_all", "FSna_all", 20, 60, strata_Percent_var)

tab_MacDl_cons <- table(Mac_Dl_train_cons$status)
tab_MacDl_cons

tab_MacDl_all <- table(Mac_Dl_train_all$status)
tab_MacDl_all

tab_Mac_Dl_all <- table(Mac_Dl_train_all$status)
tab_Mac_Dl_all

temp <- resample_graph(enhancer_class_status_scale, Mac_Dl_train_all, tab_MacDl_all, "MacDl_all1", "MacDl_all1", "MacDl_all1", 20, 60, strata_Percent_var)
temp2 <- resample_graph(enhancer_class_status_scale, Mac_Dl_train_all, tab_Mac_Dl_all, "Mac_Dl_all1", "Mac_Dl_all1", "Mac_Dl_all1", 20, 60, strata_Percent_var)


all_results_balanced <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_balanced", "Ectoderm_balanced", "Ectoderm balanced", 20, 60, strata_2var)
all_results_balanced 

all_results_unbalanced <- resample_graph(enhancer_class_Ectoderm, express, tab_ect, "Ectoderm_unbalanced", "Ectoderm_unbalanced", "Ectoderm unbalanced", 20, 60, strata_Percent_var)
all_results_unbalanced

all_results_unbalanced_mes <- resample_graph(enhancer_class_Mesoderm, express, tab_mes, "Mesoderm_unbalanced", "Mesoderm_unbalanced", "Mesoderm unbalanced", 20, 60, strata_Percent_var)
all_results_unbalanced_mes

my_features <- drop_Zeit_Twist
drop_Zeit_results <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropZeitTwi", "Ectoderm_dropZeitTwi", "Ectoderm Drop ZeitTwi", 20, 60, strata_2var)
drop_Zeit_results
my_features <- drop_White_H3K27ac
drop_White_ac <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropWhite", "Ectoderm_dropWhite", 60, "Ectoderm Drop White H3K27ac")
drop_White_ac 
my_features <- drop_Rushlow_Dorsal
drop_Rushlow <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropRushlow", "Ectoderm_dropRushlow", 60, "Ectoderm Drop Rushlow Dorsal")
drop_Rushlow
my_features <- all_features
drop_nothing <-  resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropNothing", "Ectoderm_dropNothing", 60, "Ectoderm Drop Nothing")
drop_nothing
my_features <- drop_MacDorsal
drop_MacDl <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropMacDl", "Ectoderm_dropMacDl", 60, "Ectoderm Drop MacDl")
drop_MacDl
my_features <- drop_MacTwist
drop_MacTwi <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropMacTwi", "Ectoderm_dropMacTwi", 60, "Ectoderm Drop MacTwi")
drop_MacTwi

my_features <- drop_Dl
drop_Dl_on <- resample_graph(enhancer_class_status_scale,Mac_Dl_train_cons, tab_Mac_Dl_cons, "on_off_dropDl", "on_off_dropDl","On/Off Drop Dorsal resample 20", 20, 104)
drop_Dl_on
drop_Dl_ect <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropDl", "Ectoderm_dropDl", "Ectoderm Drop Dorsal, resample 20", 20, 104)
drop_Dl_ect
my_features <- all_features
drop_nothing_ect <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropNothing", "Ectoderm_dropNothing", "Ectoderm Drop Nothing, resample 20", 20, 104)
drop_nothing_ect
my_features <- drop_Bcd
drop_Bicoid_ect <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropBcd", "Ectoderm_dropBcd", "Ectoderm Drop Bicoid, resample 20", 20, 104)
drop_Bicoid_ect
my_features <- drop_Dl
drop_Dl_mes <- resample_graph(enhancer_class_Mesoderm, balanced_mes, balanced_mes_table, "Mesoderm_dropDl", "Mesoderm_dropDl", "Mesoderm Drop Dorsal, resample 20", 20, 104)
drop_Dl_mes
my_features <- all_features
drop_nothing_mes <- resample_graph(enhancer_class_Mesoderm, balanced_mes, balanced_mes_table, "Mesoderm_dropNothing", "Mesoderm_dropNothing", "Mesoderm Drop Nothing, resample 20", 20, 104)
drop_nothing_mes
my_features <- drop_Dl
drop_Dl_end <- resample_graph(enhancer_class_Endoderm, balanced_end, balanced_end_table, "Endoderm_dropDl", "Endoderm_dropDl", "Endoderm Drop Dorsal, resample 20", 20, 104)
drop_Dl_end
my_features <- all_features
drop_nothing_end <- resample_graph(enhancer_class_Endoderm, balanced_end, balanced_end_table, "Endoderm_dropNothing", "Endoderm_dropNothing", "Endoderm Drop Nothing, resample 20", 20, 104)
drop_nothing_end
my_features <- drop_Dl
drop_Dl_A <- resample_graph(enhancer_class_A, balanced_A, balanced_A_table, "A_dropDl", "A_dropDl","Anterior Drop Dorsal, resample 20", 20, 104)
drop_Dl_A
my_features <- all_features
drop_nothing_A <- resample_graph(enhancer_class_A, balanced_A, balanced_A_table, "A_dropnothing", "A_dropnothing", "Anterior Drop Nothing, resample 20", 20, 104)
drop_nothing_A
my_features <- drop_Bcd
drop_Bcd_A <- resample_graph(enhancer_class_A, balanced_A, balanced_A_table, "A_dropBcd", "A_dropBcd", "Anterior Drop Bicoid, resample 20", 5, 104)
drop_Bcd_A
drop_Dl_A <- resample_graph(enhancer_class_A, balanced_A, balanced_A_table, "A_dropDl", "A_dropDl", "Anterior Drop Dorsal, resample 20", 5, 104)
drop_Dl_A
drop_Dl_C <- resample_graph(enhancer_class_C, balanced_C, balanced_C_table, "C_dropDl", "C_dropDl", "Center Drop Dorsal, resample 20", 5, 104)
drop_Dl_C
drop_Dl_P <- resample_graph(enhancer_class_P, balanced_P, balanced_P_table, "P_dropDl", "P_dropDl", "Posterior Drop Dorsal, resample 20", 5, 104)
drop_Dl_P
drop_Bcd_C <- resample_graph(enhancer_class_C, balanced_C, balanced_C_table, "C_dropBcd", "C_dropBcd", "Central Drop Bicoid, resample 20", 5, 104)
drop_Bcd_C
drop_Bcd_P <- resample_graph(enhancer_class_P, balanced_P, balanced_P_table, "P_dropBcd", "P_dropBcd", "Posterior Drop Bicoid, resample 20", 5, 104)
drop_Bcd_P

drop_Dl_A
drop_Dl_C
drop_Dl_P

drop_Bcd_A
drop_Bcd_C
drop_Bcd_P

#on/off
my_features <- all_features
crm_set <- Mac_Dl_train_cons
crm_tab <- tab_MacDl_cons 

my_features <- all_features
on_all <- resample_graph(enhancer_class_status_scale,Mac_Dl_train_cons, tab_Mac_Dl_cons, "on_off_all", "on_off_all","On/Off all features resample 100", 100, 60, strata_2var)
on_all

on_all_unbalanced <- resample_graph(enhancer_class_status_scale,Mac_Dl_train_all, tab_Mac_Dl_all, "on_off_all_unbalanced", "on_off_all_unbalanced","On/Off all features unbalanced resample 100", 100, 60, strata_Percent_var)
on_all_unbalanced

my_features <- all_features
on_all_Mac <- resample_graph(enhancer_class_status_scale,Mac_Dl_train_cons, tab_MacDl_cons, "on_off_all_MacDl", "on_off_all_MacDl","MacDl On/Off all features resample 100", 100, 60, strata_2var)
on_all_Mac

on_all_Mac_unbalanced <- resample_graph(enhancer_class_status_scale,Mac_Dl_train_all, tab_MacDl_all, "on_off_all_MacDl_unbalanced", "on_off_all_MacDl_unbalanced","MacDl On/Off all features resample 100 unbalanced", 100, 60, strata_Percent_var)
on_all_unbalanced

on_all_unbalanced <- resample_graph(enhancer_class_status_scale,Mac_Dl_train_all, tab_Mac_Dl_all, "on_off_all_unbalanced", "on_off_all_unbalanced","On/Off all features unbalanced resample 100", 100, 60, strata_Percent_var)
on_all_unbalanced
#Ectoderm
my_features <- all_features
Ectoderm_all <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_all", "Ectoderm_all","Ectoderm all features resample 100", 100, 104, strata_2var)
Ectoderm_all

Ectoderm_all_unbalanced <- resample_graph(enhancer_class_Ectoderm,express, tab_ect, "Ectoderm_all_unbalanced", "Ectoderm_all_unbalanced","Ectoderm all features resample 100 unbalanced", 100, 104, strata_Percent_var)
Ectoderm_all_unbalanced

#Mesoderm
my_features <- all_features
Mesoderm_all <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_all", "Mesoderm_all","Mesoderm all features resample 100", 100, 104, strata_2var)
Mesoderm_all

Mesoderm_all_unbalanced <- resample_graph(enhancer_class_Mesoderm,express, tab_mes, "Mesoderm_all_unbalanced", "Mesoderm_all_unbalanced","Mesoderm all features resample 100 unbalanced", 100, 104, strata_Percent_var)
Mesoderm_all_unbalanced

#Endoderm
my_features <- all_features
Endoderm_all <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_all", "Endoderm_all","Endoderm all features resample 100", 100, 104,strata_2var)
Endoderm_all

Endoderm_all_unbalanced <- resample_graph(enhancer_class_Endoderm,express, tab_end, "Endoderm_all_unbalanced", "Endoderm_all_unbalanced","Endoderm all features resample 100 unbalanced", 100, 104, strata_Percent_var)
Endoderm_all_unbalanced

#Anterior
my_features <- all_features
A_all <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_all", "A_all","A all features resample 100", 100, 104, strata_2var)
A_all

A_all_unbalanced <- resample_graph(enhancer_class_A,express, tab_A, "A_all_unbalanced", "A_all_unbalanced","Anterior all features resample 100 unbalanced", 100, 104, strata_Percent_var)
A_all_unbalanced

#Center
my_features <- all_features
C_all <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_all", "C_all","C all features resample 100", 100, 104, strata_2var)
C_all

C_all_unbalanced <- resample_graph(enhancer_class_C,express, tab_C, "C_all_unbalanced", "C_all_unbalanced","Center all features resample 100 unbalanced", 100, 104, strata_Percent_var)
C_all_unbalanced

#Posterior
my_features <- all_features
P_all <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_all", "P_all","P all features resample 100", 100, 104, strata_2var)
P_all

P_all_unbalanced <- resample_graph(enhancer_class_P,express, tab_P, "P_all_unbalanced", "P_all_unbalanced","Posterior all features resample 100 unbalanced", 100, 104, strata_Percent_var)
P_all_unbalanced

combo_data_on_balance <- as.data.frame(rbind(on_all, on_all_unbalanced, on_all_Mac, on_all_Mac_unbalanced))
combo_data_DV_balance <- as.data.frame(rbind(Ectoderm_all, Ectoderm_all_unbalanced, Mesoderm_all, Mesoderm_all_unbalanced, Endoderm_all, Endoderm_all_unbalanced))
combo_data_AP_balance <- as.data.frame(rbind(A_all, A_all_unbalanced, C_all, C_all_unbalanced, P_all, P_all_unbalanced))



combo_data_only_on <- as.data.frame(rbind(on_only_Zelda ,on_only_Eisen10Cad ,on_only_Eisen10Kni,on_only_Eisen10Gt,on_only_Eisen13Hb,on_only_Eisen13Bcd ,on_only_Eisen13Kr ,on_only_FurlongSnail,on_only_FurlongTwist,on_only_KK,on_only_MacBicoid,on_only_MacCad ,on_only_MacGiant,on_only_MacHb,on_only_MacHry ,on_only_MacSna ,on_only_MacTwi ,on_only_MacDl,on_only_White,on_only_Zeit ))

combo_data_drop_on <- as.data.frame(rbind(on_all,on_drop_ZeitTwi,on_drop_White_H3K27ac,on_drop_MacDorsal,on_drop_MacTwist,on_drop_MacSnail,on_drop_MacHairy,on_drop_MacHunchback,on_drop_MacGiant,on_drop_MacCaudal,on_drop_MacBicoid ,on_drop_KK_H3K4me1,on_drop_KK_H3K27ac ,on_drop_FurlongTwist,on_drop_FurlongSnail ,on_drop_Eisen13Kr  ,on_drop_Eisen13Bcd,on_drop_Eisen13Hb,on_drop_Eisen13Gt,on_drop_Eisen10Kni,on_drop_Eisen10Cad,on_drop_Zelda))
 
combo_data_only_Ectoderm <- as.data.frame(rbind(Ectoderm_only_Zelda ,Ectoderm_only_Eisen10Cad ,Ectoderm_only_Eisen10Kni,Ectoderm_only_Eisen10Gt,Ectoderm_only_Eisen13Hb,Ectoderm_only_Eisen13Bcd ,Ectoderm_only_Eisen13Kr ,Ectoderm_only_FurlongSnail,Ectoderm_only_FurlongTwist,Ectoderm_only_KK,Ectoderm_only_MacBicoid,Ectoderm_only_MacCad ,Ectoderm_only_MacGiant,Ectoderm_only_MacHb,Ectoderm_only_MacHry ,Ectoderm_only_MacSna ,Ectoderm_only_MacTwi ,Ectoderm_only_MacDl, Ectoderm_only_White,Ectoderm_only_Zeit ))

combo_data_drop_Ectoderm <- as.data.frame(rbind(Ectoderm_all,Ectoderm_drop_ZeitTwi,Ectoderm_drop_White_H3K27ac,Ectoderm_drop_MacDorsal,Ectoderm_drop_MacTwist,Ectoderm_drop_MacSnail,Ectoderm_drop_MacHairy,Ectoderm_drop_MacHunchback,Ectoderm_drop_MacGiant,Ectoderm_drop_MacCaudal,Ectoderm_drop_MacBicoid ,Ectoderm_drop_KK_H3K4me1,Ectoderm_drop_KK_H3K27ac ,Ectoderm_drop_FurlongTwist,Ectoderm_drop_FurlongSnail ,Ectoderm_drop_Eisen13Kr  ,Ectoderm_drop_Eisen13Bcd,Ectoderm_drop_Eisen13Hb,Ectoderm_drop_Eisen13Gt,Ectoderm_drop_Eisen10Kni,Ectoderm_drop_Eisen10Cad,Ectoderm_drop_Zelda))

combo_data_only_Mesoderm <- as.data.frame(rbind(Mesoderm_only_Zelda ,Mesoderm_only_Eisen10Cad ,Mesoderm_only_Eisen10Kni,Mesoderm_only_Eisen10Gt,Mesoderm_only_Eisen13Hb,Mesoderm_only_Eisen13Bcd ,Mesoderm_only_Eisen13Kr ,Mesoderm_only_FurlongSnail,Mesoderm_only_FurlongTwist,Mesoderm_only_KK,Mesoderm_only_MacBicoid,Mesoderm_only_MacCad ,Mesoderm_only_MacGiant,Mesoderm_only_MacHb,Mesoderm_only_MacHry ,Mesoderm_only_MacSna ,Mesoderm_only_MacTwi ,Mesoderm_only_MacDl, Mesoderm_only_White,Mesoderm_only_Zeit ))

combo_data_drop_Mesoderm <- as.data.frame(rbind(Mesoderm_all,Mesoderm_drop_ZeitTwi,Mesoderm_drop_White_H3K27ac,Mesoderm_drop_MacDorsal,Mesoderm_drop_MacTwist,Mesoderm_drop_MacSnail,Mesoderm_drop_MacHairy,Mesoderm_drop_MacHunchback,Mesoderm_drop_MacGiant,Mesoderm_drop_MacCaudal,Mesoderm_drop_MacBicoid ,Mesoderm_drop_KK_H3K4me1,Mesoderm_drop_KK_H3K27ac ,Mesoderm_drop_FurlongTwist, Mesoderm_drop_FurlongSnail ,Mesoderm_drop_Eisen13Kr  ,Mesoderm_drop_Eisen13Bcd,Mesoderm_drop_Eisen13Hb,Mesoderm_drop_Eisen13Gt,Mesoderm_drop_Eisen10Kni,Mesoderm_drop_Eisen10Cad,Mesoderm_drop_Zelda))

combo_data_only_Endoderm <- as.data.frame(rbind(Endoderm_only_Zelda ,Endoderm_only_Eisen10Cad ,Endoderm_only_Eisen10Kni,Endoderm_only_Eisen10Gt,Endoderm_only_Eisen13Hb,Endoderm_only_Eisen13Bcd ,Endoderm_only_Eisen13Kr ,Endoderm_only_FurlongSnail,Endoderm_only_FurlongTwist,Endoderm_only_KK,Endoderm_only_MacBicoid,Endoderm_only_MacCad ,Endoderm_only_MacGiant,Endoderm_only_MacHb,Endoderm_only_MacHry ,Endoderm_only_MacSna ,Endoderm_only_MacTwi ,Endoderm_only_MacDl,Endoderm_only_White,Endoderm_only_Zeit ))

combo_data_drop_Endoderm <- as.data.frame(rbind(Endoderm_all,Endoderm_drop_ZeitTwi,Endoderm_drop_White_H3K27ac,Endoderm_drop_MacDorsal,Endoderm_drop_MacTwist,Endoderm_drop_MacSnail,Endoderm_drop_MacHairy,Endoderm_drop_MacHunchback,Endoderm_drop_MacGiant,Endoderm_drop_MacCaudal,Endoderm_drop_MacBicoid ,Endoderm_drop_KK_H3K4me1,Endoderm_drop_KK_H3K27ac ,Endoderm_drop_FurlongTwist,Endoderm_drop_FurlongSnail ,Endoderm_drop_Eisen13Kr  ,Endoderm_drop_Eisen13Bcd,Endoderm_drop_Eisen13Hb,Endoderm_drop_Eisen13Gt,Endoderm_drop_Eisen10Kni,Endoderm_drop_Eisen10Cad,Endoderm_drop_Zelda))

combo_data_only_A <- as.data.frame(rbind(A_only_Zelda ,A_only_Eisen10Cad ,A_only_Eisen10Kni,A_only_Eisen10Gt,A_only_Eisen13Hb,A_only_Eisen13Bcd ,A_only_Eisen13Kr ,A_only_FurlongSnail, A_only_FurlongTwist,A_only_KK,A_only_MacBicoid,A_only_MacCad ,A_only_MacGiant,A_only_MacHb,A_only_MacHry ,A_only_MacSna ,A_only_MacTwi ,A_only_MacDl,A_only_White,A_only_Zeit ))
 
combo_data_drop_A <- as.data.frame(rbind(A_all,A_drop_ZeitTwi,A_drop_White_H3K27ac,A_drop_MacDorsal,A_drop_MacTwist,A_drop_MacSnail,A_drop_MacHairy,A_drop_MacHunchback,A_drop_MacGiant,A_drop_MacCaudal,A_drop_MacBicoid ,A_drop_KK_H3K4me1,A_drop_KK_H3K27ac ,A_drop_FurlongTwist,A_drop_FurlongSnail ,A_drop_Eisen13Kr  ,A_drop_Eisen13Bcd,A_drop_Eisen13Hb,A_drop_Eisen13Gt,A_drop_Eisen10Kni,A_drop_Eisen10Cad,A_drop_Zelda))

combo_data_only_C <- as.data.frame(rbind(C_only_Zelda ,C_only_Eisen10Cad ,C_only_Eisen10Kni,C_only_Eisen10Gt,C_only_Eisen13Hb,C_only_Eisen13Bcd ,C_only_Eisen13Kr ,C_only_FurlongSnail,C_only_FurlongTwist,C_only_KK,C_only_MacBicoid,C_only_MacCad ,C_only_MacGiant,C_only_MacHb,C_only_MacHry ,C_only_MacSna ,C_only_MacTwi ,C_only_MacDl,C_only_White,C_only_Zeit ))

combo_data_drop_C <- as.data.frame(rbind(C_all,C_drop_ZeitTwi,C_drop_White_H3K27ac,C_drop_MacDorsal,C_drop_MacTwist,C_drop_MacSnail,C_drop_MacHairy,C_drop_MacHunchback,C_drop_MacGiant,C_drop_MacCaudal,C_drop_MacBicoid ,C_drop_KK_H3K4me1,C_drop_KK_H3K27ac ,C_drop_FurlongTwist,C_drop_FurlongSnail ,C_drop_Eisen13Kr  ,C_drop_Eisen13Bcd,C_drop_Eisen13Hb,C_drop_Eisen13Gt,C_drop_Eisen10Kni,C_drop_Eisen10Cad,C_drop_Zelda))

combo_data_only_P <- as.data.frame(rbind(P_only_Zelda ,P_only_Eisen10Cad ,P_only_Eisen10Kni,P_only_Eisen10Gt,P_only_Eisen13Hb,P_only_Eisen13Bcd ,P_only_Eisen13Kr ,P_only_FurlongSnail,P_only_FurlongTwist,P_only_KK,P_only_MacBicoid,P_only_MacCad ,P_only_MacGiant,P_only_MacHb,P_only_MacHry ,P_only_MacSna ,P_only_MacTwi ,P_only_MacDl,P_only_White,P_only_Zeit ))

combo_data_drop_P <- as.data.frame(rbind(P_all,P_drop_ZeitTwi,P_drop_White_H3K27ac,P_drop_MacDorsal,P_drop_MacTwist,P_drop_MacSnail,P_drop_MacHairy,P_drop_MacHunchback,P_drop_MacGiant,P_drop_MacCaudal,P_drop_MacBicoid ,P_drop_KK_H3K4me1,P_drop_KK_H3K27ac ,P_drop_FurlongTwist, P_drop_FurlongSnail ,P_drop_Eisen13Kr  ,P_drop_Eisen13Bcd,P_drop_Eisen13Hb,P_drop_Eisen13Gt,P_drop_Eisen10Kni,P_drop_Eisen10Cad,P_drop_Zelda))




plot_combo_only <- function(datacombo, datatitle){
	colnames(datacombo) <- c("rf500_tpos", "rf500_tneg", "rf500_fpos", "rf500_fneg", "sd_tpos", "sd_tneg", "sd_fpos", "sd_fneg")
	par(mar=c(10,4.1,4.1,2.1))
	plot(range(1:nrow(datacombo)), range(0,100), type='n', ylab="percentage", xlab = "", main = datatitle, xaxt="n")
	lines(1:nrow(datacombo), datacombo$rf500_tpos, type='b', col="darkgreen", lwd = 3)
	lines(1:nrow(datacombo), datacombo$rf500_tneg, type='b', col="lightgreen")
	lines(1:nrow(datacombo), datacombo$rf500_fpos, type='b', col="orange", lwd = 3)
	lines(1:nrow(datacombo), datacombo$rf500_fneg, type='b', col="darkred")
	arrows(x0=1:nrow(datacombo), y0=(datacombo$rf500_tpos - 2*datacombo$sd_tpos), y1=(datacombo$rf500_tpos + 2*datacombo$sd_tpos), length=.05, angle=90, code=3)
	axis(1, at=1:nrow(datacombo), labels=rownames(datacombo),las=2)
}

plot_combo_only(combo_data_on_balance, "On/Off Dl balance")
plot_combo_only(combo_data_DV_balance, "Dorsal-ventral balance")
plot_combo_only(combo_data_AP_balance, "Anterior-Posterior balance")

plot_combo_only(combo_data_only_on, "On/Off Rushlow Dl + one factor")
plot_combo_only(combo_data_only_Ectoderm, "Ectoderm Rushlow Dl + one factor")
plot_combo_only(combo_data_only_Mesoderm, "Mesoderm Rushlow Dl + one factor")
plot_combo_only(combo_data_only_Endoderm, "Endoderm Rushlow Dl + one factor")
plot_combo_only(combo_data_only_A, "Anterior Rushlow Dl + one factor")
plot_combo_only(combo_data_only_C, "Center Rushlow Dl + one factor")
plot_combo_only(combo_data_only_P, "Posterior Rushlow Dl + one factor")

plot_combo_drop <- function(datacombo, datatitle){
	colnames(datacombo) <- c("rf500_tpos", "rf500_tneg", "rf500_fpos", "rf500_fneg", "sd_tpos", "sd_tneg", "sd_fpos", "sd_fneg")
	par(mar=c(12.5,4.1,4.1,2.1))
	plot(range(1:nrow(datacombo)), range(0,100), type='n', ylab="percentage", xlab = "", main = datatitle, xaxt="n")
	lines(1:nrow(datacombo), datacombo$rf500_tpos, type='b', col="darkgreen", lwd = 3)
	lines(1:nrow(datacombo), datacombo$rf500_tneg, type='b', col="lightgreen", lwd = 3)
	lines(1:nrow(datacombo), datacombo$rf500_fpos, type='b', col="orange", lwd = 3)
	lines(1:nrow(datacombo), datacombo$rf500_fneg, type='b', col="darkred", lwd = 3)
	arrows(x0=1:nrow(datacombo), y0=(datacombo$rf500_tpos - 2*datacombo$sd_tpos), y1=(datacombo$rf500_tpos + 2*datacombo$sd_tpos), length=.05, angle=90, code=3)
	arrows(x0=1:nrow(datacombo), y0=(datacombo$rf500_fpos - 2*datacombo$sd_fpos), y1=(datacombo$rf500_fpos + 2*datacombo$sd_fpos), length=.05, angle=90, code=3)
	#arrows(x0=1:nrow(datacombo), y0=(datacombo$rf500_tneg - 2*datacombo$sd_tneg), y1=(datacombo$rf500_tneg + 2*datacombo$sd_tneg), length=.05, angle=90, code=3)
	#arrows(x0=1:nrow(datacombo), y0=(datacombo$rf500_fneg - 2*datacombo$sd_fneg), y1=(datacombo$rf500_fneg + 2*datacombo$sd_fneg), length=.05, angle=90, code=3)
	axis(1, at=1:nrow(datacombo), labels=rownames(datacombo),las=2)
}

plot_combo_drop(combo_data_drop_on, "On/Off drop one factor")
plot_combo_drop(combo_data_drop_Ectoderm, "Ectoderm drop one factor")
plot_combo_drop(combo_data_drop_Mesoderm, "Mesoderm drop one factor")
plot_combo_drop(combo_data_drop_Endoderm, "Endoderm drop one factor")
plot_combo_drop(combo_data_drop_A, "Anterior drop one factor")
plot_combo_drop(combo_data_drop_C, "Center drop one factor")
plot_combo_drop(combo_data_drop_P, "Posterior drop one factor")