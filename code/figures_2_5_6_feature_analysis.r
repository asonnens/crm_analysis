require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 
require(devtools)
library(ggbiplot)
require(kernlab)
require(ROCR)
library(lattice)
library(plyr)
require(ROSE)
library(ggplot2)
library(reshape2)
require(caTools)
require(matrixStats)
library(caret)
library(ggfortify)
library("FactoMineR")
library("factoextra")

require(GGally)
require(scales)
require(gridExtra)
require(CCA)


if(!require("gplots")) {
     install.packages("gplots", dependencies = TRUE)
     library(gplots)
 }
if (!require("RColorBrewer")) {
     install.packages("RColorBrewer", dependencies = TRUE)
     library(RColorBrewer)
 }

my_palette <- colorRampPalette(c("coral", "seashell", "cadetblue"))(n = 299)
 

setwd("shadow_enhancers/machine_learning_input/February")

#Activity in Stage 4-6 embryos
all_data <- read.table("all_data.tsv", header = TRUE)
all_data <- all_data[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]


# Figure 2
#heatmap correlations dataset figure
#heatmap of how high scoring features overlap
high_overlap <- read.table("high_overlap.tsv", header = TRUE)
h_o <- as.matrix(high_overlap)
heatmap.2(h_o, trace = "none", margins =c(7.5,7.5), dendrogram = "none", density.info = "none", Rowv = FALSE, Colv = FALSE, col = my_palette, scale = "none", keysize = 1.2)
tmp <- h_o
tmp[upper.tri(tmp)] <- NA
heatmap.2(tmp,  trace = "none", margins =c(7.5,7.5), dendrogram = "none", density.info = "none", Rowv = FALSE, Colv = FALSE, col = my_palette, scale = "none", keysize = 1.2)

#heatmap of how all features overlap
all_overlap <- read.table("all_overlap.tsv", header = TRUE)
a_o <- as.matrix(all_overlap)
heatmap.2(a_o, trace = "none", margins =c(7.5,7.5), dendrogram = "none", density.info = "none", Rowv = FALSE, Colv = FALSE, col = my_palette, scale = "none", keysize = 1.2)
tmp <- a_o
tmp[upper.tri(tmp)] <- NA
heatmap.2(tmp,  trace = "none", margins =c(7.5,7.5), dendrogram = "none", density.info = "none", Rowv = FALSE, Colv = FALSE, col = my_palette, scale = "none", keysize = 1.2)

#heatmap of how all features correlate
high_corr <- read.table("high_correlation.tsv", header = TRUE)
h_c <- as.matrix(high_corr)/100
heatmap.2(h_c, trace = "none", margins =c(7.5,7.5), dendrogram = "none", density.info = "none", Rowv = FALSE, Colv = FALSE, col = my_palette, scale = "none", keysize = 1.2)
tmp <- h_c
tmp[upper.tri(tmp)] <- NA
heatmap.2(tmp,  trace = "none", margins =c(7.5,7.5), dendrogram = "none", density.info = "none", Rowv = FALSE, Colv = FALSE, col = my_palette, scale = "none", keysize = 1.2)


#heatmap of how all features correlate
all_corr <- read.table("all_correlation.tsv", header = TRUE)
a_c <- as.matrix(all_corr)/100
heatmap.2(a_c, trace = "none", margins =c(7.5,7.5), dendrogram = "none", density.info = "none", Rowv = FALSE, Colv = FALSE, col = my_palette, scale = "none", keysize = 1.2)
tmp <- a_c
tmp[upper.tri(tmp)] <- NA
heatmap.2(tmp,  trace = "none", margins =c(7.5,7.5), dendrogram = "none", density.info = "none", Rowv = FALSE, Colv = FALSE, col = my_palette, scale = "none", keysize = 1.2)

#figure 5
#plots of feature correlations with reporters
#plots of how features associate with bound regions
#plots including unbound regions
#to do-- make supplementary figure that includes all features
#Dorsal 2015 seq
dev.off()
temp_data <- subset(all_data, all_data$Dorsal_2015_seq > (min(all_data$Dorsal_2015_seq)))
ggplot(temp_data, aes(Dorsal_2015_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
ggplot(all_data, aes(Dorsal_2015_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())




#Dorsal 2015 seq
temp_data <- subset(all_data, all_data$Dorsal_2015_seq > (min(all_data$Dorsal_2015_seq)))
ggplot(temp_data, aes(Dorsal_2015_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dorsal 2009 chip
temp_data <- subset(all_data, all_data$Dorsal_2009_chip > (min(all_data$Dorsal_2009_chip)))
ggplot(temp_data, aes(Dorsal_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Snail 2009 chip
temp_data <- subset(all_data, all_data$Snail_2009_chip > (min(all_data$Snail_2009_chip)))
ggplot(temp_data, aes(Snail_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Twist 2009 chip
temp_data <- subset(all_data, all_data$Twist_2009_chip > (min(all_data$Twist_2009_chip)))
ggplot(temp_data, aes(Twist_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())


#Snail 2014 chip
temp_data <- subset(all_data, all_data$Snail_2014_chip > (min(all_data$Snail_2014_chip)))
ggplot(temp_data, aes(Snail_2014_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Snail_2014_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Twist 2014 chip
temp_data <- subset(all_data, all_data$Twist_2014_chip > (min(all_data$Twist_2014_chip)))
ggplot(temp_data, aes(Twist_2014_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Twist_2014_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Twist 2011 seq
temp_data <- subset(all_data, all_data$Twist_2011_seq > (min(all_data$Twist_2011_seq)))
ggplot(temp_data, aes(Twist_2011_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())


#Zelda 2011 seq
temp_data <- subset(all_data, all_data$Zelda_2011_seq > (min(all_data$Zelda_2011_seq)))
ggplot(temp_data, aes(Zelda_2011_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Zelda_2011_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Bicoid 2013 seq
temp_data <- subset(all_data, all_data$Bicoid_2013_seq > (min(all_data$Bicoid_2013_seq)))
ggplot(temp_data, aes(Bicoid_2013_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Bicoid_2013_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Bicoid 2009 chip
temp_data <- subset(all_data, all_data$Bicoid_2009_chip > (min(all_data$Bicoid_2009_chip)))
ggplot(temp_data, aes(Bicoid_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Caudal 2010 seq
temp_data <- subset(all_data, all_data$Caudal_2010_seq > (min(all_data$Caudal_2010_seq)))
ggplot(temp_data, aes(Caudal_2010_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Caudal_2010_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Caudal 2009 chip
temp_data <- subset(all_data, all_data$Caudal_2009_chip > (min(all_data$Caudal_2009_chip)))
ggplot(temp_data, aes(Caudal_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())


#Hunchback 2013 seq
temp_data <- subset(all_data, all_data$Hunchback_2013_seq > (min(all_data$Hunchback_2013_seq)))
ggplot(temp_data, aes(Hunchback_2013_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Hunchback_2013_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Hunchback 2009 chip
temp_data <- subset(all_data, all_data$Hunchback_2009_chip > (min(all_data$Hunchback_2009_chip)))
ggplot(temp_data, aes(Hunchback_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Hunchback_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Giant 2013 seq
temp_data <- subset(all_data, all_data$Giant_2013_seq > (min(all_data$Giant_2013_seq)))
ggplot(temp_data, aes(Giant_2013_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Giant_2013_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Giant 2009 chip
temp_data <- subset(all_data, all_data$Giant_2009_chip > (min(all_data$Giant_2009_chip)))
ggplot(temp_data, aes(Giant_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Giant_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Kruppel 2013 seq
temp_data <- subset(all_data, all_data$Kruppel_2013_seq > (min(all_data$Kruppel_2013_seq)))
ggplot(temp_data, aes(Kruppel_2013_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Kruppel_2013_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Knirps 2010 seq
temp_data <- subset(all_data, all_data$Knirps_2010_seq > (min(all_data$Knirps_2010_seq)))
ggplot(temp_data, aes(Knirps_2010_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Knirps_2010_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Hairy 2009 chip
temp_data <- subset(all_data, all_data$Hairy_2009_chip > (min(all_data$Hairy_2009_chip)))
ggplot(temp_data, aes(Hairy_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Hairy_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#H3K27ac 2015 seq
temp_data <- subset(all_data, all_data$H3K27ac_2015_seq > (min(all_data$H3K27ac_2015_seq)))
ggplot(temp_data, aes(H3K27ac_2015_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(H3K27ac_2015_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#H3K27ac 2010 seq
temp_data <- subset(all_data, all_data$H3K27ac_2010_seq > (min(all_data$H3K27ac_2010_seq)))
ggplot(temp_data, aes(H3K27ac_2010_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(H3K27ac_2010_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#H3K4me1 2015 seq
temp_data <- subset(all_data, all_data$H3K4me1_2015_seq > (min(all_data$H3K4me1_2015_seq)))
ggplot(temp_data, aes(H3K4me1_2015_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(H3K4me1_2015_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#p300 2010 seq
temp_data <- subset(all_data, all_data$p300_2010_seq > (min(all_data$p300_2010_seq)))
ggplot(temp_data, aes(p300_2010_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(p300_2010_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())


#Zelda solexa motif
temp_data <- subset(all_data, all_data$Zld_motif_solexa > (min(all_data$Zld_motif_solexa)))
ggplot(temp_data, aes(Zld_motif_solexa, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Zld_motif_solexa, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Zelda sanger motif
temp_data <- subset(all_data, all_data$Zld_motif_sanger > (min(all_data$Zld_motif_sanger)))
ggplot(temp_data, aes(Zld_motif_sanger, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Zld_motif_sanger, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Snail solexa motif
temp_data <- subset(all_data, all_data$Snail_motif_solexa > (min(all_data$Snail_motif_solexa)))
ggplot(temp_data, aes(Snail_motif_solexa, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Snail_motif_solexa, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Snail Sanger motif
temp_data <- subset(all_data, all_data$Snail_motif_Sanger > (min(all_data$Snail_motif_Sanger)))
ggplot(temp_data, aes(Snail_motif_Sanger, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Snail_motif_Sanger, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Snail FlyReg motif
temp_data <- subset(all_data, all_data$Snail_motif_FlyReg > (min(all_data$Snail_motif_FlyReg)))
ggplot(temp_data, aes(Snail_motif_FlyReg, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Snail_motif_FlyReg, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dorsal FlyReg motif
temp_data <- subset(all_data, all_data$Dorsal_motif_FlyReg > (min(all_data$Dorsal_motif_FlyReg)))
ggplot(temp_data, aes(Dorsal_motif_FlyReg, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Dorsal_motif_FlyReg, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dorsal NBT motif
temp_data <- subset(all_data, all_data$Dorsal_motif_NBT > (min(all_data$Dorsal_motif_NBT)))
ggplot(temp_data, aes(Dorsal_motif_NBT, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Dorsal_motif_NBT, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Twist FlyReg motif
temp_data <- subset(all_data, all_data$Twist_motif_FlyReg > (min(all_data$Twist_motif_FlyReg)))
ggplot(temp_data, aes(Twist_motif_FlyReg, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Twist_motif_FlyReg, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Twist da motif
temp_data <- subset(all_data, all_data$Twist_motif_da > (min(all_data$Twist_motif_da)))
ggplot(temp_data, aes(Twist_motif_da, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
#ggplot(all_data, aes(Twist_motif_da, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dsim
ggplot(all_data, aes(Dsim, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dsec
ggplot(all_data, aes(Dsec, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dyak
ggplot(all_data, aes(Dyak, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dere
ggplot(all_data, aes(Dere, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dana
ggplot(all_data, aes(Dana, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dpse
ggplot(all_data, aes(Dpse, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dwil
ggplot(all_data, aes(Dwil, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dvir
ggplot(all_data, aes(Dvir, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#Dgri
ggplot(all_data, aes(Dgri, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

  
   "Dsim"                "Dsec"               
[37] "Dyak"                "Dere"                "Dana"                "Dpse"               
[41] "Dwil"                "Dvir"                "Dgri"

# Figure 6
#correlation plots
short_name <- c("enhancer","status","Zld_11_s","Dl_15_s","Dl_09_c","Sna_14_c","Sna_09_c","Twi_11_s","Twi_14_c","Twi_09_c","Bcd_13_s","Bcd_09_c","Cad_10_s","Cad_09_c","Hb_13_s","Hb_09_c","Gt_13_s","Gt_09_c","Kr_13_s","Kni_10_s","Hry_09_c","H3K27ac_15","H3K27ac_10","H3K4me1_15","p300_10", "Zld_m1", "Zld_m2", "Dl_m1","Dl_m2","Sna_m1","Sna_m2","Sna_m3","Twi_m1","Twi_m2","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")
all_data_short <- all_data
colnames(all_data_short) <- short_name
current_all <- all_data_short
my_features <- current_all[,3:43]
current_all_on <- subset(current_all, status == "on")
current_all_off <- subset(current_all, status == "off")
#scaled data
data1 <- current_all_on[,3:43]
data_scale_on <- scale(data1)
data2 <- current_all_off[,3:43]
data_scale_off <- scale(data2)
#dorsal snail twist features
Dl_Sna_Twi_on <- data_scale_on[,1:8]
Dl_Sna_Twi_off <- data_scale_off[,1:8]
#AP features
AP_corr_on <- data_scale_on[,c(1,9:19)]
AP_corr_off  <- data_scale_off[,c(1,9:19)]
#Species blastz comparisons
evolution_data_on <- data_scale_on[,31:39]
evolution_data_off <- data_scale_off[,31:39]

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


#ALL DATA 
all_corr_data <- cor(as.matrix(data_scale_on))
upper_tri <- get_upper_tri(all_corr_data)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 5, hjust = 1))+
 coord_fixed()

#Dl Sna Twi on 
Dl_Sna_Twi_data <- cor(as.matrix(Dl_Sna_Twi_on))
upper_tri <- get_upper_tri(Dl_Sna_Twi_data)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
#AP on
AP_corr  <- cor(as.matrix(data_scale_on[,c(1,9:19)]))
upper_tri <- get_upper_tri(AP_corr)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
 #Dl Sna Twi off
Dl_Sna_Twi_data <- cor(as.matrix(Dl_Sna_Twi_off))
upper_tri <- get_upper_tri(Dl_Sna_Twi_data)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
#AP off
AP_corr  <- cor(as.matrix(data_scale_off[,c(1,9:19)]))
upper_tri <- get_upper_tri(AP_corr)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
#Dorsals_on
Dorsals <- cor(as.matrix(data_scale_on[,c(1:3,26,27)]))
upper_tri <- get_upper_tri(Dorsals)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 25, hjust = 1, face = "bold"))+
	 theme(axis.text.y = element_text(face = "bold",
    size = 25, hjust = 1)) + theme(legend.position="none")
 coord_fixed()
 
#Dorsals_off
Dorsals <- cor(as.matrix(data_scale_off[,c(1:3,26,27)]))
upper_tri <- get_upper_tri(Dorsals)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 25, hjust = 1, face = "bold"))+
	 theme(axis.text.y = element_text(face = "bold",
    size = 25, hjust = 1)) + theme(legend.position="none")
 coord_fixed()



#Snails_on
Snails <- cor(as.matrix(data_scale_on[,c(1,4:5,26:28)]))
upper_tri <- get_upper_tri(Snails)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()

#Twists_on 
Twists <- cor(as.matrix(data_scale_on[,c(1,6:8,31:32)]))
upper_tri <- get_upper_tri(Twists)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
 
#Twists_off
Twists <- cor(as.matrix(data_scale_off[,c(1,6:8,31:32)]))
upper_tri <- get_upper_tri(Twists)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()


#Chromatin modifiers on
chromatin <- cor(as.matrix(data_scale_on[,c(1,20:25)]))
upper_tri <- get_upper_tri(chromatin)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 25, hjust = 1, face = "bold"))+
	 theme(axis.text.y = element_text(face = "bold",
    size = 25, hjust = 1)) + theme(legend.position="none")
 coord_fixed()
 
#Chromatin modifiers off
chromatin <- cor(as.matrix(data_scale_off[,c(1,20:25)]))
upper_tri <- get_upper_tri(chromatin)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 25, hjust = 1, face = "bold"))+
	 theme(axis.text.y = element_text(face = "bold",
    size = 25, hjust = 1)) + theme(legend.position="none")
 coord_fixed()
#evolution_data_on
evolution_data <- cor(as.matrix(data_scale_on[,c(34:41)]))
upper_tri <- get_upper_tri(evolution_data)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
#evolution_data_off 
evolution_data <- cor(as.matrix(data_scale_off[,c(34:41)]))
upper_tri <- get_upper_tri(evolution_data)
melted_data <- melt(upper_tri)
melted_data <- na.omit(melted_data)
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "coral", high = "cadetblue", mid = "seashell", 
   midpoint = 0, limit = c(-0.5,1), name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
