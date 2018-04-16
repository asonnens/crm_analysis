###########################################################################
#Script for plotting cumulative blastZ scores for DNA elements, 
#from pairwise comparisons of D.mel and related species
#Produces Supplementary figure 1 and some additional figures


#packages
require(MASS)     
require(sampling)
require(devtools)
library(ggbiplot)
require(kernlab)
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
 
#personalization
my_palette <- colorRampPalette(c("coral", "seashell", "cadetblue"))(n = 299)
setwd("shadow_enhancers/machine_learning_input/February")

#loading data
#Activity in Stage 4-6 embryos, crosslisted with 41 features
all_data <- read.table("all_data.tsv", header = TRUE)
all_data <- all_data[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]

#sor
all_data_sorted <- with(all_data, all_data[order(status,enhancer),])

#stages 4-6 comparison of all highly active DNA elements with all inactive DNA elements-- active elements in hashmarks
temp <- rbind(all_data_sorted[1:1228,], all_data_sorted[6859:6925,], all_data_sorted[1229:2461,], all_data_sorted[6926:7008,], all_data_sorted[2462:3734,], all_data_sorted[7009:7079,], all_data_sorted[3735:5708,], all_data_sorted[7080:7198,], all_data_sorted[5758:6868,], all_data_sorted[7202:7250,])
plot(temp$Dsim, pch = 20, col = "darkred", ylim = c(0,10000),ylab = "BlastZ Conservation", xaxt='n', xlab = "Chromosome", xaxs = "i", main = "Conservation vs. Activity, Stage 4-6")
points(temp$Dyak, pch = 20, col = "darkgoldenrod")
points(temp$Dpse, pch = 20, col = "darkgreen")
points(temp$Dgri, pch = 20, col = "cadetblue")
rect(66,  ybottom = -1000, ytop = 11000, xleft = 1228, xright = 1296, border = "darkgray", col = "gray", density = 20)
rect(82,  ybottom = -1000, ytop = 11000, xleft = 2528, xright = 2612, border = "darkgray", col = "gray", density = 20)
rect(70,  ybottom = -1000, ytop = 11000, xleft = 3884, xright = 3956, border = "darkgray", col = "gray", density = 20)
rect(118,  ybottom = -1000, ytop = 11000, xleft = 5929, xright = 6049, border = "darkgray", col = "gray", density = 20)
rect(51,  ybottom = -1000, ytop = 11000, xleft = 7149, xright = 7208, border = "darkgray", col = "gray", density = 20)

#stages 4-6 highly active vs inactive, just chromosome 3R
chr3R <- rbind(all_data_sorted[3735:5708,], all_data_sorted[7080:7198,])
plot(chr3R$Dsim, pch = 20, col = "darkred", ylim = c(0,10000),ylab = "BlastZ Conservation", xaxt='n', xlab = "Chromosome 3R", xaxs = "i", main = "Conservation vs. Activity, Chromosome 3R, Stage 4-6")
points(chr3R$Dyak, pch = 20, col = "darkgoldenrod")
points(chr3R$Dpse, pch = 20, col = "darkgreen")
points(chr3R$Dgri, pch = 20, col = "cadetblue")
rect(118,  ybottom = -1000, ytop = 11000, xleft = 1973, xright = 2091, border = "darkgray", col = "gray", density = 20)


chr3R <- rbind(all_data_sorted[3735:5708,], all_data_sorted[7080:7198,])
chr3R_ordered <- with(chr3R, chr3R[order(enhancer),])
attach(chr3R_ordered); plot(chr3R_ordered$Dsim, ylim = c(0,10000), ylab = "BlastZ Conservation", xaxt = 'n', xlab = "Chromosome 3R", xaxs = "i", main = "Conservation vs. Activity, Chromosome 3R, Stage 4-6",col=c("coral3","darkred")[status], pch = c(20, 19)[status], cex = c(1, 1)[status]); detach(chr3R_ordered)
attach(chr3R_ordered); points(chr3R_ordered$Dyak ,col=c("burlywood","darkgoldenrod4")[status], pch = c(20, 19)[status], cex = c(1, 1)[status]); detach(chr3R_ordered)
attach(chr3R_ordered); points(chr3R_ordered$Dpse ,col=c("darkolivegreen2","darkgreen")[status], pch = c(20, 19)[status], cex = c(1, 1)[status]); detach(chr3R_ordered)
attach(chr3R_ordered); points(chr3R_ordered$Dgri ,col=c("cadetblue2","cadetblue4")[status], pch = c(20, 19)[status], cex = c(1, 1)[status]); detach(chr3R_ordered)


plot(chr3R$Dsim, pch = 20, col = "darkred", ylim = c(0,10000),ylab = "BlastZ Conservation", xaxt='n', xlab = "Chromosome 3R", xaxs = "i", main = "Conservation vs. Activity, Chromosome 3R, Stage 4-6")
points(chr3R$Dyak, pch = 20, col = "darkgoldenrod")
points(chr3R$Dpse, pch = 20, col = "darkgreen")
points(chr3R$Dgri, pch = 20, col = "cadetblue")
rect(118,  ybottom = -1000, ytop = 11000, xleft = 1973, xright = 2091, border = "darkgray", col = "gray", density = 20)

#conservation vs activity, always active vs never active
#Activity across developmental stages
stage_compare <- read.table("stage_compare.tsv",header = TRUE)
stage_compare <- stage_compare[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]

stage_compare_sorted <- with(stage_compare, stage_compare[order(status,enhancer),])
always <- subset(stage_compare_sorted, stage_compare_sorted$status == "always")
never <- subset(stage_compare_sorted,stage_compare_sorted$status == "never")
always_never <- rbind(always, never)
always_never_sorted <- with(always_never,always_never[order(enhancer,status),])

temp <- rbind(always_never_sorted[281:822,], always_never_sorted[1:49,], always_never_sorted[823:1378,], always_never_sorted[50:109,], always_never_sorted[1379:1970,], always_never_sorted[110:151,], always_never_sorted[1971:2851,], always_never_sorted[152:235,], always_never_sorted[2888:3407,], always_never_sorted[237:280,])
plot(temp$Dsim, pch = 20, col = "darkred", ylim = c(0,10000),ylab = "BlastZ Conservation", xaxt='n', xlab = "Chromosome", xaxs = "i", main = "Conservation vs. Activity, Always active vs. Never active")
points(temp$Dyak, pch = 20, col = "darkgoldenrod")
points(temp$Dpse, pch = 20, col = "darkgreen")
points(temp$Dgri, pch = 20, col = "cadetblue")
rect(49,  ybottom = -1000, ytop = 11000, xleft = 542, xright = 592, border = "darkgray", col = "gray", density = 20)
rect(59,  ybottom = -1000, ytop = 11000, xleft = 1147, xright = 1207, border = "darkgray", col = "gray", density = 20)
rect(41,  ybottom = -1000, ytop = 11000, xleft = 1798, xright = 1840, border = "darkgray", col = "gray", density = 20)
rect(83,  ybottom = -1000, ytop = 11000, xleft = 2720, xright = 2804, border = "darkgray", col = "gray", density = 20)
rect(43,  ybottom = -1000, ytop = 11000, xleft = 3323, xright = 3367, border = "darkgray", col = "gray", density = 20)

always_never_sorted_enhancer <- with(always_never,always_never[order(enhancer),])
chr3R_always_never <- always_never_sorted_enhancer[1842:2806,]
chr3R_always_never_ordered <- with(chr3R_always_never, chr3R_always_never[order(enhancer),])
chr3R_always_never_ordered$status <- factor(chr3R_always_never_ordered$status)
attach(chr3R_always_never_ordered); plot(chr3R_always_never_ordered$Dsim, ylim = c(0,10000), ylab = "BlastZ Conservation", xaxt = 'n', xlab = "Chromosome 3R", xaxs = "i", main = "Chromosome 3R, Always active vs. Never Active",col=c("darkred","coral3")[status], pch = c(19, 20)[status], cex = c(1, 1)[status])
attach(chr3R_always_never_ordered); points(chr3R_always_never_ordered$Dyak ,col=c("darkgoldenrod4", "burlywood")[status], pch = c(19, 20)[status], cex = c(1, 1)[status]); detach(chr3R_always_never_ordered)
attach(chr3R_always_never_ordered); points(chr3R_always_never_ordered$Dpse ,col=c("darkgreen","darkolivegreen2")[status], pch = c(19, 20)[status], cex = c(1, 1)[status]); detach(chr3R_ordered)
attach(chr3R_always_never_ordered); points(chr3R_always_never_ordered$Dgri ,col=c("cadetblue4", "cadetblue2")[status], pch = c(19, 20)[status], cex = c(1, 1)[status]); detach(chr3R_ordered)


stage_compare$status=factor(stage_compare$status , levels=levels(stage_compare$status)[c(4,2,3,1,6,5)])

boxplot(stage_compare$Dsim ~ stage_compare$status, ylim = c(0,10000), col = "darkred", pch = 19, outcol = "darkred")
boxplot(stage_compare$Dyak ~ stage_compare$status, ylim = c(0,10000), col = "darkgoldenrod", pch = 19, outcol = "darkgoldenrod")
boxplot(stage_compare$Dpse ~ stage_compare$status, ylim = c(0,10000), col = "darkgreen", pch = 19, outcol = "darkgreen")
boxplot(stage_compare$Dgri ~ stage_compare$status, ylim = c(0,10000), col = "cadetblue", pch = 19, outcol = "cadetblue")