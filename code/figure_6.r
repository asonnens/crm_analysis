require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 
require(devtools)
library(ggplot2)
library(ggbiplot)
require(kernlab)
require(ROCR)
library(lattice)
library(plyr)
require(ROSE)
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

tempplot <- function (a, b, object, data = NULL, scale = 1, ...) 
{
    plot.data <- ggplot2::fortify(object, data = data)
    plot.data$rownames <- rownames(plot.data)
        x.column <- a
        y.column <- b
        loadings.column <- "rotation"
        lam <- object$sdev[1L:2L]
        lam <- lam * sqrt(nrow(plot.data))
    if (scale != 0) {
        lam <- lam^scale
        plot.data[, c(x.column, y.column)] <- t(t(plot.data[, 
            c(x.column, y.column)])/lam)
    }
    plot.columns <- unique(c(x.column, y.column, colnames(plot.data)))
    plot.data <- plot.data[, plot.columns]
    if (!is.null(loadings.column)) {
        loadings.data <- as.data.frame(object[[loadings.column]][, 
            ])
        loadings.data$rownames <- rownames(loadings.data)
        loadings.columns <- unique(c(x.column, y.column, colnames(loadings.data)))
        loadings.data <- loadings.data[, loadings.columns]
    }
    else {
        loadings.data <- NULL
    }
    p <- ggbiplot(plot.data = plot.data, loadings.data = loadings.data, 
        ...)
    return(p)
}

#Activity in Stage 4-6 embryos
all_data <- read.table("all_data.tsv", header = TRUE)
all_data <- all_data[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]

#Expression patterns
ACP <- read.table("ACP.tsv", header = TRUE)
AP_DV <- read.table("AP_DV.tsv", header = TRUE)
current_target <- read.table("target_file.tsv", header = TRUE)

#stage comparison
stage_compare <- read.table("stage_compare.tsv",header = TRUE)
stage_compare <- stage_compare[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]

#PCA all data, and weights 
 
ir.pca <- prcomp(all_data[,3:43],center = TRUE,scale. = TRUE) 
screeplot(ir.pca, type = "l", main = "Scree plot PCA all reporters, 41 features", col = "cadetblue")
a <- autoplot(ir.pca, data = all_data, choices = c(1,2), colour = "status")
a + scale_color_manual(values = c("gray","cadetblue4"))
on_data <- subset(all_data, status == "on")
on_data$status <- factor(on_data$status)
off_data <- subset(all_data, status == "off")
off_data$status <- factor(off_data$status)
a <- autoplot(ir.pca, data = on_data, choices = c(1), colour = "cadetblue4")
a + scale_color_manual(values = "cadetblue4")
b <- autoplot(ir.pca, data = off_data, choices = c(1), colour = "darkgray")
b + scale_color_manual(values = "darkgray")
ggplot (ir.pca, aes (x = PC1, y = PC2,colour = "status")) + stat_density2d ()+ scale_color_manual(values = c("darkgray","cadetblue4"))
ggplot (ir.pca, aes (x = PC1, y = PC2,fill = "status")) + stat_binhex (bins=5, aes (alpha = ..count..)) + facet_grid (. ~ "status")
xyplot(ir.pca$PC1 ~ ir.pca$PC2, ir.pca, groups = all_data$status, pch= 20)
xyplot(a$PC1 ~ a$PC2 | all_data$status, pch= 20)
autoplot(prcomp(all_data[,3:43], center = TRUE, scale. = TRUE), data = all_data, colour = 'status')
temp <- scale(all_data[,3:43])
autoplot(clara(temp,2))
autoplot(fanny(temp, 2), frame = TRUE)
#plots of PC 1 vs PC 2, PC 6 vs PC 7, PC 1 vs PC 6
autoplot(ir.pca, choices = c(1,2), data = all_data,colour = "status")  + scale_color_manual(values = c("darkgray","cadetblue4"))
tempplot("PC1", "PC2",ir.pca, choices = c(1,2), data = all_data,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm")  + scale_color_manual(values = c("darkgray","cadetblue4"))
tempplot("PC3", "PC4",ir.pca, choices = c(1,2), data = all_data,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm")  + scale_color_manual(values = c("darkgray","cadetblue4"))
tempplot("PC5", "PC6",ir.pca, choices = c(1,2), data = all_data,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm")  + scale_color_manual(values = c("darkgray","cadetblue4"))

#Contributions of features to given PCs
fviz_contrib(ir.pca, choice = "var", fill = "cadetblue", axes = 1, top = 10)
fviz_contrib(ir.pca, choice = "var", fill = "cadetblue", axes = 2, top = 10)
fviz_contrib(ir.pca, choice = "var", fill = "cadetblue", axes = 3, top = 10)
fviz_contrib(ir.pca, choice = "var", fill = "cadetblue", axes = 4, top = 10)
fviz_contrib(ir.pca, choice = "var", fill = "cadetblue", axes = 5, top = 10)
fviz_contrib(ir.pca, choice = "var", fill = "cadetblue", axes = 6, top = 10)

facto_summarize(ir.pca, element = "var", result = "contrib", axes = 1)
facto_summarize(ir.pca, element = "var", result = "contrib", axes = 2)
