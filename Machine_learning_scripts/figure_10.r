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
      #size = c(180,180))
      size = c(round(variable_table[2]* (2/3)), round(variable_table[1]*(2/3))))
  my_data <- dataset[, c(1,variable_pos,start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

strata_balance <- function(dataset, variable_name, variable_table, variable_pos, dataframe_length = 14, second_val = 7, samplesize){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  var_size <- variable_table[1]
  if (variable_table[2] < variable_table[1]){ var_size <- variable_table[2]}
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
      size = c(round(var_size * (2/3)),round(var_size * (2/3))))
  my_data <- dataset[, c(1,variable_pos,start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing_initial <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  var_size_2 <- var_size - round(var_size * 2/3)
  my_strata2 <- strata(my_testing_initial, stratanames = c(variable_name), method = "srswor", 
      size = c(var_size_2,var_size_2))
  my_testing <- my_data[rownames(my_data) %in% my_strata2$ID_unit, ]
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

	




resample_default <- function(function_name, dataset, datatable, variable_name, outfilename, description="README!", reps, samplesize, stratafun){
    temp_dir <- paste("../February/", outfilename, sep = "")
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

resample_rf500 <- function(function_name, dataset, datatable, variable_name, outfilename, description="README!", reps, samplesize, stratafun, graph, validated_data, info, featurenum, format, model){
    temp_dir <- paste("",outfilename, sep = "")
    dir.create(temp_dir)
	info_file_name <- paste(temp_dir, "input_data_summary.txt", sep = "/")
	input_data_name <- paste(temp_dir, "input_data.csv", sep = "/")
	readme_file_name <- paste(temp_dir, "README.txt", sep = "/")
	write(description, file = readme_file_name)
	write.csv(dataset, file = input_data_name)
	datathing <- data.frame()
	mydataframe <- data.frame()
	par(mfrow=c(1,1))
	count <- 0
    for(i in 1:reps){
		count <- count + 1
	    temp_file_name <- paste(outfilename, i, sep = "_")
		temp_file_name <- paste(temp_file_name, ".csv", sep = "")
		temp_file_name <- paste(temp_dir, temp_file_name, sep = "/")
        myresult <- function_name(dataset, datatable, temp_file_name, info_file_name, samplesize, stratafun, count, graph, mydataframe, validated_data, info, featurenum,format, model)
			datathing[i,"rf500_acc"] <- myresult[1]
			datathing[i,"rf500_tpos"] <- myresult[2]
			datathing[i,"rf500_tneg"] <- myresult[3]
			datathing[i,"rf500_fpos"] <- myresult[4]
			datathing[i,"rf500_fneg"] <- myresult[5]
			datathing[i, "rf500_AUC"] <- myresult[6]
			datathing[i, "rf500_PRAUC"] <- myresult[7]
			mydataframe <- myresult[[8]]
     } 
    return_frame <- list(datathing, mydataframe)	 
	return(return_frame)            
}   



enhancer_class_status_scale <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize, stratafun, count, graph, mydataframe, validated_data, info, featurenum, format, model){
    enhancers_set <- stratafun(dataset, "status", datatable,2, featurenum, featurenum - 3, samplesize)   
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_training_vars <- enhancers_training2[,3:featurenum]
	enhancers_data_trainingscale <- scale(enhancers_training_vars)
    enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_trainingscale)
	enhancers_training <- enhancers_trainingscale
	enhancers_test2 <- enhancers_set[[2]]
	enhancers_data_testscale <- scale(enhancers_test2[,3:featurenum])
    enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_test <- enhancers_testscale
	actual <- enhancers_test$status
	pred.validated <- validated_data
	if (format == "PCA"){pred.validated <- predict(model, validated_data)}
	#random forest 500 trees                  
	enhancer_rf_500 <- randomForest(status ~., data = enhancers_training[,2:length(enhancers_training)],  ntree = 500, importance = TRUE)
	if (pred.validated == "none"){
	prediction_validated <- "none"
	prediction_validated_on_frame <- "none"}
	else{
	prediction_validated <- predict(enhancer_rf_500, current_target = pred.validated, type = "class")
	prediction_validated_on_frame <- cbind(validated_data[,1:2], prediction_validated)}
	outfilename <- paste("prediction_validated", info, toString(count),".tsv", sep = "_")
    write.table(prediction_validated_on_frame, outfilename , sep = "\t")
    import <- importance(enhancer_rf_500, type = 1)
	imp1 <- as.vector(import[,1])
	if (count == 1){mydataframe <- import}
	else {mydataframe <- cbind(mydataframe, imp1)}
	rf_500_predicted <- predict(enhancer_rf_500, current_target=enhancers_test, type="class")
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
	return_list <- list(rf_500_prediction_success, rf_500_actual_on, rf_500_true_off, rf_500_false_on, rf_500_false_off, AUC, prec_reca, mydataframe)
	return(return_list)             
}      

stage_compare <- read.table("Stage_compare.tsv", header = TRUE)
stage_compare <- stage_compare[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]


current_all <- stage_compare
current_tab <- table(stage_compare$status)
all_features <- c(3:43)
my_features <- all_features


#Active vs. inactive
stage_active <- subset(stage_compare,  status != "semi-specific")
stage_active$status <- factor(stage_active$status)
stage.pca <- prcomp(stage_active[,3:43],center = TRUE,scale. = TRUE) 
autoplot(stage.pca, choices = c(1,2), data = stage_active,colour = "status")
tempplot("PC1", "PC2",stage.pca, choices = c(1,2), data = stage_active,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") + scale_color_manual(values = c("purple","darkgreen", "purple", "darkgreen","darkgray"))


stage_group <- stage_active
levels(stage_group$status)[levels(stage_group$status)=="always"]<-"now"
levels(stage_group$status)[levels(stage_group$status)=="4_8"]<-"now"
levels(stage_group$status)[levels(stage_group$status)=="9_12"]<-"later"
levels(stage_group$status)[levels(stage_group$status)=="13_16"]<-"later"
stage_compare <- droplevels()
stage.pca <- prcomp(stage_group[,3:43],center = TRUE,scale. = TRUE) 
autoplot(stage.pca, choices = c(1,2), data = stage_group,colour = "status")
tempplot("PC1", "PC2",stage.pca, choices = c(1,2), data = stage_group,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") + scale_color_manual(values = c("darkred","cadetblue", "black"))
g <- lda(status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = stage_group,prior = rep(1, 2)/3)
hab.lda.values <- predict(g, status)
hab.class <- predict(g)$class
plot(hab.lda.values$x[,1], ylab=c("LDA Axis 1"), col = stage_group$status, pch = 19)

ldahist(data = hab.lda.values$x[,1], g=stage_group$status,col = "cadetblue4")
#Active now vs Always active
active_active <- subset(stage_compare,  status != "never")
active_active <- subset(active_active,  status != "semi-specific")
active_active <- subset(active_active,  status != "13_16")
active_active <- subset(active_active,  status != "9_12")
active_active$status <- factor(active_active$status)
active.pca <- prcomp(active_active[,3:43],center = TRUE,scale. = TRUE) 
tempplot("PC1", "PC2",active.pca, choices = c(1,2), data = active_active,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") 

g <- lda(status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = active_active,prior = rep(1, 2)/2)
hab.lda.values <- predict(g, active_active)
hab.class <- predict(g)$class
plot(hab.lda.values$x[,1], ylab=c("LDA Axis 1"), col = active_active$status, pch = 19)

ldahist(data = hab.lda.values$x[,1], g=active_active$status,col = "cadetblue4")
#Always vs Never
always_never <- subset(stage_compare,  status != "4_8")
always_never <- subset(always_never,  status != "semi-specific")
always_never <- subset(always_never,  status != "13_16")
always_never <- subset(always_never,  status != "9_12")
always_never$status <- factor(always_never$status)
active.pca <- prcomp(always_never[,3:43],center = TRUE,scale. = TRUE) 
tempplot("PC1", "PC2",active.pca, choices = c(1,2), data = always_never,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") 



#Always vs Never
always_never <- subset(stage_compare,  status != "4_8")
always_never <- subset(always_never,  status != "semi-specific")
always_never <- subset(always_never,  status != "13_16")
always_never <- subset(always_never,  status != "9_12")
always_never$status <- factor(always_never$status)
active.pca <- prcomp(always_never[,3:43],center = TRUE,scale. = TRUE) 
tempplot("PC1", "PC2",active.pca, choices = c(1,2), data = always_never,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") 

current_target <- read.table("target_file.tsv", header = TRUE)

prosp_target <- read.table("prospective_file.tsv", header = TRUE)





#ROC for always vs never
current_all <- always_never
current_tab <- table(current_all$status)

current_all$status <- gsub('always', 'on', current_all$status)
rownames(current_all) <- NULL
current_tab <- table(current_all$status)

all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "always_never", description = "always_never", "always_never", 10, 100, strata_balance, "ROC", current_target, "always_never", 43, "NA", "NA")
all_recall_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "always_never", description = "always_never", "always_never", 10, 100, strata_balance, "recall",prosp_target, "always_never",43, "NA", "NA")

all_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

#ROC for always vs now
current_all <- active_active
rownames(current_all) <- NULL
current_tab <- table(current_all$status)


all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "always_never", description = "always_never", "always_never", 10, 100, strata_balance, "ROC", current_target, "always_never", 43, "NA", "NA")
all_recall_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "always_never", description = "always_never", "always_never", 10, 100, strata_balance, "recall",prosp_target, "always_never",43, "NA", "NA")

all_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

#random data the same size as active now vs always active

random_data_stage <- data.frame(replicate(41,sample(0:1000,length(rownames(current_all)),rep=TRUE)))

random_data_stage <-cbind(current_all[,0:2],random_data_stage)
colnames(random_data_stage) <- colnames(current_all)

all_ROC_data <- resample_rf500(enhancer_class_status_scale, random_data_stage, current_tab, "always_never", description = "always_never", "always_never", 10, 100, strata_balance, "ROC", current_target, "always_never", 43, "NA", "NA")
random_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
random_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

 myframe <- cbind(all_means[6:7], random_means[6:7] )
 myplot <- barplot(myframe, beside = TRUE, col = c("darkgoldenrod", "darkgreen"), ylim = c(0,1), ylab = "Area Under Curve", cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.5)
 legend("topright", rownames(myframe), cex = 1.5, fill = c("darkgoldenrod", "darkgreen"))
 mysd_frame <- cbind(all_sd[6:7], random_sd[6:7])
segments(myplot, myframe - mysd_frame * 2, myplot,
         myframe + mysd_frame  * 2, lwd = 1.5)

arrows(myplot, myframe - mysd_frame * 2, myplot,
         myframe + mysd_frame  * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
	   
g <- lda(status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = current_all,prior = rep(1, 2)/2)
hab.lda.values <- predict(g, current_all)
hab.class <- predict(g)$class
plot(hab.lda.values$x[,1], ylab=c("LDA Axis 1"), col = current_all$status, pch = 19)
ldahist(data = hab.lda.values$x[,1], g=current_all$status)

g_cv <- lda(current_all$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = current_all,prior = rep(1, 2)/2, CV=TRUE)
ggplot(g_cv, aes(x[,1] fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
ggplot(all_data, aes(Dsim, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

plot(current_all$Zelda_2011_seq, current_all$H3K4me1_2015_seq,
     xlab="Zelda", ylab="H3K4me1", 
     main="Cross-Validated LD Classification", 
     pch=as.numeric(current_all$status) + 19,
     col=as.numeric(g_cv$class)) 


#always vs. later
stage_compare <- read.table("Stage_compare_always_later.tsv",header = TRUE)
current_all <- stage_compare[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","H3K27ac_2015_seq","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]
current_tab <- table(current_all$status)


my_data_temp <- current_all[, c(1,2,3:43)]

my_temp <- strata(current_all, stratanames = c("status"), method = "srswor", 
                    size = c(round(current_tab[1] * (2/3)), round(current_tab[1] * 2/3)))
my_training_temp <- my_data_temp[rownames(my_data_temp) %in% my_temp$ID_unit, ]

my_testing_initial_temp <- my_data_temp[!rownames(my_data_temp) %in% my_temp$ID_unit, ] 

  my_temp2 <- strata(my_testing_initial_temp, stratanames = c("status"), method = "srswor", 
      size = c(27,27))
  my_testing_temp <- my_data_temp[rownames(my_data_temp) %in% my_temp2$ID_unit, ]
  
 enhancer_rf_temp <- randomForest(status ~., data = my_training_temp[,2:length(my_training_temp)],  ntree = 500, importance = TRUE)
 	rf_500_pred_temp <- predict(enhancer_rf_temp, current_target=my_testing_temp, type="class")
	rf_500_temp_table <- table(actual = my_testing_temp$status, predicted = rf_500_pred_temp)
					

strata_balance <- function(dataset, variable_name, variable_table, variable_pos, dataframe_length = 14, second_val = 7, samplesize){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  var_size <- variable_table[1]
  if (variable_table[2] < variable_table[1]){ var_size <- variable_table[2]}
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
      size = c(round(var_size * (2/3)),round(var_size * (2/3))))
  my_data <- dataset[, c(1,variable_pos,start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing_intial <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  var_size_2 <- var_size - round(var_size * 2/3)
  my_strata2 <- strata(my_testing_initial, stratanames = c(variable_name), method = "srswor", 
      size = c(var_size_2,var_size_2))
  my_testing <- my_training_initial[rownames(my_training_initial) %in% my_strata2$ID_unit, ]
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}


#later vs. never
stage_compare <- read.table("later_never.tsv",header = TRUE)
current_all <- stage_compare[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]
current_tab <- table(current_all$status)


all_features <- c(3:43)
my_features <- all_features
all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 100, strata_Percent_var, "ROC", current_target, "All_Data_Reduced", 43, "NA", "NA")
all_recall_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_Percent_var, "recall",current_target, "All_Data_Reduced",43, "NA", "NA")

all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "ROC", current_all, "All_Data_Reduced", 43, "NA", "NA")
all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "recall", current_target, "All_Data_Reduced", 43, "NA", "NA")

all_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)


#Activity in Stage 4-6 embryos
all_data <- read.table("all_data.tsv", header = TRUE)
all_data <- all_data[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]

#Activity across developmental stages
stage_compare <- read.table("stage_compare.tsv",header = TRUE)
stage_compare <- stage_compare[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]

g <- lda(status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = stage_compare,prior = rep(1, 2)/2)
hab.lda.values <- predict(g, status)
hab.class <- predict(g)$class
plot(hab.lda.values$x[,1], ylab=c("LDA Axis 1"), col = stage_compare$status, pch = 19)

ldahist(data = hab.lda.values$x[,1], g=stage_compare$status,col = "cadetblue4")

stage_group <- stage_compare
stage_group <- subset(stage_group,  status != "always")
stage_group <- subset(stage_group,  status != "semi-specific")
stage_group <- subset(stage_group,  status != "never")
levels(stage_group$status)[levels(stage_group$status)=="9_12"]<-"later"
levels(stage_group$status)[levels(stage_group$status)=="13_16"]<-"later"
stage_group$status <- factor(stage_group$status)
#levels(stage_group$status)[levels(stage_group$status)=="4_8"]<-"now"
#levels(stage_group$status)[levels(stage_group$status)=="9_12"]<-"later"
#levels(stage_group$status)[levels(stage_group$status)=="13_16"]<-"later"



stage.pca <- prcomp(stage_group[,3:43],center = TRUE,scale. = TRUE) 
autoplot(stage.pca, choices = c(1,2), data = stage_group,colour = "status")
tempplot("PC1", "PC2",stage.pca, choices = c(1,2), data = stage_group,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") + scale_color_manual(values = c("coral", "cadetblue4"))

g <- lda(status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = stage_group,prior = rep(1, 2)/2, CV = FALSE)
hab.lda.values <- predict(g, stage_group)
hab.class <- predict(g)$class
plot(hab.lda.values$x[,1], ylab=c("LDA Axis 1"), col = stage_group$status, pch = 19)

ab.lda <- lda(grp ~ ., data=hab_std)

hab.lda.values <- predict(hab.lda, hab_std)
hab.class <- predict(hab.lda)$class

#create a histogram of the discriminant function values
ldahist(data = hab.lda.values$x[,1], g=grp)


ldahist(data = hab.lda.values$x[,1], g=stage_group$status,col = "cadetblue4")

g <- lda(status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = stage_group,prior = rep(1, 2)/2, CV=TRUE)
hab.lda.values <- predict(g, status)
hab.class <- predict(g)$class
plot(hab.lda.values$x[,1], ylab=c("LDA Axis 1"), col = stage_group$status, pch = 19)

ldahist(data = hab.lda.values$x[,1], g=stage_group$status,col = "cadetblue4")