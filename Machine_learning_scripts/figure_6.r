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
        myresult <- function_name(dataset, datatable, temp_file_name, info_file_name, samplesize, stratafun = strata_Percent_var, count, graph, mydataframe, validated_data, info, featurenum,format, model)
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
	#enhancers_training_vars <- enhancers_training2[,3:featurenum]
	#enhancers_data_trainingscale <- scale(enhancers_training_vars)
    #enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_trainingscale)
	enhancers_training <- enhancers_training2
	enhancers_test2 <- enhancers_set[[2]]
	#enhancers_data_testscale <- scale(enhancers_test2[,3:featurenum])
    #enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_test <- enhancers_test2 
	actual <- enhancers_test$status
	pred.validated <- validated_data
	if (format == "PCA"){pred.validated <- predict(model, validated_data)}
	#random forest 500 trees                  
	enhancer_rf_500 <- randomForest(status ~., data = enhancers_training[,2:length(enhancers_training)],  ntree = 500, importance = TRUE)
	if (pred.validated == "none"){
	prediction_validated <- "none"
	prediction_validated_on_frame <- "none"}
	else{
	prediction_validated <- predict(enhancer_rf_500, newdata = pred.validated, type = "class")
	prediction_validated_on_frame <- cbind(validated_data[,1:2], prediction_validated)}
	outfilename <- paste("prediction_validated", info, toString(count),".tsv", sep = "_")
    write.table(prediction_validated_on_frame, outfilename , sep = "\t")
    import <- importance(enhancer_rf_500, type = 1)
	imp1 <- as.vector(import[,1])
	if (count == 1){mydataframe <- import}
	else {mydataframe <- cbind(mydataframe, imp1)}
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
	return_list <- list(rf_500_prediction_success, rf_500_actual_on, rf_500_true_off, rf_500_false_on, rf_500_false_off, AUC, prec_reca, mydataframe)
	return(return_list)             
}  

enhancer_class_status_balance <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize, stratafun, count, graph, mydataframe, validated_data, info, featurenum, format, model, balance){
    balance_data <- dataset %>%
    group_by(status) %>%
    sample_n(size = samplesize)
    balance_data <- dataset %>%  na.omit()
    balance_data <- balance_data  %>%  droplevels()
    balance_data <- data.frame(balance_data)
	balancetable <- table(balance_data$status)
    enhancers_set <- stratafun(balance_data, "status", balancetable,2, featurenum, featurenum - 3, samplesize)   
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
	datasize <- nrow(balance_data)
	random_data <- data.frame(replicate(41,sample(0:1000,datasize,rep=TRUE)))
    random_data <-cbind(balance_data[,0:2],random_data)
	colnames(random_data) <- long_name
	random_set <- stratafun(random_data, "status", balancetable,2, featurenum, featurenum - 3, samplesize)   
	random_training2 <- random_set[[1]]
	random_training_vars <- random_training2[,3:featurenum]
	random_data_trainingscale <- scale(random_training_vars)
    random_trainingscale <- cbind(random_training2[,1:2],random_data_trainingscale)
	random_training <- random_trainingscale
	random_test2 <- random_set[[2]]
	random_data_testscale <- scale(random_test2[,3:featurenum])
    random_testscale <- cbind(random_test2[,1:2],random_data_testscale)
	random_test <- random_testscale
	actual_random <- random_test$status
	pred.validated <- validated_data
	if (format == "PCA"){pred.validated <- predict(model, validated_data)}
	#random forest 500 trees                  
	enhancer_rf_500 <- randomForest(status ~., data = enhancers_training[,2:length(enhancers_training)],  ntree = 500, importance = TRUE)
	random_rf_500 <- randomForest(status ~., data = random_training[,2:length(random_training)],  ntree = 500, importance = TRUE)
	if (pred.validated == "none"){
	prediction_validated <- "none"
	prediction_validated_on_frame <- "none"}
	else{
	prediction_validated <- predict(enhancer_rf_500, newdata = pred.validated, type = "class")
	prediction_validated_on_frame <- cbind(validated_data[,1:2], prediction_validated)}
	outfilename <- paste("prediction_validated", info, toString(count),".tsv", sep = "_")
    write.table(prediction_validated_on_frame, outfilename , sep = "\t")
    import <- importance(enhancer_rf_500, type = 1)
	imp1 <- as.vector(import[,1])
	if (count == 1){mydataframe <- import}
	else {mydataframe <- cbind(mydataframe, imp1)}
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
	random_rf_500_predicted <- predict(random_rf_500, newdata=random_test, type="class")
	random_rf_500_test_table <- table(actual = random_test$status, predicted = random_rf_500_predicted)
	random_rf_500_prediction_success <- 100 * sum(diag(random_rf_500_test_table)/sum(random_rf_500_test_table))
	random_rf_500_actual_on <- 100 * random_rf_500_test_table[2,2]/sum(random_rf_500_test_table[2,])
	random_rf_500_false_on <- 100 * random_rf_500_test_table[1,2]/sum(random_rf_500_test_table[,2])
	random_rf_500_true_off <- 100 *random_rf_500_test_table[1,1]/sum(random_rf_500_test_table[1,])
	random_rf_500_false_off <- 100 * random_rf_500_test_table[2,1]/sum(random_rf_500_test_table[,1])
	random_rf_500_false_pos <- 100 * random_rf_500_test_table[1,2]/sum(random_rf_500_test_table[1,])
	random_predictions=as.vector(random_rf_500$votes[,2])
	random_pred1=prediction(random_predictions,random_training$status)
	random_perf1 <- performance(random_pred1, "prec", "rec")
	random_perf_ROC=performance(random_pred1,"tpr","fpr")
	random_perf_AUC=performance(random_pred1,"auc")
	random_AUC <- random_perf_AUC@y.values[[1]]
	random_precision <- performance(random_pred1, "prec")@y.values[[1]]
	random_precision[1] <- 0
	random_recall <- performance(random_pred1, "rec")@y.values[[1]]
	random_prec_reca <- trapz(random_recall, random_precision)
	if (count == 1 && graph == "recall") {plot(perf1, col = "darkgreen", main = "Precision Recall", ylim = c(0,1))} else if (count !=1 && graph == "recall") {plot(perf1, add = TRUE, col = "darkgreen")}
	if (count == 1 && graph == "ROC") {plot(perf_ROC, col = "darkgoldenrod", main = "ROC plot")} else if (count !=1 && graph == "ROC") {plot(perf_ROC, add = TRUE, col = "darkgoldenrod")}
	if (graph == "recall") {plot(random_perf1, add = TRUE, lty = 3, col = "darkgray")}
	if (graph == "ROC") {plot(random_perf_ROC, add = TRUE, lty = 3, col = "darkgray")}
	write.csv(enhancers_test, file = summaryfile)
	#close(filename)
	return_list <- list(rf_500_prediction_success, rf_500_actual_on- random_rf_500_actual_on, rf_500_true_off-random_rf_500_true_off, rf_500_false_on-random_rf_500_false_on, rf_500_false_off-random_rf_500_false_off, AUC, prec_reca, AUC - random_AUC, prec_reca - random_prec_reca, rf_500_prediction_success-random_rf_500_prediction_success, mydataframe)
	return(return_list)             
}      
    


#Activity in Stage 4-6 embryos
all_data <- read.table("all_data.tsv", header = TRUE)
all_data <- all_data[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]


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

#plots of PC 1 vs PC 2, PC 6 vs PC 7, PC 1 vs PC 6
autoplot(ir.pca, choices = c(1,2), data = all_data,colour = "status")  + scale_color_manual(values = c("darkgray","cadetblue4"))
tempplot("PC1", "PC2",ir.pca, choices = c(1,2), data = all_data,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm")  + scale_color_manual(values = c("coral","cadetblue4"))
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

#PCA all data bound by H3K4me1 2015
H3K4 <- subset(all_data, all_data$H3K4me1_2015_seq> (min(all_data$H3K4me1_2015_seq)))
H3K4$status <- factor(H3K4$status)
row.names(H3K4) <- NULL

#unbalanced
#all data
current_all <- all_data
current_tab <- table(all_data$status)
current_target <- read.table("target_file.tsv", header = TRUE)
prosp_target <- read.table("prospective_file.tsv", header = TRUE)
newdata <- current_target
all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData", description = "allData", "allData", 10, 250, strata_Percent_var, "ROC", current_target, "All_Data_Test",43, "NA", "NA")
all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_prosp", description = "allData_prosp", "allData_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "All_Data_Test_prosp",43, "NA", "NA")
all_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

#H3K4me1
current_all <- H3K4 
current_tab <- table(current_all$status)
current_target <- read.table("target_file.tsv", header = TRUE)
prosp_target <- read.table("prospective_file.tsv", header = TRUE)

H3K4.pca <- prcomp(current_all[,3:43],center = TRUE,scale. = TRUE) 
screeplot(H3K4.pca, type = "l", main = "Scree plot PCA H3K4me1 bound reporters, 41 features", col = "cadetblue")
tempplot("PC1", "PC2",H3K4.pca, choices = c(1,2), data = H3K4,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm")  + scale_color_manual(values = c("coral","cadetblue4"))

fviz_contrib(H3K4.pca, choice = "var", fill = "cadetblue", axes = 1, top = 10)
fviz_contrib(H3K4.pca, choice = "var", fill = "cadetblue", axes = 2, top = 10)

H3K4_ROC_data <- resample_rf500(enhancer_class_status_balance, current_all, current_tab, "H3K4 ", description = "H3K4 ", "H3K4", 20, 150, strata_Percent_var, "ROC", current_target, "H3K4_Test",43, "NA", "NA")
H3K4_means <- apply(H3K4_ROC_data[[1]], MARGIN = 2, FUN = mean)
H3K4_sd <- apply(H3K4_ROC_data[[1]], MARGIN=2, FUN=sd)
H3K4_PRROC_data <- resample_rf500(enhancer_class_status_balance, current_all, current_tab, "H3K4_prosp", description = "H3K4_prosp", "H3K4_prosp", 10, 150, strata_Percent_var, "recall", current_target, "H3K4_Test_prosp",43, "NA", "NA")
all_imp <- data.frame(H3K4_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

screeplot(ir.pca, type = "l", main = "Scree plot PCA, 41 features")
lines(H3K4.pca$sdev^2[1:10], type = "b", col = "cadetblue")

all_data_lda <- lda(train$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = train,prior = rep(1, 2)/2, CV=TRUE)
all.data.lda.values <- predict(all_data_lda, test)
all.data.class <- predict(all_data_lda)$class
plot(all.data.lda.values$x[,1], ylab=c("LDA Axis 1"), col = test$status, pch = 19)

H3K4_data_lda <- lda(H3K4$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = H3K4,prior = rep(1, 2)/2, CV=TRUE)
H3K4.data.lda.values <- predict(H3K4_data_lda, H3K4)
H3K4.data.class <- predict(H3K4_data_lda)$class
plot(H3K4.data.lda.values$x[,1], ylab=c("LDA Axis 1"), col = H3K4$status, pch = 19)

H3K4_data_lda2 <- lda(H3K4$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = H3K4,prior = rep(1, 2)/2)


smp_size <- floor(0.66 * nrow(all_data))

## set the seed to make your partition reproductible
set.seed(123)
train_ind <- sample(seq_len(nrow(all_data)), size = smp_size)

train <- all_data[train_ind, ]
test <- all_data[-train_ind, ]

H3K4_data_lda <- lda(H3K4$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = H3K4,prior = rep(1, 2)/2, CV=TRUE)
H3K4.data.lda.values <- predict(H3K4_data_lda, test)
conCV1 <- rbind(tab[1, ]/sum(tab[1, ]), tab[2, ]/sum(tab[2, ]))
dimnames(conCV1) <- list(Actual = c("No", "Yes"), "Predicted (cv)" = c("No",
+     "Yes"))
> print(round(conCV1, 3))
 Pima.lda <- lda(type ~ ., data = Pima.tr)
> Pima.hat <- predict(Pima.lda)
> tabtrain <- table(Pima.tr$type, Pima.hat$class)

Rev_1 <- predict(H3K4_data_lda, all_data)
Rev1.class <- Rev_1$class
plot(Rev_1$x[,1], ylab=c("LDA Axis 1"), col = all_data$status, pch = 19)

plot(H3K4$Zelda_2011_seq, H3K4$Hunchback_2009_chip,
     xlab="Zelda", ylab="Hunchback", 
     main="Cross-Validated LD Classification", 
     pch=as.numeric(H3K4$status)+ 19,
     col=as.numeric(H3K4_data_lda$class)) 
	
all_data_lda <- lda(all_data$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = all_data,prior = rep(1, 2)/2, CV=TRUE)

plot(all_data$Zelda_2011_seq, all_data$H3K4me1_2015_seq,
     xlab="Zelda", ylab="H3K4me1", 
     main="Cross-Validated LD Classification", 
     pch=as.numeric(all_data$status) + 19,
     col=as.numeric(all_data_lda$class)) 
	 
hist(all_data$H3K4me1_2015_seq)
	 
ggplot(all.data.lda.values, aes(x[,1] fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
ggplot(all_data, aes(Dsim, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())


Rev_2 <- predict(all_data_lda, H3K4)
Rev2.class <- Rev_2$class
plot(Rev_2$x[,1], ylab=c("LDA Axis 1"), col = H3K4$status, pch = 19)

ggplot(all_data, aes(Dsim, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

all_data_lda <- lda(all_data$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = all_data,prior = rep(1, 2)/2, CV=TRUE)
all_data_lda2 <- lda(all_data$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = all_data,prior = rep(1, 2)/2)

plot(all_data$Zelda_2011_seq, all_data$H3K4me1_2015_seq,
     xlab="Zelda", ylab="H3K4me1", 
     main="Cross-Validated LD Classification", 
     pch=as.numeric(all_data_lda2$status),
     col=as.numeric(all_data_lda$class)) 
	 
g <- lda(status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = all_data,prior = rep(1, 2)/2)
hab.lda.values <- predict(g, all_data)
hab.class <- predict(g)$class
plot(hab.lda.values$x[,1], ylab=c("LDA Axis 1"), col = all_data$status, pch = 19)

ldahist(data = hab.lda.values$x[,1], g=all_data$status,col = "cadetblue4")

	 
gg <- lda(status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = H3K4,prior = rep(1, 2)/2)
gg.hab.lda.values <- predict(gg, H3K4)
gg.hab.class <- predict(gg)$class
plot(gg.hab.lda.values$x[,1], ylab=c("LDA Axis 1"), col = H3K4$status, pch = 19)

ldahist(data = gg.hab.lda.values$x[,1], g= H3K4$status,col = "seagreen")

all_data_lda <- lda(all_data$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = all_data,prior = rep(1, 2)/2, CV=TRUE)

plot(all_data$Zelda_2011_seq, all_data$H3K4me1_2015_seq,
     xlab="Zelda", ylab="H3K4me1", 
     main="Cross-Validated LD Classification", 
     pch=as.numeric(all_data$status) + 19,
     col=as.numeric(all_data_lda$class)) 

hist(as.numeric(all_data$status), col = all_data_lda$class)

dat <- rbind(as.numeric(all_data_lda$class), as.numeric(all_data$status))


library(reshape2)

dat$row <- seq_len(nrow(dat))
dat2 <- melt(dat, id.vars = "row")

library(ggplot2)

ggplot(dat2, aes(x=variable, y=value, fill=row)) + 
  geom_bar(stat="identity") +
  xlab("\nType") +
  ylab("Time\n") +
  guides(fill=FALSE) +
  theme_bw()

  
  ggboxplot(all_data_lda$class, x = "LDA classification", y = "frequency",
    color = all_data$status,
    add = "jitter", shape = all_data$status)
	
dat <- read.table(text = "AllDataOn	AllDataOff	H3K4me1on	H3K4me1off
236	6588	118	361
156	270	109	49", header = TRUE)

library(reshape2)

dat$row <- seq_len(nrow(dat))
dat2 <- melt(dat, id.vars = "row")

library(ggplot2)

ggplot(dat2, aes(x=variable, y=value, fill=row)) + 
  geom_bar(stat="identity") +
  xlab("\nLDA") +
  ylab("Frequency\n") +
  guides(fill=FALSE) +
  theme_bw()