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


#Expression patterns
ACP <- read.table("ACP.tsv", header = TRUE)
AP_DV <- read.table("AP_DV.tsv", header = TRUE)
current_target <- read.table("target_file.tsv", header = TRUE)

all_data <- read.table("all_data.tsv", header = TRUE)
all_data <- all_data[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]
all_tab <- table(current_all$status)


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
   
ACP <- ACP[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")] 
current_all <- ACP
ACP <- ACP[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]
current_all <- current_all[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","H3K27ac_2015_seq","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]
current_tab <- table(current_all$status)
current_target <-  current_target[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]


random_data <- data.frame(replicate(41,sample(0:1000,192,rep=TRUE)))
names(random_data) <- names(ACP)[3:43]
random_data <-cbind(ACP[,0:2],random_data)


#All data, ACP
all_features <- c(3:43)
current_all <- ACP
current_tab <- table(current_all$status)
my_features <- all_features
all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "ROC", current_target, "All_Data_Reduced", 43, "NA", "NA")
all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "recall", current_target, "All_Data_Reduced", 43, "NA", "NA")

all_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]

dotchart(alldata[,1], xlab = "Importance-ACP", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

#All data Reduced, ACP
current_all <- ACP
current_tab <- table(current_all$status)
current_target <- read.table("target_file.tsv", header = TRUE)
current_target <-  current_target[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","H3K27ac_2015_seq","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]
my_features <- all_features

df2 = cor(current_all[,3:43])
hc = findCorrelation(df2, cutoff=0.7) # put any value as a "cutoff" 
hc = sort(hc)
reduced_Data = current_all[,-c(hc)]
#search here
current_all <- reduced_Data
current_tab <- table(reduced_Data$status)
reduced_target <- current_target[,-c(hc)]

all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "ROC", reduced_target, "All_Data_Reduced", 39, "NA", "NA")
all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "recall", reduced_target, "All_Data_Reduced", 39, "NA", "NA")

Reduced_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
Reduced_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]

dotchart(alldata[,1], xlab = "Importance-ACP", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

#All data Reduced More, ACP
current_all <- ACP
current_tab <- table(current_all$status)
current_target <- read.table("target_file.tsv", header = TRUE)
current_target <-  current_target[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","H3K27ac_2015_seq","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]
my_features <- all_features

df2 = cor(current_all[,3:43])
hc = findCorrelation(df2, cutoff=0.5) # put any value as a "cutoff" 
hc = sort(hc)
reduced_Data_more = current_all[,-c(hc)]
#search here
current_all <- reduced_Data_more
current_all <- cbind(ACP[,1], current_all)
current_tab <- table(reduced_Data$status)
reduced_target_more <- current_target[,-c(hc)]
reduced_target_more <- cbind(current_target[,1], reduced_target_more)
all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "ROC", reduced_target_more, "All_Data_Reduced", 26, "NA", "NA")
all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "recall", reduced_target_more, "All_Data_Reduced", 26, "NA", "NA")

Reduced_means_more <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
Reduced_sd_more <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]

dotchart(alldata[,1], xlab = "Importance-ACP", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)


#Random data, ACP
current_all <- random_data
current_tab <- table(current_all$status)
colnames(current_target) <- colnames(current_all)
all_features <- c(3:43)
my_features <- all_features
all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "ROC", current_target, "All_Data_Reduced", 43, "NA", "NA")
#all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "recall", current_target, "All_Data_Reduced", 43, "NA", "NA")

random_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
random_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]

dotchart(alldata[,1], xlab = "Importance-ACP", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

ACP_means <- cbind(all_means, Reduced_means, Reduced_means_more, PCA25_means, PCA10_means, random_means)

#PCA
ir.pca <- prcomp(ACP[,3:43],center = TRUE,scale. = TRUE) 
screeplot(ir.pca, type = "l", main = "Scree plot ACP, 41 features")
tempplot("PC1", "PC2",ir.pca, choices = c(1,2), data = ACP,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") 

target_PCA <- prcomp(current_target[,3:43], center = TRUE, scale = TRUE)
PCA_target25 <- data.frame(cbind(current_target[,1:2], target_PCA$x[,0:25]))
PCA_target10 <- data.frame(cbind(current_target[,1:2], target_PCA$x[,0:10]))

PCA_ACP <- data.frame(cbind(ACP[,1:2],ir.pca$x[,0:10]))
PCA_ACP_25 <- data.frame(cbind(ACP[,1:2],ir.pca$x[,0:25]))

current_all <- PCA_ACP_25
current_tab <- table(PCA_ACP_25$status)
current_target <- PCA_target25
all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "ROC", current_target, "All_Data_Reduced", 27, "NOT", "NA")

PCA25_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
PCA25_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)



current_target <-PCA_target10
current_all <- PCA_ACP
current_tab <- table(PCA_ACP$status)

all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allDataReduced", description = "allDataReduced", "allDataReduced", 10, 250, strata_balance, "ROC", current_target, "All_Data_Reduced", 12, "NOT", "NA")

PCA10_means <- apply(all_ROC_data[[1]], MARGIN = 2, FUN = mean)
PCA10_sd <- apply(all_ROC_data[[1]], MARGIN=2, FUN=sd)

all_imp <- data.frame(all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)


ACP_means <- cbind(all_means, Reduced_means, Reduced_means_more, PCA25_means, PCA10_means, random_means)
ACP_sd <- cbind(all_sd, Reduced_sd, Reduced_sd_more, PCA25_sd, PCA10_sd, random_sd)
colnames(ACP_means) <- c("All data", "Reduced_1", "Reduced_2", "25 PC", "10 PC", "Random")


ACP_plot <- barplot(ACP_means[6:7,], beside = TRUE, col = c("darkgoldenrod", "darkgreen"), ylim = c(0,1), cex.lab = 1.5, cex.axis = 1.5, cex.names = 1)

segments(ACP_plot, ACP_means[6:7,] - ACP_sd[6:7,] * 2, ACP_plot,
         ACP_means[6:7,] + ACP_sd[6:7,] * 2, lwd = 1.5)

arrows(ACP_plot, ACP_means[6:7,] - ACP_sd[6:7,] * 2, ACP_plot,
         ACP_means[6:7,] + ACP_sd[6:7,] * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

 legend("topright", c("ROC", "Precision-Recall"), cex = 1.5, fill = c("darkgoldenrod", "darkgreen"))


ACP_models_balanced <- c(all_means[1]/100, all_means[6] + 1,all_means[7] + 2, 
		   all_reduced_means[1]/100, all_reduced_means[6]   + 1,all_reduced_means[7]   + 2,all_PCA_25[1]/100,all_PCA_25[6]   + 1,all_PCA_25[7]   + 2,
		   all_PCA_10[1]/100,all_PCA_10[6]   + 1,all_PCA_10[7])

   + 2,random_all[1]/100,random_all[6]   + 1,random_all[7]   + 2
xnames <- c(1,1,1,3,3,3,5,5,5,7,7,7,9,9,9,11,11,11,13,13,13,15,15,15,17,17,17,19,19,19,21,21,21,23,23,23,25,25,25)
plot(xnames, ACP_models_balanced, xaxt = "none", col = "black", pch = 19, ylim = c(0.1, 2.9),yaxt = "none", xlab = "")
abline(h=1,col=1,lty=1)
abline(h=2,col=1,lty=1)
abline(h=0.5,col="gray",lty=2)
abline(h=1.5,col="gray",lty=2)
abline(h=2.5,col="gray",lty=2)
abline(v=2,col="darkgray",lty=3)
abline(v=4,col="darkgray",lty=3)
abline(v=6,col="darkgray",lty=3)
abline(v=8,col="darkgray",lty=3)
abline(v=10,col="darkgray",lty=3)
abline(v=0,col="darkgray",lty=3)
abline(v=12,col="darkgray",lty=3)
abline(v=14,col="darkgray",lty=3)
abline(v=16,col="darkgray",lty=3)
abline(v=18,col="darkgray",lty=3)
abline(v=20,col="darkgray",lty=3)
abline(v=22,col="darkgray",lty=3)
abline(v=24,col="darkgray",lty=3)
abline(v=26,col="darkgray",lty=3)

lab = c("all_data_activity", "all_data_Expression", "Caudal_2009_chip", "Bicoid_2009_chip", "Bicoid_2013_seq","Caudal_2010_seq", "Knirps_2010_seq", "Hunchback_2013_seq", "Kruppel_2013_seq", "Reduced_data", "PCA_25", "PCA_10", "Random_data")
text(x=c(1,3,5,7,9,11,13,15,17,19,21,23,25), y=par()$usr[3] - 0.01*(par()$usr[4]-par()$usr[3]),labels = lab,srt=65, adj=1, xpd=TRUE, cex=0.6)

#PCA plots
 



fviz_contrib(ir.pca, choice = "var", fill = "cadetblue", axes = 5, top = 10)
fviz_contrib(ir.pca, choice = "var", fill = "cadetblue", axes = 13, top = 10)

ggplot(ACP, aes(Bicoid_2013_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
ggplot(ACP, aes(Caudal_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
ggplot(ACP, aes(Giant_2009_chip, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

AP_DV <- AP_DV[c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq","Zld_motif_sanger", "Zld_motif_solexa","Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","H3K27ac_2015_seq","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")]

AP_DV.pca <- prcomp(AP_DV[,3:43],center = TRUE,scale. = TRUE) 
screeplot(ir.pca, type = "l", main = "Scree plot AP_DV, 41 features")

tempplot("PC1", "PC2",ir.pca, choices = c(1,2), data = AP_DV,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") 

ggplot(AP_DV, aes(Twist_2011_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

ir.pca <- prcomp(ACP[,3:43],center = TRUE,scale. = TRUE) 
screeplot(ir.pca, type = "l", main = "Scree plot AP_DV, 41 features")

tempplot("PC1", "PC2",ir.pca, choices = c(1,2), data = AP_DV,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") 

ggplot(AP_DV, aes(Twist_2011_seq, fill = status)) + geom_density(alpha = 0.5)+ theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

#unbalanced
#Activity in Stage 4-6 embryos
all_data <- read.table("all_data.tsv", header = TRUE)


Dl_chip <- subset(all_data, all_data$H3K4me1_2015_seq > (min(all_data$H3K4me1_2015_seq)))

current_all <- Dl_chip 

current_tab <- table(current_all$status)

ir.pca <- prcomp(current_all[,3:43],center = TRUE,scale. = TRUE) 
screeplot(ir.pca, type = "l", main = "Scree plot filtered, 41 features")

tempplot("PC1", "PC2",ir.pca, choices = c(1,2), data = current_all,colour = "status", alpha = 1, frame = TRUE, frame.type = "norm") 



Dlchip_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_chip", description = "Dl_chip", "Dl_chip", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_chip_Test",43, "NA", "NA")
Dlchip_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_chip_prosp", description = "Dl_chip_prosp", "Dl_chip_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Dl_chip_Test_prosp",43, "NA", "NA")
Dlchip_means <- apply(Dlchip_ROC_data[[1]], MARGIN = 2, FUN = mean)
Dlchip_sd <- apply(Dlchip_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(Dlchip_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)


#linear discriminant analysis
g <- lda(status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = ACP,prior = rep(1, 2)/2)
hab.lda.values <- predict(g, ACP)
hab.class <- predict(g)$class
plot(hab.lda.values$x[,1], ylab=c("LDA Axis 1"), col = ACP$status, pch = 19)
ldahist(data = hab.lda.values$x[,1], g=ACP$status)

g_CV <- lda(ACP$status ~ Zelda_2011_seq + Dorsal_2015_seq + Dorsal_2009_chip + Snail_2014_chip + Snail_2009_chip + Twist_2011_seq + Twist_2014_chip + Twist_2009_chip + Bicoid_2013_seq + Bicoid_2009_chip + Caudal_2010_seq + Caudal_2009_chip + Hunchback_2013_seq + Hunchback_2009_chip + Giant_2013_seq + Giant_2009_chip + Kruppel_2013_seq + Knirps_2010_seq + Hairy_2009_chip + H3K27ac_2015_seq + H3K27ac_2010_seq + H3K4me1_2015_seq + p300_2010_seq + Zld_motif_sanger + Zld_motif_solexa + Dorsal_motif_FlyReg + Dorsal_motif_NBT + Snail_motif_FlyReg + Snail_motif_Sanger + Snail_motif_solexa + Twist_motif_FlyReg + Twist_motif_da + Dsim + Dsec + Dyak + Dere + Dana + Dpse + Dwil + Dvir + Dgri, data = ACP,prior = rep(1, 2)/2, CV=TRUE)

plot(ACP$Zelda_2011_seq,
     xlab="Index", ylab="Zelda ChIP-seq", 
     main="Cross-Validated LD Classification", 
     pch=as.numeric(ACP$status) + 19,
     col=as.numeric(g_cv$class)) 

