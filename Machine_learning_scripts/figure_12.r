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
long_name <- c("enhancer","status","Zelda_2011_seq","Dorsal_2015_seq","Dorsal_2009_chip","Snail_2014_chip","Snail_2009_chip","Twist_2011_seq","Twist_2014_chip","Twist_2009_chip","Bicoid_2013_seq","Bicoid_2009_chip","Caudal_2010_seq","Caudal_2009_chip","Hunchback_2013_seq","Hunchback_2009_chip","Giant_2013_seq","Giant_2009_chip","Kruppel_2013_seq","Knirps_2010_seq","Hairy_2009_chip","H3K27ac_2015_seq","H3K27ac_2010_seq","H3K4me1_2015_seq","p300_2010_seq", "Zld_motif_sanger", "Zld_motif_solexa", "Dorsal_motif_FlyReg","Dorsal_motif_NBT","Snail_motif_FlyReg","Snail_motif_Sanger","Snail_motif_solexa","Twist_motif_FlyReg","Twist_motif_da","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")
short_name <- c("enhancer","status","Zld_11_s","Dl_15_s","Dl_09_c","Sna_14_c","Sna_09_c","Twi_11_s","Twi_14_c","Twi_09_c","Bcd_13_s","Bcd_09_c","Cad_10_s","Cad_09_c","Hb_13_s","Hb_09_c","Gt_13_s","Gt_09_c","Kr_13_s","Kni_10_s","Hry_09_c","H3K27ac_15","H3K27ac_10","H3K4me1_15","p300_10", "Zld_m1", "Zld_m2", "Dl_m1","Dl_m2","Sna_m1","Sna_m2","Sna_m3","Twi_m1","Twi_m2","Dsim","Dsec","Dyak","Dere","Dana","Dpse","Dwil","Dvir","Dgri")

#Filtered reporters
Bcd13_train_cons <- read.table("cons_train_Bicoid_2013_seq.tsv", header = TRUE)
Bcd13_train_all <- read.table("all_train_Bicoid_2013_seq.tsv", header = TRUE)
Zld_train_cons <- read.table("cons_train_Zelda_2011_seq.tsv", header = TRUE)
Zld_train_all <- read.table("all_train_Zelda_2011_seq.tsv", header = TRUE)


#Random data
random_data <- data.frame(replicate(41,sample(0:1000,7250,rep=TRUE)))

random_data <-cbind(all_data[,0:2],random_data)
random_balance <- data.frame(replicate(41,sample(0:1000,600,rep=TRUE)))
name_column <- data.frame(replicate(1,sample(c("on", "off"),600, rep = TRUE)))
random_balance <-cbind(all_data[1:600,0:1], name_column, random_balance)
colnames(random_data) <- long_name
colnames(random_balance) <- long_name

random_bcd <- data.frame(replicate(41, sample(0:1000, 371, rep = TRUE)))
random_bcd <- cbind(Bcd13_train_all[,0:2], random_bcd)
colnames(random_bcd) <- long_name

random_bcd_bal <- data.frame(replicate(41, sample(0:1000, 200, rep = TRUE)))
random_bcd_bal <- cbind(Bcd13_train_cons[,0:2], random_bcd_bal)
colnames(random_bcd_bal) <- long_name

random_zld <- data.frame(replicate(41, sample(0:1000, 1232, rep = TRUE)))
random_zld <- cbind(Zld_train_all[,0:2], random_zld)
colnames(random_zld) <- long_name

random_zld_bal <- data.frame(replicate(41, sample(0:1000, 500, rep = TRUE)))
random_zld_bal <- cbind(Zld_train_cons[,0:2], random_zld_bal)
colnames(random_zld_bal) <- long_name

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
current_target <- read.table("target_file2.tsv", header = TRUE)
prosp_target <- read.table("prospective_file2.tsv", header = TRUE)
newdata <- current_target


#Random forest predictions based on all data, scaled
#Figure 7, Supplementary figure 9
#unbalanced

current_all <- all_data
current_tab <- table(all_data$status)


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

dotchart(alldata[31:41,1], xlab = "Importance", xlim = c(10,30), cex = 0.9, pch = 19)
segments(alldata[31:41,1]-alldata[31:41,2], 1:11, alldata[31:41,1]+alldata[31:41,2], 1:11)





#balanced
library(dplyr)
all_data_balance <- all_data %>%
group_by(status) %>%
sample_n(size = 300)
current_all <- all_data_balance
current_tab <- table(all_data_balance$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",43, "NA", "NA")
all_bal_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance_prosp", description = "allData_balance_prosp", "allData_balance_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "allData_balance_Test_prosp",43, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)


current_all <- Bcd13_train_all
current_tab <- table(current_all$status)
Bcd13_all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Bcd13_train_all", description = "Bcd13_train_all", "Bcd13_train_all", 10, 250, strata_Percent_var, "ROC", current_target, "Bcd13_train_all_Test",43, "NA", "NA")
Bcd13_all_means <- apply(Bcd13_all_ROC_data[[1]], MARGIN = 2, FUN = mean)
Bcd13_all_sd <- apply(Bcd13_all_ROC_data[[1]], MARGIN=2, FUN=sd)
Bcd13_all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Bcd13_all_prosp", description = "Bcd13_all_prosp", "Bcd13_all_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Bcd13_all_prosp",43, "NA", "NA")
all_imp <- data.frame(Bcd13_all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all <- Bcd13_train_cons
current_tab <- table(current_all$status)
Bcd13_cons_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Bcd13_train_cons", description = "Bcd13_train_cons", "Bcd13_train_cons", 10, 250, strata_Percent_var, "ROC", current_target, "Bcd13_train_cons_Test",43, "NA", "NA")
Bcd13_cons_means <- apply(Bcd13_cons_ROC_data[[1]], MARGIN = 2, FUN = mean)
Bcd13_cons_sd <- apply(Bcd13_cons_ROC_data[[1]], MARGIN=2, FUN=sd)
Bcd13_cons_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Bcd13_cons_prosp", description = "Bcd13_cons_prosp", "Bcd13_cons_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Bcd13_cons_prosp",43, "NA", "NA")
all_imp <- data.frame(Bcd13_cons_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)





current_all <- random_data
current_tab <- table(current_all$status)
Random_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "random_data", description = "random_data", "random_data", 10, 250, strata_Percent_var, "ROC", current_target, "random_data_Test",43, "NA", "NA")
Random_means <- apply(Random_ROC_data[[1]], MARGIN = 2, FUN = mean)
Random_sd <- apply(Random_ROC_data[[1]], MARGIN=2, FUN=sd)
Random_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Random_prosp", description = "Random_prosp", "Random_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Random_prosp",43, "NA", "NA")
all_imp <- data.frame(Random_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)



current_all <- random_balance
current_tab <- table(current_all$status)
Random_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "random_balance", description = "random_balance", "random_balance", 10, 250, strata_Percent_var, "ROC", current_target, "random_balance_Test",43, "NA", "NA")
Random_bal_means <- apply(Random_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
Random_bal_sd <- apply(Random_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
Random_bal_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Random_bal_prosp", description = "Random_bal_prosp", "Random_bal_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Random_bal_prosp",43, "NA", "NA")
all_imp <- data.frame(Random_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all <- random_bcd
current_tab <- table(current_all$status)
Random_bcd_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "random_bcd", description = "random_bcd", "random_bcd", 10, 250, strata_Percent_var, "ROC", current_target, "random_balance_Test",43, "NA", "NA")
Random_bcd_means <- apply(Random_bcd_ROC_data[[1]], MARGIN = 2, FUN = mean)
Random_bcd_sd <- apply(Random_bcd_ROC_data[[1]], MARGIN=2, FUN=sd)
Random_bcd_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Random_bcd_prosp", description = "Random_bcd_prosp", "Random_bcd_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Random_bcd_prosp",43, "NA", "NA")
all_imp <- data.frame(Random_bcd_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all <- random_bcd_bal
current_tab <- table(current_all$status)
Random_bcd_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "random_bcdbalance", description = "random_bcdbalance", "random_bcdbalance", 10, 250, strata_Percent_var, "ROC", current_target, "random_bcdbalance_Test",43, "NA", "NA")
Random_bcd_means <- apply(Random_bcd_ROC_data[[1]], MARGIN = 2, FUN = mean)
Random_bcd_sd <- apply(Random_bcd_ROC_data[[1]], MARGIN=2, FUN=sd)
Random_bcd_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Random_bcdbal_prosp", description = "Random_bcdbal_prosp", "Random_bcdbal_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Random_bcdbal_prosp",43, "NA", "NA")
all_imp <- data.frame(Random_bcd_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all <- Zld_train_all
current_tab <- table(current_all$status)
Zld_all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Zld_train_all", description = "Zld_train_all", "Zld_train_all", 10, 250, strata_Percent_var, "ROC", current_target, "Zld_train_all_Test",43, "NA", "NA")
Zld_all_means <- apply(Zld_all_ROC_data[[1]], MARGIN = 2, FUN = mean)
Zld_all_sd <- apply(Zld_all_ROC_data[[1]], MARGIN=2, FUN=sd)
Zld_all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Zld_all_prosp", description = "Zld_all_prosp", "Zld_all_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Zld_all_prosp",43, "NA", "NA")
all_imp <- data.frame(Zld_all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all <- Zld_train_cons
current_tab <- table(current_all$status)
Zld_cons_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Zld_train_cons", description = "Zld_train_cons", "Zld_train_cons", 10, 250, strata_Percent_var, "ROC", current_target, "Zld_train_cons_Test",43, "NA", "NA")
Zld_cons_means <- apply(Zld_cons_ROC_data[[1]], MARGIN = 2, FUN = mean)
Zld_cons_sd <- apply(Zld_cons_ROC_data[[1]], MARGIN=2, FUN=sd)
Zld_cons_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Zld_cons_prosp", description = "Zld_cons_prosp", "Zld_cons_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Zld_cons_prosp",43, "NA", "NA")
all_imp <- data.frame(Zld_cons_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all <- random_zld_bal
current_tab <- table(current_all$status)
Random_zld_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "random_zldbalance", description = "random_zldbalance", "random_zldbalance", 10, 250, strata_Percent_var, "ROC", current_target, "random_zldbalance_Test",43, "NA", "NA")
Random_zld_means <- apply(Random_zld_ROC_data[[1]], MARGIN = 2, FUN = mean)
Random_zld_sd <- apply(Random_zld_ROC_data[[1]], MARGIN=2, FUN=sd)
Random_zld_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Random_zldbal_prosp", description = "Random_zldbal_prosp", "Random_zldbal_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Random_zldbal_prosp",43, "NA", "NA")
all_imp <- data.frame(Random_zld_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)