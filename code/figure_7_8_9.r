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
colnames(all_data) <- short_name
#Filtered reporters
H3K4me1_train_cons <- read.table("cons_train_H3K4me1_2015_seq.tsv", header = TRUE)
H3K4me1_train_all <- read.table("all_train_H3K4me1_2015_seq.tsv", header = TRUE)
FTwi_train_cons <- read.table("cons_train_Twist_2014_chip.tsv", header = TRUE)
FTwi_train_all <- read.table("all_train_Twist_2014_chip.tsv", header = TRUE)
Zld_train_cons <- read.table("cons_train_Zelda_2011_seq.tsv", header = TRUE)
Zld_train_all <- read.table("all_train_Zelda_2011_seq.tsv", header = TRUE)
Dl15_train_cons <- read.table("cons_train_Dorsal_2015_seq.tsv", header = TRUE)
Dl15_train_all <- read.table("all_train_Dorsal_2015_seq.tsv", header = TRUE)

#Random data
random_data <- data.frame(replicate(41,sample(0:1000,7250,rep=TRUE)))

random_data <-cbind(all_data[,0:2],random_data)
random_balance <- data.frame(replicate(41,sample(0:1000,300,rep=TRUE)))
name_column <- data.frame(replicate(1,sample(c("on", "off"),300, rep = TRUE)))
random_balance <-cbind(all_data[1:300,0:1], name_column, random_balance)
colnames(random_data) <- long_name
colnames(random_balance) <- long_name

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



#Random forest predictions based on all data, scaled
#Figure 6, Supplementary figure 9
#unbalanced
current_all <- all_data
colnames(current_all) <- short_name
current_tab <- table(all_data$status)
current_target <- read.table("target_file.tsv", header = TRUE)
prosp_target <- read.table("prospective_file.tsv", header = TRUE)
colnames(prosp_target) <- short_name
newdata <- current_target
colnames(current_target) <- short_name
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
colnames(current_all) <- short_name
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

#Random forest predictions based on all data, using reporters bound by Dorsal chip-chip
colnames(all_data) <- long_name
Dl_chip <- subset(all_data, all_data$Dorsal_2009_chip > (min(all_data$Dorsal_2009_chip)))
Dl_chip_balance <- Dl_chip %>%
group_by(status) %>%
sample_n(size = 300)

#unbalanced
current_all <- Dl_chip 
colnames(current_all) <- long_name
colnames(current_target) <- long_name
colnames(prosp_target) <- long_name
current_tab <- table(current_all$status)
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

#balanced
current_all <- Dl_chip_balance 
current_tab <- table(current_all$status)
Dlchip_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_chip_balance", description = "Dl_chip_balance", "Dl_chip_balance", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_chip_balance_Test",43, "NA", "NA")
Dlchip_bal_means <- apply(Dlchip_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
Dlchip_bal_sd <- apply(Dlchip_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
Dlchip_bal_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_chip_balance_prosp", description = "Dl_chip_balance_prosp", "Dl_chip_balance_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Dl_chip_balance_Test_prosp",43, "NA", "NA")
all_imp <- data.frame(Dlchip_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)


#Random forest predictions based on all data, using reporters bound by Dorsal chip-chip AND Twist chip-chip
Dl_Twi <- subset(Dl_chip, Dl_chip$Twist_2009_chip > (min(Dl_chip$Twist_2009_chip)))
Dl_Twi_chip_balance <- Dl_Twi %>%
group_by(status) %>%
sample_n(size = 300)

#unbalanced
current_all <- Dl_Twi 
current_tab <- table(current_all$status)
DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",43, "NA", "NA")
DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
DlTwi_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi_prosp", description = "Dl_Twi_prosp", "Dl_Twi_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Dl_Twi_Test_prosp",43, "NA", "NA")
all_imp <- data.frame(DlTwi_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

#balanced
current_all <- Dl_Twi_chip_balance
current_tab <- table(current_all$status)
DlTwi_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi_chip_balance", description = "Dl_Twi_chip_balance", "Dl_Twi_chip_balance", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_chip_balance",43, "NA", "NA")
DlTwi_bal_means <- apply(DlTwi_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
DlTwi_bal_sd <- apply(DlTwi_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
DlTwi_bal_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi_balance_prosp", description = "Dl_Twi_balance_prosp", "Dl_Twi_balance_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Dl_Twi_balance_Test_prosp",43, "NA", "NA")
all_imp <- data.frame(DlTwi_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)





#Random forest predictions based on reporters bound by Dl chip-seq
current_all <- Dl15_train_all
current_tab <- table(current_all$status)
Dl15_all_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl15_train_all", description = "Dl15_train_all", "Dl15_train_all", 10, 250, strata_Percent_var, "ROC", current_target, "Dl15_train_all_Test",43, "NA", "NA")
Dl15_all_means <- apply(Dl15_all_ROC_data[[1]], MARGIN = 2, FUN = mean)
Dl15_all_sd <- apply(Dl15_all_ROC_data[[1]], MARGIN=2, FUN=sd)
Dl15_train_all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl15_train_all_prosp", description = "Dl15_train_all_prosp", "Dl15_train_all_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Dl15_train_all_prosp",43, "NA", "NA")
all_imp <- data.frame(Dl15_all_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all <- Dl15_train_cons
current_tab <- table(current_all$status)
Dl15_cons_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl15_train_cons", description = "Dl15_train_cons", "Dl15_train_cons", 10, 250, strata_Percent_var, "ROC", current_target, "Dl15_train_cons_Test",43, "NA", "NA")
Dl15_cons_means <- apply(Dl15_cons_ROC_data[[1]], MARGIN = 2, FUN = mean)
Dl15_cons_sd <- apply(Dl15_cons_ROC_data[[1]], MARGIN=2, FUN=sd)
Dl15_cons_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl15_cons_prosp", description = "Dl15_cons_prosp", "Dl15_cons_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Dl15_cons_prosp",43, "NA", "NA")
all_imp <- data.frame(Dl15_cons_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)


current_all <- H3K4me1_train_all
current_tab <- table(current_all$status)
H3K4me1_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "H3K4me1_train_all", description = "H3K4me1_train_all", "H3K4me1_train_all", 10, 250, strata_Percent_var, "ROC", current_target, "H3K4me1_train_all_Test",43, "NA", "NA")
H3K4me1_all_means <- apply(H3K4me1_ROC_data[[1]], MARGIN = 2, FUN = mean)
H3K4me1_all_sd <- apply(H3K4me1_ROC_data[[1]], MARGIN=2, FUN=sd)
H3K4me1_all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "H3K4me1_all_prosp", description = "H3K4me1_all_prosp", "H3K4me1_all_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "H3K4me1_all_prosp",43, "NA", "NA")
all_imp <- data.frame(H3K4me1_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)


current_all <- H3K4me1_train_cons
current_tab <- table(current_all$status)
H3K4me1_cons_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "H3K4me1_train_cons", description = "H3K4me1_train_cons", "H3K4me1_train_cons", 10, 250, strata_Percent_var, "ROC", current_target, "H3K4me1_train_cons_Test",43, "NA", "NA")
H3K4me1_cons_means <- apply(H3K4me1_cons_ROC_data[[1]], MARGIN = 2, FUN = mean)
H3K4me1_cons_sd <- apply(H3K4me1_cons_ROC_data[[1]], MARGIN=2, FUN=sd)
H3K4me1_cons_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "H3K4me1_cons_prosp", description = "H3K4me1_cons_prosp", "H3K4me1_cons_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "HH3K4me1_cons_prosp",43, "NA", "NA")
all_imp <- data.frame(H3K4me1_cons_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all <- FTwi_train_all
current_tab <- table(current_all$status)
Twi14_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "FTwi14_train_all", description = "FTwi14_train_all", "FTwi14_train_all", 10, 250, strata_Percent_var, "ROC", current_target, "FTwi14_train_all_Test",43, "NA", "NA")
Twi14_all_means <- apply(Twi14_ROC_data[[1]], MARGIN = 2, FUN = mean)
Twi14_all_sd <- apply(Twi14_ROC_data[[1]], MARGIN=2, FUN=sd
Twi14_all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Twi14_all_prosp", description = "Twi14_all_prosp", "Twi14_all_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Twi14_all_prosp",43, "NA", "NA")
all_imp <- data.frame(Twi14_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all <- FTwi_train_cons 
current_tab <- table(current_all$status)
Twi14_cons_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "FTwi14_train_cons", description = "FTwi14_train_cons", "FTwi14_train_cons", 10, 250, strata_Percent_var, "ROC", current_target, "FTwi14_train_cons_Test",43, "NA", "NA")
Twi14_cons_means <- apply(Twi14_cons_ROC_data[[1]], MARGIN = 2, FUN = mean)
Twi14_cons_sd <- apply(Twi14_cons_ROC_data[[1]], MARGIN=2, FUN=sd)
Twi14_all_PRROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Twi14_cons_prosp", description = "Twi14_cons_prosp", "Twi14_cons_prosp", 10, 250, strata_Percent_var, "recall", prosp_target, "Twi14_cons_prosp",43, "NA", "NA")
all_imp <- data.frame(Twi14_cons_ROC_data[2])
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

all_mean_random <- Random_means

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

#Figure 8

unbal <- c(all_means[1]/100,all_means[6] + 1,all_means[7] + 2, Dlchip_means[1]/100,Dlchip_means[6] + 1,Dlchip_means[7] + 2, 
          DlTwi_means[1]/100, DlTwi_means[6] + 1, DlTwi_means[7] + 2, Dl15_all_means[1]/100, Dl15_all_means[6] + 1,Dl15_all_means[7] + 2 , 
          Twi14_all_means[1]/100,Twi14_all_means[6] + 1, Twi14_all_means[7] + 2 ,Zld_all_means[1]/100, Zld_all_means[6] +1, Zld_all_means[7] + 2, 
          H3K4me1_all_means[1]/100, H3K4me1_all_means[6] + 1, H3K4me1_all_means[7] + 2 , Random_means[1]/100, Random_means[6] + 1, Random_means[7] + 2)

all_models <- c( Random_means[1]/100, Random_means[6]   + 1, Random_means[7]   + 2, Random_bal_means[1]/100, Random_bal_means[6]   + 1, Random_bal_means[7]   + 2,
          all_means[1]/100,all_means[6]   + 1,all_means[7]   + 2, all_bal_means[1]/100,all_bal_means[6]   + 1,all_bal_means[7]   + 2, 
		  Dlchip_means[1]/100,Dlchip_means[6]   + 1,Dlchip_means[7]   + 2, Dlchip_bal_means[1]/100,Dlchip_bal_means[6]   + 1,Dlchip_bal_means[7]   + 2,
          DlTwi_means[1]/100, DlTwi_means[6]   + 1, DlTwi_means[7]   + 2, DlTwi_bal_means[1]/100, DlTwi_bal_means[6]   + 1, DlTwi_bal_means[7]   + 2, 
          Twi14_all_means[1]/100,Twi14_all_means[6]   + 1, Twi14_all_means[7]   + 2 ,Twi14_cons_means[1]/100,Twi14_cons_means[6]   + 1, Twi14_cons_means[7]   + 2 ,
		  Zld_all_means[1]/100, Zld_all_means[6] +1, Zld_all_means[7]   + 2, Zld_cons_means[1]/100, Zld_cons_means[6] +1, Zld_cons_means[7]   + 2,
          H3K4me1_all_means[1]/100, H3K4me1_all_means[6]   + 1, H3K4me1_all_means[7]   + 2 , H3K4me1_cons_means[1]/100, H3K4me1_cons_means[6]   + 1, H3K4me1_cons_means[7]   + 2,
		  Dl15_all_means[1]/100, Dl15_all_means[6]   + 1,Dl15_all_means[7]   + 2 , Dl15_cons_means[1]/100, Dl15_cons_means[6]   + 1,Dl15_cons_means[7]   + 2 )

xnames <- c(1,1,1,3,3,3,7,7,7,9,9,9,13,13,13,15,15,15,19,19,19,21,21,21,25,25,25,27,27,27,31,31,31,33,33,33,37,37,37,39,39,39,43,43,43,45,45,45)
plot(xnames, all_models, xaxt = "none", pch = 19, col = c("black","black", "black","cyan4","cyan4","cyan4"), ylim = c(0.1, 2.9),xlim = c(1,45),yaxt = "none", xlab = "")
abline(h=1,col=1,lty=1)
abline(h=2,col=1,lty=1)
abline(h=0.5,col="gray",lty=2)
abline(h=1.5,col="gray",lty=2)
abline(h=2.5,col="gray",lty=2)
abline(v=5,col="darkgray",lty=3)
abline(v=11,col="darkgray",lty=3)
abline(v=17,col="darkgray",lty=3)
abline(v=23,col="darkgray",lty=3)
abline(v=29,col="darkgray",lty=3)
abline(v=35,col="darkgray",lty=3)
abline(v=41,col="darkgray",lty=3)

lab = c("random", "all_data", "Dl_2009", "Dl_Twi_2009", "Twi_2014", "Zld_2011", "H3K4me1_2015", "Dl_2015")

text(x=c(2,8,14,20,26,32,38,44), y=par()$usr[3] - 0.01*(par()$usr[4]-par()$usr[3]),labels = lab,srt=45, adj=1, xpd=TRUE)

all_models_sd <- c(all_sd[1]/100,all_sd[6] ,all_sd[7] , all_bal_sd[1]/100,all_bal_sd[6] ,all_bal_sd[7] , 
          Dlchip_sd[1]/100,Dlchip_sd[6] ,Dlchip_sd[7] , Dlchip_bal_sd[1]/100,Dlchip_bal_sd[6] ,Dlchip_bal_sd[7] ,
          DlTwi_sd[1]/100, DlTwi_sd[6] , DlTwi_sd[7] , DlTwi_bal_sd[1]/100, DlTwi_bal_sd[6] , DlTwi_bal_sd[7] , 
		  Dl15_all_sd[1]/100, Dl15_all_sd[6] ,Dl15_all_sd[7]  , Dl15_cons_sd[1]/100, Dl15_cons_sd[6] ,Dl15_cons_sd[7]  , 
          Twi14_all_sd[1]/100,Twi14_all_sd[6] , Twi14_all_sd[7]  ,Twi14_cons_sd[1]/100,Twi14_cons_sd[6] , Twi14_cons_sd[7]  ,
		  Zld_all_sd[1]/100, Zld_all_sd[6], Zld_all_sd[7] , Zld_cons_sd[1]/100, Zld_cons_sd[6], Zld_cons_sd[7] ,
          H3K4me1_all_sd[1]/100, H3K4me1_all_sd[6] , H3K4me1_all_sd[7]  , H3K4me1_cons_sd[1]/100, H3K4me1_cons_sd[6] , H3K4me1_cons_sd[7]  , 
		  Random_sd[1]/100, Random_sd[6] , Random_sd[7] , Random_bal_sd[1]/100, Random_bal_sd[6] , Random_bal_sd[7] )
		  
	
current_all <- Dl_Twi 
current_tab <- table(current_all$status)
DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",43, "NA", "NA")
DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(DlTwi_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)


current_all_drop <- current_all
current_all_drop$Zelda_2011_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",42, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop1 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Snail_2014_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",41, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop2 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Giant_2009_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",40, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop3 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Bicoid_2013_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",39, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop4 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Kruppel_2013_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",38, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop5 <- all_bal_means


current_all_drop <- current_all
current_all_drop$Bicoid_2009_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",37, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop6 <- all_bal_means


current_all_drop <- current_all
current_all_drop$Twist_2014_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",36, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop7 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Twist_2011_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",35, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop8 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Hairy_2009_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",34, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop9 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Hunchback_2009_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",33, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop10 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Knirps_2010_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",32, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop11 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Caudal_2010_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",31, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop12 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Twist_2009_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",30, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop13 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Caudal_2009_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",29, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop14 <- all_bal_means

current_all_drop <- current_all
current_all_drop$p300_2010_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",28, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop15 <- all_bal_means


current_all_drop <- current_all
current_all_drop$Giant_2013_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",27, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop16 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Hunchback_2013_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",26, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop17 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dorsal_2015_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",25, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop18 <- all_bal_means


current_all_drop <- current_all
current_all_drop$Zld_motif_sanger <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",24, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop19 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Snail_2009_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",23, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop20 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Zld_motif_solexa <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",22, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop21 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dere <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",21, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop22 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dsim <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",20, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop23 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dyak <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",19, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop24 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dorsal_2009_chip <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",18, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop24_5 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dana <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",17, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop25 <- all_bal_means

current_all_drop <- current_all
current_all_drop$H3K4me1_2015_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",16, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop26 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dana <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",15, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop27 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dpse <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",14, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop28 <- all_bal_means

current_all_drop <- current_all
current_all_drop$H3K27ac_2015_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",13, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop29 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Snail_motif_FlyReg <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",12, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop30 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dorsal_motif_NBT <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",11, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop31 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dsec <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",10, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop32 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dorsal_motif_FlyReg <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",9, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop33 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Twist_motif_da <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",8, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop34 <- all_bal_means


current_all_drop <- current_all
current_all_drop$Snail_motif_Sanger <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",7, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop35 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dgri <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",6, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop36 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Snail_motif_solexa <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",5, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop37 <- all_bal_means

current_all_drop <- current_all
current_all_drop$H3K27ac_2010_seq <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",4, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop38 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Twist_motif_FlyReg <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",3, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop39 <- all_bal_means

current_all_drop <- current_all
current_all_drop$Dvir <- NULL
current_all <- current_all_drop
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",2, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop40 <- all_bal_means
#top performing unbalanced model


current_all <- random_data
current_tab <- table(current_all$status)
Random_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "random_data", description = "random_data", "random_data", 10, 250, strata_Percent_var, "ROC", current_target, "random_data_Test",43, "NA", "NA")
Random_means <- apply(Random_ROC_data[[1]], MARGIN = 2, FUN = mean)
Random_sd <- apply(Random_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(Random_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_random <- Random_means


Dl_Twi_unbalanced <- c(DlTwi_means[1]/100,DlTwi_means[6]   + 1,DlTwi_means[7]   + 2, all_mean_drop1[1]/100,all_mean_drop1[6]   + 1,all_mean_drop1[7]   + 2,
		   all_mean_drop2[1]/100,all_mean_drop2[6]   + 1,all_mean_drop2[7]   + 2,all_mean_drop3[1]/100,all_mean_drop3[6]   + 1,all_mean_drop3[7]   + 2,
		   all_mean_drop4[1]/100,all_mean_drop4[6]   + 1,all_mean_drop4[7]   + 2,all_mean_drop5[1]/100,all_mean_drop5[6]   + 1,all_mean_drop5[7]   + 2,
		   all_mean_drop6[1]/100,all_mean_drop6[6]   + 1,all_mean_drop6[7]   + 2,all_mean_drop7[1]/100,all_mean_drop7[6]   + 1,all_mean_drop7[7]   + 2,
		   all_mean_drop8[1]/100,all_mean_drop8[6]   + 1,all_mean_drop8[7]   + 2,all_mean_drop9[1]/100,all_mean_drop9[6]   + 1,all_mean_drop9[7]   + 2,
		   all_mean_drop10[1]/100,all_mean_drop10[6]   + 1,all_mean_drop10[7]   + 2,all_mean_drop11[1]/100,all_mean_drop11[6]   + 1,all_mean_drop11[7] + 2,
		   all_mean_drop12[1]/100,all_mean_drop12[6]   + 1,all_mean_drop12[7]   + 2,
		   all_mean_drop13[1]/100,all_mean_drop13[6]   + 1,all_mean_drop13[7]   + 2,all_mean_drop14[1]/100,all_mean_drop14[6]   + 1,all_mean_drop14[7]   + 2,
		   all_mean_drop15[1]/100,all_mean_drop15[6]   + 1,all_mean_drop15[7]   + 2,all_mean_drop16[1]/100,all_mean_drop16[6]   + 1,all_mean_drop16[7]   + 2,
		   all_mean_drop17[1]/100,all_mean_drop17[6]   + 1,all_mean_drop17[7]   + 2,all_mean_drop18[1]/100,all_mean_drop18[6]   + 1,all_mean_drop18[7]   + 2,
		   all_mean_drop19[1]/100,all_mean_drop19[6]   + 1,all_mean_drop19[7]   + 2,all_mean_drop20[1]/100,all_mean_drop20[6]   + 1,all_mean_drop20[7]   + 2,
		   all_mean_drop21[1]/100,all_mean_drop21[6]   + 1,all_mean_drop21[7]   + 2,all_mean_drop22[1]/100,all_mean_drop22[6]   + 1,all_mean_drop22[7]   + 2,
		   all_mean_drop23[1]/100,all_mean_drop23[6]   + 1,all_mean_drop23[7]   + 2,all_mean_drop24[1]/100,all_mean_drop24[6]   + 1,all_mean_drop24[7]   + 2,all_mean_drop24_5[1]/100,all_mean_drop24_5[6]   + 1,all_mean_drop24_5[7]   + 2,
		   all_mean_drop25[1]/100,all_mean_drop25[6]   + 1,all_mean_drop25[7]   + 2,all_mean_drop26[1]/100,all_mean_drop26[6]   + 1,all_mean_drop26[7]   + 2,
		   all_mean_drop27[1]/100,all_mean_drop27[6]   + 1,all_mean_drop27[7]   + 2,all_mean_drop28[1]/100,all_mean_drop28[6]   + 1,all_mean_drop28[7]   + 2,
		   all_mean_drop29[1]/100,all_mean_drop29[6]   + 1,all_mean_drop29[7]   + 2,all_mean_drop30[1]/100,all_mean_drop30[6]   + 1,all_mean_drop30[7]   + 2,
		   all_mean_drop31[1]/100,all_mean_drop31[6]   + 1,all_mean_drop31[7]   + 2,all_mean_drop32[1]/100,all_mean_drop32[6]   + 1,all_mean_drop32[7]   + 2,
		   all_mean_drop33[1]/100,all_mean_drop33[6]   + 1,all_mean_drop33[7]   + 2,all_mean_drop34[1]/100,all_mean_drop34[6]   + 1,all_mean_drop34[7]   + 2,
		   all_mean_drop35[1]/100,all_mean_drop35[6]   + 1,all_mean_drop35[7]   + 2,all_mean_drop36[1]/100,all_mean_drop36[6]   + 1,all_mean_drop36[7]   + 2,
		   all_mean_drop37[1]/100,all_mean_drop37[6]   + 1,all_mean_drop37[7]   + 2,all_mean_drop38[1]/100,all_mean_drop38[6]   + 1,all_mean_drop38[7]   + 2,
		   all_mean_drop39[1]/100,all_mean_drop39[6]   + 1,all_mean_drop39[7] +2,
          all_mean_random[1]/100, all_mean_random[6]   + 1, all_mean_random[7]   + 2)


xnames <- c(1,1,1,3,3,3,5,5,5,7,7,7,9,9,9,11,11,11,13,13,13,15,15,15,17,17,17,19,19,19,21,21,21,23,23,23,25,25,25,27,27,27,29,29,29,31,31,31,33,33,33,35,35,35,37,37,37,39,39,39,41,41,41,43,43,43,45,45,45,47,47,47,49,49,49,51,51,51,53,53,53,55,55,55,57,57,57,59,59,59,61,61,61,63,63,63,65,65,65, 67,67,67,69,69,69,71,71,71,73,73,73,75,75,75,77,77,77,79,79,79,81,81,81,83,83,83)
plot(xnames, Dl_Twi_unbalanced, xaxt = "none", col = "black", pch = 19, ylim = c(0.1, 2.9),yaxt = "none", xlab = "")
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
abline(v=28,col="darkgray",lty=3)
abline(v=30,col="darkgray",lty=3)
abline(v=32,col="darkgray",lty=3)
abline(v=34,col="darkgray",lty=3)
abline(v=36,col="darkgray",lty=3)
abline(v=38,col="darkgray",lty=3)
abline(v=40,col="darkgray",lty=3)
abline(v=42,col="darkgray",lty=3)
abline(v=44,col="darkgray",lty=3)
abline(v=46,col="darkgray",lty=3)
abline(v=48,col="darkgray",lty=3)
abline(v=50,col="darkgray",lty=3)
abline(v=52,col="darkgray",lty=3)
abline(v=54,col="darkgray",lty=3)
abline(v=56,col="darkgray",lty=3)
abline(v=58,col="darkgray",lty=3)
abline(v=60,col="darkgray",lty=3)
abline(v=62,col="darkgray",lty=3)
abline(v=64,col="darkgray",lty=3)
abline(v=66,col="darkgray",lty=3)
abline(v=68,col="darkgray",lty=3)
abline(v=70,col="darkgray",lty=3)
abline(v=72,col="darkgray",lty=3)
abline(v=74,col="darkgray",lty=3)
abline(v=76,col="darkgray",lty=3)
abline(v=78,col="darkgray",lty=3)
abline(v=80,col="darkgray",lty=3)
abline(v=82,col="darkgray",lty=3)
abline(v=84,col="darkgray",lty=3)

lab = c("all data", "Zelda_2011_seq", "Snail_2014_chip", "Giant_2009_chip", "Bicoid_2013_seq", "Kruppel_2013_seq", "Bicoid_2009_chip", "Twist_2014_chip", "Twist_2011_seq", "Hairy_2009_chip", "Hunchback_2009_chip", "Knirps_2010_seq", "Caudal_2010_seq", "Twist_2009_chip", "Caudal_2009_chip", "p300_2010_seq", "Giant_2013_seq", "Hunchback_2013_seq", "Dorsal_2015_seq", "Zld_motif_sanger", "Snail_2009_chip", "Zld_motif_solexa", "Dere", "Dsim", "Dyak", "Dorsal_2009_chip", "Dana", "H3K4me1_2015_seq", "Dsim", "Dpse", "H3K27ac_2015_seq", "Dorsal_motif_NBT", "Snail_motif_FlyReg", "Dsec", "Dorsal_motif_FlyReg","Twist_motif_da", "Snail_motif_Sanger", "Dgri", "Snail_motif_solexa", "H3K27ac_2010_seq","Twist_motif_FlyReg","random_data")
text(x=c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83), y=par()$usr[3] - 0.01*(par()$usr[4]-par()$usr[3]),labels = lab,srt=65, adj=1, xpd=TRUE, cex=0.6)


current_all <- Dl_Twi 
current_tab <- table(current_all$status)
DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",43, "NA", "NA")
DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(DlTwi_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

current_all_drop_end <- current_all
current_all_drop_end$Twist_motif_da <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",42, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end1 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Snail_motif_Sanger <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",41, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end2 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Snail_motif_solexa <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",40, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end3 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Dorsal_motif_NBT <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",39, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end4 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$H3K4me1_2015_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",38, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end5 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Snail_motif_FlyReg <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",37, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end6 <- all_bal_means

current_all_drop_end$Dorsal_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",36, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end7 <- all_bal_means

current_all_drop_end$Dwil <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",35, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end8 <- all_bal_means

current_all_drop_end$Dorsal_motif_FlyReg <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",34, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end9 <- all_bal_means

current_all_drop_end$Twist_motif_FlyReg <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",33, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end10 <- all_bal_means

current_all_drop_end$H3K27ac_2015_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",32, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end11 <- all_bal_means

current_all_drop_end$H3K27ac_2010_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",31, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end12 <- all_bal_means

current_all_drop_end$Snail_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",30, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end13 <- all_bal_means

current_all_drop_end$Dgri <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",29, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end14 <- all_bal_means

current_all_drop_end$Dana <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",28, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end15 <- all_bal_means

current_all_drop_end$Dpse <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",27, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end16 <- all_bal_means

current_all_drop_end$Zld_motif_solexa <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",26, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end17 <- all_bal_means

current_all_drop_end$Dsec <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",25, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end18 <- all_bal_means

current_all_drop_end$Dsim <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",24, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end19 <- all_bal_means

current_all_drop_end$Dvir <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",23, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end20 <- all_bal_means

current_all_drop_end$Caudal_2010_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",22, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end21 <- all_bal_means

current_all_drop_end$Hairy_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",21, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end22 <- all_bal_means

current_all_drop_end$Caudal_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",20, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end23 <- all_bal_means

current_all_drop_end$Giant_2013_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",19, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end24 <- all_bal_means


current_all_drop_end$Dorsal_2015_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",18, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end25 <- all_bal_means

current_all_drop_end$p300_2010_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",17, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end26 <- all_bal_means

current_all_drop_end$Dyak <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",16, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end27 <- all_bal_means

current_all_drop_end$Twist_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",15, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end28 <- all_bal_means


current_all_drop_end$Dere <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",14, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end29 <- all_bal_means

current_all_drop_end$Knirps_2010_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",13, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end30 <- all_bal_means

current_all_drop_end$Zld_motif_sanger <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",12, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end31 <- all_bal_means

current_all_drop_end$Twist_2014_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",11, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end32 <- all_bal_means

current_all_drop_end$Bicoid_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",10, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end33 <- all_bal_means


current_all_drop_end$Hunchback_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",9, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end34 <- all_bal_means

current_all_drop_end$Hunchback_2013_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",8, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end35 <- all_bal_means

current_all_drop_end$Twist_2011_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",7, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end36 <- all_bal_means

current_all_drop_end$Giant_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",6, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end37 <- all_bal_means

current_all_drop_end$Kruppel_2013_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",5, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end38 <- all_bal_means

current_all_drop_end$Bicoid_2013_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",4, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end39 <- all_bal_means

current_all_drop_end$Snail_2014_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",3, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_mean_drop_end40 <- all_bal_means




Dl_Twi_Reverse <- c(DlTwi_means[1]/100,DlTwi_means[6]   + 1,DlTwi_means[7]   + 2, all_mean_drop_end1[1]/100,all_mean_drop_end1[6]   + 1,all_mean_drop_end1[7]   + 2,
		   all_mean_drop_end2[1]/100,all_mean_drop_end2[6]   + 1,all_mean_drop_end2[7]   + 2,all_mean_drop_end3[1]/100,all_mean_drop_end3[6]   + 1,all_mean_drop_end3[7]   + 2,
		   all_mean_drop_end4[1]/100,all_mean_drop_end4[6]   + 1,all_mean_drop_end4[7]   + 2,all_mean_drop_end5[1]/100,all_mean_drop_end5[6]   + 1,all_mean_drop_end5[7]   + 2,
		   all_mean_drop_end6[1]/100,all_mean_drop_end6[6]   + 1,all_mean_drop_end6[7]   + 2,all_mean_drop_end7[1]/100,all_mean_drop_end7[6]   + 1,all_mean_drop_end7[7]   + 2,
		   all_mean_drop_end8[1]/100,all_mean_drop_end8[6]   + 1,all_mean_drop_end8[7]   + 2,all_mean_drop_end9[1]/100,all_mean_drop_end9[6]   + 1,all_mean_drop_end9[7]   + 2,
		   all_mean_drop_end10[1]/100,all_mean_drop_end10[6]   + 1,all_mean_drop_end10[7]   + 2,all_mean_drop_end11[1]/100,all_mean_drop_end11[6]   + 1,all_mean_drop_end11[7] + 2,
		   all_mean_drop_end12[1]/100,all_mean_drop_end12[6]   + 1,all_mean_drop_end12[7]   + 2,
		   all_mean_drop_end13[1]/100,all_mean_drop_end13[6]   + 1,all_mean_drop_end13[7]   + 2,all_mean_drop_end14[1]/100,all_mean_drop_end14[6]   + 1,all_mean_drop_end14[7]   + 2,
		   all_mean_drop_end15[1]/100,all_mean_drop_end15[6]   + 1,all_mean_drop_end15[7]   + 2,all_mean_drop_end16[1]/100,all_mean_drop_end16[6]   + 1,all_mean_drop_end16[7]   + 2,
		   all_mean_drop_end17[1]/100,all_mean_drop_end17[6]   + 1,all_mean_drop_end17[7]   + 2,all_mean_drop_end18[1]/100,all_mean_drop_end18[6]   + 1,all_mean_drop_end18[7]   + 2,
		   all_mean_drop_end19[1]/100,all_mean_drop_end19[6]   + 1,all_mean_drop_end19[7]   + 2,all_mean_drop_end20[1]/100,all_mean_drop_end20[6]   + 1,all_mean_drop_end20[7]   + 2,
		   all_mean_drop_end21[1]/100,all_mean_drop_end21[6]   + 1,all_mean_drop_end21[7]   + 2,all_mean_drop_end22[1]/100,all_mean_drop_end22[6]   + 1,all_mean_drop_end22[7]   + 2,
		   all_mean_drop_end23[1]/100,all_mean_drop_end23[6]   + 1,all_mean_drop_end23[7]   + 2,all_mean_drop_end24[1]/100,all_mean_drop_end24[6]   + 1,all_mean_drop_end24[7]   + 2,
		   all_mean_drop_end25[1]/100,all_mean_drop_end25[6]   + 1,all_mean_drop_end25[7]   + 2,all_mean_drop_end26[1]/100,all_mean_drop_end26[6]   + 1,all_mean_drop_end26[7]   + 2,
		   all_mean_drop_end27[1]/100,all_mean_drop_end27[6]   + 1,all_mean_drop_end27[7]   + 2,all_mean_drop_end28[1]/100,all_mean_drop_end28[6]   + 1,all_mean_drop_end28[7]   + 2,
		   all_mean_drop_end29[1]/100,all_mean_drop_end29[6]   + 1,all_mean_drop_end29[7]   + 2,all_mean_drop_end30[1]/100,all_mean_drop_end30[6]   + 1,all_mean_drop_end30[7]   + 2,
		   all_mean_drop_end31[1]/100,all_mean_drop_end31[6]   + 1,all_mean_drop_end31[7]   + 2,all_mean_drop_end32[1]/100,all_mean_drop_end32[6]   + 1,all_mean_drop_end32[7]   + 2,
		   all_mean_drop_end33[1]/100,all_mean_drop_end33[6]   + 1,all_mean_drop_end33[7]   + 2,all_mean_drop_end34[1]/100,all_mean_drop_end34[6]   + 1,all_mean_drop_end34[7]   + 2,
		   all_mean_drop_end35[1]/100,all_mean_drop_end35[6]   + 1,all_mean_drop_end35[7]   + 2,all_mean_drop_end36[1]/100,all_mean_drop_end36[6]   + 1,all_mean_drop_end36[7]   + 2,
		   all_mean_drop_end37[1]/100,all_mean_drop_end37[6]   + 1,all_mean_drop_end37[7]   + 2,all_mean_drop_end38[1]/100,all_mean_drop_end38[6]   + 1,all_mean_drop_end38[7]   + 2,
		   all_mean_drop_end39[1]/100,all_mean_drop_end39[6]   + 1,all_mean_drop_end39[7] +2,all_mean_drop_end40[1]/100,all_mean_drop_end40[6]   + 1,all_mean_drop_end40[7]   + 2,
          all_mean_random[1]/100, all_mean_random[6]   + 1, all_mean_random[7]   + 2)


xnames <- c(1,1,1,3,3,3,5,5,5,7,7,7,9,9,9,11,11,11,13,13,13,15,15,15,17,17,17,19,19,19,21,21,21,23,23,23,25,25,25,27,27,27,29,29,29,31,31,31,33,33,33,35,35,35,37,37,37,39,39,39,41,41,41,43,43,43,45,45,45,47,47,47,49,49,49,51,51,51,53,53,53,55,55,55,57,57,57,59,59,59,61,61,61,63,63,63,65,65,65, 67,67,67,69,69,69,71,71,71,73,73,73,75,75,75,77,77,77,79,79,79,81,81,81,83,83,83)
plot(xnames, Dl_Twi_Reverse, xaxt = "none", col = "black", pch = 19, ylim = c(0.1, 2.9),yaxt = "none", xlab = "")
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
abline(v=28,col="darkgray",lty=3)
abline(v=30,col="darkgray",lty=3)
abline(v=32,col="darkgray",lty=3)
abline(v=34,col="darkgray",lty=3)
abline(v=36,col="darkgray",lty=3)
abline(v=38,col="darkgray",lty=3)
abline(v=40,col="darkgray",lty=3)
abline(v=42,col="darkgray",lty=3)
abline(v=44,col="darkgray",lty=3)
abline(v=46,col="darkgray",lty=3)
abline(v=48,col="darkgray",lty=3)
abline(v=50,col="darkgray",lty=3)
abline(v=52,col="darkgray",lty=3)
abline(v=54,col="darkgray",lty=3)
abline(v=56,col="darkgray",lty=3)
abline(v=58,col="darkgray",lty=3)
abline(v=60,col="darkgray",lty=3)
abline(v=62,col="darkgray",lty=3)
abline(v=64,col="darkgray",lty=3)
abline(v=66,col="darkgray",lty=3)
abline(v=68,col="darkgray",lty=3)
abline(v=70,col="darkgray",lty=3)
abline(v=72,col="darkgray",lty=3)
abline(v=74,col="darkgray",lty=3)
abline(v=76,col="darkgray",lty=3)
abline(v=78,col="darkgray",lty=3)
abline(v=80,col="darkgray",lty=3)
abline(v=82,col="darkgray",lty=3)
abline(v=84,col="darkgray",lty=3)

lab = c("all data", "Twist_motif_da", "Snail_motif_Sanger", "Snail_motif_solexa", "Dorsal_motif_NBT", "H3K4me1_2015_seq", "Snail_motif_FlyReg", "Dorsal_2009_chip", "Dwil", "Dorsal_motif_FlyReg", "Twist_motif_FlyReg", "H3K27ac_2015_seq", "H3K27ac_2010_seq", "Snail_2009_chip", "Dgri", "Dana", "Dpse", "Zld_motif_solexa", 
         "Dsec", "Dsim", "Dvir", "Caudal_2010_seq", "Hairy_2009_chip", "Caudal_2009_chip", "Giant_2013_seq", "Dorsal_2015_seq", "p300_2010_seq", "Dyak", "Twist_2009_chip", "Dere", "Knirps_2010_seq", "Zld_motif_sanger", "Twist_2014_chip", "Bicoid_2009_chip", "Hunchback_2009_chip","Hunchback_2013_seq", "Twist_2011_seq", "Giant_2009_chip", "Kruppel_2013_seq", "Bicoid_2013_seq","Snail_2014_chip","random_data")
text(x=c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83), y=par()$usr[3] - 0.01*(par()$usr[4]-par()$usr[3]),labels = lab,srt=65, adj=1, xpd=TRUE, cex=0.6)


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

all_bal_means_no_drop <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Twist_motif_da <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",42, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end1 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Twist_motif_FlyReg <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",41, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end2 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Snail_motif_Sanger <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",40, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end3 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Snail_motif_solexa <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",39, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end4 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Snail_motif_FlyReg <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",38, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end5 <- all_bal_means

current_all_drop_end <- current_all
current_all_drop_end$Dorsal_motif_FlyReg <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",37, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end6 <- all_bal_means

current_all_drop_end$Zld_motif_solexa <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",36, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end7 <- all_bal_means

current_all_drop_end$H3K27ac_2015_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",35, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end8 <- all_bal_means

current_all_drop_end$Dorsal_motif_NBT <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",34, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end9 <- all_bal_means

current_all_drop_end$Knirps_2010_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",33, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end10 <- all_bal_means

current_all_drop_end$H3K27ac_2010_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",32, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end11 <- all_bal_means

current_all_drop_end$Snail_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",31, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end12 <- all_bal_means

current_all_drop_end$p300_2010_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",30, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end13 <- all_bal_means

current_all_drop_end$Dorsal_2015_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",29, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end14 <- all_bal_means

current_all_drop_end$Dere <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",28, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end15 <- all_bal_means

current_all_drop_end$Dsim <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",27, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end16 <- all_bal_means

current_all_drop_end$Dpse <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",26, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end17 <- all_bal_means

current_all_drop_end$Giant_2013_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",25, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end18 <- all_bal_means

current_all_drop_end$Zld_motif_sanger <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",24, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end19 <- all_bal_means

current_all_drop_end$Dwil <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",23, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end20 <- all_bal_means

current_all_drop_end$Dana <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",22, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end21 <- all_bal_means

current_all_drop_end$Dgri <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",21, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end22 <- all_bal_means

current_all_drop_end$Bicoid_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",20, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end23 <- all_bal_means

current_all_drop_end$Dsec <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",19, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end24 <- all_bal_means


current_all_drop_end$Dyak <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",18, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end25 <- all_bal_means

current_all_drop_end$Dvir <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",17, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end26 <- all_bal_means

current_all_drop_end$Twist_2014_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",16, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end27 <- all_bal_means

current_all_drop_end$Hairy_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",15, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end28 <- all_bal_means


current_all_drop_end$Caudal_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",14, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end29 <- all_bal_means

current_all_drop_end$Kruppel_2013_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",13, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end30 <- all_bal_means

current_all_drop_end$Caudal_2010_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",12, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end31 <- all_bal_means

current_all_drop_end$Dorsal_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",11, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end32 <- all_bal_means

current_all_drop_end$Giant_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",10, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end33 <- all_bal_means


current_all_drop_end$Hunchback_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",9, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end34 <- all_bal_means

current_all_drop_end$Bicoid_2013_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",8, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end35 <- all_bal_means

current_all_drop_end$Twist_2011_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",7, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end36 <- all_bal_means

current_all_drop_end$Twist_2009_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",6, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end37 <- all_bal_means

current_all_drop_end$H3K4me1_2015_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",5, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end38 <- all_bal_means

current_all_drop_end$Hunchback_2013_seq <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",4, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end39 <- all_bal_means

current_all_drop_end$Snail_2014_chip <- NULL
current_all <- current_all_drop_end
current_tab <- table(current_all$status)
all_bal_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "allData_balance", description = "allData_balance", "allData_balance", 10, 250, strata_Percent_var, "ROC", current_target, "allData_balance_Test",3, "NA", "NA")

all_bal_means <- apply(all_bal_ROC_data[[1]], MARGIN = 2, FUN = mean)
all_bal_sd <- apply(all_bal_ROC_data[[1]], MARGIN=2, FUN=sd)
all_imp <- data.frame(all_bal_ROC_data[2])
all_mean <- rowMeans(all_imp)
all_sd_imp <- rowSds(as.matrix(all_imp))
all_temp <- cbind(all_mean, all_sd_imp)
alldata <- all_temp[order(all_mean),]
dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)

all_balance_drop_end40 <- all_bal_means

current_all <- random_balance
colnames(current_all) <- short_name
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

all_bal_means_no_drop <- all_bal_means
all_balance_Reverse <- c(all_bal_means_no_drop[1]/100,all_bal_means_no_drop[6]   + 1,all_bal_means_no_drop[7]   + 2, all_balance_drop_end1[1]/100,all_balance_drop_end1[6]   + 1,all_balance_drop_end1[7]   + 2,
		   all_balance_drop_end2[1]/100,all_balance_drop_end2[6]   + 1,all_balance_drop_end2[7]   + 2,all_balance_drop_end3[1]/100,all_balance_drop_end3[6]   + 1,all_balance_drop_end3[7]   + 2,
		   all_balance_drop_end4[1]/100,all_balance_drop_end4[6]   + 1,all_balance_drop_end4[7]   + 2,all_balance_drop_end5[1]/100,all_balance_drop_end5[6]   + 1,all_balance_drop_end5[7]   + 2,
		   all_balance_drop_end6[1]/100,all_balance_drop_end6[6]   + 1,all_balance_drop_end6[7]   + 2,all_balance_drop_end7[1]/100,all_balance_drop_end7[6]   + 1,all_balance_drop_end7[7]   + 2,
		   all_balance_drop_end8[1]/100,all_balance_drop_end8[6]   + 1,all_balance_drop_end8[7]   + 2,all_balance_drop_end9[1]/100,all_balance_drop_end9[6]   + 1,all_balance_drop_end9[7]   + 2,
		   all_balance_drop_end10[1]/100,all_balance_drop_end10[6]   + 1,all_balance_drop_end10[7]   + 2,all_balance_drop_end11[1]/100,all_balance_drop_end11[6]   + 1,all_balance_drop_end11[7] + 2,
		   all_balance_drop_end12[1]/100,all_balance_drop_end12[6]   + 1,all_balance_drop_end12[7]   + 2,
		   all_balance_drop_end13[1]/100,all_balance_drop_end13[6]   + 1,all_balance_drop_end13[7]   + 2,all_balance_drop_end14[1]/100,all_balance_drop_end14[6]   + 1,all_balance_drop_end14[7]   + 2,
		   all_balance_drop_end15[1]/100,all_balance_drop_end15[6]   + 1,all_balance_drop_end15[7]   + 2,all_balance_drop_end16[1]/100,all_balance_drop_end16[6]   + 1,all_balance_drop_end16[7]   + 2,
		   all_balance_drop_end17[1]/100,all_balance_drop_end17[6]   + 1,all_balance_drop_end17[7]   + 2,all_balance_drop_end18[1]/100,all_balance_drop_end18[6]   + 1,all_balance_drop_end18[7]   + 2,
		   all_balance_drop_end19[1]/100,all_balance_drop_end19[6]   + 1,all_balance_drop_end19[7]   + 2,all_balance_drop_end20[1]/100,all_balance_drop_end20[6]   + 1,all_balance_drop_end20[7]   + 2,
		   all_balance_drop_end21[1]/100,all_balance_drop_end21[6]   + 1,all_balance_drop_end21[7]   + 2,all_balance_drop_end22[1]/100,all_balance_drop_end22[6]   + 1,all_balance_drop_end22[7]   + 2,
		   all_balance_drop_end23[1]/100,all_balance_drop_end23[6]   + 1,all_balance_drop_end23[7]   + 2,all_balance_drop_end24[1]/100,all_balance_drop_end24[6]   + 1,all_balance_drop_end24[7]   + 2,
		   all_balance_drop_end25[1]/100,all_balance_drop_end25[6]   + 1,all_balance_drop_end25[7]   + 2,all_balance_drop_end26[1]/100,all_balance_drop_end26[6]   + 1,all_balance_drop_end26[7]   + 2,
		   all_balance_drop_end27[1]/100,all_balance_drop_end27[6]   + 1,all_balance_drop_end27[7]   + 2,all_balance_drop_end28[1]/100,all_balance_drop_end28[6]   + 1,all_balance_drop_end28[7]   + 2,
		   all_balance_drop_end29[1]/100,all_balance_drop_end29[6]   + 1,all_balance_drop_end29[7]   + 2,all_balance_drop_end30[1]/100,all_balance_drop_end30[6]   + 1,all_balance_drop_end30[7]   + 2,
		   all_balance_drop_end31[1]/100,all_balance_drop_end31[6]   + 1,all_balance_drop_end31[7]   + 2,all_balance_drop_end32[1]/100,all_balance_drop_end32[6]   + 1,all_balance_drop_end32[7]   + 2,
		   all_balance_drop_end33[1]/100,all_balance_drop_end33[6]   + 1,all_balance_drop_end33[7]   + 2,all_balance_drop_end34[1]/100,all_balance_drop_end34[6]   + 1,all_balance_drop_end34[7]   + 2,
		   all_balance_drop_end35[1]/100,all_balance_drop_end35[6]   + 1,all_balance_drop_end35[7]   + 2,all_balance_drop_end36[1]/100,all_balance_drop_end36[6]   + 1,all_balance_drop_end36[7]   + 2,
		   all_balance_drop_end37[1]/100,all_balance_drop_end37[6]   + 1,all_balance_drop_end37[7]   + 2,all_balance_drop_end38[1]/100,all_balance_drop_end38[6]   + 1,all_balance_drop_end38[7]   + 2,
		   all_balance_drop_end39[1]/100,all_balance_drop_end39[6]   + 1,all_balance_drop_end39[7] +2,all_balance_drop_end40[1]/100,all_balance_drop_end40[6]   + 1,all_balance_drop_end40[7]   + 2,
          Random_bal_means[1]/100, Random_bal_means[6]   + 1, Random_bal_means[7]   + 2)


xnames <- c(1,1,1,3,3,3,5,5,5,7,7,7,9,9,9,11,11,11,13,13,13,15,15,15,17,17,17,19,19,19,21,21,21,23,23,23,25,25,25,27,27,27,29,29,29,31,31,31,33,33,33,35,35,35,37,37,37,39,39,39,41,41,41,43,43,43,45,45,45,47,47,47,49,49,49,51,51,51,53,53,53,55,55,55,57,57,57,59,59,59,61,61,61,63,63,63,65,65,65, 67,67,67,69,69,69,71,71,71,73,73,73,75,75,75,77,77,77,79,79,79,81,81,81,83,83,83)
plot(xnames, all_balance_Reverse, xaxt = "none", col = "cyan4", pch = 19, ylim = c(0.1, 2.9),yaxt = "none", xlab = "")
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
abline(v=28,col="darkgray",lty=3)
abline(v=30,col="darkgray",lty=3)
abline(v=32,col="darkgray",lty=3)
abline(v=34,col="darkgray",lty=3)
abline(v=36,col="darkgray",lty=3)
abline(v=38,col="darkgray",lty=3)
abline(v=40,col="darkgray",lty=3)
abline(v=42,col="darkgray",lty=3)
abline(v=44,col="darkgray",lty=3)
abline(v=46,col="darkgray",lty=3)
abline(v=48,col="darkgray",lty=3)
abline(v=50,col="darkgray",lty=3)
abline(v=52,col="darkgray",lty=3)
abline(v=54,col="darkgray",lty=3)
abline(v=56,col="darkgray",lty=3)
abline(v=58,col="darkgray",lty=3)
abline(v=60,col="darkgray",lty=3)
abline(v=62,col="darkgray",lty=3)
abline(v=64,col="darkgray",lty=3)
abline(v=66,col="darkgray",lty=3)
abline(v=68,col="darkgray",lty=3)
abline(v=70,col="darkgray",lty=3)
abline(v=72,col="darkgray",lty=3)
abline(v=74,col="darkgray",lty=3)
abline(v=76,col="darkgray",lty=3)
abline(v=78,col="darkgray",lty=3)
abline(v=80,col="darkgray",lty=3)
abline(v=82,col="darkgray",lty=3)
abline(v=84,col="darkgray",lty=3)

lab = c("Dl_Twi", "Twist_motif_da", "Twist_motif_FlyReg", "Snail_motif_Sanger", "Snail_motif_solexa", "Snail_motif_FlyReg", "Dorsal_motif_FlyReg", "Zld_motif_solexa", "H3K27ac_2015_seq", "Dorsal_motif_NBT", "Knirps_2010_seq", "H3K27ac_2010_seq", "Snail_2009_chip", "p300_2010_seq", "Dorsal_2015_seq", "Dere", "Dsim","Dpse", "Giant_2013_seq", "Zld_motif_sanger", "Dwil", "Dana", "Dgri", "Bicoid_2009_chip", "Dsec", "Dyak", "Dvir","Twist_2014_chip", "Hairy_2009_chip", "Caudal_2009_chip", "Kruppel_2013_seq", "Caudal_2010_seq", "Dorsal_2009_chip","Giant_2009_chip", "Hunchback_2009_chip", "Bicoid_2013_seq","Twist_2011_seq", "Twist_2009_chip", "H3K4me1_2015_seq", "Hunchback_2013_seq", "Snail_2014_chip","random_data")

text(x=c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83), y=par()$usr[3] - 0.01*(par()$usr[4]-par()$usr[3]),labels = lab,srt=65, adj=1, xpd=TRUE, cex=0.6)

