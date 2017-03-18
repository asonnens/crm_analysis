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
require(dplyr)


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

all_data_balance <- all_data %>%
group_by(status) %>%
sample_n(size = 300)

#Filtered reporters
H3K4me1_train_cons <- read.table("cons_train_H3K4me1_2015_seq.tsv", header = TRUE)
H3K4me1_train_all <- read.table("all_train_H3K4me1_2015_seq.tsv", header = TRUE)
FTwi_train_cons <- read.table("cons_train_Twist_2014_chip.tsv", header = TRUE)
FTwi_train_all <- read.table("all_train_Twist_2014_chip.tsv", header = TRUE)
Zld_train_cons <- read.table("cons_train_Zelda_2011_seq.tsv", header = TRUE)
Zld_train_all <- read.table("all_train_Zelda_2011_seq.tsv", header = TRUE)
Dl15_train_cons <- read.table("cons_train_Dorsal_2015_seq.tsv", header = TRUE)
Dl15_train_all <- read.table("all_train_Dorsal_2015_seq.tsv", header = TRUE)
Bcd13_train_cons <- read.table("cons_train_Bicoid_2013_seq.tsv", header = TRUE)
Bcd13_train_all <- read.table("all_train_Bicoid_2013_seq.tsv", header = TRUE)
Bcd09_train_cons <- read.table("cons_train_Bicoid_2009_chip.tsv", header = TRUE)
Bcd09_train_all <- read.table("all_train_Bicoid_2009_chip.tsv", header = TRUE)
Cad09_train_cons <- read.table("cons_train_Caudal_2009_chip.tsv", header = TRUE)
Cad09_train_all <- read.table("all_train_Caudal_2009_chip.tsv", header = TRUE)
Gt09_train_cons <- read.table("cons_train_Giant_2009_chip.tsv", header = TRUE)
Gt09_train_all <- read.table("all_train_Giant_2009_chip.tsv", header = TRUE)
Gt13_train_cons <- read.table("cons_train_Giant_2013_seq.tsv", header = TRUE)
Gt13_train_all <- read.table("all_train_Giant_2013_seq.tsv", header = TRUE)
Kr13_train_cons <- read.table("cons_train_Kruppel_2013_seq.tsv", header = TRUE)
Kr13_train_all <- read.table("all_train_Kruppel_2013_seq.tsv", header = TRUE)
Kni10_train_cons <- read.table("cons_train_Knirps_2010_seq.tsv", header = TRUE)
Kni10_train_all <- read.table("all_train_Knirps_2010_seq.tsv", header = TRUE)
Hry09_train_cons <- read.table("cons_train_Hairy_2009_chip.tsv", header = TRUE)
Hry09_train_all <- read.table("all_train_Hairy_2009_chip.tsv", header = TRUE)

current_target <- read.table("target_file.tsv", header = TRUE)
prosp_target <- read.table("prospective_file.tsv", header = TRUE)
newdata <- current_target

#Random data
random_data <- data.frame(replicate(41,sample(0:1000,7250,rep=TRUE)))

random_data <-cbind(all_data[,0:2],random_data)
random_balance <- data.frame(replicate(41,sample(0:1000,300,rep=TRUE)))
name_column <- data.frame(replicate(1,sample(c("on", "off"),300, rep = TRUE)))
random_balance <-cbind(all_data[1:300,0:1], name_column, random_balance)
colnames(random_data) <- long_name
colnames(random_balance) <- long_name

random_bcd <- data.frame(replicate(41, sample(0:1000, 371, rep = TRUE)))
random_bcd <- cbind(Bcd13_train_all[,0:2], random_bcd)
colnames(random_bcd) <- long_name

random_bcd_bal <- data.frame(replicate(41, sample(0:1000, 100, rep = TRUE)))
random_bcd_bal <- cbind(Bcd13_train_cons[,0:2], random_bcd_bal)
colnames(random_bcd_bal) <- long_name

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

Leave_one_out <- function(mydataset, mydatabalance){
    current_all <- mydataset
	current_tab <- table(mydataset$status)
	datalength1 <- length(rownames(mydataset))
	datalength2 <- length(rownames(mydatabalance))
	dataval <- length(current_all)
	random_data <- data.frame(replicate(41,sample(0:1000,datalength1,rep=TRUE)))
	random_data <-cbind(mydataset[,0:2],random_data)
	name_column <- data.frame(replicate(1,sample(c("on", "off"),datalength2, rep = TRUE)))
	random_balance <- data.frame(replicate(41,sample(0:1000,datalength2,rep=TRUE)))
	random_balance <-cbind(mydatabalance[,1], name_column, random_balance)
	colnames(random_data) <- long_name
	colnames(random_balance) <- long_name
	l <- vector(mode="numeric", length=126) #unbalanced, dropping most important
	l2 <- vector(mode="character", length=42)  #most important to least important features
	l3 <- vector(mode="numeric", length=126)#unbalanced, dropping least important
	l4 <- vector(mode="character", length=42) #least important to most important features
	l5 <- vector(mode="numeric", length=126) #balanced, dropping most important
	l6 <- vector(mode="character", length=42)  #most important to least important features
	l7 <- vector(mode="numeric", length=126)#balanced, dropping least important
	l8 <- vector(mode="character", length=42) #least important to most important features
	DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",43, "NA", "NA")
    DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
    DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
	l[1] <- DlTwi_means[1]/100
	l[2] <- DlTwi_means[6]   + 1
	l[3] <- DlTwi_means[7]   + 2
	l3[1] <- DlTwi_means[1]/100
	l3[2] <- DlTwi_means[6]   + 1
	l3[3] <- DlTwi_means[7]   + 2
    all_imp <- data.frame(DlTwi_ROC_data[2])
    all_mean <- rowMeans(all_imp)
    all_sd_imp <- rowSds(as.matrix(all_imp))
    all_temp <- cbind(all_mean, all_sd_imp)
	alldata <- all_temp[order(all_mean),]
	dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
    segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)
	firstname <- rownames(alldata)
	firstname2 <- firstname
	count = 4
	count2 = 1
	for (i in 1:39){
		temptemp <- firstname
	    temptemptemp <- tail(temptemp, n = 1)
		l2[count2] <- temptemptemp
		drop_val <- grep(temptemptemp, colnames(current_all))
		current_all <- current_all[, -drop_val]
	    DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",dataval - i, "NA", "NA")
        DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
        DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
        all_imp <- data.frame(DlTwi_ROC_data[2])
        all_mean <- rowMeans(all_imp)
        all_sd_imp <- rowSds(as.matrix(all_imp))
        all_temp <- cbind(all_mean, all_sd_imp)
	    alldata <- all_temp[order(all_mean),]
		if (length(alldata) > 1){
		dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
        segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)}
		firstname <- rownames(alldata)
		l[count] <- DlTwi_means[1]/100
		l[count + 1] <- DlTwi_means[6]   + 1
		l[count + 2] <- DlTwi_means[7]   + 2
		count <- count + 3
		count2 <- count2 + 1
		l2[40] <- tail(firstname, n = 2)
		l2[41] <- tail(firstname, n = 3)
		}
	count <- 4
	count2 <- 1
	current_all <- mydataset
	colnames(current_all) <- long_name
	for (j in 1:39){
		temptemp <- firstname2
	    temptemptemp <- head(temptemp, n = 1)
		l4[count2] <- temptemptemp
		drop_val <- grep(temptemptemp, colnames(current_all))
		current_all <- current_all[, -drop_val]
	    DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",dataval - j, "NA", "NA")
        DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
        DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
        all_imp <- data.frame(DlTwi_ROC_data[2])
        all_mean <- rowMeans(all_imp)
        all_sd_imp <- rowSds(as.matrix(all_imp))
        all_temp <- cbind(all_mean, all_sd_imp)
	    alldata <- all_temp[order(all_mean),]
		if (length(alldata) > 1){
		dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
        segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)}
		firstname2 <- rownames(alldata)
		l3[count] <- DlTwi_means[1]/100
		l3[count + 1] <- DlTwi_means[6]   + 1
		l3[count + 2] <- DlTwi_means[7]   + 2
		count <- count + 3
		count2 <- count2 + 1
		l4[40] <- tail(firstname, n = 2)
		l4[41] <- tail(firstname, n = 3)
		}
    colnames(mydatabalance) <- long_name
 	current_all <- mydatabalance
	current_tab <- table(current_all$status)
	DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",43, "NA", "NA")
    DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
    DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
	l5[1] <- DlTwi_means[1]/100
	l5[2] <- DlTwi_means[6]   + 1
	l5[3] <- DlTwi_means[7]   + 2
	l7[1] <- DlTwi_means[1]/100
	l7[2] <- DlTwi_means[6]   + 1
	l7[3] <- DlTwi_means[7]   + 2
    all_imp <- data.frame(DlTwi_ROC_data[2])
    all_mean <- rowMeans(all_imp)
    all_sd_imp <- rowSds(as.matrix(all_imp))
    all_temp <- cbind(all_mean, all_sd_imp)
	alldata <- all_temp[order(all_mean),]
	if (length(alldata) > 1){
		dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
        segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)}
	firstname <- rownames(alldata)
	firstname2 <- firstname
	count = 4
	count2 = 1
	for (i in 1:39){
		temptemp <- firstname
	    temptemptemp <- tail(temptemp, n = 1)
		l6[count2] <- temptemptemp
		drop_val <- grep(temptemptemp, colnames(current_all))
		current_all <- current_all[, -drop_val]
	    DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",dataval - i, "NA", "NA")
        DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
        DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
        all_imp <- data.frame(DlTwi_ROC_data[2])
        all_mean <- rowMeans(all_imp)
        all_sd_imp <- rowSds(as.matrix(all_imp))
        all_temp <- cbind(all_mean, all_sd_imp)
	    alldata <- all_temp[order(all_mean),]
		if (length(alldata) > 1){
		dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
        segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)}
		firstname <- rownames(alldata)
		l5[count] <- DlTwi_means[1]/100
		l5[count + 1] <- DlTwi_means[6]   + 1
		l5[count + 2] <- DlTwi_means[7]   + 2
		count <- count + 3
		count2 <- count2 + 1
		l6[40] <- tail(firstname, n = 2)
		l6[41] <- tail(firstname, n = 3)
		}
	count <- 4
	count2 <- 1
	colnames(mydatabalance) <- long_name
	current_all <- mydatabalance
	for (j in 1:39){
		temptemp <- firstname2
	    temptemptemp <- head(temptemp, n = 1)
		l8[count2] <- temptemptemp
		l8[]
		drop_val <- grep(temptemptemp, colnames(current_all))
		current_all <- current_all[, -drop_val]
	    DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",dataval - j, "NA", "NA")
        DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
        DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
        all_imp <- data.frame(DlTwi_ROC_data[2])
        all_mean <- rowMeans(all_imp)
        all_sd_imp <- rowSds(as.matrix(all_imp))
        all_temp <- cbind(all_mean, all_sd_imp)
	    alldata <- all_temp[order(all_mean),]
		if (length(alldata) > 1){
		dotchart(alldata[,1], xlab = "Importance", xlim = c(-5,30), cex = 0.7, pch = 19)
        segments(alldata[,1]-alldata[,2], 1:43, alldata[,1]+alldata[,2], 1:43)}
		firstname2 <- rownames(alldata)
		l7[count] <- DlTwi_means[1]/100
		l7[count + 1] <- DlTwi_means[6]   + 1
		l7[count + 2] <- DlTwi_means[7]   + 2
		count <- count + 3
		count2 <- count2 + 1
		l8[40] <- tail(firstname, n = 2)
		l8[41] <- tail(firstname, n = 3)
		}
	current_all <- random_data
	colnames(current_all) <- long_name
	current_tab <- table(random_data$status)
    DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",43, "NA", "NA")
    DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
    DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
    all_imp <- data.frame(DlTwi_ROC_data[2])
    all_mean <- rowMeans(all_imp)
    all_sd_imp <- rowSds(as.matrix(all_imp))
    all_temp <- cbind(all_mean, all_sd_imp)
	alldata <- all_temp[order(all_mean),]
	l[124] <- DlTwi_means[1]/100
	l[125] <- DlTwi_means[6]   + 1
	l[126] <- DlTwi_means[7]   + 2
    l2[41] <- "random"
	l3[124] <- DlTwi_means[1]/100
	l3[125] <- DlTwi_means[6]   + 1
	l3[126] <- DlTwi_means[7]   + 2
    l4[41] <- "random"
	current_all <- random_balance
	colnames(current_all) <- long_name	
	current_tab <- table(random_balance$status)
    DlTwi_ROC_data <- resample_rf500(enhancer_class_status_scale, current_all, current_tab, "Dl_Twi ", description = "Dl_Twi ", "Dl_Twi", 10, 250, strata_Percent_var, "ROC", current_target, "Dl_Twi_Test",43, "NA", "NA")
    DlTwi_means <- apply(DlTwi_ROC_data[[1]], MARGIN = 2, FUN = mean)
    DlTwi_sd <- apply(DlTwi_ROC_data[[1]], MARGIN=2, FUN=sd)
    all_imp <- data.frame(DlTwi_ROC_data[2])
    all_mean <- rowMeans(all_imp)
    all_sd_imp <- rowSds(as.matrix(all_imp))
    all_temp <- cbind(all_mean, all_sd_imp)
	alldata <- all_temp[order(all_mean),]
	l5[124] <- DlTwi_means[1]/100
	l5[125] <- DlTwi_means[6]   + 1
	l5[126] <- DlTwi_means[7]   + 2
    l6[41] <- "random"
	l7[124] <- DlTwi_means[1]/100
	l7[125] <- DlTwi_means[6]   + 1
	l7[126] <- DlTwi_means[7]   + 2
    l8[41] <- "random"
    return_frame <- list(l, l2, l3, l4, l5, l6, l7, l8)	 
	return(return_frame)

}

myframe <- Leave_one_out(Bcd13_train_all, Bcd13_train_cons)
xnames <- c(1,1,1,3,3,3,5,5,5,7,7,7,9,9,9,11,11,11,13,13,13,15,15,15,17,17,17,19,19,19,21,21,21,23,23,23,25,25,25,27,27,27,29,29,29,31,31,31,33,33,33,35,35,35,37,37,37,39,39,39,41,41,41,43,43,43,45,45,45,47,47,47,49,49,49,51,51,51,53,53,53,55,55,55,57,57,57,59,59,59,61,61,61,63,63,63,65,65,65, 67,67,67,69,69,69,71,71,71,73,73,73,75,75,75,77,77,77,79,79,79,81,81,81,83,83,83)
plot(xnames, myframe[[1]], xaxt = "none", col = "black", pch = "+", ylim = c(0.05, 3.0),yaxt = "none", xlab = "", main = "Leave one out- Bcd_13_s")
points(xnames, myframe[[3]], col = "black", pch = "-")
points(xnames, myframe[[5]], col = "cadetblue4", pch = "+")
points(xnames, myframe[[7]], col = "cadetblue4", pch = "-")
abline(h=1,col=1,lty=1)
abline(h=2,col=1,lty=1)
abline(h=0.5,col="gray",lty=2)
abline(h=1.5,col="gray",lty=2)
abline(h=2.5,col="gray",lty=2)
abline(h=3,col="gray",lty=1)
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

myframe2 <- Leave_one_out(all_data, all_data_balance)
xnames <- c(1,1,1,3,3,3,5,5,5,7,7,7,9,9,9,11,11,11,13,13,13,15,15,15,17,17,17,19,19,19,21,21,21,23,23,23,25,25,25,27,27,27,29,29,29,31,31,31,33,33,33,35,35,35,37,37,37,39,39,39,41,41,41,43,43,43,45,45,45,47,47,47,49,49,49,51,51,51,53,53,53,55,55,55,57,57,57,59,59,59,61,61,61,63,63,63,65,65,65, 67,67,67,69,69,69,71,71,71,73,73,73,75,75,75,77,77,77,79,79,79,81,81,81,83,83,83)
plot(xnames, myframe2[[1]], xaxt = "none", col = "black", pch = "-", ylim = c(0.5, 3.0),yaxt = "none", xlab = "", main = "Leave one out- all data")
points(xnames, myframe2[[3]], col = "black", pch = "+")
points(xnames, myframe2[[5]], col = "cadetblue4", pch = "+")
points(xnames, myframe2[[7]], col = "cadetblue4", pch = "-")
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