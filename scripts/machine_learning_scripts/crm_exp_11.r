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

#maybe two functions
#is this an enhancer, if yes, what's its expression pattern

setwd("shadow_enhancers/machine_learning_input/September")


Mac_Dl_train_cons <- read.table("cons_train_MacArthur_Dorsal.tsv", header = TRUE)
Mac_Dl_train_lax <- read.table("lax_train_MacArthur_Dorsal.tsv", header = TRUE)
Mac_Dl_train_rand <- read.table("rand_train_MacArthur_Dorsal.tsv", header = TRUE)
Mac_Dl_put <- read.table("Dorsal_putative_Sept1.tsv", header = TRUE)

Mac_Twi_train_cons <- read.table("cons_train_MacArthur_Twist.tsv", header = TRUE)
Mac_Twi_train_lax <- read.table("lax_train_MacArthur_Twist.tsv", header = TRUE)
Mac_Twi_train_rand <- read.table("rand_train_MacArthur_Twist.tsv", header = TRUE)
Mac_Twi_put <- read.table("Twist_putative_Sept1.tsv", header = TRUE)

RDl_train_cons <- read.table("cons_train_Rushlow_Dorsal.tsv", header = TRUE)
RDl_train_lax <- read.table("lax_train_Rushlow_Dorsal.tsv", header = TRUE)
RDl_train_rand <- read.table("rand_train_Rushlow_Dorsal.tsv", header = TRUE)
RDl_put <- read.table("DorsalR_putative_Sept1.tsv", header = TRUE)

Zld_train_cons <- read.table("cons_train_Zelda180.tsv", header = TRUE)
Zld_train_lax <- read.table("lax_train_Zelda180.tsv", header = TRUE)
Zld_train_rand <- read.table("rand_train_Zelda180.tsv", header = TRUE)
Zld_put <- read.table("Zelda180_putative_Sept1.tsv", header = TRUE)

express <- read.table("expression_training.tsv", header = TRUE)

tab_MacDl_cons <- table(Mac_Dl_train_cons$status)
tab_MacDl_cons
tab_MacDl_lax <- table(Mac_Dl_train_lax$status)
tab_MacDl_lax
tab_MacDl_rand <- table(Mac_Dl_train_rand$status)
tab_MacDl_rand 

tab_MacTwi_cons <- table(Mac_Twi_train_cons$status)
tab_MacTwi_cons
tab_MacTwi_lax <- table(Mac_Twi_train_lax$status)
tab_MacTwi_lax
tab_MacTwi_rand <- table(Mac_Twi_train_rand$status)
tab_MacTwi_rand 

tab_RDl_cons <- table(RDl_train_cons$status)
tab_RDl_cons
tab_RDl_lax <- table(RDl_train_lax$status)
tab_RDl_lax
tab_RDl_rand <- table(RDl_train_rand$status)
tab_RDl_rand 

tab_Zld_cons <- table(Zld_train_cons$status)
tab_Zld_cons
tab_Zld_lax <- table(Zld_train_lax$status)
tab_Zld_lax
tab_Zld_rand <- table(Zld_train_rand$status)
tab_Zld_rand 


my_features <- all_features
all_features <- c(3:24)

Dorsal_snail_twist_allsources <- c(10,11,19:22,24)
Dorsal_snail_twist_allsources_Zelda <- c(3,10,11,19:22,24)
Dorsal_snail_twist_noMacArthur <- c(10,11,22,24)

drop_Zeit_Twist <- c(3:23)
drop_White_H3K27ac <- c(3:22,24)
drop_Rushlow_Dorsal <- c(3:21,23:24)
drop_MacDorsal <- c(3:20,22:24)
drop_MacTwist <- c(3:19,21:24)
drop_MacSnail <- c(3:18,20:24)
drop_MacHairy <- c(3:17:19:24)
drop_MacHunchback <- c(3:16,18:24)
drop_MacGiant <- c(3:15,17:24)
drop_MacCaudal <- c(3:14,16:24)
drop_MacBicoid <- c(3:13,15:24)
drop_KK_H3K4me1 <- c(3:12,14:24)
drop_KK_H3K27ac <- c(3:11,13:24)
drop_FurlongTwist <- c(3:10,12:24)
drop_FurlongSnail <- c(3:9,11:24)
drop_Eisen13Kr <- c(3:8,10:24)
drop_Eisen13Bcd <- c(3:7,9:24)
drop_Eisen13Hb <- c(3:6,8:24)
drop_Eisen13Gt <- c(3:5,7:24)
drop_Eisen10Kni <- c(3:4,6:24)
drop_Eisen10Cad <- c(3,5:24)
drop_Zelda <- c(4:24)
drop_Dl <- c(3:20,23:24)
drop_Bcd <- c(3:7,9:13,15:24)
drop_Cad <- c(3,5:14,16:24)
only_Zelda <- c(3,22)
only_Eisen10Cad <- c(4,22)
only_Eisen10Kni <- c(5,22)
only_Eisen13Gt <- c(6,22)
only_Eisen13Hb <- c(7,22)
only_Eisen13Bcd <- c(8,22)
only_Eisen13Kr <- c(9,22)
only_FurlongSnail <- c(10,22)
only_FurlongTwist <- c(11,22)
only_KK <- c(12:13,22)
only_MacBicoid <- c(14,22)
only_MacCad <- c(15,22)
only_MacGiant <- c(16,22)
only_MacHb <- c(17,22)
only_MacHry <- c(18,22)
only_MacSna <- c(19,22)
only_MacTwi <- c(20,22)
only_MacDl <- c(21,22)
only_White <- c(22,23)
only_Zeit <- c(22,24)



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



resample_default <- function(function_name, dataset, datatable, variable_name, outfilename, description="README!", reps, samplesize){
    temp_dir <- paste("../September/", outfilename, sep = "")
    dir.create(temp_dir)
	info_file_name <- paste(temp_dir, "input_data_summary.txt", sep = "/")
	input_data_name <- paste(temp_dir, "input_data.csv", sep = "/")
	readme_file_name <- paste(temp_dir, "README.txt", sep = "/")
	write(description, file = readme_file_name)
	write.csv(crm_set, file = input_data_name)
	#summary_info <- str(dataset)
	#write.table(datatable, informational_file)
	#writeLines(summary_info, informational_file)
	#close(informational_file)
	datathing <- data.frame()
    for(i in 1:reps){
	    temp_file_name <- paste(outfilename, i, sep = "_")
		temp_file_name <- paste(temp_file_name, ".csv", sep = "")
		temp_file_name <- paste(temp_dir, temp_file_name, sep = "/")
        myresult <- suppressWarnings(function_name(dataset, datatable, temp_file_name, info_file_name, samplesize))
			datathing[i,"svmrad_acc"] <- myresult[[1]][1]
			datathing[i,"svmrad_tpos"] <- myresult[[2]][1]
			datathing[i,"svmrad_tneg"] <- myresult[[3]][1]
			datathing[i,"svmrad_fpos"] <- myresult[[4]][1]
			datathing[i,"svmrad_fneg"] <- myresult[[5]][1]
			datathing[i,"svmsig_acc"] <- myresult[[1]][2]
			datathing[i,"svmsig_tpos"] <- myresult[[2]][2]
			datathing[i,"svmsig_tneg"] <- myresult[[3]][2]
			datathing[i,"svmsig_fpos"] <- myresult[[4]][2]
			datathing[i,"svmsig_fneg"] <- myresult[[5]][2]
			datathing[i,"rf500_acc"] <- myresult[[1]][3]
			datathing[i,"rf500_tpos"] <- myresult[[2]][3]
			datathing[i,"rf500_tneg"] <- myresult[[3]][3]
			datathing[i,"rf500_fpos"] <- myresult[[4]][3]
			datathing[i,"rf500_fneg"] <- myresult[[5]][3]
			
     }      
	return(datathing)            
}   


enhancer_class_status_scale <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize){
    enhancers_set <- strata_2var(crm_set, "status", crm_tab,2, 24, 21, samplesize)   
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
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(status ~., data = enhancers_training[,2:length(enhancers_training)], kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$status, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_false_off <- 100 * svm_radial_test_table[2,1]/sum(svm_radial_test_table[,1])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(status ~., data = enhancers_training[,2:length(enhancers_training)],kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$status, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_false_off <- 100 * svm_sigmoid_test_table[2,1]/sum(svm_sigmoid_test_table[,1])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
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
	prediction_list <- c(svm_radial_prediction_success, svm_sigmoid_prediction_success, rf_500_prediction_success)
	true_positives_list <- c(svm_rad_actual_on, svm_sig_actual_on, rf_500_actual_on)
	true_negatives_list <- c(svm_rad_true_off, svm_sig_true_off, rf_500_true_off)
	false_positives_list <- c(svm_rad_false_on, svm_sig_false_on, rf_500_false_on)
	false_negatives_list <- c(svm_rad_false_off, svm_sig_false_off, rf_500_false_off)
	name_list <- as.vector(enhancers_test$name)
	output_list <- cbind(name_list, actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(enhancers_test, file = summaryfile)
	write.csv(output_list, file = filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list, false_negatives_list)
	return(return_list)             
}       




Mac_dorsal_rf_500 <- randomForest(status ~., data = Mac_Dl_train_cons[,2:length(Mac_Dl_train_cons)], ntree = 500)

prediction_MDl <- predict(Mac_dorsal_rf_500, newdata = Mac_Dl_put, type = "class")
Macdl_on_frame <- cbind(Mac_Dl_put[,1:2], prediction_MDl)


write.table(Macdl_on_frame, "Mac_Dl_on_predict.tsv", sep="\t") 

Rdorsal_rf_500 <- randomForest(status ~., data = RDl_train_cons[,2:length(RDl_train_cons)], ntree = 500)

prediction_RDl <- predict(Rdorsal_rf_500, newdata = RDl_put, type = "class")
Rdl_on_frame <- cbind(RDl_put[,1:2], prediction_RDl)


write.table(Rdl_on_frame, "RDl_on_predict.tsv", sep="\t")


enhancer_class_Ectoderm <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize){
    enhancers_set <- strata_2var(dataset, "Ectoderm", datatable, 8, 32, 21, samplesize) 
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
	actual <- enhancers_test$Ectoderm
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(Ectoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$Ectoderm, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_false_off <- 100 * svm_radial_test_table[2,1]/sum(svm_radial_test_table[,1])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(Ectoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$Ectoderm, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_false_off <- 100 * svm_sigmoid_test_table[2,1]/sum(svm_sigmoid_test_table[,1])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(Ectoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$Ectoderm, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	prediction_list <- c(svm_radial_prediction_success, svm_sigmoid_prediction_success, rf_500_prediction_success)
	true_positives_list <- c(svm_rad_actual_on, svm_sig_actual_on, rf_500_actual_on)
	true_negatives_list <- c(svm_rad_true_off, svm_sig_true_off, rf_500_true_off)
	false_positives_list <- c(svm_rad_false_on, svm_sig_false_on, rf_500_false_on)
	false_negatives_list <- c(svm_rad_false_off, svm_sig_false_off, rf_500_false_off)
	output_list <- cbind(actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(output_list, filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list, false_negatives_list)
	return(return_list)             
}        

enhancer_class_Mesoderm <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize){
    enhancers_set <- strata_2var(dataset, "Mesoderm", datatable, 10, 32, 21, samplesize) 
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
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(Mesoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$Mesoderm, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_false_off <- 100 * svm_radial_test_table[2,1]/sum(svm_radial_test_table[,1])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(Mesoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$Mesoderm, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_false_off <- 100 * svm_sigmoid_test_table[2,1]/sum(svm_sigmoid_test_table[,1])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(Mesoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$Mesoderm, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	prediction_list <- c(svm_radial_prediction_success, svm_sigmoid_prediction_success, rf_500_prediction_success)
	true_positives_list <- c(svm_rad_actual_on, svm_sig_actual_on, rf_500_actual_on)
	true_negatives_list <- c(svm_rad_true_off, svm_sig_true_off, rf_500_true_off)
	false_positives_list <- c(svm_rad_false_on, svm_sig_false_on, rf_500_false_on)
	false_negatives_list <- c(svm_rad_false_off, svm_sig_false_off, rf_500_false_off)
	output_list <- cbind(actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(output_list, filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list, false_negatives_list)
	return(return_list)             
}  

enhancer_class_Endoderm <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize,variable_list){
    enhancers_set <- strata_2var(dataset, "Endoderm", datatable, 9, 32, 21, samplesize) 
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
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(Endoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$Endoderm, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_false_off <- 100 * svm_radial_test_table[2,1]/sum(svm_radial_test_table[,1])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(Endoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$Endoderm, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_false_off <- 100 * svm_sigmoid_test_table[2,1]/sum(svm_sigmoid_test_table[,1])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(Endoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$Endoderm, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	prediction_list <- c(svm_radial_prediction_success, svm_sigmoid_prediction_success, rf_500_prediction_success)
	true_positives_list <- c(svm_rad_actual_on, svm_sig_actual_on, rf_500_actual_on)
	true_negatives_list <- c(svm_rad_true_off, svm_sig_true_off, rf_500_true_off)
	false_positives_list <- c(svm_rad_false_on, svm_sig_false_on, rf_500_false_on)
	false_negatives_list <- c(svm_rad_false_off, svm_sig_false_off, rf_500_false_off)
	output_list <- cbind(actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(output_list, filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list, false_negatives_list)
	return(return_list)             
} 

enhancer_class_A <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize){
    enhancers_set <- strata_2var(dataset, "A", datatable, 5, 32, 21, samplesize) 
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
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(A ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$A, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_false_off <- 100 * svm_radial_test_table[2,1]/sum(svm_radial_test_table[,1])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(A ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$A, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_false_off <- 100 * svm_sigmoid_test_table[2,1]/sum(svm_sigmoid_test_table[,1])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(A ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$A, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	prediction_list <- c(svm_radial_prediction_success, svm_sigmoid_prediction_success, rf_500_prediction_success)
	true_positives_list <- c(svm_rad_actual_on, svm_sig_actual_on, rf_500_actual_on)
	true_negatives_list <- c(svm_rad_true_off, svm_sig_true_off, rf_500_true_off)
	false_positives_list <- c(svm_rad_false_on, svm_sig_false_on, rf_500_false_on)
	false_negatives_list <- c(svm_rad_false_off, svm_sig_false_off, rf_500_false_off)
	output_list <- cbind(actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(output_list, filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list, false_negatives_list)
	return(return_list)             
}        

enhancer_class_C <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize){
    enhancers_set <- strata_2var(dataset, "C", datatable, 6, 32, 21, samplesize) 
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
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(C ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$C, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_false_off <- 100 * svm_radial_test_table[2,1]/sum(svm_radial_test_table[,1])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(C ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$C, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_false_off <- 100 * svm_sigmoid_test_table[2,1]/sum(svm_sigmoid_test_table[,1])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(C ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$C, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	prediction_list <- c(svm_radial_prediction_success, svm_sigmoid_prediction_success, rf_500_prediction_success)
	true_positives_list <- c(svm_rad_actual_on, svm_sig_actual_on, rf_500_actual_on)
	true_negatives_list <- c(svm_rad_true_off, svm_sig_true_off, rf_500_true_off)
	false_positives_list <- c(svm_rad_false_on, svm_sig_false_on, rf_500_false_on)
	false_negatives_list <- c(svm_rad_false_off, svm_sig_false_off, rf_500_false_off)
	output_list <- cbind(actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(output_list, filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list, false_negatives_list)
	return(return_list)             
} 

enhancer_class_P <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize){
    enhancers_set <- strata_2var(dataset, "P", datatable, 7, 32, 21, samplesize) 
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
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(P ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$P, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_false_off <- 100 * svm_radial_test_table[2,1]/sum(svm_radial_test_table[,1])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(P ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$P, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_false_off <- 100 * svm_sigmoid_test_table[2,1]/sum(svm_sigmoid_test_table[,1])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(P ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$P, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_false_off <- 100 * rf_500_test_table[2,1]/sum(rf_500_test_table[,1])
	prediction_list <- c(svm_radial_prediction_success, svm_sigmoid_prediction_success, rf_500_prediction_success)
	true_positives_list <- c(svm_rad_actual_on, svm_sig_actual_on, rf_500_actual_on)
	true_negatives_list <- c(svm_rad_true_off, svm_sig_true_off, rf_500_true_off)
	false_positives_list <- c(svm_rad_false_on, svm_sig_false_on, rf_500_false_on)
	false_negatives_list <- c(svm_rad_false_off, svm_sig_false_off, rf_500_false_off)
	output_list <- cbind(actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(output_list, filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list, false_negatives_list)
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

 
colnames(RDl_put)[2] <- "Ectoderm"
balanced_ect_comp <- balanced_ect[,c(8,11:length(balanced_ect))]
Rdorsal_Ectoderm_rf_500 <- randomForest(Ectoderm ~., data = balanced_ect_comp, ntree = 500)
prediction_RDl_Ectoderm <- predict(Rdorsal_Ectoderm_rf_500, newdata = RDl_put, type = "class")
Rdl_Ectoderm_frame <- cbind(RDl_put[,1:2], prediction_RDl_Ectoderm)


colnames(RDl_put)[2] <- "Endoderm"
balanced_end_comp <- balanced_end[,c(9,11:length(balanced_end))]
Rdorsal_Endoderm_rf_500 <- randomForest(Endoderm ~., data = balanced_end_comp, ntree = 500)
prediction_RDl_Endoderm <- predict(Rdorsal_Endoderm_rf_500, newdata = RDl_put, type = "class")
Rdl_Endoderm_frame <- cbind(RDl_put[,1:2], prediction_RDl_Endoderm)

colnames(RDl_put)[2] <- "Mesoderm"
balanced_mes_comp <- balanced_mes[,c(10:length(balanced_mes))]
Rdorsal_Mesoderm_rf_500 <- randomForest(Mesoderm ~., data = balanced_mes_comp, ntree = 500)
prediction_RDl_Mesoderm <- predict(Rdorsal_Mesoderm_rf_500, newdata = RDl_put, type = "class")
Rdl_Mesoderm_frame <- cbind(RDl_put[,1:2], prediction_RDl_Mesoderm)

write.table(Rdl_Ectoderm_frame, "RDl_Ectoderm_predict.tsv", sep="\t")
write.table(Rdl_Endoderm_frame, "RDl_Endoderm_predict.tsv", sep="\t")
write.table(Rdl_Mesoderm_frame, "RDl_Mesoderm_predict.tsv", sep="\t")


resample_graph <- function(function_choice, dataset, tab_choice, run_name, readme_info, graph_title, reps, samplenum){
    mymodel <- resample_default(function_choice, dataset, tab_choice, run_name, readme_info, description = "README!", reps, samplenum)
	my_means <- apply(mymodel, MARGIN = 2, FUN = mean)
	my_sd <- apply(mymodel, MARGIN=2, FUN=sd)
	par(mgp = c(0,1,0))
	barplot(my_means, ylim=c(-6,100), col=c(33,"darkgreen","lightgreen","orange", "darkred"), , xaxt="n", main = graph_title)
	arrows(seq(.7, 18, 1.2), (my_means - my_sd), seq(.7, 18, 1.2), (my_means + my_sd), length=.05, angle=90, code=3)
	text(y=-4, x=3, "SVM Radial")
	text(y=-4, x=9, "SVM Sigmoid")
	text(y=-4, x=15, "RandomForest500")
	return_list <- list(my_means[12:15], my_sd[12:15])
	return_vals <- c(return_list[[1]], return_list[[2]])
	return(return_vals)
}



my_features <- drop_Zeit_Twist
drop_Zeit_results <- resample_graph(enhancer_class_Ectoderm, balanced_ect, balanced_ect_table, "Ectoderm_dropZeitTwi", "Ectoderm_dropZeitTwi", "Ectoderm Drop ZeitTwi", 20, 60)
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
drop_Dl_on <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_dropDl", "on_off_dropDl","On/Off Drop Dorsal resample 20", 20, 104)
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
on_all <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_all", "on_off_all","On/Off all features resample 100", 100, 60)
on_all
my_features <- drop_Zeit_Twist
on_drop_ZeitTwi <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_ZeitTwi", "on_off_drop_ZeitTwi","On/Off drop Zeitlinger Twi resample 100", 100, 60)
on_drop_ZeitTwi
my_features <- drop_White_H3K27ac
on_drop_White_H3K27ac <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_White_H3K27ac", "on_off_drop_White_H3K27ac","On/Off drop White H3K27ac resample 100", 100, 60)
on_drop_White_H3K27ac
my_features <- drop_MacDorsal 
on_drop_MacDorsal <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_MacDorsal", "on_off_drop_MacDorsal","On/Off drop MacArthur Dorsal resample 100", 100, 60)
on_drop_MacDorsal
my_features <- drop_MacTwist
on_drop_MacTwist <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_MacTwist", "on_off_drop_MacTwist","On/Off drop MacArthur Twist resample 100", 100, 60)
on_drop_MacTwist
my_features <- drop_MacSnail
on_drop_MacSnail <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_MacSnail", "on_off_drop_MacSnail","On/Off drop MacArthur Snail resample 100", 100, 60)
on_drop_MacSnail
my_features <- drop_MacHairy
on_drop_MacHairy <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_MacHairy", "on_off_drop_MacHairy","On/Off drop MacArthur Hairy resample 100", 100, 60)
on_drop_MacHairy
my_features <- drop_MacHunchback
on_drop_MacHunchback <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_MacHunchback", "on_off_drop_MacHunchback","On/Off drop MacArthur Hunchback resample 100", 100, 60)
on_drop_MacHunchback
my_features <- drop_MacGiant
on_drop_MacGiant <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_MacGiant", "on_off_drop_MacGiant","On/Off drop MacArthur Giant resample 100", 100, 60)
on_drop_MacGiant
my_features <- drop_MacCaudal
on_drop_MacCaudal <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_MacCaudal", "on_off_drop_MacCaudal","On/Off drop MacArthur Caudal resample 100", 100, 60)
on_drop_MacCaudal
my_features <- drop_MacBicoid 
on_drop_MacBicoid  <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_MacBicoid", "on_off_drop_MacBicoid","On/Off drop MacArthur Bicoid resample 100", 100, 60)
on_drop_MacBicoid 
my_features <- drop_KK_H3K4me1 
on_drop_KK_H3K4me1   <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_KK_H3K4me1", "on_off_drop_KK_H3K4me1", "On/Off drop KK H3K4me1 resample 100", 100, 60)
on_drop_KK_H3K4me1  
my_features <- drop_KK_H3K27ac
on_drop_KK_H3K27ac   <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_KK_H3K27ac", "on_off_drop_KK_H3K27ac", "On/Off drop KK H3K27ac resample 100", 100, 60)
on_drop_KK_H3K27ac 
my_features <- drop_FurlongTwist
on_drop_FurlongTwist  <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_FurlongTwist", "on_off_drop_FurlongTwist", "On/Off drop Furlong Twist resample 100", 100, 60)
on_drop_FurlongTwist 
my_features <- drop_FurlongSnail
on_drop_FurlongSnail  <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_FurlongSnail", "on_off_drop_FurlongSnail", "On/Off drop Furlong Snail resample 100", 100, 60)
on_drop_FurlongSnail 
my_features <- drop_Eisen13Kr 
on_drop_Eisen13Kr <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_Eisen13Kr", "on_off_drop_Eisen13Kr", "On/Off drop Eisen 2013 Kr resample 100", 100, 60)
on_drop_Eisen13Kr  
my_features <- drop_Eisen13Bcd
on_drop_Eisen13Bcd <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_Eisen13Bcd", "on_off_drop_Eisen13Bcd", "On/Off drop Eisen 2013 Bcd resample 100", 100, 60)
on_drop_Eisen13Bcd
my_features <- drop_Eisen13Hb
on_drop_Eisen13Hb <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_Eisen13Hb", "on_off_drop_Eisen13Hb", "On/Off drop Eisen 2013 Hb resample 100", 100, 60)
on_drop_Eisen13Hb
my_features <- drop_Eisen13Gt
on_drop_Eisen13Gt <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_Eisen13Gt", "on_off_drop_Eisen13Gt", "On/Off drop Eisen 2013 Gt resample 100", 100, 60)
on_drop_Eisen13Gt
my_features <- drop_Eisen10Kni
on_drop_Eisen10Kni <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_Eisen10Kni", "on_off_drop_Eisen10Kni", "On/Off drop Eisen 2010 Kni resample 100", 100, 60)
on_drop_Eisen10Kni
my_features <- drop_Eisen10Cad
on_drop_Eisen10Cad <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_Eisen10Cad", "on_off_drop_Eisen10Cad", "On/Off drop Eisen 2010 Cad resample 100", 100, 60)
on_drop_Eisen10Cad
my_features <- drop_Zelda
on_drop_Zelda <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_drop_Zelda", "on_off_drop_Zelda", "On/Off drop Zelda resample 100", 100, 60)
on_drop_Zelda
my_features <- drop_Dl 
my_features <- drop_Bcd 
my_features <- drop_Cad 
my_features <- only_Zelda
on_only_Zelda <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_Zelda", "on_off_only_Zelda", "On/Off only Zelda resample 100", 100, 60)
on_only_Zelda 
my_features <- only_Eisen10Cad 
on_only_Eisen10Cad <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_Eisen10Cad", "on_off_only_Eisen10Cad", "On/Off only Eisen10Cad resample 100", 100, 60)
on_only_Eisen10Cad 
my_features <- only_Eisen10Kni 
on_only_Eisen10Kni <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_Eisen10Kni", "on_off_only_Eisen10Kni", "On/Off only Eisen10Kni resample 100", 100, 60)
on_only_Eisen10Kni
my_features <- only_Eisen13Gt
on_only_Eisen10Gt <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_Eisen10Gt", "on_off_only_Eisen10Gt", "On/Off only Eisen10Gt resample 100", 100, 60)
on_only_Eisen10Gt
my_features <- only_Eisen13Hb
on_only_Eisen13Hb <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_Eisen13Hb", "on_off_only_Eisen13Hb", "On/Off only Eisen13Hb resample 100", 100, 60)
on_only_Eisen13Hb
my_features <- only_Eisen13Bcd
on_only_Eisen13Bcd <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_Eisen13Bcd", "on_off_only_Eisen13Bcd", "On/Off only Eisen13Bcd resample 100", 100, 60)
on_only_Eisen13Bcd 
my_features <- only_Eisen13Kr
on_only_Eisen13Kr <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_Eisen13Kr", "on_off_only_Eisen13Kr", "On/Off only Eisen13Kr resample 100", 100, 60)
on_only_Eisen13Kr 
my_features <- only_FurlongSnail
on_only_FurlongSnail <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_FurlongSnail", "on_off_only_FurlongSnail", "On/Off only FurlongSnail resample 100", 100, 60)
on_only_FurlongSnail 
my_features <- only_FurlongTwist
on_only_FurlongTwist <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_FurlongTwist", "on_off_only_FurlongTwist", "On/Off only FurlongTwist resample 100", 100, 60)
on_only_FurlongTwist
my_features <- only_KK
on_only_KK <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_KK", "on_off_only_KK", "On/Off only Kurtulus resample 100", 100, 60)
on_only_KK
my_features <- only_MacBicoid
on_only_MacBicoid <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_MacBicoid", "on_off_only_MacBicoid", "On/Off only MacBicoid resample 100", 100, 60)
on_only_MacBicoid
my_features <- only_MacCad
on_only_MacCad <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_MacCad", "on_off_only_MacCad", "On/Off only MacCad resample 100", 100, 60)
on_only_MacCad 
my_features <- only_MacGiant
on_only_MacGiant <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_MacGiant", "on_off_only_MacGiant", "On/Off only MacGiant resample 100", 100, 60)
on_only_MacGiant  
my_features <- only_MacHb
on_only_MacHb <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_MacHb", "on_off_only_MacHb", "On/Off only MacHb resample 100", 100, 60)
on_only_MacHb
my_features <- only_MacHry
on_only_MacHry <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_MacHry", "on_off_only_MacHry", "On/Off only MacHry resample 100", 100, 60)
on_only_MacHry 
my_features <- only_MacSna
on_only_MacSna <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_MacSna", "on_off_only_MacSna", "On/Off only MacSna resample 100", 100, 60)
on_only_MacSna 
my_features <- only_MacTwi
on_only_MacTwi <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_MacTwi", "on_off_only_MacTwi", "On/Off only MacTwi resample 100", 100, 60)
on_only_MacTwi 
my_features <- only_MacDl
on_only_MacDl <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_MacDl", "on_off_only_MacDl", "On/Off only MacDl resample 100", 100, 60)
on_only_MacDl  
my_features <- only_White
on_only_White <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_White", "on_off_only_White", "On/Off only White resample 100", 100, 60)
on_only_White 
my_features <- only_Zeit
on_only_Zeit <- resample_graph(enhancer_class_status_scale,RDl_train_cons, tab_RDl_cons, "on_off_only_Zeit", "on_off_only_Zeit", "On/Off only Zeit resample 100", 100, 60)
on_only_Zeit 


#Ectoderm
my_features <- all_features
Ectoderm_all <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_all", "Ectoderm_all","Ectoderm all features resample 100", 100, 104)
Ectoderm_all
my_features <- drop_Zeit_Twist
Ectoderm_drop_ZeitTwi <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_ZeitTwi", "Ectoderm_drop_ZeitTwi","Ectoderm drop Zeitlinger Twi resample 100", 100, 104)
Ectoderm_drop_ZeitTwi
my_features <- drop_White_H3K27ac
Ectoderm_drop_White_H3K27ac <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_White_H3K27ac", "Ectoderm_drop_White_H3K27ac","Ectoderm drop White H3K27ac resample 100", 100, 104)
Ectoderm_drop_White_H3K27ac
my_features <- drop_MacDorsal 
Ectoderm_drop_MacDorsal <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_MacDorsal", "Ectoderm_drop_MacDorsal","Ectoderm drop MacArthur Dorsal resample 100", 100, 104)
Ectoderm_drop_MacDorsal
my_features <- drop_MacTwist
Ectoderm_drop_MacTwist <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_MacTwist", "Ectoderm_drop_MacTwist","Ectoderm drop MacArthur Twist resample 100", 100, 104)
Ectoderm_drop_MacTwist
my_features <- drop_MacSnail
Ectoderm_drop_MacSnail <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_MacSnail", "Ectoderm_drop_MacSnail","Ectoderm drop MacArthur Snail resample 100", 100, 104)
Ectoderm_drop_MacSnail
my_features <- drop_MacHairy
Ectoderm_drop_MacHairy <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_MacHairy", "Ectoderm_drop_MacHairy","Ectoderm drop MacArthur Hairy resample 100", 100, 104)
Ectoderm_drop_MacHairy
my_features <- drop_MacHunchback
Ectoderm_drop_MacHunchback <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_MacHunchback", "Ectoderm_drop_MacHunchback","Ectoderm drop MacArthur Hunchback resample 100", 100, 104)
Ectoderm_drop_MacHunchback
my_features <- drop_MacGiant
Ectoderm_drop_MacGiant <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_MacGiant", "Ectoderm_drop_MacGiant","Ectoderm drop MacArthur Giant resample 100", 100, 104)
Ectoderm_drop_MacGiant
my_features <- drop_MacCaudal
Ectoderm_drop_MacCaudal <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_MacCaudal", "Ectoderm_drop_MacCaudal","Ectoderm drop MacArthur Caudal resample 100", 100, 104)
Ectoderm_drop_MacCaudal
my_features <- drop_MacBicoid 
Ectoderm_drop_MacBicoid  <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_MacBicoid", "Ectoderm_drop_MacBicoid","Ectoderm drop MacArthur Bicoid resample 100", 100, 104)
Ectoderm_drop_MacBicoid 
my_features <- drop_KK_H3K4me1 
Ectoderm_drop_KK_H3K4me1   <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_KK_H3K4me1", "Ectoderm_drop_KK_H3K4me1", "Ectoderm drop KK H3K4me1 resample 100", 100, 104)
Ectoderm_drop_KK_H3K4me1  
my_features <- drop_KK_H3K27ac
Ectoderm_drop_KK_H3K27ac   <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_KK_H3K27ac", "Ectoderm_drop_KK_H3K27ac", "Ectoderm drop KK H3K27ac resample 100", 100, 104)
Ectoderm_drop_KK_H3K27ac 
my_features <- drop_FurlongTwist
Ectoderm_drop_FurlongTwist  <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_FurlongTwist", "Ectoderm_drop_FurlongTwist", "Ectoderm drop Furlong Twist resample 100", 100, 104)
Ectoderm_drop_FurlongTwist 
my_features <- drop_FurlongSnail
Ectoderm_drop_FurlongSnail  <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_FurlongSnail", "Ectoderm_drop_FurlongSnail", "Ectoderm drop Furlong Snail resample 100", 100, 104)
Ectoderm_drop_FurlongSnail 
my_features <- drop_Eisen13Kr 
Ectoderm_drop_Eisen13Kr <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_Eisen13Kr", "Ectoderm_drop_Eisen13Kr", "Ectoderm drop Eisen 2013 Kr resample 100", 100, 104)
Ectoderm_drop_Eisen13Kr  
my_features <- drop_Eisen13Bcd
Ectoderm_drop_Eisen13Bcd <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_Eisen13Bcd", "Ectoderm_drop_Eisen13Bcd", "Ectoderm drop Eisen 2013 Bcd resample 100", 100, 104)
Ectoderm_drop_Eisen13Bcd
my_features <- drop_Eisen13Hb
Ectoderm_drop_Eisen13Hb <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_Eisen13Hb", "Ectoderm_drop_Eisen13Hb", "Ectoderm drop Eisen 2013 Hb resample 100", 100, 104)
Ectoderm_drop_Eisen13Hb
my_features <- drop_Eisen13Gt
Ectoderm_drop_Eisen13Gt <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_Eisen13Gt", "Ectoderm_drop_Eisen13Gt", "Ectoderm drop Eisen 2013 Gt resample 100", 100, 104)
Ectoderm_drop_Eisen13Gt
my_features <- drop_Eisen10Kni
Ectoderm_drop_Eisen10Kni <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_Eisen10Kni", "Ectoderm_drop_Eisen10Kni", "Ectoderm drop Eisen 2010 Kni resample 100", 100, 104)
Ectoderm_drop_Eisen10Kni
my_features <- drop_Eisen10Cad
Ectoderm_drop_Eisen10Cad <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_Eisen10Cad", "Ectoderm_drop_Eisen10Cad", "Ectoderm drop Eisen 2010 Cad resample 100", 100, 104)
Ectoderm_drop_Eisen10Cad
my_features <- drop_Zelda
Ectoderm_drop_Zelda <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_drop_Zelda", "Ectoderm_drop_Zelda", "Ectoderm drop Zelda resample 100", 100, 104)
Ectoderm_drop_Zelda
my_features <- drop_Dl 
my_features <- drop_Bcd 
my_features <- drop_Cad 
my_features <- only_Zelda
Ectoderm_only_Zelda <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_Zelda", "Ectoderm_only_Zelda", "Ectoderm only Zelda resample 100", 100, 104)
Ectoderm_only_Zelda 
my_features <- only_Eisen10Cad 
Ectoderm_only_Eisen10Cad <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_Eisen10Cad", "Ectoderm_only_Eisen10Cad", "Ectoderm only Eisen10Cad resample 100", 100, 104)
Ectoderm_only_Eisen10Cad 
my_features <- only_Eisen10Kni 
Ectoderm_only_Eisen10Kni <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_Eisen10Kni", "Ectoderm_only_Eisen10Kni", "Ectoderm only Eisen10Kni resample 100", 100, 104)
Ectoderm_only_Eisen10Kni
my_features <- only_Eisen13Gt
Ectoderm_only_Eisen10Gt <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_Eisen10Gt", "Ectoderm_only_Eisen10Gt", "Ectoderm only Eisen10Gt resample 100", 100, 104)
Ectoderm_only_Eisen10Gt
my_features <- only_Eisen13Hb
Ectoderm_only_Eisen13Hb <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_Eisen13Hb", "Ectoderm_only_Eisen13Hb", "Ectoderm only Eisen13Hb resample 100", 100, 104)
Ectoderm_only_Eisen13Hb
my_features <- only_Eisen13Bcd
Ectoderm_only_Eisen13Bcd <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_Eisen13Bcd", "Ectoderm_only_Eisen13Bcd", "Ectoderm only Eisen13Bcd resample 100", 100, 104)
Ectoderm_only_Eisen13Bcd 
my_features <- only_Eisen13Kr
Ectoderm_only_Eisen13Kr <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_Eisen13Kr", "Ectoderm_only_Eisen13Kr", "Ectoderm only Eisen13Kr resample 100", 100, 104)
Ectoderm_only_Eisen13Kr 
my_features <- only_FurlongSnail
Ectoderm_only_FurlongSnail <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_FurlongSnail", "Ectoderm_only_FurlongSnail", "Ectoderm only FurlongSnail resample 100", 100, 104)
Ectoderm_only_FurlongSnail 
my_features <- only_FurlongTwist
Ectoderm_only_FurlongTwist <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_FurlongTwist", "Ectoderm_only_FurlongTwist", "Ectoderm only FurlongTwist resample 100", 100, 104)
Ectoderm_only_FurlongTwist
my_features <- only_KK
Ectoderm_only_KK <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_KK", "Ectoderm_only_KK", "Ectoderm only Kurtulus resample 100", 100, 104)
Ectoderm_only_KK
my_features <- only_MacBicoid
Ectoderm_only_MacBicoid <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_MacBicoid", "Ectoderm_only_MacBicoid", "Ectoderm only MacBicoid resample 100", 100, 104)
Ectoderm_only_MacBicoid
my_features <- only_MacCad
Ectoderm_only_MacCad <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_MacCad", "Ectoderm_only_MacCad", "Ectoderm only MacCad resample 100", 100, 104)
Ectoderm_only_MacCad 
my_features <- only_MacGiant
Ectoderm_only_MacGiant <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_MacGiant", "Ectoderm_only_MacGiant", "Ectoderm only MacGiant resample 100", 100, 104)
Ectoderm_only_MacGiant  
my_features <- only_MacHb
Ectoderm_only_MacHb <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_MacHb", "Ectoderm_only_MacHb", "Ectoderm only MacHb resample 100", 100, 104)
Ectoderm_only_MacHb
my_features <- only_MacHry
Ectoderm_only_MacHry <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_MacHry", "Ectoderm_only_MacHry", "Ectoderm only MacHry resample 100", 100, 104)
Ectoderm_only_MacHry 
my_features <- only_MacSna
Ectoderm_only_MacSna <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_MacSna", "Ectoderm_only_MacSna", "Ectoderm only MacSna resample 100", 100, 104)
Ectoderm_only_MacSna 
my_features <- only_MacTwi
Ectoderm_only_MacTwi <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_MacTwi", "Ectoderm_only_MacTwi", "Ectoderm only MacTwi resample 100", 100, 104)
Ectoderm_only_MacTwi 
my_features <- only_MacDl
Ectoderm_only_MacDl <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_MacDl", "Ectoderm_only_MacDl", "Ectoderm only MacDl resample 100", 100, 104)
Ectoderm_only_MacDl  
my_features <- only_White
Ectoderm_only_White <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_White", "Ectoderm_only_White", "Ectoderm only White resample 100", 100, 104)
Ectoderm_only_White 
my_features <- only_Zeit
Ectoderm_only_Zeit <- resample_graph(enhancer_class_Ectoderm,balanced_ect, balanced_ect_table, "Ectoderm_only_Zeit", "Ectoderm_only_Zeit", "Ectoderm only Zeit resample 100", 100, 104)
Ectoderm_only_Zeit 

#Mesoderm
my_features <- all_features
Mesoderm_all <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_all", "Mesoderm_all","Mesoderm all features resample 100", 100, 104)
Mesoderm_all
my_features <- drop_Zeit_Twist
Mesoderm_drop_ZeitTwi <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_ZeitTwi", "Mesoderm_drop_ZeitTwi","Mesoderm drop Zeitlinger Twi resample 100", 100, 104)
Mesoderm_drop_ZeitTwi
my_features <- drop_White_H3K27ac
Mesoderm_drop_White_H3K27ac <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_White_H3K27ac", "Mesoderm_drop_White_H3K27ac","Mesoderm drop White H3K27ac resample 100", 100, 104)
Mesoderm_drop_White_H3K27ac
my_features <- drop_MacDorsal 
Mesoderm_drop_MacDorsal <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_MacDorsal", "Mesoderm_drop_MacDorsal","Mesoderm drop MacArthur Dorsal resample 100", 100, 104)
Mesoderm_drop_MacDorsal
my_features <- drop_MacTwist
Mesoderm_drop_MacTwist <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_MacTwist", "Mesoderm_drop_MacTwist","Mesoderm drop MacArthur Twist resample 100", 100, 104)
Mesoderm_drop_MacTwist
my_features <- drop_MacSnail
Mesoderm_drop_MacSnail <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_MacSnail", "Mesoderm_drop_MacSnail","Mesoderm drop MacArthur Snail resample 100", 100, 104)
Mesoderm_drop_MacSnail
my_features <- drop_MacHairy
Mesoderm_drop_MacHairy <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_MacHairy", "Mesoderm_drop_MacHairy","Mesoderm drop MacArthur Hairy resample 100", 100, 104)
Mesoderm_drop_MacHairy
my_features <- drop_MacHunchback
Mesoderm_drop_MacHunchback <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_MacHunchback", "Mesoderm_drop_MacHunchback","Mesoderm drop MacArthur Hunchback resample 100", 100, 104)
Mesoderm_drop_MacHunchback
my_features <- drop_MacGiant
Mesoderm_drop_MacGiant <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_MacGiant", "Mesoderm_drop_MacGiant","Mesoderm drop MacArthur Giant resample 100", 100, 104)
Mesoderm_drop_MacGiant
my_features <- drop_MacCaudal
Mesoderm_drop_MacCaudal <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_MacCaudal", "Mesoderm_drop_MacCaudal","Mesoderm drop MacArthur Caudal resample 100", 100, 104)
Mesoderm_drop_MacCaudal
my_features <- drop_MacBicoid 
Mesoderm_drop_MacBicoid  <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_MacBicoid", "Mesoderm_drop_MacBicoid","Mesoderm drop MacArthur Bicoid resample 100", 100, 104)
Mesoderm_drop_MacBicoid 
my_features <- drop_KK_H3K4me1 
Mesoderm_drop_KK_H3K4me1   <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_KK_H3K4me1", "Mesoderm_drop_KK_H3K4me1", "Mesoderm drop KK H3K4me1 resample 100", 100, 104)
Mesoderm_drop_KK_H3K4me1  
my_features <- drop_KK_H3K27ac
Mesoderm_drop_KK_H3K27ac   <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_KK_H3K27ac", "Mesoderm_drop_KK_H3K27ac", "Mesoderm drop KK H3K27ac resample 100", 100, 104)
Mesoderm_drop_KK_H3K27ac 
my_features <- drop_FurlongTwist
Mesoderm_drop_FurlongTwist  <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_FurlongTwist", "Mesoderm_drop_FurlongTwist", "Mesoderm drop Furlong Twist resample 100", 100, 104)
Mesoderm_drop_FurlongTwist 
my_features <- drop_FurlongSnail
Mesoderm_drop_FurlongSnail  <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_FurlongSnail", "Mesoderm_drop_FurlongSnail", "Mesoderm drop Furlong Snail resample 100", 100, 104)
Mesoderm_drop_FurlongSnail 
my_features <- drop_Eisen13Kr 
Mesoderm_drop_Eisen13Kr <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_Eisen13Kr", "Mesoderm_drop_Eisen13Kr", "Mesoderm drop Eisen 2013 Kr resample 100", 100, 104)
Mesoderm_drop_Eisen13Kr  
my_features <- drop_Eisen13Bcd
Mesoderm_drop_Eisen13Bcd <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_Eisen13Bcd", "Mesoderm_drop_Eisen13Bcd", "Mesoderm drop Eisen 2013 Bcd resample 100", 100, 104)
Mesoderm_drop_Eisen13Bcd
my_features <- drop_Eisen13Hb
Mesoderm_drop_Eisen13Hb <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_Eisen13Hb", "Mesoderm_drop_Eisen13Hb", "Mesoderm drop Eisen 2013 Hb resample 100", 100, 104)
Mesoderm_drop_Eisen13Hb
my_features <- drop_Eisen13Gt
Mesoderm_drop_Eisen13Gt <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_Eisen13Gt", "Mesoderm_drop_Eisen13Gt", "Mesoderm drop Eisen 2013 Gt resample 100", 100, 104)
Mesoderm_drop_Eisen13Gt
my_features <- drop_Eisen10Kni
Mesoderm_drop_Eisen10Kni <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_Eisen10Kni", "Mesoderm_drop_Eisen10Kni", "Mesoderm drop Eisen 2010 Kni resample 100", 100, 104)
Mesoderm_drop_Eisen10Kni
my_features <- drop_Eisen10Cad
Mesoderm_drop_Eisen10Cad <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_Eisen10Cad", "Mesoderm_drop_Eisen10Cad", "Mesoderm drop Eisen 2010 Cad resample 100", 100, 104)
Mesoderm_drop_Eisen10Cad
my_features <- drop_Zelda
Mesoderm_drop_Zelda <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_drop_Zelda", "Mesoderm_drop_Zelda", "Mesoderm drop Zelda resample 100", 100, 104)
Mesoderm_drop_Zelda
my_features <- drop_Dl 
my_features <- drop_Bcd 
my_features <- drop_Cad 
my_features <- only_Zelda
Mesoderm_only_Zelda <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_Zelda", "Mesoderm_only_Zelda", "Mesoderm only Zelda resample 100", 100, 104)
Mesoderm_only_Zelda 
my_features <- only_Eisen10Cad 
Mesoderm_only_Eisen10Cad <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_Eisen10Cad", "Mesoderm_only_Eisen10Cad", "Mesoderm only Eisen10Cad resample 100", 100, 104)
Mesoderm_only_Eisen10Cad 
my_features <- only_Eisen10Kni 
Mesoderm_only_Eisen10Kni <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_Eisen10Kni", "Mesoderm_only_Eisen10Kni", "Mesoderm only Eisen10Kni resample 100", 100, 104)
Mesoderm_only_Eisen10Kni
my_features <- only_Eisen13Gt
Mesoderm_only_Eisen10Gt <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_Eisen10Gt", "Mesoderm_only_Eisen10Gt", "Mesoderm only Eisen10Gt resample 100", 100, 104)
Mesoderm_only_Eisen10Gt
my_features <- only_Eisen13Hb
Mesoderm_only_Eisen13Hb <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_Eisen13Hb", "Mesoderm_only_Eisen13Hb", "Mesoderm only Eisen13Hb resample 100", 100, 104)
Mesoderm_only_Eisen13Hb
my_features <- only_Eisen13Bcd
Mesoderm_only_Eisen13Bcd <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_Eisen13Bcd", "Mesoderm_only_Eisen13Bcd", "Mesoderm only Eisen13Bcd resample 100", 100, 104)
Mesoderm_only_Eisen13Bcd 
my_features <- only_Eisen13Kr
Mesoderm_only_Eisen13Kr <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_Eisen13Kr", "Mesoderm_only_Eisen13Kr", "Mesoderm only Eisen13Kr resample 100", 100, 104)
Mesoderm_only_Eisen13Kr 
my_features <- only_FurlongSnail
Mesoderm_only_FurlongSnail <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_FurlongSnail", "Mesoderm_only_FurlongSnail", "Mesoderm only FurlongSnail resample 100", 100, 104)
Mesoderm_only_FurlongSnail 
my_features <- only_FurlongTwist
Mesoderm_only_FurlongTwist <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_FurlongTwist", "Mesoderm_only_FurlongTwist", "Mesoderm only FurlongTwist resample 100", 100, 104)
Mesoderm_only_FurlongTwist
my_features <- only_KK
Mesoderm_only_KK <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_KK", "Mesoderm_only_KK", "Mesoderm only Kurtulus resample 100", 100, 104)
Mesoderm_only_KK
my_features <- only_MacBicoid
Mesoderm_only_MacBicoid <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_MacBicoid", "Mesoderm_only_MacBicoid", "Mesoderm only MacBicoid resample 100", 100, 104)
Mesoderm_only_MacBicoid
my_features <- only_MacCad
Mesoderm_only_MacCad <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_MacCad", "Mesoderm_only_MacCad", "Mesoderm only MacCad resample 100", 100, 104)
Mesoderm_only_MacCad 
my_features <- only_MacGiant
Mesoderm_only_MacGiant <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_MacGiant", "Mesoderm_only_MacGiant", "Mesoderm only MacGiant resample 100", 100, 104)
Mesoderm_only_MacGiant  
my_features <- only_MacHb
Mesoderm_only_MacHb <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_MacHb", "Mesoderm_only_MacHb", "Mesoderm only MacHb resample 100", 100, 104)
Mesoderm_only_MacHb
my_features <- only_MacHry
Mesoderm_only_MacHry <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_MacHry", "Mesoderm_only_MacHry", "Mesoderm only MacHry resample 100", 100, 104)
Mesoderm_only_MacHry 
my_features <- only_MacSna
Mesoderm_only_MacSna <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_MacSna", "Mesoderm_only_MacSna", "Mesoderm only MacSna resample 100", 100, 104)
Mesoderm_only_MacSna 
my_features <- only_MacTwi
Mesoderm_only_MacTwi <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_MacTwi", "Mesoderm_only_MacTwi", "Mesoderm only MacTwi resample 100", 100, 104)
Mesoderm_only_MacTwi 
my_features <- only_MacDl
Mesoderm_only_MacDl <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_MacDl", "Mesoderm_only_MacDl", "Mesoderm only MacDl resample 100", 100, 104)
Mesoderm_only_MacDl  
my_features <- only_White
Mesoderm_only_White <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_White", "Mesoderm_only_White", "Mesoderm only White resample 100", 100, 104)
Mesoderm_only_White 
my_features <- only_Zeit
Mesoderm_only_Zeit <- resample_graph(enhancer_class_Mesoderm,balanced_mes, balanced_mes_table, "Mesoderm_only_Zeit", "Mesoderm_only_Zeit", "Mesoderm only Zeit resample 100", 100, 104)
Mesoderm_only_Zeit

#Endoderm
my_features <- all_features
Endoderm_all <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_all", "Endoderm_all","Endoderm all features resample 100", 100, 104)
Endoderm_all
my_features <- drop_Zeit_Twist
Endoderm_drop_ZeitTwi <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_ZeitTwi", "Endoderm_drop_ZeitTwi","Endoderm drop Zeitlinger Twi resample 100", 100, 104)
Endoderm_drop_ZeitTwi
my_features <- drop_White_H3K27ac
Endoderm_drop_White_H3K27ac <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_White_H3K27ac", "Endoderm_drop_White_H3K27ac","Endoderm drop White H3K27ac resample 100", 100, 104)
Endoderm_drop_White_H3K27ac
my_features <- drop_MacDorsal 
Endoderm_drop_MacDorsal <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_MacDorsal", "Endoderm_drop_MacDorsal","Endoderm drop MacArthur Dorsal resample 100", 100, 104)
Endoderm_drop_MacDorsal
my_features <- drop_MacTwist
Endoderm_drop_MacTwist <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_MacTwist", "Endoderm_drop_MacTwist","Endoderm drop MacArthur Twist resample 100", 100, 104)
Endoderm_drop_MacTwist
my_features <- drop_MacSnail
Endoderm_drop_MacSnail <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_MacSnail", "Endoderm_drop_MacSnail","Endoderm drop MacArthur Snail resample 100", 100, 104)
Endoderm_drop_MacSnail
my_features <- drop_MacHairy
Endoderm_drop_MacHairy <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_MacHairy", "Endoderm_drop_MacHairy","Endoderm drop MacArthur Hairy resample 100", 100, 104)
Endoderm_drop_MacHairy
my_features <- drop_MacHunchback
Endoderm_drop_MacHunchback <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_MacHunchback", "Endoderm_drop_MacHunchback","Endoderm drop MacArthur Hunchback resample 100", 100, 104)
Endoderm_drop_MacHunchback
my_features <- drop_MacGiant
Endoderm_drop_MacGiant <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_MacGiant", "Endoderm_drop_MacGiant","Endoderm drop MacArthur Giant resample 100", 100, 104)
Endoderm_drop_MacGiant
my_features <- drop_MacCaudal
Endoderm_drop_MacCaudal <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_MacCaudal", "Endoderm_drop_MacCaudal","Endoderm drop MacArthur Caudal resample 100", 100, 104)
Endoderm_drop_MacCaudal
my_features <- drop_MacBicoid 
Endoderm_drop_MacBicoid  <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_MacBicoid", "Endoderm_drop_MacBicoid","Endoderm drop MacArthur Bicoid resample 100", 100, 104)
Endoderm_drop_MacBicoid 
my_features <- drop_KK_H3K4me1 
Endoderm_drop_KK_H3K4me1   <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_KK_H3K4me1", "Endoderm_drop_KK_H3K4me1", "Endoderm drop KK H3K4me1 resample 100", 100, 104)
Endoderm_drop_KK_H3K4me1  
my_features <- drop_KK_H3K27ac
Endoderm_drop_KK_H3K27ac   <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_KK_H3K27ac", "Endoderm_drop_KK_H3K27ac", "Endoderm drop KK H3K27ac resample 100", 100, 104)
Endoderm_drop_KK_H3K27ac 
my_features <- drop_FurlongTwist
Endoderm_drop_FurlongTwist  <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_FurlongTwist", "Endoderm_drop_FurlongTwist", "Endoderm drop Furlong Twist resample 100", 100, 104)
Endoderm_drop_FurlongTwist 
my_features <- drop_FurlongSnail
Endoderm_drop_FurlongSnail  <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_FurlongSnail", "Endoderm_drop_FurlongSnail", "Endoderm drop Furlong Snail resample 100", 100, 104)
Endoderm_drop_FurlongSnail 
my_features <- drop_Eisen13Kr 
Endoderm_drop_Eisen13Kr <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_Eisen13Kr", "Endoderm_drop_Eisen13Kr", "Endoderm drop Eisen 2013 Kr resample 100", 100, 104)
Endoderm_drop_Eisen13Kr  
my_features <- drop_Eisen13Bcd
Endoderm_drop_Eisen13Bcd <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_Eisen13Bcd", "Endoderm_drop_Eisen13Bcd", "Endoderm drop Eisen 2013 Bcd resample 100", 100, 104)
Endoderm_drop_Eisen13Bcd
my_features <- drop_Eisen13Hb
Endoderm_drop_Eisen13Hb <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_Eisen13Hb", "Endoderm_drop_Eisen13Hb", "Endoderm drop Eisen 2013 Hb resample 100", 100, 104)
Endoderm_drop_Eisen13Hb
my_features <- drop_Eisen13Gt
Endoderm_drop_Eisen13Gt <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_Eisen13Gt", "Endoderm_drop_Eisen13Gt", "Endoderm drop Eisen 2013 Gt resample 100", 100, 104)
Endoderm_drop_Eisen13Gt
my_features <- drop_Eisen10Kni
Endoderm_drop_Eisen10Kni <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_Eisen10Kni", "Endoderm_drop_Eisen10Kni", "Endoderm drop Eisen 2010 Kni resample 100", 100, 104)
Endoderm_drop_Eisen10Kni
my_features <- drop_Eisen10Cad
Endoderm_drop_Eisen10Cad <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_Eisen10Cad", "Endoderm_drop_Eisen10Cad", "Endoderm drop Eisen 2010 Cad resample 100", 100, 104)
Endoderm_drop_Eisen10Cad
my_features <- drop_Zelda
Endoderm_drop_Zelda <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_drop_Zelda", "Endoderm_drop_Zelda", "Endoderm drop Zelda resample 100", 100, 104)
Endoderm_drop_Zelda
my_features <- drop_Dl 
my_features <- drop_Bcd 
my_features <- drop_Cad 
my_features <- only_Zelda
Endoderm_only_Zelda <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_Zelda", "Endoderm_only_Zelda", "Endoderm only Zelda resample 100", 100, 104)
Endoderm_only_Zelda 
my_features <- only_Eisen10Cad 
Endoderm_only_Eisen10Cad <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_Eisen10Cad", "Endoderm_only_Eisen10Cad", "Endoderm only Eisen10Cad resample 100", 100, 104)
Endoderm_only_Eisen10Cad 
my_features <- only_Eisen10Kni 
Endoderm_only_Eisen10Kni <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_Eisen10Kni", "Endoderm_only_Eisen10Kni", "Endoderm only Eisen10Kni resample 100", 100, 104)
Endoderm_only_Eisen10Kni
my_features <- only_Eisen13Gt
Endoderm_only_Eisen10Gt <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_Eisen10Gt", "Endoderm_only_Eisen10Gt", "Endoderm only Eisen10Gt resample 100", 100, 104)
Endoderm_only_Eisen10Gt
my_features <- only_Eisen13Hb
Endoderm_only_Eisen13Hb <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_Eisen13Hb", "Endoderm_only_Eisen13Hb", "Endoderm only Eisen13Hb resample 100", 100, 104)
Endoderm_only_Eisen13Hb
my_features <- only_Eisen13Bcd
Endoderm_only_Eisen13Bcd <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_Eisen13Bcd", "Endoderm_only_Eisen13Bcd", "Endoderm only Eisen13Bcd resample 100", 100, 104)
Endoderm_only_Eisen13Bcd 
my_features <- only_Eisen13Kr
Endoderm_only_Eisen13Kr <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_Eisen13Kr", "Endoderm_only_Eisen13Kr", "Endoderm only Eisen13Kr resample 100", 100, 104)
Endoderm_only_Eisen13Kr 
my_features <- only_FurlongSnail
Endoderm_only_FurlongSnail <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_FurlongSnail", "Endoderm_only_FurlongSnail", "Endoderm only FurlongSnail resample 100", 100, 104)
Endoderm_only_FurlongSnail 
my_features <- only_FurlongTwist
Endoderm_only_FurlongTwist <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_FurlongTwist", "Endoderm_only_FurlongTwist", "Endoderm only FurlongTwist resample 100", 100, 104)
Endoderm_only_FurlongTwist
my_features <- only_KK
Endoderm_only_KK <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_KK", "Endoderm_only_KK", "Endoderm only Kurtulus resample 100", 100, 104)
Endoderm_only_KK
my_features <- only_MacBicoid
Endoderm_only_MacBicoid <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_MacBicoid", "Endoderm_only_MacBicoid", "Endoderm only MacBicoid resample 100", 100, 104)
Endoderm_only_MacBicoid
my_features <- only_MacCad
Endoderm_only_MacCad <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_MacCad", "Endoderm_only_MacCad", "Endoderm only MacCad resample 100", 100, 104)
Endoderm_only_MacCad 
my_features <- only_MacGiant
Endoderm_only_MacGiant <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_MacGiant", "Endoderm_only_MacGiant", "Endoderm only MacGiant resample 100", 100, 104)
Endoderm_only_MacGiant  
my_features <- only_MacHb
Endoderm_only_MacHb <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_MacHb", "Endoderm_only_MacHb", "Endoderm only MacHb resample 100", 100, 104)
Endoderm_only_MacHb
my_features <- only_MacHry
Endoderm_only_MacHry <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_MacHry", "Endoderm_only_MacHry", "Endoderm only MacHry resample 100", 100, 104)
Endoderm_only_MacHry 
my_features <- only_MacSna
Endoderm_only_MacSna <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_MacSna", "Endoderm_only_MacSna", "Endoderm only MacSna resample 100", 100, 104)
Endoderm_only_MacSna 
my_features <- only_MacTwi
Endoderm_only_MacTwi <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_MacTwi", "Endoderm_only_MacTwi", "Endoderm only MacTwi resample 100", 100, 104)
Endoderm_only_MacTwi 
my_features <- only_MacDl
Endoderm_only_MacDl <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_MacDl", "Endoderm_only_MacDl", "Endoderm only MacDl resample 100", 100, 104)
Endoderm_only_MacDl  
my_features <- only_White
Endoderm_only_White <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_White", "Endoderm_only_White", "Endoderm only White resample 100", 100, 104)
Endoderm_only_White 
my_features <- only_Zeit
Endoderm_only_Zeit <- resample_graph(enhancer_class_Endoderm,balanced_end, balanced_end_table, "Endoderm_only_Zeit", "Endoderm_only_Zeit", "Endoderm only Zeit resample 100", 100, 104)
Endoderm_only_Zeit

#Anterior
my_features <- all_features
A_all <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_all", "A_all","A all features resample 100", 100, 104)
A_all
my_features <- drop_Zeit_Twist
A_drop_ZeitTwi <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_ZeitTwi", "A_drop_ZeitTwi","A drop Zeitlinger Twi resample 100", 100, 104)
A_drop_ZeitTwi
my_features <- drop_White_H3K27ac
A_drop_White_H3K27ac <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_White_H3K27ac", "A_drop_White_H3K27ac","A drop White H3K27ac resample 100", 100, 104)
A_drop_White_H3K27ac
my_features <- drop_MacDorsal 
A_drop_MacDorsal <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_MacDorsal", "A_drop_MacDorsal","A drop MacArthur Dorsal resample 100", 100, 104)
A_drop_MacDorsal
my_features <- drop_MacTwist
A_drop_MacTwist <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_MacTwist", "A_drop_MacTwist","A drop MacArthur Twist resample 100", 100, 104)
A_drop_MacTwist
my_features <- drop_MacSnail
A_drop_MacSnail <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_MacSnail", "A_drop_MacSnail","A drop MacArthur Snail resample 100", 100, 104)
A_drop_MacSnail
my_features <- drop_MacHairy
A_drop_MacHairy <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_MacHairy", "A_drop_MacHairy","A drop MacArthur Hairy resample 100", 100, 104)
A_drop_MacHairy
my_features <- drop_MacHunchback
A_drop_MacHunchback <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_MacHunchback", "A_drop_MacHunchback","A drop MacArthur Hunchback resample 100", 100, 104)
A_drop_MacHunchback
my_features <- drop_MacGiant
A_drop_MacGiant <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_MacGiant", "A_drop_MacGiant","A drop MacArthur Giant resample 100", 100, 104)
A_drop_MacGiant
my_features <- drop_MacCaudal
A_drop_MacCaudal <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_MacCaudal", "A_drop_MacCaudal","A drop MacArthur Caudal resample 100", 100, 104)
A_drop_MacCaudal
my_features <- drop_MacBicoid 
A_drop_MacBicoid  <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_MacBicoid", "A_drop_MacBicoid","A drop MacArthur Bicoid resample 100", 100, 104)
A_drop_MacBicoid 
my_features <- drop_KK_H3K4me1 
A_drop_KK_H3K4me1   <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_KK_H3K4me1", "A_drop_KK_H3K4me1", "A drop KK H3K4me1 resample 100", 100, 104)
A_drop_KK_H3K4me1  
my_features <- drop_KK_H3K27ac
A_drop_KK_H3K27ac   <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_KK_H3K27ac", "A_drop_KK_H3K27ac", "A drop KK H3K27ac resample 100", 100, 104)
A_drop_KK_H3K27ac 
my_features <- drop_FurlongTwist
A_drop_FurlongTwist  <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_FurlongTwist", "A_drop_FurlongTwist", "A drop Furlong Twist resample 100", 100, 104)
A_drop_FurlongTwist 
my_features <- drop_FurlongSnail
A_drop_FurlongSnail  <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_FurlongSnail", "A_drop_FurlongSnail", "A drop Furlong Snail resample 100", 100, 104)
A_drop_FurlongSnail 
my_features <- drop_Eisen13Kr 
A_drop_Eisen13Kr <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_Eisen13Kr", "A_drop_Eisen13Kr", "A drop Eisen 2013 Kr resample 100", 100, 104)
A_drop_Eisen13Kr  
my_features <- drop_Eisen13Bcd
A_drop_Eisen13Bcd <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_Eisen13Bcd", "A_drop_Eisen13Bcd", "A drop Eisen 2013 Bcd resample 100", 100, 104)
A_drop_Eisen13Bcd
my_features <- drop_Eisen13Hb
A_drop_Eisen13Hb <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_Eisen13Hb", "A_drop_Eisen13Hb", "A drop Eisen 2013 Hb resample 100", 100, 104)
A_drop_Eisen13Hb
my_features <- drop_Eisen13Gt
A_drop_Eisen13Gt <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_Eisen13Gt", "A_drop_Eisen13Gt", "A drop Eisen 2013 Gt resample 100", 100, 104)
A_drop_Eisen13Gt
my_features <- drop_Eisen10Kni
A_drop_Eisen10Kni <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_Eisen10Kni", "A_drop_Eisen10Kni", "A drop Eisen 2010 Kni resample 100", 100, 104)
A_drop_Eisen10Kni
my_features <- drop_Eisen10Cad
A_drop_Eisen10Cad <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_Eisen10Cad", "A_drop_Eisen10Cad", "A drop Eisen 2010 Cad resample 100", 100, 104)
A_drop_Eisen10Cad
my_features <- drop_Zelda
A_drop_Zelda <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_drop_Zelda", "A_drop_Zelda", "A drop Zelda resample 100", 100, 104)
A_drop_Zelda
my_features <- drop_Dl 
my_features <- drop_Bcd 
my_features <- drop_Cad 
my_features <- only_Zelda
A_only_Zelda <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_Zelda", "A_only_Zelda", "A only Zelda resample 100", 100, 104)
A_only_Zelda 
my_features <- only_Eisen10Cad 
A_only_Eisen10Cad <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_Eisen10Cad", "A_only_Eisen10Cad", "A only Eisen10Cad resample 100", 100, 104)
A_only_Eisen10Cad 
my_features <- only_Eisen10Kni 
A_only_Eisen10Kni <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_Eisen10Kni", "A_only_Eisen10Kni", "A only Eisen10Kni resample 100", 100, 104)
A_only_Eisen10Kni
my_features <- only_Eisen13Gt
A_only_Eisen10Gt <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_Eisen10Gt", "A_only_Eisen10Gt", "A only Eisen10Gt resample 100", 100, 104)
A_only_Eisen10Gt
my_features <- only_Eisen13Hb
A_only_Eisen13Hb <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_Eisen13Hb", "A_only_Eisen13Hb", "A only Eisen13Hb resample 100", 100, 104)
A_only_Eisen13Hb
my_features <- only_Eisen13Bcd
A_only_Eisen13Bcd <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_Eisen13Bcd", "A_only_Eisen13Bcd", "A only Eisen13Bcd resample 100", 100, 104)
A_only_Eisen13Bcd 
my_features <- only_Eisen13Kr
A_only_Eisen13Kr <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_Eisen13Kr", "A_only_Eisen13Kr", "A only Eisen13Kr resample 100", 100, 104)
A_only_Eisen13Kr 
my_features <- only_FurlongSnail
A_only_FurlongSnail <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_FurlongSnail", "A_only_FurlongSnail", "A only FurlongSnail resample 100", 100, 104)
A_only_FurlongSnail 
my_features <- only_FurlongTwist
A_only_FurlongTwist <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_FurlongTwist", "A_only_FurlongTwist", "A only FurlongTwist resample 100", 100, 104)
A_only_FurlongTwist
my_features <- only_KK
A_only_KK <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_KK", "A_only_KK", "A only Kurtulus resample 100", 100, 104)
A_only_KK
my_features <- only_MacBicoid
A_only_MacBicoid <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_MacBicoid", "A_only_MacBicoid", "A only MacBicoid resample 100", 100, 104)
A_only_MacBicoid
my_features <- only_MacCad
A_only_MacCad <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_MacCad", "A_only_MacCad", "A only MacCad resample 100", 100, 104)
A_only_MacCad 
my_features <- only_MacGiant
A_only_MacGiant <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_MacGiant", "A_only_MacGiant", "A only MacGiant resample 100", 100, 104)
A_only_MacGiant  
my_features <- only_MacHb
A_only_MacHb <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_MacHb", "A_only_MacHb", "A only MacHb resample 100", 100, 104)
A_only_MacHb
my_features <- only_MacHry
A_only_MacHry <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_MacHry", "A_only_MacHry", "A only MacHry resample 100", 100, 104)
A_only_MacHry 
my_features <- only_MacSna
A_only_MacSna <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_MacSna", "A_only_MacSna", "A only MacSna resample 100", 100, 104)
A_only_MacSna 
my_features <- only_MacTwi
A_only_MacTwi <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_MacTwi", "A_only_MacTwi", "A only MacTwi resample 100", 100, 104)
A_only_MacTwi 
my_features <- only_MacDl
A_only_MacDl <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_MacDl", "A_only_MacDl", "A only MacDl resample 100", 100, 104)
A_only_MacDl  
my_features <- only_White
A_only_White <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_White", "A_only_White", "A only White resample 100", 100, 104)
A_only_White 
my_features <- only_Zeit
A_only_Zeit <- resample_graph(enhancer_class_A,balanced_A, balanced_A_table, "A_only_Zeit", "A_only_Zeit", "A only Zeit resample 100", 100, 104)
A_only_Zeit

#Center
my_features <- all_features
C_all <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_all", "C_all","C all features resample 100", 100, 104)
C_all
my_features <- drop_Zeit_Twist
C_drop_ZeitTwi <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_ZeitTwi", "C_drop_ZeitTwi","C drop Zeitlinger Twi resample 100", 100, 104)
C_drop_ZeitTwi
my_features <- drop_White_H3K27ac
C_drop_White_H3K27ac <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_White_H3K27ac", "C_drop_White_H3K27ac","C drop White H3K27ac resample 100", 100, 104)
C_drop_White_H3K27ac
my_features <- drop_MacDorsal 
C_drop_MacDorsal <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_MacDorsal", "C_drop_MacDorsal","C drop MacArthur Dorsal resample 100", 100, 104)
C_drop_MacDorsal
my_features <- drop_MacTwist
C_drop_MacTwist <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_MacTwist", "C_drop_MacTwist","C drop MacArthur Twist resample 100", 100, 104)
C_drop_MacTwist
my_features <- drop_MacSnail
C_drop_MacSnail <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_MacSnail", "C_drop_MacSnail","C drop MacArthur Snail resample 100", 100, 104)
C_drop_MacSnail
my_features <- drop_MacHairy
C_drop_MacHairy <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_MacHairy", "C_drop_MacHairy","C drop MacArthur Hairy resample 100", 100, 104)
C_drop_MacHairy
my_features <- drop_MacHunchback
C_drop_MacHunchback <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_MacHunchback", "C_drop_MacHunchback","C drop MacArthur Hunchback resample 100", 100, 104)
C_drop_MacHunchback
my_features <- drop_MacGiant
C_drop_MacGiant <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_MacGiant", "C_drop_MacGiant","C drop MacArthur Giant resample 100", 100, 104)
C_drop_MacGiant
my_features <- drop_MacCaudal
C_drop_MacCaudal <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_MacCaudal", "C_drop_MacCaudal","C drop MacArthur Caudal resample 100", 100, 104)
C_drop_MacCaudal
my_features <- drop_MacBicoid 
C_drop_MacBicoid  <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_MacBicoid", "C_drop_MacBicoid","C drop MacArthur Bicoid resample 100", 100, 104)
C_drop_MacBicoid 
my_features <- drop_KK_H3K4me1 
C_drop_KK_H3K4me1   <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_KK_H3K4me1", "C_drop_KK_H3K4me1", "C drop KK H3K4me1 resample 100", 100, 104)
C_drop_KK_H3K4me1  
my_features <- drop_KK_H3K27ac
C_drop_KK_H3K27ac   <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_KK_H3K27ac", "C_drop_KK_H3K27ac", "C drop KK H3K27ac resample 100", 100, 104)
C_drop_KK_H3K27ac 
my_features <- drop_FurlongTwist
C_drop_FurlongTwist  <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_FurlongTwist", "C_drop_FurlongTwist", "C drop Furlong Twist resample 100", 100, 104)
C_drop_FurlongTwist 
my_features <- drop_FurlongSnail
C_drop_FurlongSnail  <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_FurlongSnail", "C_drop_FurlongSnail", "C drop Furlong Snail resample 100", 100, 104)
C_drop_FurlongSnail 
my_features <- drop_Eisen13Kr 
C_drop_Eisen13Kr <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_Eisen13Kr", "C_drop_Eisen13Kr", "C drop Eisen 2013 Kr resample 100", 100, 104)
C_drop_Eisen13Kr  
my_features <- drop_Eisen13Bcd
C_drop_Eisen13Bcd <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_Eisen13Bcd", "C_drop_Eisen13Bcd", "C drop Eisen 2013 Bcd resample 100", 100, 104)
C_drop_Eisen13Bcd
my_features <- drop_Eisen13Hb
C_drop_Eisen13Hb <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_Eisen13Hb", "C_drop_Eisen13Hb", "C drop Eisen 2013 Hb resample 100", 100, 104)
C_drop_Eisen13Hb
my_features <- drop_Eisen13Gt
C_drop_Eisen13Gt <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_Eisen13Gt", "C_drop_Eisen13Gt", "C drop Eisen 2013 Gt resample 100", 100, 104)
C_drop_Eisen13Gt
my_features <- drop_Eisen10Kni
C_drop_Eisen10Kni <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_Eisen10Kni", "C_drop_Eisen10Kni", "C drop Eisen 2010 Kni resample 100", 100, 104)
C_drop_Eisen10Kni
my_features <- drop_Eisen10Cad
C_drop_Eisen10Cad <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_Eisen10Cad", "C_drop_Eisen10Cad", "C drop Eisen 2010 Cad resample 100", 100, 104)
C_drop_Eisen10Cad
my_features <- drop_Zelda
C_drop_Zelda <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_drop_Zelda", "C_drop_Zelda", "C drop Zelda resample 100", 100, 104)
C_drop_Zelda
my_features <- drop_Dl 
my_features <- drop_Bcd 
my_features <- drop_Cad 
my_features <- only_Zelda
C_only_Zelda <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_Zelda", "C_only_Zelda", "C only Zelda resample 100", 100, 104)
C_only_Zelda 
my_features <- only_Eisen10Cad 
C_only_Eisen10Cad <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_Eisen10Cad", "C_only_Eisen10Cad", "C only Eisen10Cad resample 100", 100, 104)
C_only_Eisen10Cad 
my_features <- only_Eisen10Kni 
C_only_Eisen10Kni <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_Eisen10Kni", "C_only_Eisen10Kni", "C only Eisen10Kni resample 100", 100, 104)
C_only_Eisen10Kni
my_features <- only_Eisen13Gt
C_only_Eisen10Gt <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_Eisen10Gt", "C_only_Eisen10Gt", "C only Eisen10Gt resample 100", 100, 104)
C_only_Eisen10Gt
my_features <- only_Eisen13Hb
C_only_Eisen13Hb <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_Eisen13Hb", "C_only_Eisen13Hb", "C only Eisen13Hb resample 100", 100, 104)
C_only_Eisen13Hb
my_features <- only_Eisen13Bcd
C_only_Eisen13Bcd <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_Eisen13Bcd", "C_only_Eisen13Bcd", "C only Eisen13Bcd resample 100", 100, 104)
C_only_Eisen13Bcd 
my_features <- only_Eisen13Kr
C_only_Eisen13Kr <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_Eisen13Kr", "C_only_Eisen13Kr", "C only Eisen13Kr resample 100", 100, 104)
C_only_Eisen13Kr 
my_features <- only_FurlongSnail
C_only_FurlongSnail <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_FurlongSnail", "C_only_FurlongSnail", "C only FurlongSnail resample 100", 100, 104)
C_only_FurlongSnail 
my_features <- only_FurlongTwist
C_only_FurlongTwist <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_FurlongTwist", "C_only_FurlongTwist", "C only FurlongTwist resample 100", 100, 104)
C_only_FurlongTwist
my_features <- only_KK
C_only_KK <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_KK", "C_only_KK", "C only Kurtulus resample 100", 100, 104)
C_only_KK
my_features <- only_MacBicoid
C_only_MacBicoid <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_MacBicoid", "C_only_MacBicoid", "C only MacBicoid resample 100", 100, 104)
C_only_MacBicoid
my_features <- only_MacCad
C_only_MacCad <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_MacCad", "C_only_MacCad", "C only MacCad resample 100", 100, 104)
C_only_MacCad 
my_features <- only_MacGiant
C_only_MacGiant <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_MacGiant", "C_only_MacGiant", "C only MacGiant resample 100", 100, 104)
C_only_MacGiant  
my_features <- only_MacHb
C_only_MacHb <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_MacHb", "C_only_MacHb", "C only MacHb resample 100", 100, 104)
C_only_MacHb
my_features <- only_MacHry
C_only_MacHry <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_MacHry", "C_only_MacHry", "C only MacHry resample 100", 100, 104)
C_only_MacHry 
my_features <- only_MacSna
C_only_MacSna <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_MacSna", "C_only_MacSna", "C only MacSna resample 100", 100, 104)
C_only_MacSna 
my_features <- only_MacTwi
C_only_MacTwi <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_MacTwi", "C_only_MacTwi", "C only MacTwi resample 100", 100, 104)
C_only_MacTwi 
my_features <- only_MacDl
C_only_MacDl <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_MacDl", "C_only_MacDl", "C only MacDl resample 100", 100, 104)
C_only_MacDl  
my_features <- only_White
C_only_White <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_White", "C_only_White", "C only White resample 100", 100, 104)
C_only_White 
my_features <- only_Zeit
C_only_Zeit <- resample_graph(enhancer_class_C,balanced_C, balanced_C_table, "C_only_Zeit", "C_only_Zeit", "C only Zeit resample 100", 100, 104)
C_only_Zeit

#Posterior
my_features <- all_features
P_all <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_all", "P_all","P all features resample 100", 100, 104)
P_all
my_features <- drop_Zeit_Twist
P_drop_ZeitTwi <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_ZeitTwi", "P_drop_ZeitTwi","P drop Zeitlinger Twi resample 100", 100, 104)
P_drop_ZeitTwi
my_features <- drop_White_H3K27ac
P_drop_White_H3K27ac <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_White_H3K27ac", "P_drop_White_H3K27ac","P drop White H3K27ac resample 100", 100, 104)
P_drop_White_H3K27ac
my_features <- drop_MacDorsal 
P_drop_MacDorsal <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_MacDorsal", "P_drop_MacDorsal","P drop MacArthur Dorsal resample 100", 100, 104)
P_drop_MacDorsal
my_features <- drop_MacTwist
P_drop_MacTwist <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_MacTwist", "P_drop_MacTwist","P drop MacArthur Twist resample 100", 100, 104)
P_drop_MacTwist
my_features <- drop_MacSnail
P_drop_MacSnail <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_MacSnail", "P_drop_MacSnail","P drop MacArthur Snail resample 100", 100, 104)
P_drop_MacSnail
my_features <- drop_MacHairy
P_drop_MacHairy <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_MacHairy", "P_drop_MacHairy","P drop MacArthur Hairy resample 100", 100, 104)
P_drop_MacHairy
my_features <- drop_MacHunchback
P_drop_MacHunchback <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_MacHunchback", "P_drop_MacHunchback","P drop MacArthur Hunchback resample 100", 100, 104)
P_drop_MacHunchback
my_features <- drop_MacGiant
P_drop_MacGiant <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_MacGiant", "P_drop_MacGiant","P drop MacArthur Giant resample 100", 100, 104)
P_drop_MacGiant
my_features <- drop_MacCaudal
P_drop_MacCaudal <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_MacCaudal", "P_drop_MacCaudal","P drop MacArthur Caudal resample 100", 100, 104)
P_drop_MacCaudal
my_features <- drop_MacBicoid 
P_drop_MacBicoid  <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_MacBicoid", "P_drop_MacBicoid","P drop MacArthur Bicoid resample 100", 100, 104)
P_drop_MacBicoid 
my_features <- drop_KK_H3K4me1 
P_drop_KK_H3K4me1   <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_KK_H3K4me1", "P_drop_KK_H3K4me1", "P drop KK H3K4me1 resample 100", 100, 104)
P_drop_KK_H3K4me1  
my_features <- drop_KK_H3K27ac
P_drop_KK_H3K27ac   <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_KK_H3K27ac", "P_drop_KK_H3K27ac", "P drop KK H3K27ac resample 100", 100, 104)
P_drop_KK_H3K27ac 
my_features <- drop_FurlongTwist
P_drop_FurlongTwist  <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_FurlongTwist", "P_drop_FurlongTwist", "P drop Furlong Twist resample 100", 100, 104)
P_drop_FurlongTwist 
my_features <- drop_FurlongSnail
P_drop_FurlongSnail  <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_FurlongSnail", "P_drop_FurlongSnail", "P drop Furlong Snail resample 100", 100, 104)
P_drop_FurlongSnail 
my_features <- drop_Eisen13Kr 
P_drop_Eisen13Kr <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_Eisen13Kr", "P_drop_Eisen13Kr", "P drop Eisen 2013 Kr resample 100", 100, 104)
P_drop_Eisen13Kr  
my_features <- drop_Eisen13Bcd
P_drop_Eisen13Bcd <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_Eisen13Bcd", "P_drop_Eisen13Bcd", "P drop Eisen 2013 Bcd resample 100", 100, 104)
P_drop_Eisen13Bcd
my_features <- drop_Eisen13Hb
P_drop_Eisen13Hb <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_Eisen13Hb", "P_drop_Eisen13Hb", "P drop Eisen 2013 Hb resample 100", 100, 104)
P_drop_Eisen13Hb
my_features <- drop_Eisen13Gt
P_drop_Eisen13Gt <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_Eisen13Gt", "P_drop_Eisen13Gt", "P drop Eisen 2013 Gt resample 100", 100, 104)
P_drop_Eisen13Gt
my_features <- drop_Eisen10Kni
P_drop_Eisen10Kni <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_Eisen10Kni", "P_drop_Eisen10Kni", "P drop Eisen 2010 Kni resample 100", 100, 104)
P_drop_Eisen10Kni
my_features <- drop_Eisen10Cad
P_drop_Eisen10Cad <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_Eisen10Cad", "P_drop_Eisen10Cad", "P drop Eisen 2010 Cad resample 100", 100, 104)
P_drop_Eisen10Cad
my_features <- drop_Zelda
P_drop_Zelda <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_drop_Zelda", "P_drop_Zelda", "P drop Zelda resample 100", 100, 104)
P_drop_Zelda

#Rushlow Dl + one other
my_features <- only_Zelda
P_only_Zelda <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_Zelda", "P_only_Zelda", "P only Zelda resample 100", 100, 104)
P_only_Zelda 
my_features <- only_Eisen10Cad 
P_only_Eisen10Cad <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_Eisen10Cad", "P_only_Eisen10Cad", "P only Eisen10Cad resample 100", 100, 104)
P_only_Eisen10Cad 
my_features <- only_Eisen10Kni 
P_only_Eisen10Kni <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_Eisen10Kni", "P_only_Eisen10Kni", "P only Eisen10Kni resample 100", 100, 104)
P_only_Eisen10Kni
my_features <- only_Eisen13Gt
P_only_Eisen10Gt <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_Eisen10Gt", "P_only_Eisen10Gt", "P only Eisen10Gt resample 100", 100, 104)
P_only_Eisen10Gt
my_features <- only_Eisen13Hb
P_only_Eisen13Hb <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_Eisen13Hb", "P_only_Eisen13Hb", "P only Eisen13Hb resample 100", 100, 104)
P_only_Eisen13Hb
my_features <- only_Eisen13Bcd
P_only_Eisen13Bcd <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_Eisen13Bcd", "P_only_Eisen13Bcd", "P only Eisen13Bcd resample 100", 100, 104)
P_only_Eisen13Bcd 
my_features <- only_Eisen13Kr
P_only_Eisen13Kr <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_Eisen13Kr", "P_only_Eisen13Kr", "P only Eisen13Kr resample 100", 100, 104)
P_only_Eisen13Kr 
my_features <- only_FurlongSnail
P_only_FurlongSnail <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_FurlongSnail", "P_only_FurlongSnail", "P only FurlongSnail resample 100", 100, 104)
P_only_FurlongSnail 
my_features <- only_FurlongTwist
P_only_FurlongTwist <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_FurlongTwist", "P_only_FurlongTwist", "P only FurlongTwist resample 100", 100, 104)
P_only_FurlongTwist
my_features <- only_KK
P_only_KK <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_KK", "P_only_KK", "P only Kurtulus resample 100", 100, 104)
P_only_KK
my_features <- only_MacBicoid
P_only_MacBicoid <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_MacBicoid", "P_only_MacBicoid", "P only MacBicoid resample 100", 100, 104)
P_only_MacBicoid
my_features <- only_MacCad
P_only_MacCad <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_MacCad", "P_only_MacCad", "P only MacCad resample 100", 100, 104)
P_only_MacCad 
my_features <- only_MacGiant
P_only_MacGiant <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_MacGiant", "P_only_MacGiant", "P only MacGiant resample 100", 100, 104)
P_only_MacGiant  
my_features <- only_MacHb
P_only_MacHb <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_MacHb", "P_only_MacHb", "P only MacHb resample 100", 100, 104)
P_only_MacHb
my_features <- only_MacHry
P_only_MacHry <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_MacHry", "P_only_MacHry", "P only MacHry resample 100", 100, 104)
P_only_MacHry 
my_features <- only_MacSna
P_only_MacSna <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_MacSna", "P_only_MacSna", "P only MacSna resample 100", 100, 104)
P_only_MacSna 
my_features <- only_MacTwi
P_only_MacTwi <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_MacTwi", "P_only_MacTwi", "P only MacTwi resample 100", 100, 104)
P_only_MacTwi 
my_features <- only_MacDl
P_only_MacDl <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_MacDl", "P_only_MacDl", "P only MacDl resample 100", 100, 104)
P_only_MacDl  
my_features <- only_White
P_only_White <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_White", "P_only_White", "P only White resample 100", 100, 104)
P_only_White 
my_features <- only_Zeit
P_only_Zeit <- resample_graph(enhancer_class_P,balanced_P, balanced_P_table, "P_only_Zeit", "P_only_Zeit", "P only Zeit resample 100", 100, 104)
P_only_Zeit

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
	#lines(1:nrow(datacombo), datacombo$rf500_tneg, type='b', col="lightgreen")
	lines(1:nrow(datacombo), datacombo$rf500_fpos, type='b', col="orange", lwd = 3)
	#lines(1:nrow(datacombo), datacombo$rf500_fneg, type='b', col="darkred")
	arrows(x0=1:nrow(datacombo), y0=(datacombo$rf500_tpos - 2*datacombo$sd_tpos), y1=(datacombo$rf500_tpos + 2*datacombo$sd_tpos), length=.05, angle=90, code=3)
	axis(1, at=1:nrow(datacombo), labels=rownames(datacombo),las=2)
}

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