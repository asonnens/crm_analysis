require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 
require(kernlab)
require(ROCR)
library(lattice)
library(plyr)

#maybe two functions
#is this an enhancer, if yes, what's its expression pattern

setwd("shadow_enhancers/machine_learning_input")

#crm_off <- read.table("on_vs_off4.tsv", header = TRUE) #training vs testing for on/off
crm_off_zelda <- read.table("on_vs_off5.tsv", header = TRUE)
crm_exp <- read.table("expression_300.tsv", header = TRUE)
test_test <- read.table("temp_test2.tsv", header = TRUE) #putative enhancers (based on Dl binding) around brk, vnd, zen, sog
crm_exp2 <- read.table("expression_301.tsv", header = TRUE)
crm_off <- read.csv("on_off_histone1.tsv", header = TRUE)

tab_on <- table(crm_off$status)
tab_on
length(tab_on)
tab_on2 <- table(crm_off_zelda$status)
tab_on2
length(tab_on)

tab_ectoderm <- table(crm_exp$Ectoderm)
tab_endoderm <- table(crm_exp$Endoderm)
tab_mesoderm <- table(crm_exp$Mesoderm)
tab_DV <- table(crm_exp2$DV)
tab_A <- table(crm_exp$A)
tab_C <- table(crm_exp$C)
tab_P <- table(crm_exp$P)
tab_ectoderm
tab_endoderm
tab_mesoderm
tab_A
tab_C
tab_P

tab_test <- table(test_test$status)
tab_test
length(tab_test)

strata_2var <- function(dataset, variable_name, variable_table, variable_pos, dataframe_length = 14, second_val = 7, samples = 300){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
      size = c(samples,samples))
  my_data <- dataset[, c(1,variable_pos,start_position:dataframe_length)]
  #my_data <- dataset[, c(1,variable_pos,9,10,5,6,7, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

strata_5var <- function(dataset, variable_name, variable_table, variable_pos, dataframe_length = 14, second_val = 7, samples = 300){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
      size = c(samples,samples, samples, samples, samples))
  my_data <- dataset[, c(1,variable_pos,start_position:dataframe_length)]
  #my_data <- dataset[, c(1,variable_pos,9,10,5,6,7, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

resample_default <- function(function_name, dataset, datatable, reps, variable_name, outfilename, description="README!", samplesize = 100){
    temp_dir <- paste("../machine_learning_output/", outfilename, sep = "")
    dir.create(temp_dir)
	info_file_name <- paste(temp_dir, "input_data_summary.txt", sep = "/")
	input_data_name <- paste(temp_dir, "input_data.csv", sep = "/")
	readme_file_name <- paste(temp_dir, "README.txt", sep = "/")
	write(description, file = readme_file_name)
	write.csv(crm_off, file = input_data_name)
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



resample_default2 <- function(function_name, dataset, datatable, reps, variable_name, outfilename){
    temp_dir <- paste("../machine_learning_output/", outfilename, sep = "")
    dir.create(temp_dir)
	datathing <- data.frame()
    for(i in 1:reps){
	    temp_file_name <- paste(outfilename, i, sep = "_")
		temp_file_name <- paste(temp_file_name, ".csv", sep = "")
		temp_file_name <- paste(temp_dir, temp_file_name, sep = "/")
        myresult <- suppressWarnings(function_name(dataset, datatable, temp_file_name))
			datathing[i,"svmrad_Ect"] <- myresult[[1]][1]
			datathing[i,"svmrad_EctEndMes"] <- myresult[[2]][1]
			datathing[i,"svmrad_EctMes"] <- myresult[[3]][1]
			datathing[i,"svmrad_End"] <- myresult[[4]][1]
			datathing[i,"svmrad_Mes"] <- myresult[[5]][1] 
		    datathing[i,"svmsig_Ect"] <- myresult[[1]][2]			
			datathing[i,"svmsig_EctEndMes"] <- myresult[[2]][2]
			datathing[i,"svmsig_EctMes"] <- myresult[[3]][2]
			datathing[i,"svmsig_End"] <- myresult[[4]][2]
			datathing[i,"svmsig_Mes"] <- myresult[[5]][2]    
		    datathing[i,"rf500_Ect"] <- myresult[[1]][3]
			datathing[i,"rf500_EctEndMes"] <- myresult[[2]][3]
			datathing[i,"rf500_EctMes"] <- myresult[[3]][3]
			datathing[i,"rf500_End"] <- myresult[[4]][3]
			datathing[i,"rf500_Mes"] <- myresult[[5]][3] 
			
     }      
	return(datathing)            
}                                         

enhancers_set <- strata_2var(crm_off, "status", tab_on, 2, 27, 25, 342)   
enhancers_training2 <- enhancers_set[[1]]
enhancers_data_scale <- scale(enhancers_training2[,c(4:28)])                       
enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_scale)
enhancers_test2 <- enhancers_set[[2]]
enhancers_data_testscale <- scale(enhancers_test2[,c(4:28)])
enhancers_testcale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)

enhancer_class_status_scale <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize){
    enhancers_set <- strata_2var(crm_off, "status", tab_on, 2, 27, 25, 342)   
    enhancers_training2 <- enhancers_set[[1]]
    enhancers_data_scale <- scale(enhancers_training2[,c(4:28)])                       
    enhancers_trainingscale <- cbind(enhancers_training2[,1:2],enhancers_data_scale)
    enhancers_test2 <- enhancers_set[[2]]
    enhancers_data_testscale <- scale(enhancers_test2[,c(4:28)])
    enhancers_testscale <- cbind(enhancers_test2[,1:2],enhancers_data_testscale)
	enhancers_train_original <- enhancers_training2[,c(1:2, 4,6,9,12,17,19,21,23,25)]
	enhancers_test_original <- enhancers_test2[,c(1:2, 4,6,9,12,17,19,21,23,25)]
	enhancers_train_scale <- enhancers_trainingscale[,c(1:2,3,5,8,11,16,18,20,22,24)]
	enhancers_test_scale <- enhancers_testscale[,c(1:2,3,5,8,11,16,18,20,22,24)]
	enhancers_train_25 <- enhancers_training2[,c(1:2, 4,7,10,13,18,20,22,24,26)]
	enhancers_test_25 <- enhancers_test2[,c(1:2, 4,7,10,13,18,20,22,24,26)]
	enhancers_train_original_Zelda180 <- enhancers_training2[,c(1:2, 5,6,9,12,17,19,21,23,25)]
	enhancers_test_original_Zelda180 <- enhancers_test2[,c(1:2, 5,6,9,12,17,19,21,23,25)]
	enhancers_train_histone <- enhancers_training2[,c(1:2, 4,6,9,12,17,19,21,23,25,27)]
	enhancers_test_histone <- enhancers_test2[,c(1:2, 4,6,9,12,17,19,21,23,25,27)]
	enhancers_train_Rushlow <- enhancers_training2[,c(1:2, 4,6,8,9,12,17,19,21,23,25)]
	enhancers_test_Rushlow <- enhancers_test2[,c(1:2, 4,6,8,9,12,17,19,21,23,25)]
	enhancers_training <- enhancers_train_histone
	enhancers_test <- enhancers_test_histone
	enhancers_weights <- table(enhancers_test$status) 
	actual <- enhancers_test$status
		#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(status ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$status, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_false_off <- 100 * svm_radial_test_table[2,1]/sum(svm_radial_test_table[,1])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(status ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$status, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_false_off <- 100 * svm_sigmoid_test_table[2,1]/sum(svm_sigmoid_test_table[,1])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(status ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$status, predicted = rf_500_predicted)
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
	name_list <- as.vector(enhancers_test$name)
	output_list <- cbind(name_list, actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(enhancers_test, file = summaryfile)
	write.csv(output_list, file = filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list, false_negatives_list)
	return(return_list)             
}       

                                    
enhancer_class_status <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize){
    enhancers_set <- strata_2var(dataset, "status", datatable, 2, 13, 10, samplesize)
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_data_scale <- scale(enhancers_training2[,c(4:11)])
	enhancers_training_scale <- cbind(enhancers_training2[,c(1:2)], enhancers_data_scale)
	enhancers_training <- enhancers_training2[,c(1:4,7,9,10)] #dl, sna, tw plus zelda
	#enhancers_training <- enhancers_training2[,c(1:2,6,8,9)] #dl, sna, tw
	#enhancers_training <- enhancers_training2[,c(1:3,6:length(enhancers_training2))] #plus zelda
	enhancers_test2 <- enhancers_set[[2]]
	enhancers_test <- enhancers_test2[,c(1:3,6,8,9)] #dl sna tw zelda
	#enhancers_test <- enhancers_test2[,c(1:2,6,8,9)] #dl sna tw
	#enhancers_test <- enhancers_test2[,c(1:3,6:length(enhancers_test2))] #plus zelda
	enhancers_weights <- table(enhancers_test$status)
	actual <- enhancers_test$status
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(status ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$status, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_false_off <- 100 * svm_radial_test_table[2,1]/sum(svm_radial_test_table[,1])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(status ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$status, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_false_off <- 100 * svm_sigmoid_test_table[2,1]/sum(svm_sigmoid_test_table[,1])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(status ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$status, predicted = rf_500_predicted)
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
	name_list <- as.vector(enhancers_test$name)
	output_list <- cbind(name_list, actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(enhancers_test, file = summaryfile)
	write.csv(output_list, file = filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list, false_negatives_list)
	return(return_list)             
}       

enhancer_class_Ectoderm <- function(dataset, datatable, filename = "temp", summaryfile = "temp_summary", samplesize){
    enhancers_set <- strata_2var(dataset, "Ectoderm", datatable, 8, 21, 10, samplesize) 
	enhancers_training2 <- enhancers_set[[1]]
	enhancers_training <- enhancers_training2[,c(1:3,6,8:9)] #dl sna tw zld
	#enhancers_training <- enhancers_training2[,c(1:2,6,8:9)] #dl sna tw
	#enhancers_training <- enhancers_training2[,c(1:2,6:length(enhancers_training2))]  #8 + zelda
	enhancers_test2 <- enhancers_set[[2]]
	#enhancers_test <- enhancers_test2
	#enhancers_test <- enhancers_test2[,c(1:3, 6:length(enhancers_test2))]  #8 + zelda
	#enhancers_test <- enhancers_test2[,c(1:2,6,8:9)] #dl sna tw
	enhancers_test <- enhancers_test2[,c(1:3,6,8:9)] #dl sna  zld
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

enhancer_class_Endoderm <- function(dataset, datatable, filename = "temp"){
    enhancers_set <- strata_2var(dataset, "Endoderm", datatable, 9, 21, 10, 60) 
	enhancers_training <- enhancers_set[[1]]
	enhancers_test <- enhancers_set[[2]]
	enhancers_weights <- table(enhancers_test$Endoderm)
	actual <- enhancers_test$Endoderm
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(Endoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$Endoderm, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(Endoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$Endoderm, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(Endoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$Endoderm, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	prediction_list <- c(svm_radial_prediction_success, svm_sigmoid_prediction_success, rf_500_prediction_success)
	true_positives_list <- c(svm_rad_actual_on, svm_sig_actual_on, rf_500_actual_on)
	true_negatives_list <- c(svm_rad_true_off, svm_sig_true_off, rf_500_true_off)
	false_positives_list <- c(svm_rad_false_on, svm_sig_false_on, rf_500_false_on)
	output_list <- cbind(actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(output_list, filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list)
	return(return_list)             
}     

enhancer_class_Mesoderm <- function(dataset, datatable, filename = "temp"){
    enhancers_set <- strata_2var(dataset, "Mesoderm", datatable, 10, 21, 10, 60) 
	enhancers_training <- enhancers_set[[1]]
	enhancers_test <- enhancers_set[[2]]
	enhancers_weights <- table(enhancers_test$Mesoderm)
	actual <- enhancers_test$Mesoderm
	Endoderm_1 <- resample_default(enhancer_class_Endoderm, crm_exp, tab_endoderm, 100, "End", "Endoderm_model_comparison1")        
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(Mesoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$Mesoderm, predicted = radial_svm_predicted)
	svm_rad_actual_on <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_false_on <- 100 * svm_radial_test_table[1,2]/sum(svm_radial_test_table[,2])
	svm_rad_true_off <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(Mesoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$Mesoderm, predicted = sigmoid_svm_predicted)
	svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
	svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(Mesoderm ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$Mesoderm, predicted = rf_500_predicted)
	rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
	rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
	rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	prediction_list <- c(svm_radial_prediction_success, svm_sigmoid_prediction_success, rf_500_prediction_success)
	true_positives_list <- c(svm_rad_actual_on, svm_sig_actual_on, rf_500_actual_on)
	true_negatives_list <- c(svm_rad_true_off, svm_sig_true_off, rf_500_true_off)
	false_positives_list <- c(svm_rad_false_on, svm_sig_false_on, rf_500_false_on)
	output_list <- cbind(actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(output_list, filename)
	#close(filename)
	return_list <- list(prediction_list, true_positives_list, true_negatives_list, false_positives_list)
	return(return_list)             
}                    

enhancer_class_DV <- function(dataset, datatable, filename = "temp"){
    enhancers_set <- strata_5var(dataset, "DV", datatable, 4, 21, 10, 25) 
	enhancers_training <- enhancers_set[[1]]
	enhancers_test <- enhancers_set[[2]]
	enhancers_weights <- table(enhancers_test$DV)
	actual <- enhancers_test$DV
	#radial SVM                     
	enhancer_svm_radial <- suppressWarnings(svm(DV ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "radial"))
	radial_svm_predicted <- predict(enhancer_svm_radial, newdata=enhancers_test, type="class")
	svm_radial_test_table <- table(actual = enhancers_test$DV, predicted = radial_svm_predicted)
	svm_rad_actual_Ectoderm <- 100 * svm_radial_test_table[1,1]/sum(svm_radial_test_table[1,])
	svm_rad_actual_EctEndMes <- 100 * svm_radial_test_table[2,2]/sum(svm_radial_test_table[2,])
	svm_rad_actual_EctMes <- 100 * svm_radial_test_table[3,3]/sum(svm_radial_test_table[3,])
	svm_rad_actual_End <- 100 * svm_radial_test_table[4,4]/sum(svm_radial_test_table[4,])
	svm_rad_actual_Mes <- 100 * svm_radial_test_table[5,5]/sum(svm_radial_test_table[5,])
	svm_radial_prediction_success <- 100 * sum(diag(svm_radial_test_table)/sum(svm_radial_test_table))
	#sigmoidal SVM                  
	enhancer_svm_sigmoid <- svm(DV ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, kernel = "sigmoid")
	sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=enhancers_test, type="class")
	svm_sigmoid_test_table <- table(actual = enhancers_test$DV, predicted = sigmoid_svm_predicted)
	svm_sig_actual_Ectoderm <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
	svm_sig_actual_EctEndMes <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
	svm_sig_actual_EctMes <- 100 * svm_sigmoid_test_table[3,3]/sum(svm_sigmoid_test_table[3,])
	svm_sig_actual_End <- 100 * svm_sigmoid_test_table[4,4]/sum(svm_sigmoid_test_table[4,])
	svm_sig_actual_Mes <- 100 * svm_sigmoid_test_table[5,5]/sum(svm_sigmoid_test_table[5,])
	svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))
	#random forest 500 trees        
	enhancer_rf_500 <- randomForest(DV ~., data = enhancers_training[,2:length(enhancers_training)], class.weights = enhancers_weights, ntree = 500)
	rf_500_predicted <- predict(enhancer_rf_500, newdata=enhancers_test, type="class")
	rf_500_test_table <- table(actual = enhancers_test$DV, predicted = rf_500_predicted)
	rf_500_actual_Ectoderm <- 100 * rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
	rf_500_actual_EctEndMes <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
	rf_500_actual_EctMes <- 100 * rf_500_test_table[3,3]/sum(rf_500_test_table[3,])
	rf_500_actual_End <- 100 * rf_500_test_table[4,4]/sum(rf_500_test_table[4,])
	rf_500_actual_Mes <- 100 * rf_500_test_table[5,5]/sum(rf_500_test_table[5,])
	Ectoderm_list <- c(svm_rad_actual_Ectoderm, svm_sig_actual_Ectoderm, rf_500_actual_Ectoderm)
	EctEndMes_list <- c(svm_rad_actual_EctEndMes, svm_sig_actual_EctEndMes, rf_500_actual_EctEndMes)
	EctMes_list <- c(svm_rad_actual_EctMes, svm_sig_actual_EctMes, rf_500_actual_EctMes)
	End_list <- c(svm_rad_actual_End, svm_sig_actual_End, rf_500_actual_End)
    Mes_list <- c(svm_rad_actual_Mes, svm_sig_actual_Mes, rf_500_actual_Mes)
	output_list <- cbind(actual, radial_svm_predicted, sigmoid_svm_predicted, rf_500_predicted)
	write.csv(output_list, filename)
	close(filename)
	return_list <- list(Ectoderm_list, EctEndMes_list, EctMes_list, End_list, Mes_list)
	return(return_list)             
}                                     
                                    
cats <- resample_default(enhancer_class_status_scale, crm_off, tab_on, 100, "on", "on_off_model_histonesweird", "This is a temporary file", 342)
dogmeans <- apply(cats, MARGIN=2, FUN=mean)
barplot(dogs)
hamsters <- apply(cats, MARGIN=2, FUN=sd)
stanD <- hamsters
par(las=2)
par(las=1)
par(mgp = c(0, 1, 0))
#par(mar = c(1, 5 ,8 ,2.1))
#barplot(dogmeans, ylim=c(-6,100), col=c(33,"darkgreen","lightgreen","orange"), , xaxt="n")

par(mfrow=c(1,1))
barplot(dogmeans, ylim=c(-8,100), col=c(33,"darkgreen","lightgreen","orange", "darkred"), xaxt="n", main = "Model comparison- On vs. Off weird histones, 100 reps")
#barplot(dogmeans, ylim=c(-6,100), col=c(33,33,33,33,"cornsilk4","cornsilk4","cornsilk4","cornsilk4",3,3,3,3))

arrows(seq(.7, 18, 1.2), (dogmeans - stanD), seq(.7, 18, 1.2), (dogmeans + stanD), length=.05, angle=90, code=3)
text(y=-4, x=3, "SVM Radial")
text(y=-4, x=9, "SVM Sigmoid")
text(y=-4, x=15, "RandomForest500")
par(xpd=TRUE)
legend(9,-7, c("overall accuracy","true positive", "true negative", "false positive", "false negative"), text.col = c(33, "darkgreen", "lightgreen", "orange", "darkred"))  # puts text in the legend 

on_off_zelda <- resample_default(enhancer_class_status, crm_off_zelda, tab_on2, 100, "on", "on_off_model_comparison_Dl_sna_tw_Zelda1", "This is a temporary file", 384)
zelda_off_means <- apply(on_off_zelda, MARGIN=2, FUN=mean)
barplot(zelda_off_means)
zelda_off_sd <- apply(on_off_zelda, MARGIN=2, FUN=sd)
stanD_zelda <- zelda_off_sd
par(las=2)
par(las=1)
par(mgp = c(0, 1, 0))
#par(mar = c(1, 5 ,8 ,2.1))
#barplot(dogmeans, ylim=c(-6,100), col=c(33,"darkgreen","lightgreen","orange"), , xaxt="n")

par(mfrow=c(1,1))
barplot(zelda_off_means, ylim=c(-8,100), col=c(33,"darkgreen","lightgreen","orange", "darkred"), xaxt="n", main = "Model comparison- On vs. Off, 100 reps, dl, sna, tw,+ zld")
#barplot(dogmeans, ylim=c(-6,100), col=c(33,33,33,33,"cornsilk4","cornsilk4","cornsilk4","cornsilk4",3,3,3,3))

arrows(seq(.7, 18, 1.2), (zelda_off_means - stanD_zelda), seq(.7, 18, 1.2), (zelda_off_means + stanD_zelda), length=.05, angle=90, code=3)
text(y=-4, x=3, "SVM Radial")
text(y=-4, x=9, "SVM Sigmoid")
text(y=-4, x=15, "RandomForest500")
par(xpd=TRUE)


lty=c(1,1), # gives the legend appropriate symbols (lines)

lwd=c(2.5,2.5),col=c("blue","red")) # gives the legend lines the correct color and width



Ect_1 <- resample_default(enhancer_class_Ectoderm, crm_exp, tab_ectoderm, 100, "Ect", "Ect_model_comparison2_dl_sna_tw_zld", "this is a temporary file", 60)        
Ect_means <- apply(Ect_1, MARGIN=2, FUN=mean)
barplot(Ect_means)
Ect_sd <- apply(Ect_1, MARGIN=2, FUN=sd)
par(mgp = c(0, 1, 0))
barplot(Ect_means, ylim=c(-6,100), col=c(33,"darkgreen","lightgreen","orange", "darkred"), , xaxt="n", main = "Model comparison- Ectoderm dl, sna, tw, zld, 100 reps") 
arrows(seq(.7, 18, 1.2), (Ect_means - Ect_sd), seq(.7, 18, 1.2), (Ect_means + Ect_sd), length=.05, angle=90, code=3)
text(y=-4, x=3, "SVM Radial")
text(y=-4, x=9, "SVM Sigmoid")
text(y=-4, x=15, "RandomForest500")                 


Mesoderm_1 <- resample_default(enhancer_class_Mesoderm, crm_exp, tab_mesoderm, 100, "Mes", "Mesoderm_model_comparison1")        
Mes_means <- apply(Mesoderm_1, MARGIN=2, FUN=mean)
barplot(Mes_means)
Mes_sd <- apply(Mesoderm_1, MARGIN=2, FUN=sd)
par(mgp = c(0, 1, 0))
barplot(Mes_means, ylim=c(-6,100), col=c(33,"darkgreen","lightgreen","darkred"), , xaxt="n", main = "Model comparison- Mesoderm, 100 reps") 
arrows(seq(.7, 18, 1.2), (Mes_means - Mes_sd), seq(.7, 18, 1.2), (Mes_means + Mes_sd), length=.05, angle=90, code=3)
text(y=-4, x=2, "SVM Radial")
text(y=-4, x=7, "SVM Sigmoid")
text(y=-4, x=12, "RandomForest500")    

Endoderm_1 <- resample_default(enhancer_class_Endoderm, crm_exp, tab_mesoderm, 100, "End", "Endoderm_model_comparison1")        
End_means <- apply(Endoderm_1, MARGIN=2, FUN=mean)
barplot(End_means)
End_sd <- apply(Endoderm_1, MARGIN=2, FUN=sd)
par(mgp = c(0, 1, 0))
barplot(End_means, ylim=c(-6,100), col=c(33,"darkgreen","lightgreen","darkred"), , xaxt="n", main = "Model comparison- Endoderm, 100 reps") 
arrows(seq(.7, 18, 1.2), (End_means - End_sd), seq(.7, 18, 1.2), (End_means + End_sd), length=.05, angle=90, code=3)
text(y=-4, x=2, "SVM Radial")
text(y=-4, x=7, "SVM Sigmoid")
text(y=-4, x=12, "RandomForest500")

DV_1 <- resample_default2(enhancer_class_DV, crm_exp2, tab_DV, 100, "DV", "DV_model_comparison1")        
DV_means <- apply(DV_1, MARGIN=2, FUN=mean)
DV_means_plot <- c(DV_means[1], DV_means[6], DV_means[11], DV_means[2], DV_means[7], DV_means[12], DV_means[3],
 DV_means[8], DV_means[13], DV_means[4], DV_means[9], DV_means[14], DV_means[5], DV_means[10], DV_means[15])

DV_sd <- apply(DV_1, MARGIN=2, FUN=sd)
DV_sd_plot <- c(DV_sd[1], DV_sd[6], DV_sd[11], DV_sd[2], DV_sd[7], DV_sd[12], DV_sd[3],
 DV_sd[8], DV_sd[13], DV_sd[4], DV_sd[9], DV_sd[14], DV_sd[5], DV_sd[10], DV_sd[15])

par(mgp = c(0, 1, 0))D
barplot(DV_means_plot, ylim=c(-6,100), col=c(1,2,3), , xaxt="n", main = "Model comparison- DV, 100 reps") 
arrows(seq(.7, 18, 1.2), (DV_means_plot - DV_sd_plot), seq(.7, 18, 1.2), (DV_means_plot + DV_sd_plot), length=.05, angle=90, code=3)
text(y=-4, x=2, "Ect")
text(y=-4, x=5.5, "EctEndMes")
text(y=-4, x=9.3, "EctMes")
text(y=-4, x=12.8, "End")   
text(y=-4, x=16, "Mes") 


crm_on_train_test <- strata_2var(crm_off, "status", tab_on, 2, 10, 7, 300)
on_training <- crm_on_train_test[[1]]
on_test <- crm_on_train_test[[2]]
active_training <- on_training[,2:length(on_training)]
active_test <- on_test[,2:length(on_test)]
active_weights <- table(active_test$status)
#active_rf_50 <- randomForest(status ~., data = active_training, class.weights = active_weights, ntree = 50)
active_rf_500 <- randomForest(status ~., data = active_training, class.weights = active_weights, ntree = 500)
active_svm_lin <- svm(status ~., data = active_training, class.weights = active_weights, kernel = "linear")
active_svm_rad <- svm(status ~., data = active_training, class.weights = active_weights, kernel = "radial")
active_svm_sig <- svm(status ~., data = active_training, class.weights = active_weights, kernel = "sigmoid")
active_rf_50_results <- predict(active_rf_50, newdata = active_test, type = "class")
active_rf_500_results <- predict(active_rf_500, newdata = active_test, type = "class")
active_svm_linear_results <- predict(active_svm_lin, newdata=active_test, type="class")
active_svm_radial_results <- predict(active_svm_rad, newdata=active_test, type="class")
active_svm_sigmoid_results <- predict(active_svm_sig, newdata=active_test, type="class")
rf500_test_table <- table(actual = actual_status, predicted = active_rf_500_results)
sum(diag(rf500_test_table)/sum(rf500_test_table))
svmlin_test_table <- table(actual = actual_status, predicted = active_svm_linear_results)
sum(diag(svmlin_test_table)/sum(svmlin_test_table))
svmrad_test_table <- table(actual = actual_status, predicted = active_svm_radial_results)
sum(diag(svmrad_test_table)/sum(svmrad_test_table))
svmsig_test_table <- table(actual = actual_status, predicted = active_svm_sigmoid_results)
sum(diag(svmsig_test_table)/sum(svmsig_test_table))


actual_status <- active_test$status
on_frame <- cbind(active_test[,1:2], active_rf_50_results, active_rf_500_results, active_svm_linear_results, active_svm_radial_results, active_svm_sigmoid_results)

library(ggplot2); library(reshape2)
dat3 <- melt(on_frame, id.var = 'status')
ggplot(dat3, aes(variable, status)) + geom_tile(aes(fill = value),
   colour = "white") + scale_fill_manual(values=c("red", "blue", "black"))

rf_50 = c()
rf_500 = c()
svm_lin = c()	
svm_rad = c()
svm_sig = c()

for(i in (1:length(actual_status))){if ((actual_status[i] == "on") && (active_rf_50_results[i] == "on")){rf_500[i] <- 2} 
                                      else if ((actual_status[i] == "on") && (active_rf_50_results[i] == "off")) {rf_50[i] <- -2}
									  else if ((actual_status[i] == "off") && (active_rf_50_results[i] == "on")) {rf_50[i] <- -1}
									  else if ((actual_status[i] == "off") && (active_rf_50_results[i] == "off")){rf_50[i] <- 1}
									  else {print(active_rf_50_results[i])}}
for(i in (1:length(actual_status))){if ((actual_status[i] == "on") && (active_rf_500_results[i] == "on")){rf_500[i] <- 2} 
                                      else if ((actual_status[i] == "on") && (active_rf_500_results[i] == "off")) {rf_500[i] <- -2}
									  else if ((actual_status[i] == "off") && (active_rf_500_results[i] == "on")) {rf_500[i] <- -1}
									  else {rf_500[i] <- 1}}		   
for(i in (1:length(actual_status))){if ((actual_status[i] == "on") && (active_svm_linear_results[i] == "on")){svm_lin[i] <- 2} 
                                      else if ((actual_status[i] == "on") && (active_svm_linear_results[i] == "off")) {svm_lin[i] <- -2}
									  else if ((actual_status[i] == "off") && (active_svm_linear_results[i] == "on")) {svm_lin[i] <- -1}
									  else {svm_lin[i] <- 1}}
for(i in (1:length(actual_status))){if ((actual_status[i] == "on") && (active_svm_radial_results[i] == "on")){svm_rad[i] <- 2} 
                                      else if ((actual_status[i] == "on") && (active_svm_radial_results[i] == "off")) {svm_rad[i] <- -2}
									  else if ((actual_status[i] == "off") && (active_svm_radial_results[i] == "on")) {svm_rad[i] <- -1}
									  else {svm_rad[i] <- 1}}
for(i in (1:length(actual_status))){if ((actual_status[i] == "on") && (active_svm_sigmoid_results[i] == "on")){svm_sig[i] <- 2} 
                                      else if ((actual_status[i] == "on") && (active_svm_sigmoid_results[i] == "off")) {svm_sig[i] <- -2}
									  else if ((actual_status[i] == "off") && (active_svm_sigmoid_results[i] == "on")) {svm_sig[i] <- -1}
									  else {svm_sig[i] <- 1}}

on_frame <- cbind(active_test[,1], rf_500, svm_lin, svm_rad, svm_sig)
temp_frame <- on_frame[,2:5]
heatmap.2(as.matrix(temp_frame),symm=FALSE,trace="none", Colv = FALSE, Rowv = FALSE, margins=c(5,1) , xaxt="n", main = "Activity plot", cexCol=1, col = c("darkred", "orange", "white", "lightgreen", "darkgreen"))

temp_frame <- new_frame[,2:5]


on_set <- strata_2var(crm_off, "status", tab_on, 2, 10, 7, 373)
on_train <- on_set[[1]]
on_test <- on_set[[2]]
on_weights <- table(on_test$status)
actual_on <- on_test$status
enhancer_rf_500 <- randomForest(status ~., data = on_train[,2:length(on_train)], class.weights = on_weights, ntree = 500)
rf_500_predicted <- predict(enhancer_rf_500, newdata=on_test, type="class")
rf_500_test_table <- table(actual = on_test$status, predicted = rf_500_predicted)
rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
rf_500_prediction_success
rf_500_actual_on
rf_500_false_on

ect_set <- strata_2var(crm_exp, "Ectoderm", tab_ectoderm, 8, 21, 7, 60) 
ect_training <- ect_set[[1]]
ect_test <- ect_set[[2]]
ect_weights <- table(ect_test$Ectoderm)
actual_ect <- ect_test$Ectoderm
enhancer_svm_sigmoid <- svm(Ectoderm ~., data = ect_training[,2:length(ect_training)], class.weights = ect_weights, kernel = "sigmoid")
sigmoid_svm_predicted <- predict(enhancer_svm_sigmoid, newdata=ect_test, type="class")
svm_sigmoid_test_table <- table(actual = ect_test$Ectoderm, predicted = sigmoid_svm_predicted)
svm_sig_actual_on <- 100 * svm_sigmoid_test_table[2,2]/sum(svm_sigmoid_test_table[2,])
svm_sig_false_on <- 100 * svm_sigmoid_test_table[1,2]/sum(svm_sigmoid_test_table[,2])
svm_sig_true_off <- 100 * svm_sigmoid_test_table[1,1]/sum(svm_sigmoid_test_table[1,])
svm_sigmoid_prediction_success <- 100 * sum(diag(svm_sigmoid_test_table)/sum(svm_sigmoid_test_table))

on_train <- on_set[[1]]
on_test <- on_set[[2]]
on_weights <- table(on_test$status)
actual_on <- on_test$status
enhancer_rf_500 <- randomForest(status ~., data = on_train[,2:length(on_train)], class.weights = on_weights, ntree = 500)
rf_500_predicted <- predict(enhancer_rf_500, newdata=on_test, type="class")
rf_500_test_table <- table(actual = on_test$status, predicted = rf_500_predicted)
rf_500_prediction_success <- 100 * sum(diag(rf_500_test_table)/sum(rf_500_test_table))
rf_500_actual_on <- 100 * rf_500_test_table[2,2]/sum(rf_500_test_table[2,])
rf_500_false_on <- 100 * rf_500_test_table[1,2]/sum(rf_500_test_table[,2])
rf_500_true_off <- 100 *rf_500_test_table[1,1]/sum(rf_500_test_table[1,])
rf_500_prediction_success
rf_500_actual_on
rf_500_false_on

real_data <- read.table("predict_16.tsv", header = TRUE)
temp <- predict(enhancer_rf_500, newdata = real_data, type = "class")
real_on_frame <- cbind(real_data[,1], temp, real_data[,2:length(real_data)])
off_frame <- subset(real_on_frame, temp)

predicted_on <- subset(real_on_frame, temp == "on")

svm_sig_predicted_ext <- predict(enhancer_svm_sigmoid, newdata = predicted_on, type = "class")
ect_on_frame <- cbind(predicted_on[,1:2], svm_sig_predicted_ext)
#16 genes
#88 putative enhancers
#55 predicted on
#51 predicted to be ectoderm