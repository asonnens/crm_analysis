require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 
require(kernlab)
require(ROCR)

#maybe two functions
#is this an enhancer, if yes, what's its expression pattern

setwd("shadow_enhancers/machine_learning_input")

crm_off <- read.table("on_vs_off4.tsv", header = TRUE) 
crm_exp <- read.table("expression_300.tsv", header = TRUE)
test_test <- read.table("temp_test2.tsv", header = TRUE) #putative enhancers (based on Dl binding) around brk, vnd, zen, sog

tab_on <- table(crm_off$status)
tab_on
length(tab_on)

tab_ectoderm <- table(crm_exp$Ectoderm)
tab_endoderm <- table(crm_exp$Endoderm)
tab_mesoderm <- table(crm_exp$Mesoderm)
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

svm_enhancers_linear <- function(dataset, datatable, variable_pos, variable_name){
    svm_enhancers <- strata_5var(dataset, variable_name, datatable, 2/3, variable_pos, 10, 7, 150)
	#svm_enhancers <- strata_5var(dataset, variable_name, datatable, 2/3, variable_pos, 21, 10, 15)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_weights <- table(svm_test$status)
	svm_model <- svm(status ~., data = svm_training, class.weights = svm_weights, kernel = "linear")
	svm_training_table <- table(actual = svm_training$status,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$status,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

strata_5var <- function(dataset, variable_name, variable_table, percentage, variable_pos, dataframe_length = 14, second_val = 7, training_num = 30){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
      size = c(300,300))
      #size = c(round(variable_table[1] *(percentage)), round(variable_table[2] * (percentage))))
  my_data <- dataset[, c(1,variable_pos,9,10, start_position:dataframe_length)]
  #my_data <- dataset[, c(1,variable_pos,9,10,5,6,7, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

resample_default <- function(function_name, dataset, variable_name, datatable, reps, variable_pos){
	out_list <- list()
    for(i in 1:reps){
        myresult <- function_name(dataset, datatable, variable_pos, variable_name)
		out_list <- unlist(c(out_list, myresult))
     }
	mean_score <- mean(out_list)
	error_score <- sd(out_list)
	return_list <- list(mean_score,error_score,out_list) 
	return(return_list)
}

svm_enhancer_comparison_linear <- resample_default(svm_enhancers_linear, crm_off, "status", tab_on, 100, 2)
svm_enhancer_comparison_linear 