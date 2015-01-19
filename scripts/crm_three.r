require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 

#maybe two functions
#is this an enhancer, if yes, what's its expression pattern

setwd("shadow_enhancers")
crm_off <- read.table("on_vs_off3.tsv", header = TRUE) #training vs testing for on/off
test_test <- read.table("temp_test2.tsv", header = TRUE) #putative enhancers (based on Dl binding) around brk, vnd, zen, sog

tab_expr <- table(crm_one$AP)
tab_expr
length(tab_expr)

tab_on <- table(crm_off$status)
tab_on
length(tab_on)

tab_test <- table(test_test$status)
tab_test
length(tab_test)

strata_5var <- function(dataset, variable_name, variable_table, percentage, variable_pos, dataframe_length = 14, second_val = 7, training_num = 30){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
      size = c(training_num,training_num,training_num))
      size = c(round(variable_table[1] *(percentage)), 
	  round(variable_table[2] * (percentage)),
	  round(variable_table[3] * (percentage)))
	#  round(variable_table[4] * (percentage)))
  my_data <- dataset[, c(variable_pos, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

resample_default <- function(function_name, dataset, datatable, reps){
	out_list <- list()
    for(i in 1:reps){
        myresult <- function_name(dataset, datatable)
		out_list <- unlist(c(out_list, myresult))
     }
	mean_score <- mean(out_list)
	error_score <- sd(out_list)
	return_list <- list(mean_score,error_score,out_list) 
	return(return_list)
}

svm_enhancers_fun_linear <- function(dataset, datatable, kernel_type = "linear"){
    svm_enhancers <- strata_5var(dataset, "status", datatable, 2/3, 2, 10, 7, 150)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(status ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$status,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$status,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_enhancers_fun_radial <- function(dataset, datatable, kernel_type = "radial"){
    svm_enhancers <- strata_5var(dataset, "status", datatable, 2/3, 2, 10, 7, 150)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(status ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$status,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$status,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}
svm_enhancers_fun_polynomial <- function(dataset, datatable, kernel_type = "polynomial"){
    svm_enhancers <- strata_5var(dataset, "status", datatable, 2/3, 2, 10, 7, 150)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(status ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$status,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$status,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_enhancers_fun_sigmoid <- function(dataset, datatable, kernel_type = "sigmoid"){
    svm_enhancers <- strata_5var(dataset, "status", datatable, 2/3, 2, 10, 7, 150)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(status ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$status,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$status,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}


svm_enhancer_comparison_linear <- resample_default(svm_enhancers_fun_linear, crm_off, tab_on, 100)
svm_enhancer_comparison_linear
svm_enhancer_comparison_radial <- resample_default(svm_enhancers_fun_radial, crm_off, tab_on, 100)
svm_enhancer_comparison_radial
svm_enhancer_comparison_polynomial <- resample_default(svm_enhancers_fun_radial, crm_off, tab_on, 100)
svm_enhancer_comparison_polynomial
svm_enhancer_comparison_sigmoid <- resample_default(svm_enhancers_fun_radial, crm_off, tab_on, 100)
svm_enhancer_comparison_sigmoid

rf_enhancer_fun_100 <- function(dataset, datatable, trees = 100){
    rf_enhancer <- strata_5var(dataset, "status", datatable, 2/3, 2, 10, 7, 150)
	rf_training <- rf_enhancer[[1]]
	rf_test <- rf_enhancer[[2]]
	random_forest_model <- randomForest(status ~., data = rf_training, ntree = trees)
	rf_training_table <- table(actual = rf_training$status,
                       predicted = predict(random_forest_model, type="class"))
	rf_test_table <- table(actual = rf_test$status,
                       predicted = predict(random_forest_model, newdata=rf_test, type="class"))
	prediction_success <- 100 * sum(diag(rf_test_table)/sum(rf_test_table))
	return(prediction_success)
}

rf_enhancer_fun_500 <- function(dataset, datatable, trees = 500){
    rf_enhancer <- strata_5var(dataset, "status", datatable, 2/3, 2, 10, 7, 150)
	rf_training <- rf_enhancer[[1]]
	rf_test <- rf_enhancer[[2]]
	random_forest_model <- randomForest(status ~., data = rf_training, ntree = trees)
	rf_training_table <- table(actual = rf_training$status,
                       predicted = predict(random_forest_model, type="class"))
	rf_test_table <- table(actual = rf_test$status,
                       predicted = predict(random_forest_model, newdata=rf_test, type="class"))
	prediction_success <- 100 * sum(diag(rf_test_table)/sum(rf_test_table))
	return(prediction_success)
}

rf_enhancer_fun_1000 <- function(dataset, datatable, trees = 1000){
    rf_enhancer <- strata_5var(dataset, "status", datatable, 2/3, 2, 10, 7, 150)
	rf_training <- rf_enhancer[[1]]
	rf_test <- rf_enhancer[[2]]
	random_forest_model <- randomForest(status ~., data = rf_training, ntree = trees)
	rf_training_table <- table(actual = rf_training$status,
                       predicted = predict(random_forest_model, type="class"))
	rf_test_table <- table(actual = rf_test$status,
                       predicted = predict(random_forest_model, newdata=rf_test, type="class"))
	prediction_success <- 100 * sum(diag(rf_test_table)/sum(rf_test_table))
	return(prediction_success)
}


rf_enhancer_comparison_100 <- resample_default(rf_enhancer_fun_100, crm_off, tab_on, 100)
rf_enhancer_comparison_100 #about 84%
rf_enhancer_comparison_500 <- resample_default(rf_enhancer_fun_100, crm_off, tab_on, 100)
rf_enhancer_comparison_500 #about 84%
rf_enhancer_comparison_1000 <- resample_default(rf_enhancer_fun_1000, crm_off, tab_on, 100)
rf_enhancer_comparison_1000 #about 85%

crm_on_svm <- svm(status ~ ., data = on_training, kernel = polynomial)
crm_on_svm_training_table <- table(actual = on_training$status,
                          predicted = predict(crm_on_svm, type="class"))
						  
crm_on_svm_test_table <- table(actual = on_test$status,
                          predicted = predict(crm_on_svm, newdata=on_test, type="class"))

crm_on_svm_test_table  
100*sum(diag(crm_on_svm_test_table)/sum(crm_on_svm_test_table))


enhancer_random_forest_model <- randomForest(status ~., data = on_training, ntree = 1000)
rf_training_table <- table(actual = on_training$status,
                          predicted = predict(enhancer_random_forest_model, type="class"))
rf_test_table <- table(actual = on_test$status,
                          predicted = predict(enhancer_random_forest_model, newdata=on_test, type="class"))
100 * sum(diag(rf_test_table)/sum(rf_test_table))



crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 10, 7, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]
nrow(on_training)
nrow(on_test)


crm_test_test <- strata_5var(test_test, "status", tab_test, 1, 4, 13, 7, 55)
test_test_set <- crm_test_test[[1]]
nrow(test_test_set)

crm_test_test_svm <- (predicted = predict(crm_on_svm, newdata = test_test_set, type = "class"))
crm_test_test_svm
crm_test_test_rf <- (predicted = predict(enhancer_random_forest_model, newdata = test_test_set, type = "class"))


new_data_frame <- cbind(test_test_set, crm_test_test_svm, crm_test_test_rf)

write.table(new_data_frame, "R_export2.txt", sep="\t")
