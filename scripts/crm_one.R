require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(mda)      # flexible discriminant analysis and mixture discrimant analysis
require(class)    # k nearest neighbours (knn)
require(adabag)   # bagging()
require(ada)      # boosting function, ada()
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 
require(plyr)

#maybe two functions
#is this an enhancer, if yes, what's its expression pattern

setwd("shadow_enhancers")
crm_one <- read.table("stark_unchecked_AP_DV.tsv", header = TRUE)
crm_off <- read.table("on_vs_off2.tsv", header = TRUE)
test_test <- read.table("temp_test.tsv", header = TRUE)

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

crm_one_train_test <- strata_5var(crm_one, "AP", tab_expr, 0.5, 2, 11) 
crm_training <- crm_one_train_test[[1]]
crm_test <- crm_one_train_test[[2]]
nrow(crm_training)
nrow(crm_test)
nrow(crm_one)

crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 0.5, 2, 10, 7, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]
nrow(on_training)
nrow(on_test)

crm_test_test <- strata_5var(test_test, "status", tab_test, 0.5, 4, 13, 7, 29)
test_test_set <- crm_test_test[[1]]
nrow(test_test_set)

linDiscrim_crm <- lda(formula=AP ~.,data = crm_training, tol = 1.0e-8, CV=FALSE) #modified tolerance
crm_lda_training_table <- table(actual = crm_training$AP,
                          predicted = predict(linDiscrim_crm, newdata=crm_training)$class)
crm_lda_training_table
100*sum(diag(crm_lda_training_table)/sum(crm_lda_training_table))
crm_lda_test_table <- table(actual = crm_test$AP,
                          predicted = predict(linDiscrim_crm, newdata=crm_test)$class)
crm_lda_test_table    
100*sum(diag(crm_lda_test_table)/sum(crm_lda_test_table))

crm_svm <- svm(AP ~., data = crm_training)
crm_svm_training_table <- table(actual = crm_training$AP,
                          predicted = predict(crm_svm, type="class"))
crm_svm_training_table           
100*sum(diag(crm_svm_training_table)/sum(crm_svm_training_table))               

crm_svm_test_table <- table(actual = crm_test$AP,
                          predicted = predict(crm_svm, newdata=crm_test, type="class"))
crm_svm_test_table  
100*sum(diag(crm_svm_test_table)/sum(crm_svm_test_table))

crm_on_svm <- svm(status ~ ., data = on_training)
crm_on_svm_training_table <- table(actual = on_training$status,
                          predicted = predict(crm_on_svm, type="class"))
						  
crm_on_svm_test_table <- table(actual = on_test$status,
                          predicted = predict(crm_on_svm, newdata=on_test, type="class"))

crm_on_svm_test_table  
100*sum(diag(crm_on_svm_test_table)/sum(crm_on_svm_test_table))

crm_test_test_svm <- (predicted = predict(crm_on_svm, newdata = test_test_set, type = "class"))