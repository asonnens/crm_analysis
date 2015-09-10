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

crm_off <- read.table("on_vs_off4.tsv", header = TRUE) #training vs testing for on/off
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

#function for svm linear on enhancer status (on vs. off)

svm_enhancers_linear <- function(dataset, datatable, variable_pos, variable_name){
    svm_enhancers <- strata_5var(dataset, variable_name, datatable, 2/3, variable_pos, 10, 7, 300)
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
      size = c(training_num,training_num))
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

svm_enhancer_comparison_linear <- resample_default(svm_enhancers_linear, crm_off, "status", tab_on, 10, 2)
svm_enhancer_comparison_linear 

#predicting on and off, proceeding with predicted list

svm_on_off_predict <- svm(status ~., data = crm_off, class.weights = tab_on, kernel = "linear")

strata_var <- function(dataset, variable_name, variable_table, variable_pos, dataframe_length = 14, second_val = 7, training_num = 30){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", size = c(50,50))
      #size = c(round(variable_table[1] *(percentage)), round(variable_table[2] * (percentage))))
  my_data <- dataset[, c(1,variable_pos, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}


svm_enhancers_linear <- function(dataset, datatable, variable_pos, variable_name){
    #svm_enhancers <- strata_5var(dataset, variable_name, datatable, variable_pos, 10, 7, 150)
	svm_enhancers <- strata_var(dataset, variable_name, datatable, variable_pos, 21, 10, 60)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_weights <- table(svm_test$Ectoderm)
	svm_model <- svm(Ectoderm ~., data = svm_training, class.weights = svm_weights, kernel = "linear")
	svm_training_table <- table(actual = svm_training$Ectoderm,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$Ectoderm,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}


svm_ectoderm_linear <- resample_default(svm_enhancers_linear, crm_exp, "Ectoderm", tab_ectoderm, 10, 8)
svm_ectoderm_linear

crm_exp_train_test <- strata_var(crm_exp, "Ectoderm", tab_ectoderm, 1/2, 8, 21, 10, 20)
crm_exp_train_test <- strata_var(crm_exp, "Ectoderm", tab_ectoderm, 2/3, 8, 21, 10, 15)



crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 13, 7, 150)
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


temp_svm <- ksvm(status ~., data =on_training,type = "C-svc", kernel = "vanilladot")
temp_predict <- predict(temp_svm, on_test)
temp_svm_results <- table(on_test$status, temp_predict)
sum(diag(temp_svm_results))/sum(temp_svm_results)
ypredscore = predict(temp_svm,on_test,type="decision")
pred <- prediction(ypredscore, on_test$status)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, main = "ROC linear SVM 'on vs off'")



crm_exp <- read.table("expression_300.tsv", header = TRUE) 

tab_on <- table(crm_exp$Ectoderm)
tab_on
length(tab_on)

strata_5var <- function(dataset, variable_name, variable_table, percentage, variable_pos, dataframe_length = 14, second_val = 7, training_num = 30){
#this function stratifies a dataset by sixteen variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - second_val
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
      size = c(50,50))
      #size = c(round(variable_table[1] *(percentage)), round(variable_table[2] * (percentage))))
  my_data <- dataset[, c(1,variable_pos,9,10,5,6,7, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

variable_temp<- colnames(svm_training[2])
which(colnames(svm_training)=="Ectoderm")
run_svm_linear <- function(variable_name, variable_table, train_test){
    svm_training <- train_test[[1]]
    svm_test <- train_test[[2]]
	variable_pos <- which(colnames(svm_training)==variable_name)
	variable_info_train <- svm_training[ ,names(svm_training) %in% variable_name]
	variable_info_test <- svm_test[ ,names(svm_test) %in% variable_name]
	variable_training <- svm_training[,c(variable_pos,8:length(svm_training))]
	variable_test <- svm_test[,c(variable_pos,8:length(svm_test))]
	variable_ID <- colnames(variable_training[variable_pos])
	variable_weights <- table(variable_info_test)
	#variable_formula <- as.formula(variable_ID ~., variable_training, na.action = na.omit)
	svm_linear <- svm(Ectoderm ~., data = variable_training, class.weights = variable_weights, kernel = "linear", na.action = na.omit)
	svm_linear_results <- predict(svm_linear, newdata= variable_test, type="class")
	svm_test_table <- table(actual = variable_test[ ,names(variable_test) %in% variable_name], predicted =svm_linear_results) 
	accuracy <- sum(diag(svm_test_table)/sum(svm_test_table))
	return(accuracy)
}


table(svm_training[ ,names(svm_training) %in% "Ectoderm"])

crm_exp_train_test <- strata_5var(crm_exp, "Ectoderm", tab_ectoderm, 1/2, 8, 21, 10, 20)
run_svm_linear <- ("Ectoderm", tab_ectoderm, crm_exp_train_test)

svm_training <- crm_exp_train_test[[1]]
svm_test <- crm_exp_train_test[[2]]
ectoderm_training <- svm_training[,8:length(svm_training)]
ectoderm_weights <- table(ectoderm_test$Ectoderm)
ectoderm_svm_linear <- svm(Ectoderm ~., data = ectoderm_training, class.weights = ectoderm_weights, kernel = "linear")
ectoderm_svm_linear_results <- predict(ectoderm_svm_linear, newdata=ectoderm_test, type="class")
actual_ectoderm <- svm_test$Ectoderm
svm_test_table <- table(actual = svm_test$Ectoderm,
              predicted = ectoderm_svm_linear_results)
ectoderm_test_table <- table(actual = ectoderm_test$Ectoderm, predicted =ectoderm_svm_linear_results)  
ectoderm_test_table
sum(diag(ectoderm_test_table)/sum(ectoderm_test_table))
ect_svm_lin = c()
for(i in (1:length(actual_ectoderm))){if ((actual_ectoderm[i] == "on") && (ectoderm_svm_linear_results[i] == "on")){ect_svm_lin[i] <- 2} 
                                      else if ((actual_ectoderm[i] == "on") && (ectoderm_svm_linear_results[i] == "off")) {ect_svm_lin[i] <- -2}
									  else if ((actual_ectoderm[i] == "off") && (ectoderm_svm_linear_results[i] == "on")) {ect_svm_lin[i] <- -1}
									  else {ect_svm_lin[i] <- 1}}

									

endoderm_training <- cbind(svm_training$Endoderm, svm_training[,9:length(svm_training)])
mesoderm_training <- cbind(svm_training$Mesoderm, svm_training[,10:length(svm_training)])
ectoderm_test <- cbind(svm_test$Ectoderm, svm_test[,8:length(svm_training)])
endoderm_test <- cbind(svm_test$Endoderm, svm_test[,8:length(svm_training)])
mesoderm_test <- cbind(svm_test$Mesoderm, svm_test[,8:length(svm_training)])
colnames(ectoderm_training)[1] <- "Ectoderm"
colnames(endoderm_training)[1] <- "Endoderm"
colnames(mesoderm_training)[1] <- "Mesoderm"
colnames(ectoderm_test)[1] <- "Ectoderm"
colnames(endoderm_test)[1] <- "Endoderm"
colnames(mesoderm_test)[1] <- "Mesoderm"
ectoderm_weights <- table(ectoderm_test$Ectoderm)
endoderm_weights <- table(endoderm_test$Endoderm)
mesoderm_weights <- table(mesoderm_test$Mesoderm)
ectoderm_rf_500 <- randomForest(Ectoderm ~., data = ectoderm_training, class.weights = ectoderm_weights, ntree = 500)
ectoderm_svm_linear <- svm(Ectoderm ~., data = ectoderm_training, class.weights = ectoderm_weights, kernel = "linear")
ectoderm_svm_radial <- svm(Ectoderm ~., data = ectoderm_training, class.weights = ectoderm_weights, kernel = "radial")
ectoderm_svm_sigmoid <- svm(Ectoderm ~., data = ectoderm_training, class.weights = ectoderm_weights, kernel = "sigmoid")
endoderm_rf_500 <- randomForest(Endoderm ~., data = endoderm_training, class.weights = endoderm_weights, ntree = 500)
endoderm_svm_linear <- svm(Endoderm ~., data = endoderm_training, class.weights = endoderm_weights, kernel = "linear")
endoderm_svm_radial <- svm(Endoderm ~., data = endoderm_training, class.weights = endoderm_weights, kernel = "radial")
endoderm_svm_sigmoid <- svm(Endoderm ~., data = endoderm_training, class.weights = endoderm_weights, kernel = "sigmoid")
mesoderm_rf_500 <- randomForest(Mesoderm ~., data = mesoderm_training, class.weights = mesoderm_weights, ntree = 500)
mesoderm_svm_linear <- svm(Mesoderm ~., data = mesoderm_training, class.weights = mesoderm_weights, kernel = "linear")
mesoderm_svm_radial <- svm(Mesoderm ~., data = mesoderm_training, class.weights = mesoderm_weights, kernel = "radial")
mesoderm_svm_sigmoid <- svm(Mesoderm ~., data = mesoderm_training, class.weights = mesoderm_weights, kernel = "sigmoid")
ectoderm_training_table <- table(actual = ectoderm_training$Ectoderm,
                       predicted = predict(ectoderm_svm_linear, type="class"))
ectoderm_rf_500_results <- predict(ectoderm_rf_500, ectoderm_test)
ectoderm_svm_linear_results <- predict(ectoderm_svm_linear, newdata=ectoderm_test, type="class")
ectoderm_svm_radial_results <- predict(ectoderm_svm_radial, newdata=ectoderm_test, type="class")
ectoderm_svm_sigmoid_results <- predict(ectoderm_svm_sigmoid, newdata=ectoderm_test, type="class")
endoderm_training_table <- table(actual = endoderm_training$Endoderm,
                       predicted = predict(endoderm_svm_linear, type="class"))
endoderm_rf_500_results <- predict(endoderm_rf_500, endoderm_test)
endoderm_svm_linear_results <- predict(endoderm_svm_linear, newdata=endoderm_test, type="class")
endoderm_svm_radial_results <- predict(endoderm_svm_radial, newdata=endoderm_test, type="class")
endoderm_svm_sigmoid_results <- predict(endoderm_svm_sigmoid, newdata=endoderm_test, type="class")
mesoderm_training_table <- table(actual = mesoderm_training$Mesoderm,
                       predicted = predict(mesoderm_svm_linear, type="class"))
mesoderm_rf_500_results <- predict(mesoderm_rf_500, mesoderm_test)
mesoderm_svm_linear_results <- predict(mesoderm_svm_linear, newdata=mesoderm_test, type="class")
mesoderm_svm_radial_results <- predict(mesoderm_svm_radial, newdata=mesoderm_test, type="class")
mesoderm_svm_sigmoid_results <- predict(mesoderm_svm_sigmoid, newdata=mesoderm_test, type="class")


actual_ectoderm <- svm_test$Ectoderm
actual_endoderm <- svm_test$Endoderm
actual_mesoderm <- svm_test$Mesoderm

svm_test_table <- table(actual = svm_test$Ectoderm,
                       predicted = ectoderm_svm_linear_results)

ect_rf_500 = c()
ect_svm_lin = c()	
ect_svm_rad = c()
ect_svm_sig = c()
end_rf_500 = c()
end_svm_lin = c()
end_svm_rad = c()
end_svm_sig = c()
mes_rf_500 = c()
mes_svm_lin = c()	
mes_svm_rad = c()
mes_svm_sig = c()
b_ect_svm = c()
u_ect_svm = c()
for(i in (1:length(actual_ectoderm))){if ((actual_ectoderm[i] == "on") && (ectoderm_rf_500_results[i] == "on")){ect_rf_500[i] <- 2} 
                                      else if ((actual_ectoderm[i] == "on") && (ectoderm_rf_500_results[i] == "off")) {ect_rf_500[i] <- -2}
									  else if ((actual_ectoderm[i] == "off") && (ectoderm_rf_500_results[i] == "on")) {ect_rf_500[i] <- -1}
									  else {ect_rf_500[i] <- 1}}		   
for(i in (1:length(actual_ectoderm))){if ((actual_ectoderm[i] == "on") && (ectoderm_svm_linear_results[i] == "on")){ect_svm_lin[i] <- 2} 
                                      else if ((actual_ectoderm[i] == "on") && (ectoderm_svm_linear_results[i] == "off")) {ect_svm_lin[i] <- -2}
									  else if ((actual_ectoderm[i] == "off") && (ectoderm_svm_linear_results[i] == "on")) {ect_svm_lin[i] <- -1}
									  else {ect_svm_lin[i] <- 1}}
for(i in (1:length(actual_ectoderm))){if ((actual_ectoderm[i] == "on") && (ectoderm_svm_radial_results[i] == "on")){ect_svm_rad[i] <- 2} 
                                      else if ((actual_ectoderm[i] == "on") && (ectoderm_svm_radial_results[i] == "off")) {ect_svm_rad[i] <- -2}
									  else if ((actual_ectoderm[i] == "off") && (ectoderm_svm_radial_results[i] == "on")) {ect_svm_rad[i] <- -1}
									  else {ect_svm_rad[i] <- 1}}
for(i in (1:length(actual_ectoderm))){if ((actual_ectoderm[i] == "on") && (ectoderm_svm_sigmoid_results[i] == "on")){ect_svm_sig[i] <- 2} 
                                      else if ((actual_ectoderm[i] == "on") && (ectoderm_svm_sigmoid_results[i] == "off")) {ect_svm_sig[i] <- -2}
									  else if ((actual_ectoderm[i] == "off") && (ectoderm_svm_sigmoid_results[i] == "on")) {ect_svm_sig[i] <- -1}
									  else {ect_svm_sig[i] <- 1}}
for(i in (1:length(actual_endoderm))){if ((actual_endoderm[i] == "on") && (endoderm_rf_500_results[i] == "on")){end_rf_500[i] <- 2} 
                                      else if ((actual_endoderm[i] == "on") && (endoderm_rf_500_results[i] == "off")) {end_rf_500[i] <- -2}
									  else if ((actual_endoderm[i] == "off") && (endoderm_rf_500_results[i] == "on")) {end_rf_500[i] <- -1}
									  else {end_rf_500[i] <- 1}}
for(i in (1:length(actual_endoderm))){if ((actual_endoderm[i] == "on") && (endoderm_svm_linear_results[i] == "on")){end_svm_lin[i] <- 2} 
                                      else if ((actual_endoderm[i] == "on") && (endoderm_svm_linear_results[i] == "off")) {end_svm_lin[i] <- -2}
									  else if ((actual_endoderm[i] == "off") && (endoderm_svm_linear_results[i] == "on")) {end_svm_lin[i] <- -1}
									  else {end_svm_lin[i] <- 1}}
for(i in (1:length(actual_endoderm))){if ((actual_endoderm[i] == "on") && (endoderm_svm_radial_results[i] == "on")){end_svm_rad[i] <- 2} 
                                      else if ((actual_endoderm[i] == "on") && (endoderm_svm_radial_results[i] == "off")) {end_svm_rad[i] <- -2}
									  else if ((actual_endoderm[i] == "off") && (endoderm_svm_radial_results[i] == "on")) {end_svm_rad[i] <- -1}
									  else {end_svm_rad[i] <- 1}}
for(i in (1:length(actual_endoderm))){if ((actual_endoderm[i] == "on") && (endoderm_svm_sigmoid_results[i] == "on")){end_svm_sig[i] <- 2} 
                                      else if ((actual_endoderm[i] == "on") && (endoderm_svm_sigmoid_results[i] == "off")) {end_svm_sig[i] <- -2}
									  else if ((actual_endoderm[i] == "off") && (endoderm_svm_sigmoid_results[i] == "on")) {end_svm_sig[i] <- -1}
									  else {end_svm_sig[i] <- 1}}
for(i in (1:length(actual_mesoderm))){if ((actual_mesoderm[i] == "on") && (mesoderm_rf_500_results[i] == "on")){mes_rf_500[i] <- 2} 
                                      else if ((actual_mesoderm[i] == "on") && (mesoderm_rf_500_results[i] == "off")) {mes_rf_500[i] <- -2}
									  else if ((actual_mesoderm[i] == "off") && (mesoderm_rf_500_results[i] == "on")) {mes_rf_500[i] <- -1}
									  else {mes_rf_500[i] <- 1}}
for(i in (1:length(actual_mesoderm))){if ((actual_mesoderm[i] == "on") && (mesoderm_svm_linear_results[i] == "on")){mes_svm_lin[i] <- 2} 
                                      else if ((actual_mesoderm[i] == "on") && (mesoderm_svm_linear_results[i] == "off")) {mes_svm_lin[i] <- -2}
									  else if ((actual_mesoderm[i] == "off") && (mesoderm_svm_linear_results[i] == "on")) {mes_svm_lin[i] <- -1}
									  else {mes_svm_lin[i] <- 1}}
for(i in (1:length(actual_mesoderm))){if ((actual_mesoderm[i] == "on") && (mesoderm_svm_radial_results[i] == "on")){mes_svm_rad[i] <- 2} 
                                      else if ((actual_mesoderm[i] == "on") && (mesoderm_svm_radial_results[i] == "off")) {mes_svm_rad[i] <- -2}
									  else if ((actual_mesoderm[i] == "off") && (mesoderm_svm_radial_results[i] == "on")) {mes_svm_rad[i] <- -1}
									  else {mes_svm_rad[i] <- 1}}
for(i in (1:length(actual_mesoderm))){if ((actual_mesoderm[i] == "on") && (mesoderm_svm_sigmoid_results[i] == "on")){mes_svm_sig[i] <- 2} 
                                      else if ((actual_mesoderm[i] == "on") && (mesoderm_svm_sigmoid_results[i] == "off")) {mes_svm_sig[i] <- -2}
									  else if ((actual_mesoderm[i] == "off") && (mesoderm_svm_sigmoid_results[i] == "on")) {mes_svm_sig[i] <- -1}
									  else {mes_svm_sig[i] <- 1}}

new_frame <- cbind(svm_test[,1:7], ect_rf_500, ect_svm_lin, ect_svm_rad, ect_svm_sig, 
                                   end_rf_500, end_svm_lin, end_svm_rad, end_svm_sig, 
								   mes_rf_500, mes_svm_lin, mes_svm_rad, mes_svm_sig)

temp_frame <- new_frame[,8:19]
heatmap.2(as.matrix(temp_frame),symm=FALSE,trace="none", col = c("darkred", "orange", "lightgreen", "darkgreen"))
heatmap.2(as.matrix(temp_frame),symm=FALSE,trace="none", Colv = FALSE, Rowv = FALSE, margins=c(8,1), col = c("darkred", "orange", "white", "lightgreen", "darkgreen"))

ectoderm_test_table <- table(actual = ectoderm_test$Ectoderm, predicted =ectoderm_svm_linear_results)  
ectoderm_test_table
sum(diag(ectoderm_test_table)/sum(ectoderm_test_table))


write.table(new_frame, "R_test_unbalanced2.txt", sep="\t")



svm_enhancers_ectoderm_linear <- function(dataset, datatable, kernel_type = "linear"){
    svm_enhancers <- strata_5var(dataset, "Ectoderm", datatable, 2/3, 8, 21, 10, 15)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(Ectoderm ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$Ectoderm,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$Ectoderm,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_enhancer_ectoderm_comparison_linear <- resample_default(svm_enhancers_ectoderm_linear, crm_exp, tab_on, 100)
svm_enhancer_ectoderm_comparison_linear 

rf_enhancer_fun_500 <- function(dataset, datatable, trees = 500){
    rf_enhancer <- strata_5var(dataset, "Mesoderm", datatable, 2/3, 10, 21, 10, 15)
	rf_training <- rf_enhancer[[1]]
	rf_test <- rf_enhancer[[2]]
	random_forest_model <- randomForest(Mesoderm ~., data = rf_training, ntree = trees)
	rf_training_table <- table(actual = rf_training$Mesoderm,
                       predicted = predict(random_forest_model, type="class"))
	rf_test_table <- table(actual = rf_test$Mesoderm,
                       predicted = predict(random_forest_model, newdata=rf_test, type="class"))
	prediction_success <- 100 * sum(diag(rf_test_table)/sum(rf_test_table))
	return(prediction_success)
}
tab_on <- table(crm_exp$Endoderm)
tab_on
length(tab_on)


crm_exp_train_test <- strata_5var(crm_exp, "Endoderm", tab_on, 2/3, 9, 21, 10, 15)
on_training <- crm_exp_train_test[[1]]
on_test <- crm_exp_train_test[[2]]

svm_enhancers_endoderm_linear <- function(dataset, datatable, kernel_type = "linear"){
    svm_enhancers <- strata_5var(dataset, "Endoderm", datatable, 2/3, 9, 21, 10, 15)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(Endoderm ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$Endoderm,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$Endoderm,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_enhancer_comparison_linear <- resample_default(svm_enhancers_endoderm_linear, crm_exp, tab_on, 100)
svm_enhancer_comparison_linear 

tab_on <- table(crm_exp$Mesoderm)
tab_on
length(tab_on)

crm_exp_train_test <- strata_5var(crm_exp, "Mesoderm", tab_on, 2/3, 10, 21, 10, 15)
on_training <- crm_exp_train_test[[1]]
on_test <- crm_exp_train_test[[2]]

svm_enhancers_mesoderm_linear <- function(dataset, datatable, kernel_type = "linear"){
    svm_enhancers <- strata_5var(dataset, "Mesoderm", datatable, 2/3, 10, 21, 10, 15)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(Mesoderm ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$Mesoderm,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$Mesoderm,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_enhancer_comparison_linear <- resample_default(svm_enhancers_mesoderm_linear, crm_exp, tab_on, 100)
svm_enhancer_comparison_linear 

svm_enhancers_fun_radial <- function(dataset, datatable, kernel_type = "radial"){
    svm_enhancers <- strata_5var(dataset, "P", datatable, 2/3, 7, 21, 10, 15)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(P~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$P,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$P,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_enhancers_Endoderm_radial <- function(dataset, datatable, kernel_type = "radial"){
    svm_enhancers <- strata_5var(dataset, "Endoderm", datatable, 2/3, 9, 21, 10, 15)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(Endoderm ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$Endoderm,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$Endoderm,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_enhancers_fun_polynomial <- function(dataset, datatable, kernel_type = "polynomial"){
    svm_enhancers <- strata_5var(dataset, "Mesoderm", datatable, 2/3, 10, 21, 10, 15)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(Mesoderm ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$Mesoderm,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$Mesoderm,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_enhancers_fun_sigmoid <- function(dataset, datatable, kernel_type = "sigmoid"){
    svm_enhancers <- strata_5var(dataset, "Mesoderm", datatable, 2/3, 10, 21, 10, 15)
	svm_training <- svm_enhancers[[1]]
	svm_test <- svm_enhancers[[2]]
	svm_model <- svm(Mesoderm ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$Mesoderm,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$Mesoderm,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}


svm_enhancer_comparison_linear <- resample_default(svm_enhancers_ectoderm_linear, crm_exp, tab_on, 100)
svm_enhancer_comparison_linear #66% Ectoderm, 72% endoderm, 50% mesoderm, 47% A, 49% C, 48% P
svm_enhancer_comparison_radial <- resample_default(svm_enhancers_fun_radial, crm_exp, tab_on, 100)
svm_enhancer_comparison_radial #68% Ectoderm, #72% Endoderm, 51% Mesoderm, 47% A, 51% C, 48% P
svm_enhancer_comparison_polynomial <- resample_default(svm_enhancers_fun_radial, crm_exp, tab_on, 100)
svm_enhancer_comparison_polynomial #68% Ectoderm, 68% Endoderm,50% Mesoderm
svm_enhancer_comparison_sigmoid <- resample_default(svm_enhancers_fun_radial, crm_exp, tab_on, 100)
svm_enhancer_comparison_sigmoid #67% Ectoderm, 48% endoderm, 50% Mesoderm

rf_enhancer_comparison_500 <- resample_default(rf_enhancer_fun_500, crm_exp, tab_on, 100)
rf_enhancer_comparison_500 #64% Ectoderm, 74% Endoderm, 49% mesoderm  57% A, 54% C, 48% P



temp_svm <- ksvm(A ~., data =on_training,type = "C-svc", kernel = "vanilladot")
temp_predict <- predict(temp_svm, on_test)
temp_svm_results <- table(on_test$A, temp_predict)
sum(diag(temp_svm_results))/sum(temp_svm_results)

dataset, "P", datatable, 2/3, 7, 21, 10, 15

#original MacArthur data
crm_exp_train_test <- strata_5var(crm_exp, "P", tab_on, 2/3, 7, 21, 10, 15)
on_training <- crm_exp_train_test[[1]]
on_test <- crm_exp_train_test[[2]]

crm_exp_train_test <- strata_5var(crm_exp, "P", tab_on, 2/3, 7, 20, 2, 15)
on_training <- crm_exp_train_test[[1]]
on_test <- crm_exp_train_test[[2]]
#minus Dorsal
crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 13, 6, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]
#minus Dorsal and Hairy
crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 13, 5, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]
#minus dorsal, hairy and snail
crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 13, 4, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]
#minus dorsal, hairy, snail and twist
crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 13, 3, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]

#just bicoid, cad, and gt
crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 12, 2, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]
#just cad and gt
crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 12, 1, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]
#just giant and hunchback
crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 13, 1, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]



svm_enhancers_fun_linear <- function(dataset, datatable, kernel_type = "linear"){
    svm_enhancers <- strata_5var(dataset, "status", datatable, 2/3, 2, 13, 10, 150)
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
svm_enhancer_comparison_linear #82%


crm_off_train_test <- strata_5var(crm_off, "status", tab_on, 2/3, 2, 12, 1, 150)
on_training <- crm_off_train_test[[1]]
on_test <- crm_off_train_test[[2]]

temp_svm <- ksvm(P ~., data =on_training,type = "C-svc", kernel = "vanilladot")
temp_predict <- predict(temp_svm, on_test)
temp_svm_results <- table(on_test$P, temp_predict)
sum(diag(temp_svm_results))/sum(temp_svm_results)
ypredscore = predict(temp_svm,on_test,type="decision")
pred <- prediction(ypredscore, on_test$P)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, main = "ROC linear SVM, gt & cad,  'on vs off'")

new_data_frame <- cbind(test_test[,1],test_test_set, crm_test_test_svm, crm_test_test_rf)

write.table(new_data_frame, "R_export3.txt", sep="\t")
