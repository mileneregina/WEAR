# Clear the workspace by removing all objects
rm(list = ls())
# List all objects in the environment
ls()
# Run garbage collection to free memory
gc()

# Load necessary libraries
library(glmnet)
library(corrplot)
library(randomForest)
library(formattable)
library(boot)
library(pls)
library(GGally)
library(lavaan)
library(dplyr)
library(glmnet)
library(plotmo)
library(MASS)
library(leaps)
library(rpart)
library(rpart.plot)
library(mlbench)
library(caret)
library(tidyverse)
library(readxl)
library(broom)
library(knitr)
library(ggplot2)
library(gridExtra)
library(ggfortify)
library(car)
library(readr)
library(xtable)
library(gtable)
library(caret)
library(dplyr)
library(arrow)
library(ranger)

# Define the file name to store the downloaded dataset
destfile <- "arquivo.zip"

# Define the URL for downloading the dataset
url <- "https://archive.ics.uci.edu/static/public/265/physicochemical+properties+of+protein+tertiary+structure.zip"

# Download the dataset from the URL
download.file(url, destfile)

# Unzip the downloaded file into a specific directory
unzip(destfile, exdir = "data_protein")

# List the files extracted in the directory
list.files("data_protein")

# Load the dataset into a dataframe
df <- read.csv("data_protein/CASP.csv")

# Initialize variables to store Mean Squared Errors (MSE) for different models
mse_forest_gold = NULL
mse_tree_gold  = NULL
mse_lm_gold = NULL
mse_lasso_gold = NULL

mse_forest_weight = NULL
mse.tree.weight= NULL
mse.lm.weight= NULL
mse.lasso.weight = NULL

mse_forest_mean = NULL
mse.tree.mean= NULL
mse.lm.mean= NULL
mse.lasso.mean = NULL

# Initialize variables to store Standard Errors (SE) for different models
se_forest_gold = NULL
se_tree_gold  = NULL
se_lm_gold = NULL
se_lasso_gold = NULL

se_forest_weight = NULL
se.tree.weight= NULL
se.lm.weight= NULL
se.lasso.weight = NULL

se_forest_mean = NULL
se.tree.mean= NULL
se.lm.mean= NULL
se.lasso.mean = NULL

# Rename the target variable from "RMSD" to "y"
names(df)[names(df) == "RMSD"] <- "y"

# Set a seed for reproducibility
set.seed(123)

# Define key parameters
count_row <- nrow(df) # Number of rows in dataset
count_specialist <- 4 # Number of specialists
count_varibles <- ncol(df) # Number of variables in dataset
tamanho_matrix <- (count_varibles-4) * count_row # Matrix size calculation
percentual_train <- 0.7 # Percentage for training data
percentual_validation <- 0.1 # Percentage for validation data
percentual_test <- 1 - percentual_train - percentual_validation # Remaining percentage for testing
total_mse = 1 # Placeholder for MSE calculations

# Create matrices to store specialist predictions and deviations
specialist = matrix(data = NA, nrow = count_row, ncol = count_specialist)
deviation.specialist = matrix(data = NA, nrow = count_row, ncol = count_specialist)
lambdas = matrix(data = NA, nrow = total_mse, ncol = count_specialist)
weight.specialist.reg.linear = matrix(data = NA, nrow = total_mse, ncol = count_specialist)
weight.specialist.forest = matrix(data = NA, nrow = total_mse, ncol = count_specialist)
weight.specialist.tree = matrix(data = NA, nrow = total_mse, ncol = count_specialist)
weight.specialist.lasso = matrix(data = NA, nrow = total_mse, ncol = count_specialist)
mean_var_est_specialist = matrix(data = NA, nrow =1, ncol = count_specialist)

# Define deviations for specialists
devitation_specialist =  c(1, 2, 5, 15)
deviation_t = t(devitation_specialist)

# Generate noisy simulated specialist assessments using normal distribution
for (j in 1:count_specialist) {
  specialist[,j] =  df[,1] + rnorm(count_row, 0, deviation_t[,j])
}

# Extract covariates from dataset
covariable = df[,-1]

# Split data into training, validation, and testing sets
set.seed(123)
fractionTraining   <- percentual_train
fractionValidation <- percentual_validation
fractionTest       <- percentual_test

# Compute sample sizes for each set
sampleSizeTraining   <- floor(fractionTraining   * nrow(df))
sampleSizeValidation <- floor(fractionValidation * nrow(df))
sampleSizeTest       <- floor(fractionTest       * nrow(df))

# Generate indices for splitting data into training, validation, and testing sets
indicesTraining    <- sort(sample(seq_len(nrow(df)), size=sampleSizeTraining))
indicesNotTraining <- setdiff(seq_len(nrow(df)), indicesTraining)
indicesValidation  <- sort(sample(indicesNotTraining, size=sampleSizeValidation))
indicesTest        <- setdiff(indicesNotTraining, indicesValidation)

# Assign data to training, validation, and testing sets
train   <- df[indicesTraining, ]
validation <- df[indicesValidation, ]
test      <- df[indicesTest, ]

# Assign specialist assessments to corresponding sets
specialist.train <- specialist[indicesTraining, ]
specialist.test <- specialist[indicesTest, ]
specialist.validation <- specialist[indicesValidation, ]

# Assign covariate data to corresponding sets
covariable.train <- covariable[indicesTraining, ]
covariable.test <- covariable[indicesTest, ]
covariable.validation <- covariable[indicesValidation, ]

# Display dimensions of validation, test, and training sets
dim(covariable.validation)
dim(covariable.test)
dim(covariable.train)


#forest


# Train a Random Forest model using all features
fit_forest_gold <- ranger(y ~ ., data = train)
# Predict using the trained model on the test set
pred_forest_gold <- predict(fit_forest_gold, data = test)$predictions
# Compute Mean Squared Error (MSE) for the forest model
mse_forest_gold <- sum((pred_forest_gold - test$y)^2) / length(test$y)
# Compute squared errors
w_forest_gold <- (pred_forest_gold - test$y)^2 
# Compute standard deviation of squared errors
w_forest_gold_deviation <- sd(w_forest_gold)
# Compute standard error
se_forest_gold <- w_forest_gold_deviation / sqrt(length(test$y))

# Compute the mean predictions across specialists
mean.specialist = rowMeans(specialist) 
# Create a dataset combining mean specialist predictions and covariates
banco.specialist = data.frame(mean.specialist, covariable)
# Split the dataset into training, validation, and testing sets
train.specialist <- banco.specialist[indicesTraining, ] 
test.specialist <- banco.specialist[indicesTest, ] 
validation.specialist <- banco.specialist[indicesValidation, ]

# Train a Random Forest model using mean specialist predictions
fit_forest_mean_specialist <- ranger(mean.specialist ~ ., data = train.specialist)
# Predict using the trained model on the test set
pred_forest_mean_specialist <- predict(fit_forest_mean_specialist, data = test.specialist)$predictions

# Compute MSE for the mean specialist model
mse_forest_mean_specialist <- sum((pred_forest_mean_specialist - test$y)^2) / length(test$y) 
# Compute squared errors
w_forest_mean_specialist <- (pred_forest_mean_specialist - test$y)^2 
# Compute standard deviation of squared errors
w_forest_mean_specialist_deviation <- sd(w_forest_mean_specialist) 
# Compute standard error
se_forest_mean_specialist <- w_forest_mean_specialist_deviation / sqrt(length(test$y))

#####################################################
# Compute MSE for individual specialists using Random Forest
specialist.mse.weight_forest <- NULL
for (i in 1:count_specialist) {
  specialistecialistatest <- specialist.test[, i]
  specialistecialistatrain <- specialist.train[, i]
  specialistecilistavalidation <- specialist.validation[, i]
  df.weight.test <- data.frame(specialistecialistatest, covariable.test)
  df.weight.train <- data.frame(specialistecialistatrain, covariable.train)
  df.weight.validation <- data.frame(specialistecilistavalidation, covariable.validation)
  
  # Train a Random Forest model for each specialist
  model_forest <- ranger(specialistecialistatrain ~ ., data = df.weight.train)
  # Predict using the trained model
  pred.specialist_forest <- predict(model_forest, data = df.weight.validation)$predictions
  
  # Compute MSE for each specialist
  specialist.mse.weight_forest[i] <- sum((pred.specialist_forest - specialist.validation[, i])^2) / length(validation$y)
}

# Compute weighted specialist predictions
specialist.train.weight_forest <- specialist.train 
specialist.test.weight_forest <- specialist.test

for (j in 1:count_specialist) {
  specialist.train.weight_forest[, j] <- (specialist.train[, j] * (1 / specialist.mse.weight_forest[j]))
  specialist.test.weight_forest[, j] <- (specialist.test[, j] * (1 / specialist.mse.weight_forest[j]))
}

# Store the weights of specialists in a matrix
weight.specialist.forest <- t(specialist.mse.weight_forest)

# Compute row sums for training and test specialist data
specialist.train.weight.soma_forest <- rowSums(specialist.train.weight_forest)
specialist.test.weight.soma_forest <- rowSums(specialist.test.weight_forest)

# Convert MSE specialist weights to a data frame
specialist.mse.weight_forest <- data.frame(specialist.mse.weight_forest)

# Compute denominator for weighted specialist calculations
denominador_forest <- colSums((1 / specialist.mse.weight_forest))

mean.train.weight_notas_forest <- specialist.train.weight.soma_forest / denominador_forest
mean.test.weight_notas_forest <- specialist.test.weight.soma_forest / denominador_forest

df.mean.train.weight_notas_forest <- data.frame(mean.train.weight_notas_forest, covariable.train)
df.mean.test.weight_notas_forest <- data.frame(mean.test.weight_notas_forest, covariable.test)

# Create final dataset for model training
fit10_mean.weight_forest <- ranger(mean.train.weight_notas_forest ~ ., data = df.mean.train.weight_notas_forest)
predict10_mean.weight_forest <- predict(fit10_mean.weight_forest, data = df.mean.test.weight_notas_forest)$predictions

# Compute MSE for weighted forest model
mse_forest_weight = sum((predict10_mean.weight_forest - test$y)^2) / length(test$y) 
# Compute standard deviation of squared errors
w6 = (predict10_mean.weight_forest - test$y)^2 
w6_deviation = sd(w6)
# Compute standard error
se_forest_weight = w6_deviation / sqrt(length(test$y))

######################### Decision Tree Model ##############################

# Train a Decision Tree model
fit_tree_gold = rpart(y ~ ., data = train) 
# Predict using the trained model
pred_tree_gold <- predict(fit_tree_gold, newdata = test)

# Compute MSE for the tree model
mse.tree.gold = sum((pred_tree_gold - test$y)^2) / length(test$y)
# Compute squared errors
w_tree_gold = (pred_tree_gold - test$y)^2 
# Compute standard deviation of squared errors
w_tree_gold_deviation = sd(w_tree_gold) 
# Compute standard error
se_tree_gold = w_tree_gold_deviation / sqrt(length(test$y))


#############tree 

# Compute the mean predictions across specialists
mean.specialist = rowMeans(specialist) 
# Create a dataset combining the mean specialist predictions and covariates
banco.specialist = data.frame(mean.specialist, covariable)
# Split the specialist dataset into training, validation, and testing sets
train.specialist <- banco.specialist[indicesTraining, ] 
test.specialist <- banco.specialist[indicesTest,] 
validation.specialist <- banco.specialist[indicesValidation, ]

# Train a Decision Tree model using mean specialist predictions
fit_tree_mean_specialist = rpart(mean.specialist ~ ., data = train.specialist)
# Predict using the trained model on the test data
pred_tree_mean_specialist <- predict(fit_tree_mean_specialist, newdata = test.specialist)

# Compute MSE for the mean specialist model
mse.tree.mean = sum((pred_tree_mean_specialist - test$y)^2) / length(test$y)
# Compute squared errors
w_tree_mean_specialist = (pred_tree_mean_specialist - test$y)^2 
# Compute standard deviation of squared errors
w_tree_mean_specialist_deviation = sd(w_tree_mean_specialist) 
# Compute standard error
se_tree_mean = w_tree_mean_specialist_deviation / sqrt(length(test$y))

#####################################################
# Compute MSE for individual specialists using Decision Trees
specialist.mse.weight_tree = NULL
for (i in 1:count_specialist) {
  specialistecialistatest = specialist.test[, i]
  specialistecialistatrain = specialist.train[, i]
  specialistecilistavalidation = specialist.validation[, i]
  df.weight.test = data.frame(specialistecialistatest, covariable.test)
  df.weight.train = data.frame(specialistecialistatrain, covariable.train)
  df.weight.validation = data.frame(specialistecilistavalidation, covariable.validation)
  
  # Train a Decision Tree model for each specialist
  model_tree = rpart(specialistecialistatrain ~ ., df.weight.train)
  # Predict using the trained model
  pred.specialist_tree <- predict(model_tree, newdata = df.weight.validation)
  
  # Compute MSE for each specialist
  specialist.mse.weight_tree[i] <- sum((pred.specialist_tree - specialist.validation[, i])^2) / length(validation$y)
}

# Compute weighted specialist predictions
specialist.train.weight_tree = specialist.train 
specialist.test.weight_tree = specialist.test

for (j in 1:count_specialist) {
  specialist.train.weight_tree[, j] <- (specialist.train[, j] * (1 / specialist.mse.weight_tree[j]))
  specialist.test.weight_tree[, j] <- (specialist.test[, j] * (1 / specialist.mse.weight_tree[j]))
}

# Compute final weighted specialist model
weight.specialist.tree = t(specialist.mse.weight_tree)

# Compute row sums for training and test specialist data
specialist.train.weight.soma_tree = rowSums(specialist.train.weight_tree)
specialist.test.weight.soma_tree = rowSums(specialist.test.weight_tree)

specialist.mse.weight_tree = data.frame(specialist.mse.weight_tree) 

denominador_tree = colSums((1 / specialist.mse.weight_tree))

mean.train.weight_notas_tree <- specialist.train.weight.soma_tree / denominador_tree 
mean.test.weight_notas_tree <- specialist.test.weight.soma_tree / denominador_tree

df.mean.train.weight_notas_tree <- data.frame(mean.train.weight_notas_tree,covariable.train)
df.mean.test.weight_notas_tree <- data.frame(mean.test.weight_notas_tree,covariable.test)

# Create final dataset for model training
fit10_mean.weight_tree <- rpart(mean.train.weight_notas_tree ~ ., data = df.mean.train.weight_notas_tree)
predict10_mean.weight_tree <- predict(fit10_mean.weight_tree, newdata = df.mean.test.weight_notas_tree)

# Compute MSE for weighted tree model
mse.tree.weight = sum((predict10_mean.weight_tree - test$y)^2) / length(test$y) 
# Compute standard deviation of squared errors
w6 = (predict10_mean.weight_tree - test$y)^2 
w6_deviation = sd(w6)
# Compute standard error
se.tree.weight = w6_deviation / sqrt(length(test$y))

######################### Linear Regression Model ##############################

# Train a Linear Regression model
fit_lm_gold = lm(y ~ ., data = train)
# Predict using the trained model
pred_lm_gold <- predict(fit_lm_gold, newdata = test)

# Compute MSE for the linear regression model
mse_lm_gold = sum((pred_lm_gold - test$y)^2) / length(test$y)
# Compute squared errors
w_lm_gold = (pred_lm_gold - test$y)^2 
# Compute standard deviation of squared errors
w_lm_gold_deviation = sd(w_lm_gold) 
# Compute standard error
se_lm_gold = w_lm_gold_deviation / sqrt(length(test$y))

# Compute the mean predictions across specialists for linear regression
mean.specialist = rowMeans(specialist) 
banco.specialist = data.frame(mean.specialist, covariable)
# Split the specialist dataset into training, validation, and testing sets
train.specialist <- banco.specialist[indicesTraining, ] 
test.specialist <- banco.specialist[indicesTest,] 
validation.specialist <- banco.specialist[indicesValidation, ]

# Train a Linear Regression model using mean specialist predictions
fit_lm_mean_specialist = lm(mean.specialist ~ ., data = train.specialist)
# Predict using the trained model
pred_lm_mean_specialist <- predict(fit_lm_mean_specialist, newdata = test.specialist)

# Compute MSE for the mean specialist model
mse.lm.mean = sum((pred_lm_mean_specialist - test$y)^2) / length(test$y) 
# Compute squared errors
w_lm_mean_specialist = (pred_lm_mean_specialist - test$y)^2 
# Compute standard deviation of squared errors
w_lm_mean_specialist_deviation = sd(w_lm_mean_specialist)
# Compute standard error
se_lm_mean = w_lm_mean_specialist_deviation / sqrt(length(test$y))


# Compute MSE for individual specialists using Linear Regression
specialist.mse.weight_lm = NULL
for (i in 1:count_specialist) {
  specialistecialistatest = specialist.test[, i]
  specialistecialistatrain = specialist.train[, i]
  specialistecilistavalidation = specialist.validation[, i]
  df.weight.test = data.frame(specialistecialistatest, covariable.test)
  df.weight.train = data.frame(specialistecialistatrain, covariable.train)
  df.weight.validation = data.frame(specialistecilistavalidation, covariable.validation)
  
  # Train a Linear Regression model for each specialist
  model_lm = lm(specialistecialistatrain ~ ., df.weight.train)
  # Predict using the trained model
  pred.specialist_lm <- predict(model_lm, newdata = df.weight.validation)
  
  # Compute MSE for each specialist
  specialist.mse.weight_lm[i] <- sum((pred.specialist_lm - specialist.validation[, i])^2) / length(validation$y)
}

# Compute weighted specialist predictions
specialist.train.weight_lm = specialist.train 
specialist.test.weight_lm = specialist.test

for (j in 1:count_specialist) {
  specialist.train.weight_lm[, j] <- (specialist.train[, j] * (1 / specialist.mse.weight_lm[j]))
  specialist.test.weight_lm[, j] <- (specialist.test[, j] * (1 / specialist.mse.weight_lm[j]))
}

# Compute final weighted specialist model
weight.specialist.lm = t(specialist.mse.weight_lm)

# Compute row sums for training and test specialist data
specialist.train.weight.soma_lm = rowSums(specialist.train.weight_lm) 
specialist.test.weight.soma_lm = rowSums(specialist.test.weight_lm)

specialist.mse.weight_lm = data.frame(specialist.mse.weight_lm) 

denominador_lm = colSums((1 / specialist.mse.weight_lm))

mean.train.weight_notas_lm <- specialist.train.weight.soma_lm / denominador_lm
mean.test.weight_notas_lm <- specialist.test.weight.soma_lm / denominador_lm

df.mean.train.weight_notas_lm <-data.frame(mean.train.weight_notas_lm,covariable.train)
df.mean.test.weight_notas_lm <-data.frame(mean.test.weight_notas_lm,covariable.test)

# Create final dataset for model training
fit10_mean.weight_lm <- lm(mean.train.weight_notas_lm ~ ., data = df.mean.train.weight_notas_lm)
predict10_mean.weight_lm <- predict(fit10_mean.weight_lm, newdata = df.mean.test.weight_notas_lm, type="response")

# Compute MSE for weighted linear regression model
mse.lm.weight = sum((predict10_mean.weight_lm - test$y)^2) / length(test$y) 
# Compute standard deviation of squared errors
w6 = (predict10_mean.weight_lm - test$y)^2 
w6_deviation = sd(w6)
# Compute standard error
se.lm.weight = w6_deviation / sqrt(length(test$y))

######################### LASSO Regression ##############################

# Prepare data matrices for LASSO regression
x_lasso_y_train = model.matrix(~.-1, train[, -c(1)]) 
y_lasso_y_train = train$y 
x_test = model.matrix(~.-1, test[, -c(1)]) 
y_test = test$y

# Train LASSO regression model
model_y = cv.glmnet(x_lasso_y_train, y_lasso_y_train, alpha = 1)
# Predict using the trained LASSO model
pred.specialist_y <- predict(model_y, s = model_y$lambda.min, newx = x_test)

# Compute MSE for LASSO regression
mse_lasso_gold = sum((pred.specialist_y - test$y)^2) / length(test$y) 
# Compute squared errors
w_lasso_gold = (pred.specialist_y - test$y)^2 
# Compute standard deviation of squared errors
w_lasso_gold_deviation = sd(w_lasso_gold) 
# Compute standard error
se_lasso_gold = w_lasso_gold_deviation / sqrt(length(test$y))

######################### LASSO Regression with Mean Specialist ##############################

# Prepare data matrices for LASSO regression using mean specialist predictions
x_lasso_y_mean_train = model.matrix(~.-1, train.specialist[, -c(1)])
y_lasso_y_mean_train = train.specialist$mean.specialist 
x_test_mean = model.matrix(~.-1, test.specialist[, -c(1)])

# Train LASSO regression model using mean specialist data
model_mean = cv.glmnet(x_lasso_y_mean_train, y_lasso_y_mean_train, alpha = 1)
# Predict using the trained LASSO model
pred.specialist_mean <- predict(model_mean, s = model_mean$lambda.min, newx = x_test_mean)

# Compute MSE for LASSO regression using mean specialist data
mse_lasso_mean_specialist = sum((pred.specialist_mean - test$y)^2) / length(test$y)
# Compute squared errors
w_lasso_mean_specialist = (pred.specialist_mean - test$y)^2 
# Compute standard deviation of squared errors
w_lasso_mean_specialist_deviation = sd(w_lasso_mean_specialist) 
# Compute standard error
se_lasso_mean = w_lasso_mean_specialist_deviation / sqrt(length(test$y))

# Initialize variables for LASSO regression
model = NULL 
df.weight = NULL 
specialistecialistaconsiderado = NULL 
specialist.mse.weight = NULL
specialist.mse.weight.lasso = NULL

# Compute MSE for individual specialists using LASSO regression
for (i in 1:count_specialist) {
  specialistecialistatest = specialist.test[, i]
  specialistecialistatrain = specialist.train[, i]
  specialistecilistavalidation = specialist.validation[, i]
  df.weight.test = data.frame(specialistecialistatest, covariable.test)
  df.weight.train = data.frame(specialistecialistatrain, covariable.train)
  df.weight.validation = data.frame(specialistecilistavalidation, covariable.validation)
  
  # Prepare data matrices for LASSO regression
  x_train_l = model.matrix(~.-1, df.weight.train[, -c(1)])
  y_train = df.weight.train$specialistecialistatrain
  x_test = model.matrix(~.-1, df.weight.test[, -c(1)])
  y_test = df.weight.test$specialistecialistatest
  x_validation = model.matrix(~.-1, df.weight.validation[, -c(1)])
  y_validation = df.weight.validation$specialistecilistavalidation
  
  # Train LASSO regression model
  model = cv.glmnet(x_train_l, y_train, alpha = 1)
  # Predict using the trained model
  pred.specialist <- predict(model, s = model$lambda.min, newx = x_validation)
  
  # Compute MSE for each specialist
  specialist.mse.weight.lasso[i] <- sum((pred.specialist - specialist.validation[, i])^2) / length(validation$y)
}

# Compute weighted specialist predictions
specialist.train.weight = specialist.train 
specialist.test.weight = specialist.test 
specialist.validation.weight = specialist.validation

for (j in 1:count_specialist) {
  specialist.train.weight[, j] <- (specialist.train[, j] * (1 / specialist.mse.weight.lasso[j]))
  specialist.test.weight[, j] <- (specialist.test[, j] * (1 / specialist.mse.weight.lasso[j]))
}

# Compute final weighted specialist model
weight.specialist.lasso = t(specialist.mse.weight.lasso) 
specialist.train.weight = data.frame(specialist.train.weight) 
specialist.test.weight = data.frame(specialist.test.weight)

attach(specialist.train.weight) 
attach(specialist.test.weight)

# Compute row sums for training and test specialist data
specialist.train.weight.soma = rowSums(specialist.train.weight) 
specialist.test.weight.soma = rowSums(specialist.test.weight) 
specialist.mse.weight = data.frame(specialist.mse.weight.lasso)

denominador = colSums(1 / specialist.mse.weight)

mean.train.weight_notas <- specialist.train.weight.soma / denominador
mean.test.weight_notas <- specialist.test.weight.soma / denominador

df.mean.train.weight_notas <-data.frame(mean.train.weight_notas,covariable.train)
df.mean.test.weight_notas <-data.frame(mean.test.weight_notas,covariable.test) 
x_train_weight_lasso =model.matrix(~.-1, df.mean.train.weight_notas[,-c(1)]) 
y_train_weight_lasso= df.mean.train.weight_notas$mean.train.weight_notas
x_test_weight_lasso = model.matrix(~.-1,df.mean.test.weight_notas[,-c(1)]) 
y_test_weight_lasso = df.mean.test.weight_notas$mean.test.weight_notas 
dim(x_train_weight_lasso)


# Create final dataset for model training
fit10_mean.weight_lasso <- cv.glmnet(x_train_weight_lasso, y_train_weight_lasso, alpha = 1)
predict_mean_weight_lasso <- predict(fit10_mean.weight_lasso, s = fit10_mean.weight_lasso$lambda.min, newx = x_test_weight_lasso)

# Compute MSE for weighted LASSO regression model
mse.lasso.weight = sum((predict_mean_weight_lasso - test$y)^2) / length(test$y) 
# Compute standard deviation of squared errors
w10 = (predict_mean_weight_lasso - test$y)^2 
w10_deviation = sd(w10)
# Compute standard error
se.lasso_weight = w10_deviation / sqrt(length(test$y))

# Train a Linear Regression model
modelo1 <- lm(y ~ ., data = train)
# Display model summary
summary(modelo1)
# Predict using the trained model
pred.mod1 = predict(modelo1, newdata = test)

# Compute MSE for the linear regression model
mse.mod1 = sum((pred.mod1 - test$y)^2) / length(test$y)
# Compute squared errors
w1 = (pred.mod1 - test$y)^2
# Compute standard deviation of squared errors
w1_deviation = sd(w1)
# Compute standard error
se.mod1 = w1_deviation / sqrt(length(test$y))

######################### MODEL WITH MEAN SPECIALIST ##############################

# Compute mean predictions from specialists
mean.specialist = rowMeans(specialist)
# Create a dataset combining mean specialist predictions and covariates
banco.specialist = data.frame(mean.specialist, covariable)
# Split the dataset into training and testing
train.specialist <- banco.specialist[indicesTraining, ]
test.specialist <- banco.specialist[indicesTest, ]

# Train a Linear Regression model using mean specialist predictions
modelo2 <- lm(mean.specialist ~ ., data = train.specialist)
# Predict using the trained model
pred.mod2 = predict(modelo2, newdata = test.specialist)

# Compute MSE for the mean specialist model
mse.mod2 = sum((pred.mod2 - test$y)^2) / length(test$y)
# Compute squared errors
w2 = (pred.mod2 - test$y)^2
# Compute standard deviation of squared errors
w2_deviation = sd(w2)
# Compute standard error
se.mod2 = w2_deviation / sqrt(length(test$y))

######################### WEIGHTED MODEL #################################

# Add a column of ones for bias term in regression
covariable.alpha <- data.frame(rep(1, count_row), covariable)
covariable.alpha.train <- covariable.alpha[indicesTraining, ]
covariable.alpha.test <- covariable.alpha[indicesTest, ]
covariable.alpha.validation <- covariable.alpha[indicesValidation, ]
specialist.train <- specialist[indicesTraining, ]
specialist.test <- specialist[indicesTest, ]
specialist.validation <- specialist[indicesValidation, ]

# Initialize coefficient matrix
coeficientes = data.frame(coef(modelo2))
# Replace NA values with 0 in coefficients
coeficientes[is.na(coeficientes[,1]), 1] <- 0

lambda = NULL 
lambdas = NULL 
iteracao = 20
for(i in 1:iteracao){
  #PASSO E
  passo.ezao =  (as.matrix(covariable.alpha.train) %*% as.matrix((coeficientes[,1])))
  #PASSO M
  variabilidade = matrix (NA, ncol = count_specialist, nrow = count_row*percentual_train)
  for (k in 1:count_specialist) {
    for (j in 1:count_row*percentual_train){
      variabilidade[j,k] = ((specialist.train[j,k]) - passo.ezao[j])^2
    }
  }
  sum.variabilidade = colSums(variabilidade)
  reverse_lambda = sum.variabilidade / (count_row * percentual_train)
  lambda = 1 / reverse_lambda
  
  # Ridge Regression Regularization
  lambda_ridge <- 10000
  XtX <- t(as.matrix(covariable.alpha.train)) %*% as.matrix(covariable.alpha.train)
  XtX_regularizado <- XtX + diag(lambda_ridge, ncol(XtX))
  parte1w <- solve(XtX_regularizado)
  parte3w = (as.matrix(specialist.train) %*% lambda) / sum(lambda)
  parte4w = as.matrix(t(covariable.alpha.train)) %*% parte3w
  w = parte1w %*% parte4w
  coeficientes = w
}

lambdas = matrix(c(rnorm(count_specialist*total_mse,0,1)), ncol = count_specialist)
lambdas = t(lambda)

######################### MODEL COMPARISON #################################

predict.com.em = as.matrix(covariable.alpha.test) %*% as.matrix(coeficientes[,1])
mse.mod3 = sum((predict.com.em - test$y)^2) / length(test$y)
w3 = (predict.com.em - test$y)^2
w3_deviation = sd(w3)
se.mod3 = w3_deviation / sqrt(length(test$y))

# Model names
MODELS = c(
  "RAYKAR",
  "LINEAR REGRESSION WITH WEIGHTED MEAN",
  "FOREST WITH WEIGHTED MEAN",
  "TREE WITH WEIGHTED MEAN",
  "LASSO WITH WEIGHTED MEAN",
  "LINEAR REGRESSION WITH ARITHMETIC MEAN",
  "FOREST WITH ARITHMETIC MEAN",
  "TREE WITH ARITHMETIC MEAN",
  "LASSO WITH ARITHMETIC MEAN",
  "LINEAR REGRESSION WITH TRUE Y",
  "FOREST WITH TRUE Y",
  "TREE WITH TRUE Y",
  "LASSO WITH TRUE Y"
)

# MSE and standard errors for each model
mean.MSE = round(c(
  mse.mod3,
  mse.lm.weight,
  mse_forest_weight,
  mse.tree.weight,
  mse.lasso.weight,
  mse.lm.mean,
  mse_forest_mean_specialist,
  mse.tree.mean,
  mse_lasso_mean_specialist,
  mse_lm_gold,
  mse_forest_gold,
  mse.tree.gold,
  mse_lasso_gold
), 4)

se = round(c(
  se.mod3,
  se.lm.weight,
  se_forest_weight,
  se.tree.weight,
  se.lasso_weight,
  se_lm_mean,
  se_forest_mean_specialist,
  se_tree_mean,
  se_lasso_mean,
  se_lm_gold,
  se_forest_gold,
  se_tree_gold,
  se_lasso_gold
), 4)

# Create results dataframe
result = data.frame(MODELS, mean.MSE, se)
result

# Sort results by MSE
result_1 = result %>% arrange((mean.MSE))
formattable(result_1)



lamb_mean = matrix(data = NA, nrow = 1, ncol = count_specialist)
for(count_specialist in 1:count_specialist){
  lamb_mean[1,count_specialist] = mean(abs(lambdas[,count_specialist]))
}
lamb_mean

weights_wear <- data.frame(
  weight1 = c(lamb_mean[1,1]),
  weight2 = c(lamb_mean[1,2]),
  weight3 = c(lamb_mean[1,3]),
  weight4 = c(lamb_mean[1,4])
)

round(weights_wear,4)