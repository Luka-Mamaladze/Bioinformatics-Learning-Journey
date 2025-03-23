#### Libraries ####
library(rpart)
library(rpart.plot)
library(caret)
library(randomForest)


#### Question 1 ####
setwd("C:/Users/luka_/Downloads/Machine Learning")
covid_data <- read.csv("covid_19_lab_data_2025.csv")
covid_data$deceased <- as.factor(covid_data$deceased)

set.seed(100) # Set seed for reproductive results

# Create training split at 70/30
trainSet <- sample(nrow(covid_data), 0.7*nrow(covid_data), replace = FALSE)

# Assign test and training data
train_data <- covid_data[trainSet, ]
test_data <- covid_data[-trainSet, ]

# Create decision tree from training data
model_dt <- rpart(deceased ~ ., data = train_data, method = "class")

# Plot decision tree
rpart.plot(model_dt)

# Assign vector for variable importance
importance <- summary(model_dt)$variable.importance
# Plot the variable importance with title
barplot(importance, main = "Variable Importance for the Decision Tree")

# Use the trained decision tree model to predict the class labels for the test data
predictions_dt <- predict(model_dt, newdata = test_data, type = "class")
# Confusion matrix for the decision tree
cm_dt <- confusionMatrix(predictions_dt, test_data$deceased, positive="y")
print(cm_dt) # Print the confusion matrix


#### Question 2 ####
# Impute missing data
for(i in 1:ncol(covid_data)) {
  covid_data[,i][is.na(covid_data[, i])] <- mean(covid_data[, i], na.rm = TRUE)
}

# Recreate the training data with the imputed data
trainSet <- sample(nrow(covid_data), 0.7*nrow(covid_data), replace = FALSE)
train_data <- covid_data[trainSet, ]
test_data <- covid_data[-trainSet, ]

# Random Forest model
model_rf <- randomForest(deceased ~ ., data = train_data, ntree = 100, importance = T)

# Generate Predictions
predictions_rf <- predict(model_rf, newdata = test_data, type = "class")

# Create Confusion Matrix
cm_rf <- confusionMatrix(predictions_rf, test_data$deceased, positive = "y")
print(cm_rf) # Print confusion matrix

#### Question 3 ####

accuracy_values <- vector() # Empty vector to hold accuracies

# Loop for accuracy with varying mtry
for (mtry_value in 1:10) {
  model_rf_mtry <- randomForest(deceased ~ ., data = train_data, ntree = 50, mtry = mtry_value)
  predictions_rf_mtry <- predict(model_rf_mtry, newdata = test_data)
  cm_rf_mtry <- confusionMatrix(predictions_rf_mtry, test_data$deceased, positive = "y")
  accuracy_values <- c(accuracy_values, cm_rf_mtry$overall[[1]])
} 

# Plot accuracy against mtry
plot(1:10, accuracy_values, 
     xlab = "mtry", 
     ylab = "Accuracy", 
     main = "Accuracy vs. mtry")

#### Question 4 ####

# Retrain data with a 50/50 split
trainSet <- sample(nrow(covid_data), 0.5*nrow(covid_data), replace = FALSE)
train_data <- covid_data[trainSet, ]
test_data <- covid_data[-trainSet, ]

# Create random forest 
model_rf_50 <- randomForest(deceased ~ ., data = train_data, ntree = 100, importance = T)

# Generate Predictions
predictions_rf_50 <- predict(model_rf_50, newdata = test_data, type = "class")

# Create Confusion Matrix
cm_rf <- confusionMatrix(predictions_rf_50, test_data$deceased, positive = "y")
print(cm_rf)
