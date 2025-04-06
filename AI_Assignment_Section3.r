#### Setup ####
if(T){
  library(readxl)
  library(rpart)
  library(rpart.plot)
  library(caret)
  library(randomForest)
}

df <- read_excel("Dataset_3.xlsx", sheet = 1) 
df$diagnosis <- as.factor(df$diagnosis)

train_data <- subset(df, set == "train")
test_data <- subset(df, set == "test")
train_data <- train_data[,-c(1,2)]
test_data <- test_data[,-c(1,2)]

#### DT ####

set.seed(100)
model <- rpart(diagnosis ~ ., data = train_data, method = "class")

rpart.plot(model)

predictions_dt <- predict(model, newdata = test_data, type = "class")

cm_dt <- confusionMatrix(predictions_dt, test_data$diagnosis)
print(cm_dt) 

mean(c(cm_dt$byClass["Sensitivity"], cm_dt$byClass["Specificity"]))

#### Forest ####
set.seed(100)
model_rf <- randomForest(diagnosis ~ ., data = train_data, ntree = 100, importance = T)

predictions_rf <- predict(model_rf, newdata = test_data, type = "class")

cm_rf <- confusionMatrix(predictions_rf, test_data$diagnosis)
print(cm_rf) 

mean(c(cm_rf$byClass["Sensitivity"], cm_rf$byClass["Specificity"]))

#### Q3 ####

importance <- summary(model)$variable.importance

barplot(importance, main = "Variable Importance for the Decision Tree")

# View importance
importance(model_rf)

# Optionally plot
varImpPlot(model_rf)

importance_data <- as.data.frame(importance(model_rf))
importance_data <- importance_data[order(-importance_data$MeanDecreaseGini),]
rownames(head(importance_data, 6)) %in% names(importance[1:6])

cbind("Decision Tree (X)" = sort(names(importance[1:6])), 
      "Random Forest (Y)" = sort(rownames(head(importance_data, 6))),
      "X in Y?" = rownames(head(importance_data, 6)) %in% names(importance[1:6]))
