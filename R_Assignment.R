#### Libraries ####
library(ggplot2)
library(survival)

#### Question 1 ####
# Load in the csv and clean up data frame
data_assessment <- read.csv("Assignment_regression_dataset_v2.csv", row.names = "X")
data_assessment <- data_assessment[,c(1:7)]
data_assessment <- na.omit(data_assessment)

# Plot to make sure that relationship between age and meal.cal is linear
ggplot(data_assessment, aes(x = meal.cal, y = age)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Scatter Plot of Age vs. Meal Caloric Intake",
       x = "Meal Caloric Intake",
       y = "Age")

# Basic visualisation of meal.cal using boxplots and histograms
boxplot(data_assessment$meal.cal, main = "Meal Calories - Boxplot", horizontal = FALSE)
hist(data_assessment$meal.cal, main = "Meal Calories - Histogram", breaks = 20)

# Build the linear model to identify relationship between age and meal.cal
model <- lm(age ~ meal.cal, data = data_assessment)
summary(model)

# Diagnostic plots to ensure model is adhering to assumptions
par(mfrow = c(2, 2))
plot(model)

# Build the optimal linear model to predict age
model <- lm(age ~ meal.cal + Treatment_type + ph.ecog, data = data_assessment)
summary(model)

##### Question 2 ####
# Convert treatment type to a factor (since its categorical data)
data_assessment$Treatment_type <- as.factor(data_assessment$Treatment_type)

# Optional step which can be used to figure out which variables belong in the optimal model
best <- step(model_treatment)
summary(best)

# Build the optimal glm to predict treatment type
model_treatment <- glm(Treatment_type ~ meal.cal + wt.loss + age, 
                       data = data_assessment, 
                       family = binomial)
summary(model_treatment)

#### Question 3 ####
# Create the survival object
surv_obj <- Surv(time = data_assessment$Survival_months, event = data_assessment$status)

# Create formula for kaplan meier plot
km_fit <- survfit(surv_obj ~ Treatment_type, data = data_assessment)

# Plot the Kaplan-Meier curve
par(mfrow = c(1, 1))
plot(km_fit, 
     col = c("blue", "red"), 
     xlab = "Survival Time (Months)", 
     ylab = "Survival Probability",
     main = "Kaplan-Meier Survival Curve by Treatment Type")
legend("topright", legend = c("Chemotherapy", "Radiation"), col = c("blue", "red"), lwd = 2)

# Build the cox model to see if treatment type predicts status 
cox <- coxph(surv_obj ~ Treatment_type, 
                   data = data_assessment)

# Build the optimal cox model to see if treatment type and ph.ecog predicts status 
cox <- coxph(surv_obj ~ Treatment_type + ph.ecog, 
                   data = data_assessment)

# Summaries the results of the cox model
summary(cox)
