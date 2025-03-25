#### Question 1 ####
data <- read.delim("asthma.txt")
attach(data)

head(data)

library(doBy)
summaryBy(attack ~ gender, data = data, FUN = c(mean, sd, median))
summaryBy(attack ~ res_inf, data = data, FUN = c(mean, sd, median))
summaryBy(attack ~ ghq12, data = data, FUN = c(mean, sd, median))

hist(attack, main = "Histogram of Attack Distribution", 
     xlab = "Number of attacks", 
     ylab="Count")

library(ggplot2)
ggplot(data, aes(x = gender, y = attack)) +
  geom_boxplot() +
  labs(title = "Asthma Episodes by Gender", 
       x = "Gender", 
       y = "Number of Asthma Episodes")

# Boxplot for res-inf vs. attack
ggplot(data, aes(x = res_inf, y = attack)) +
  geom_boxplot() +
  labs(title = "Asthma Episodes by Respiratory Infection", 
       x = "Respiratory Infection", 
       y = "Number of Asthma Episodes")

# Scatter plot for attack vs. ghq12
ggplot(data, aes(x = ghq12, y = attack)) +
  geom_point() +
  labs(title = "Asthma Episodes vs. GHQ12 Score", 
       x = "GHQ12 Score", 
       y = "Number of Asthma Episodes")

#### Question 2 ####
library(ggplot2)

ggplot(data, aes(ghq12, attack)) +
  geom_point() + geom_smooth()

ggplot(data, aes(ghq12, attack)) + 
  geom_point() + 
  geom_smooth(method=lm) + 
  labs(title = "Linear regression line of GHQ12 and Attack")


poisson_model <- glm(attack ~ ghq12, family = poisson, data = data)

plot(data$ghq12, residuals(poisson_model, type = "pearson"),
     main = "Residuals vs GHQ12",
     xlab = "GHQ12 Score", ylab = "Pearson Residuals",
     pch = 16, col = "blue")
abline(h=0, col="red", lty=2)  

#### Question 3 ####

#### Question 4 ####

library(AER)

# Model 1
model1 <- glm(attack ~ gender, family = "poisson", data = data)
summary(model1)
dispersiontest(model1)
confint(model1)

pchisq(model1$deviance, df = model1$df.residual, lower.tail = F)
qchisq(0.95, df.residual(model1))
deviance(model1)

# Model 2
model2 <- glm(attack ~ ghq12, family = "poisson", data = data)
summary(model2)
dispersiontest(model2)
confint(model2)

pchisq(model2$deviance, df = model2$df.residual, lower.tail = F)
qchisq(0.95, df.residual(model2))
deviance(model2)

# Model 3
model3 <- glm(attack ~ res_inf, family = "poisson", data = data)
summary(model3)
dispersiontest(model3)
confint(model3)

pchisq(model3$deviance, df = model3$df.residual, lower.tail = F)
qchisq(0.95, df.residual(model3))
deviance(model3)


#### Question 5 ####

model4 <- glm(attack ~ gender + ghq12 + res_inf, family = "poisson", data = data)
summary(model4)
dispersiontest(model4)
confint(model4)

pchisq(model4$deviance, df = model4$df.residual, lower.tail = F)
qchisq(0.95, df.residual(model4))
deviance(model4)

#### Question 6 ####

model5 <- glm(attack ~ ghq12 + res_inf, family = "poisson", data = data)
summary(model5)
confint(model5)
dispersiontest(model5)

qchisq(0.95, df.residual(model5))
deviance(model5)

#### Question 7 ####

ggplot(data, aes(x = ghq12, y = attack, color = as.factor(gender))) +
  geom_point() + 
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = T) +  
  labs(title = "Attack vs GHQ12 by Gender",
       x = "GHQ12",
       y = "Attack Count",
       color = "Gender") +
  theme_minimal()

ggplot(data, aes(x=ghq12, y=attack, color=res_inf)) +
  geom_point() + 
  geom_smooth(method = "glm", method.args=list(family=('poisson'))) +
  labs(title = "Attack vs GHQ12 by res_inf",
       x = "GHQ12",
       y = "Attack Count",
       color = "res_inf")

#### Question 8 ####

x <- glm(attack ~ ghq12 + res_inf, family = "poisson", data = data)
y <- glm(attack ~ ghq12 * res_inf, family = "poisson", data = data)

summary(x)
dispersiontest(x)
pchisq(x$deviance, df = x$df.residual, lower.tail = F)
confint(x)
qchisq(0.95, df.residual(x))
deviance(x)

summary(y)
dispersiontest(y)
pchisq(y$deviance, df = y$df.residual, lower.tail = F)
confint(y)
qchisq(0.95, df.residual(y))
deviance(y)

anova(x,y)

t.test(ghq12 ~ gender, data = data)  
chisq.test(table(data$gender, data$res_inf)) 
