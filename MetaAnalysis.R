#### Libraries ####
library(survival)
library(coxme)
library(lme4)
library(car)
library(MASS)
library(metafor)
library(meta)
library(ggplot2)

#### Setup ####
setwd("D:/Assignments/CB2/Meta Analysis")
setwd("C:/Users/luka_/Downloads/Meta Analysis")
data <- read.csv("C:/Users/luka_/Downloads/Meta Analysis/Meta_analysis_dataset_v2.csv")
data <- read.csv("D:/Assignments/CB2/Meta Analysis/Meta_analysis_dataset_v2.csv")
# Check for NA values and remove NA values
any(is.na(data)) # TRUE
data = na.omit(data)
# There was a lowercase F found, all changed to uppercase
table(data$Sex)
data$Sex <- toupper(data$Sex)
attach(data)





# Quick visualisation of data to ensure quality and understand data
table(OS_event)
table(Cohort, Hospital)
table(PBRM1_mutation)
table(BAP1_mutation)
table(PIK3CA_mutation)
table(Advanced_disease)


#### Q1 A ####

# Create random effects model
model_re=glmer(Advanced_disease ~ Age_at_diagnosis +(1|Cohort), 
                         family = binomial,
                         data = data)

# Summarise model
summary(model_re)
confint(model_re)

table(data$Advanced_disease, data$BAP1_mutation)

# Technically optimal model, but it isn't significant
model_optimal=glmer(Advanced_disease ~ 
                      Sex * Age_at_diagnosis + (1|Cohort), 
                   family = binomial,
                   data = data)
summary(model_optimal)
confint(model_optimal)

# Heterogeneity
x<-matrix(0,6,2)
a<-1
while(a<7){
  x[a,1]<-summary(glm(Advanced_disease[which(Cohort == a)] ~ 
                        Age_at_diagnosis[which(Cohort == a)] 
                      * as.factor(Sex[which(Cohort == a)]) ,family=binomial, data=data))$coefficients[2,1]
  x[a,2]<-summary(glm(Advanced_disease[which(Cohort == a)]  ~ 
                        Age_at_diagnosis[which(Cohort == a)] 
                      * as.factor(Sex[which(Cohort == a)]),family=binomial, data=data))$coefficients[2,2]
  a<-a+1
}
rma(x[,1], sei=x[,2], data=x)

##############

surv_model <- coxph(Surv(OS,OS_event) ~ PIK3CA_mutation + frailty(Cohort))
summary(surv_model)

cox.zph(surv_model)
confint(surv_model)

opt_model <- coxph(Surv(OS,OS_event) ~ as.factor(Sex) + Age_at_diagnosis + 
                     Advanced_disease + PBRM1_mutation + frailty(Cohort), data=data)
summary(opt_model)

cox.zph(opt_model)
confint(opt_model)





### test for heterogeneity ###
x<-matrix(0,6,5)
a<-1
while(a<7)
{
  x[a,1]<-summary(coxph(Surv(OS[which(Cohort == a)],OS_event[which(Cohort == a)])~PIK3CA_mutation[which(Cohort == a)],data=data))$coefficients[1]
  x[a,2]<-summary(coxph(Surv(OS[which(Cohort == a)],OS_event[which(Cohort == a)])~PIK3CA_mutation[which(Cohort == a)],data=data))$coefficients[3]
  x[a,3]<-summary(coxph(Surv(OS[which(Cohort == a)],OS_event[which(Cohort == a)])~PIK3CA_mutation[which(Cohort == a)],data=data))$conf.int[1]
  x[a,4]<-summary(coxph(Surv(OS[which(Cohort == a)],OS_event[which(Cohort == a)])~PIK3CA_mutation[which(Cohort == a)],data=data))$conf.int[3]
  x[a,5]<-summary(coxph(Surv(OS[which(Cohort == a)],OS_event[which(Cohort == a)])~PIK3CA_mutation[which(Cohort == a)],data=data))$conf.int[4]
  a<-a+1
}

hetero_test<-metagen(x[,1],x[,2],method.tau.ci="QP")
summary(hetero_test)

### generate the forest plot
colnames(x)<-c("LOG(HR)","SE(LOG(HR))","HR","CI:2.5","CI:97.5")

y<-c(0,0,as.matrix(exp(summary(opt_model)$coefficients)),exp(confint(opt_model)))
z<-rbind(x,y)
row.names(z)<-c("study 1","study 2","study 3","study 4","study 5","study 6","Total")

ggplot(data=as.data.frame(z), aes(y=1:7, x=HR,
                                  xmin=z[,4], 
                                  xmax=z[,5])) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  labs(title='PIK3CA Forest Plot', x='Hazard Ratio', y = 'Cohort') +
  geom_vline(xintercept=1, color='red', linetype='dashed', alpha=.8) +
  scale_y_continuous(breaks=1:nrow(z),labels=row.names(z))
