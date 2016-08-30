#load dependencies
library(lme4)
library(arm)

lmm.data  <- read.csv("/home/alban/Documents/internship/R/mle2.csv",                      
                      header=TRUE, sep=",", dec = ".")
head(lmm.data)


#MLexamp <- glm(amp ~ beta + gamma, data = lmm.data)
#display(MLexamp)
#AIC(MLexamp)

anova(lm(amp ~ beta + gamma, lmm.data))
#slicing to areas by ROI
#newdata <- mydata[ which(mydata$gender=='F' 
#                         & mydata$age > 65), ]
LOdata <- lmm.data[ which(lmm.data$roi=='LO'),]
head(LOdata)
anova(lm(amp ~ beta + gamma, LOdata))

#modi <- lmer(size ~ condition + (condition|subject), data = lmm.data)
#display(modi)
#modi <- lmer(f ~ i + (i|s), data=d)
#summary(modi)
#anova(modi, type = 3, test = "F")