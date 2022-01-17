

Data <- cbind(usedTestGroup, ktspStatsTestRes$statistics)


Data <- as.data.frame(Data)
colnames(Data) <- c("usedTestGroup", "Stat")
roc(usedTestGroup~Stat, as.data.frame(Data))
levels(Data$usedTestGroup) <- c(1,2)

ROCStrap <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  auc <- roc(formula, data=d)
  return(auc$auc)
}


##
set.seed(333)
bootobject <- boot(data= Data, statistic= ROCStrap, R= 2000, formula = usedTestGroup ~ Stat) 

boot.ci(bootobject, type="norm")

X <- bootobject$t


hist(X)

#mean(X)+ 1.96*(sd(X)/sqrt(2000))

#mean(X)- 1.96*(sd(X)/sqrt(2000))

quantile(X, c(0.025, 0.975))

qnorm(p=c(.025, .975), 
             mean=mean(X), 
             sd=sd(X))



#############
Data2 <- cbind(usedTestGroup, ktspStatsTestUnRes$statistics)


Data2 <- as.data.frame(Data2)
colnames(Data2) <- c("usedTestGroup", "Stat")
roc(usedTestGroup~Stat, as.data.frame(Data2))
levels(Data2$usedTestGroup) <- c(1,2)


ROCStrap <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  auc <- roc(formula, data=d)
  return(auc$auc)
}


##
set.seed(333)
bootobject2 <- boot(data= Data2, statistic= ROCStrap, R= 2000, formula = usedTestGroup ~ Stat) 

boot.ci(bootobject2, type="norm")

Y <- bootobject2$t

hist(Y)

Z <- X-Y
hist(Z)
abline(v = mean(Z))

bootobject3 <- bootobject2

bootobject3$t <- Z

boot.ci(bootobject3, type="norm")


#mean(Z)+ 1.96*(sd(Z)/sqrt(2000))

#mean(Z)- 1.96*(sd(Z)/sqrt(2000))

#library(MKmisc)
#quantileCI(Z)

quantile(Z, c(0.025, 0.975))

ggqqplot(Z)

ggdensity(Z)




Mechanistic <- data.frame(Value = X)
Agnostic <- data.frame(Value = Y)

Mechanistic$modelType <- "Mech"
Agnostic$modelType <- "Agnostic"

ModelCompare <- rbind(Mechanistic, Agnostic)

ggplot(ModelCompare, aes(Value, fill = modelType)) + geom_density(alpha = 0.5)

ggplot(ModelCompare, aes(Value, fill = modelType)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
