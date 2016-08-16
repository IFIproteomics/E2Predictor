library(pROC)


testset$predicted <- pr.nn_
r1 <- roc(has_ligands ~ predicted, testset, plot=TRUE, smooth=F)
print(r1)

testset2$predicted.nn.2 <- pr2.nn_
r2 <- roc(has_ligands ~ predicted.nn.2, testset2, plot=TRUE, smooth=F)
print(r2)


library(randomForest)
forestFit <- randomForest(x=trainset[,2:ncol(trainset)], y=as.factor(trainset[,"has_ligands"]),
                            importance=TRUE, do.trace=100, ntree = 500)
forestPredict <- predict(forestFit, newdata = testset)

#forestFit$importance

testset$predicted.rf <- as.numeric(forestPredict)
r3 <- roc(has_ligands ~ predicted.rf, testset, plot=TRUE, smooth=F, add=TRUE, percent=r1$percent)
print(r3)


library(caret)
postResample(testset$predicted, testset$has_ligands)

library(rpart)
decTreeFit <- rpart(f, data = trainset, method = "anova")
printcp(decTreeFit)
plotcp(decTreeFit)
plot(decTreeFit, uniform=TRUE)
summary(decTreeFit)

decTreePredict <- predict(decTreeFit, testset)
colors <- ifelse(omics.model.sc$has_ligands == 0, "black", "red")
plot(decTreePredict, pch=19, col=colors)
plot(density(decTreePredict[testset$has_ligands==F]))
lines(density(decTreePredict[testset$has_ligands==T]), col="red")
testset$predicted.dt <- as.numeric(decTreePredict)
r4 <- roc(has_ligands ~ predicted.dt, testset, plot=TRUE, smooth=F)
print(r4)
