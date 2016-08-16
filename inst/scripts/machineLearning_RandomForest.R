library(E2Predictor)
library(randomForest)
library(stringr)
library(readr)
library(devtools)

#relative importance function
source_url('https://gist.github.com/fawda123/6206737/raw/2e1bc9cbc48d1a56d2a79dd1d33f414213f5f1b1/gar_fun.r')
#plotting function
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')


### Configuration #########################

initConfiguration(netMHCpath = "/Users/napedro/external_sources/netMHC-4.0/netMHC",
                  working.path = "/Users/napedro/CloudStation")

cell.lines <- names(E2Predictor.Config$hla_alleles)
conditions <- c("Control", "Bortezomib", "DMOG", "IFNg", "Rapamycin")
cell.line <- cell.lines[1]

is_perturbation <- F
condition <- conditions[5]
addgrepcond <- ""  # " II "

percentage_training <- 0.6
hidden.layers <- c(4, 2)

set.seed(1234567890)
###########################################



experim.folder <- ifelse(is_perturbation, "perturbations", "cell.lines")
modeling.folder <- file.path(E2Predictor.Config$working.path, "tables.for.modeling", experim.folder, cell.line)
RF.folder <- file.path(E2Predictor.Config$working.path, "data_analyses", "machine_learning", "RandomForest", experim.folder, cell.line)

if(is_perturbation){
    modeling.folder <- file.path(modeling.folder, paste(condition, addgrepcond, sep = "_") )
    RF.folder <- file.path(RF.folder, paste(condition, addgrepcond, sep = "_") )
}
mkdir(RF.folder)


modeling.file <- file.path(modeling.folder, "omics_table_for_modeling.RData")


# load data for modeling
omics <- readRDS(modeling.file)

omics <- omics %>% mutate(has_ligands = ifelse(is.na(ligands), FALSE, TRUE))  # variable to be predicted



omics.model <- omics
#rename variables for the model
renaming.vars <- str_extract(names(omics.model), "\\[GO:.*\\]")
renaming.vars <- gsub("\\[", "", renaming.vars)
renaming.vars <- gsub("\\]", "", renaming.vars)
renaming.vars <- gsub(":", "", renaming.vars)

renaming.vars.index <- which(!is.na(renaming.vars))
renaming.vars <-renaming.vars[renaming.vars.index]
names(omics.model)[renaming.vars.index] <- renaming.vars

omics.model.variables <- c("has_ligands", "PROTEIN_QUANTITY", "RPKM.avg", grep("^GO.*", names(omics.model), value = T))

omics.model <- omics.model[, omics.model.variables]

# Normalize variables
omics.model[is.na(omics.model)] <- 0.0
maxs <- apply(omics.model, 2, max)
mins <- apply(omics.model, 2, min)

num.vars <- ncol(omics.model) - 1

omics.model.sc <- as.data.frame(scale(omics.model, center = mins, scale = maxs - mins))

# Sampling for training
trainset.index <- sample(1:nrow(omics.model.sc), size = round(nrow(omics.model.sc) * percentage_training, 0), replace = F )
trainset <- omics.model.sc[trainset.index, ]
testset <- omics.model.sc[-trainset.index, ]

omics.model.sc.ligands <- omics.model.sc %>% filter(has_ligands == T)
omics.model.sc.noligands <- omics.model.sc %>% filter(has_ligands == F)
trainset.index.ligands <- sample(1:nrow(omics.model.sc.ligands), size = round(nrow(omics.model.sc.ligands) * percentage_training,0), replace = F)
trainset.index.noligands <- sample(1:nrow(omics.model.sc.noligands), size=length(trainset.index.ligands))
trainset <- rbind(omics.model.sc.ligands[trainset.index.ligands, ], omics.model.sc.noligands[trainset.index.noligands, ])

# Train the neural network
n <- names(trainset)
f <- as.formula(paste("has_ligands ~", paste(n[!n %in% "has_ligands"], collapse = " + ")))
#nn <- neuralnet(f,data=trainset,hidden=hidden.layers,linear.output=F)
nn <- neuralnet(f,data=trainset,hidden=hidden.layers,linear.output=F, stepmax = 1e+06)

#save the neural network
nn_name <- paste(omics.model.variables, sep="_", collapse = "_")
nn_name <- paste("HL", hidden.layers, sep="_", collapse ="_")

saveRDS(nn, file=file.path(RF.folder, paste("ANN.trained", nn_name, "Rds", sep = ".") ))



#Predicting values with the neural network
pr.nn <- compute(nn, testset[, 2:ncol(testset)])

pr.nn_ <- pr.nn$net.result*(max(omics.model$has_ligands)-min(omics.model$has_ligands))+min(omics.model$has_ligands)
test.r <- (testset$has_ligands)*(max(omics.model$has_ligands)-min(omics.model$has_ligands))+min(omics.model$has_ligands)

MSE.nn <- sum((test.r - pr.nn_)^2)/nrow(testset)




# Plot the neural network
pdf(file.path(RF.folder,  paste("NeuralNetwork_map", nn_name, "pdf", sep = ".") ), width = 10, height = 10)

#relative importance of input variables for Y1
rel.imp<-gar.fun('Y1',nn,bar.plot=F)$rel.imp

#color vector based on relative importance of input values
cols<-colorRampPalette(c('green','red'))(num.vars)[rank(rel.imp)]

#plot model with new color vector
#separate colors for input vectors using a list for 'circle.col'
plot.nnet(nn, circle.col = list(cols, 'lightblue') )
dev.off()


## Compare with a linear model fit
lm.fit <- glm(f, data = trainset)
pr.lm <- predict(lm.fit, testset)
MSE.lm <- sum((pr.lm - testset$has_ligands)^2)/nrow(testset)



contingencyTable <- function(predicted, test_var, threshold = 0.6){

    predicted_boolean <- ifelse(predicted > threshold, TRUE, FALSE)

    ct <- data.frame(
        TruePos = (test_var == predicted_boolean & test_var == TRUE),
        TrueNeg = (test_var == predicted_boolean & test_var == FALSE),
        FalsePos = (!(test_var == predicted_boolean) & test_var == FALSE),
        FalseNeg = (!(test_var == predicted_boolean) & test_var == TRUE)
    ) %>% summarise_all(sum)

    return(ct)
}

ct.lm <- contingencyTable(pr.lm, testset$has_ligands, 0.6)
ct.nn <- contingencyTable(pr.nn_, testset$has_ligands, 0.6)

colors <- ifelse(testset$has_ligands == 0, "black", "red")
plot(pr.nn_,pch=19, col=colors)
plot(density(pr.nn_[testset$has_ligands==F]))
lines(density(pr.nn_[testset$has_ligands==T]), col="red")

plot(density(pr.lm[testset$has_ligands==F]))
lines(density(pr.lm[testset$has_ligands==T]), col="red")


readr::write_tsv(ct.nn, path = file.path(RF.folder,  paste("contingencyTable", nn_name, "tsv", sep = ".") ))


neuralnet::print.nn(nn)
