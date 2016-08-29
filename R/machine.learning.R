
#' prelim.ML
#' @description makes some preliminary steps useful for Machine Learning (A PCA and a correlation analysis)
#'
#' @param df data frame for analysis
#' @param ncp number of dimensions kept in the results (by default 5) of PCA analysis
#' @param dim.reduction.factor Correlation threshold over which correlated variables will be dropped
#' @param saveInFolder when not NULL, folder path where result files will be written
#' @param suffix when saveInFolder is not NULL, it adds a suffix to file names
#'
#' @return list with max four objects: PCA and correlation of variables, and if variable reduction is possible, PCA and correlation of reduced variables.
#'
#' @export
#'
prelim.ML <- function(df, ncp = 5, dim.reduction.factor = 0.7, saveInFolder = NULL, suffix = ""){

    prelim.ML.obj <- list()

    pca.summary.file <- ifelse(is.null(saveInFolder),
                               "",
                               file.path(saveInFolder, paste("PCA_summary", suffix , "txt", sep = ".") ) )

    pca.red.summary.file <- ifelse(is.null(saveInFolder),
                                   "",
                                   file.path(saveInFolder, paste("PCA_summary_reduced", suffix , "txt", sep = ".") ) )


    # Only numeric columns can be used in PCA and correlation
    df.num <- df[,  sapply(df, is.numeric)]

    # Do a Principal Component Analysis
    pca.df <- FactoMineR::PCA(df.num, scale.unit=TRUE, ncp=ncp, graph = F)
    prelim.ML.obj$pca <- pca.df
    summary(pca.df, nbelements = Inf , nbind = Inf, file = pca.summary.file)

    #dimdesc(pca.df)

    if(!is.null(saveInFolder)){
        pdf(file.path(saveInFolder, paste("PCA_graphs", suffix, "pdf", sep=".")), width = 10, height = 10)
        devNum = dev.cur()
    }
    plot(pca.df, cex=0.6, choix = "ind", unselect = 0.7, label="ind.sup" )
    plot(pca.df, cex=0.6,  choix = "var")
    if(!is.null(saveInFolder)){
        dev.off(devNum)
    }

    #Calculate variable correlations
    corr.df <- cor(df.num)
    prelim.ML.obj$correlation <- corr.df
    tl.cex <- ifelse(ncol(df.num) > 30, 0.6, 1.0)
    if(!is.null(saveInFolder)){
        pdf(file.path(saveInFolder, paste("Correlation_graphs", suffix, "pdf", sep=".")), width = 10, height = 10)
        devNum = dev.cur()
    }
    corrplot::corrplot(corr.df, order="hclust", tl.cex = tl.cex, method = "ellipse", diag = F, tl.col = "darkblue")
    if(!is.null(saveInFolder)){
        dev.off(devNum)
    }

    # Perform a similar analysis with reduced number of variables
    highlyCor <- caret::findCorrelation(corr.df, dim.reduction.factor)

    if(length(highlyCor) > 0){
        prelim.ML.obj$highly.correlated.vars <- names(df.num[, highlyCor])
        df.red <- df.num[,-highlyCor]
        corr.df.red <- cor(df.red)
        prelim.ML.obj$correlation.reduced <- corr.df.red

        if(!is.null(saveInFolder)){
            pdf(file.path(saveInFolder, paste("Correlation_Reduced_graphs", suffix, "pdf", sep=".")), width = 10, height = 10)
            devNum = dev.cur()
        }
        tl.cex.red <- ifelse(ncol(df.red) > 30, 0.6, 1.0)
        corrplot::corrplot(corr.df.red, order = "hclust", tl.cex = tl.cex.red, method = "ellipse", diag = F, tl.col = "darkblue")
        if(!is.null(saveInFolder)){
            dev.off(devNum)
        }

        pca.df.red <- FactoMineR::PCA(df.red, scale.unit = T, ncp=5, graph = F)
        prelim.ML.obj$pca.reduced <- pca.df.red
        summary(pca.df.red, nbelements = Inf , nbind = Inf, file = pca.red.summary.file)

        if(!is.null(saveInFolder)){
            pdf(file.path(saveInFolder, paste("PCA_graphs_red", suffix, "pdf", sep=".")), width = 10, height = 10)
            devNum = dev.cur()
        }

        plot(pca.df.red, cex=0.6, choix = "ind", unselect = 0.7, label="ind.sup" )
        plot(pca.df.red, cex=0.6,  choix = "var")

        if(!is.null(saveInFolder)){
            dev.off(devNum)
        }
    }

    return(prelim.ML.obj)
}



#' get.training.sample
#'
#' @description separates a data frame into two data frames for training and testing in ML classification.
#'
#' @param df data frame containing a whole data set.
#' @param class_category boolean column name used for classification.
#' @param percentage_training percentage of data used for training. It takes the portion percentage_training from the smallest class and a equivalent number of elements from the biggest class.
#'
#' @return list with two data frames: trainset and testset
#'
#' @export
#'
get_training_sample <- function(df, class_category, percentage_training){

    # Sampling for training
    df[, class_category] = as.logical(df[, class_category])

    df.class.true <- df[df[,class_category] == T,] # %>% filter_(class_category == T)
    df.class.false<- df[df[,class_category] == F,] # %>% filter_(class_category == F)

    #set the size of the training depending on the smallest dataset
    trainset.size <- round(min(nrow(df.class.true), nrow(df.class.false)) * percentage_training, 0)

    trainset.index.true  <- sample(1:nrow(df.class.true),  size = trainset.size, replace = F)
    trainset.index.false <- sample(1:nrow(df.class.false), size = trainset.size, replace = F)

    trainset <- rbind(df.class.true[ trainset.index.true, ], df.class.false[ trainset.index.false, ])
    testset <-  rbind(df.class.true[-trainset.index.true, ], df.class.false[-trainset.index.false,])

    train.test.sets <- list(trainset = trainset, testset = testset)

    return(train.test.sets)
}

#' contingencyTable
#' @description calculates the contingency table, given test and predicted values, and a threshold for the predicted values.
#'
#' @param predicted vector with predicted values (if it is a score, a threshold should be given)
#' @param test_var vector with classification
#'
#' @return data frame with the four values of the contingency table
#'
#' @export
#'
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


#WARNING: not finished yet!!
#' retrieve.data.modeling
#'
#' @description retrieves data for modeling from a defined data structure
#'
#' @param cell.line cell line from which data should be retrieved
#' @param is_perturbation
retrieve.data.modeling <- function(cell.line, is_perturbation, ML_method){

    # Define file and folder names
    experim.folder <- ifelse(is_perturbation, "perturbations", "cell.lines")
    modeling.folder <- file.path(E2Predictor.Config$working.path, "tables.for.modeling", experim.folder, cell.line)
    ML.folder <- file.path(E2Predictor.Config$working.path, "data_analyses", "machine_learning", ML_method, experim.folder, cell.line)

    if(is_perturbation){
        modeling.folder <- file.path(modeling.folder, paste(condition, addgrepcond, sep = "_") )
        ML.folder <- file.path(ML.folder, paste(condition, addgrepcond, sep = "_") )
    }
    mkdir(ML.folder)

    modeling.file <- file.path(modeling.folder, "omics_table_for_modeling.RData")

    # load data for modeling
    return(readRDS(modeling.file))
}


pacompilar <- function(){

library(E2Predictor)
library(neuralnet)
library(stringr)
library(readr)
library(devtools)
library(pROC)

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
hidden.layers <- c(5, 2)

set.seed(1234567890)
###########################################

run.ML <- function(cell.line, is_perturbation = F, ML_method = "ANN", ML_parameters = NULL, suffix = NULL){

    ML_methods = c(neural.network = "ANN", random.forest = "RForest" )

    ML_parameters_example = list(hidden.layers = c(c(5,2)), percentage_training = 0.6)

    #Sanity checks
    if(!ML_method %in% ML_methods){
        ML_methods_str = paste( ML_methods, collapse = ", ")
        stop.message = paste("Wrong selection of machine learning method. ML_method should be one of the following strings: ", ML_methods_str)
        stop(stop.message)
    }
    if(ML_method == "ANN"){
        if( !"hidden.layers" %in% names(ML_parameters) ){
            stop("When ANN is selected, a hidden.layers vector variable should be passed within the ML_parameters list")
        }
    }

    # Define file and folder names
    experim.folder <- ifelse(is_perturbation, "perturbations", "cell.lines")
    modeling.folder <- file.path(E2Predictor.Config$working.path, "tables.for.modeling", experim.folder, cell.line)
    ML.folder <- file.path(E2Predictor.Config$working.path, "data_analyses", "machine_learning", ML_method, experim.folder, cell.line)

    if(is_perturbation){
        modeling.folder <- file.path(modeling.folder, paste(condition, addgrepcond, sep = "_") )
        ML.folder <- file.path(ML.folder, paste(condition, addgrepcond, sep = "_") )
    }
    mkdir(ML.folder)

    modeling.file <- file.path(modeling.folder, "omics_table_for_modeling.RData")


    # load data for modeling
    #WARNING: take ligands that have an associated allele
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

    omics.GO <- omics[, grep("\\[GO:.*\\]", names(omics))]
    omics.GO.pos <- omics.GO[,1:15]
    omics.GO.neg <- omics.GO[,16:30]

    #omics.model$GO.positive <- rowSums(omics.GO.pos)
    #omics.model$GO.negative <- rowSums(omics.GO.neg)

    omics.model.variables <- c("has_ligands", "PROTEIN_QUANTITY", "RPKM.avg",  grep("^GO.*", names(omics.model), value = T)) #
    #omics.model.variables <- c("has_ligands", "PROTEIN_QUANTITY", "RPKM.avg", "GO.positive", "GO.negative") #

    omics.model <- omics.model[, omics.model.variables]

    # Normalize variables
    omics.model[is.na(omics.model)] <- 0.0

    maxs <- apply(omics.model, 2, max)
    mins <- apply(omics.model, 2, min)

    num.vars <- ncol(omics.model) - 1

    omics.model.sc <- as.data.frame(scale(omics.model, center = mins, scale = maxs - mins))


}





# Train the neural network
n <- names(trainset)
f <- as.formula(paste("has_ligands ~", paste(n[!n %in% "has_ligands"], collapse = " + ")))
#nn <- neuralnet(f,data=trainset,hidden=hidden.layers,linear.output=F)
nn <- neuralnet(f,data=trainset,hidden=hidden.layers,linear.output=F, stepmax = 1e+06)

#save the neural network
nn_name <- paste(omics.model.variables, sep="_", collapse = "_")
nn_name <- paste(nn_name, "HL", hidden.layers, sep="_", collapse ="_")

saveRDS(nn, file=file.path(ML.folder, paste("ANN.trained.GOandprotein", nn_name, "Rds", sep = ".") ))



#Predicting values with the neural network
pr.nn <- compute(nn, testset[, 2:ncol(testset)])

pr.nn_ <- pr.nn$net.result*(max(omics.model$has_ligands)-min(omics.model$has_ligands))+min(omics.model$has_ligands)
test.r <- (testset$has_ligands)*(max(omics.model$has_ligands)-min(omics.model$has_ligands))+min(omics.model$has_ligands)

MSE.nn <- sum((test.r - pr.nn_)^2)/nrow(testset)




# Plot the neural network
pdf(file.path(ML.folder,  paste("NeuralNetwork_map_GOandprotein", nn_name, "pdf", sep = ".") ), width = 10, height = 10)

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


ct.lm <- contingencyTable(pr.lm, testset$has_ligands, 0.6)
ct.nn <- contingencyTable(pr.nn_, testset$has_ligands, 0.6)

readr::write_tsv(ct.nn, path = file.path(ML.folder,  paste("contingencyTable_GOandprotein", nn_name, "tsv", sep = ".") ))


neuralnet::print.nn(nn)

nn <- readRDS("/Users/napedro/CloudStation/data_analyses/machine_learning/ANN/cell.lines/JY/ANN.trained.GO.HL_5_HL_2.Rds")
pr.alldataset.nn <- compute(nn, omics.model.sc[, 2:ncol(omics.model.sc)])
pr.alldataset.nn_ <- pr.alldataset.nn$net.result*(max(omics.model$has_ligands)-min(omics.model$has_ligands))+min(omics.model$has_ligands)

colors <- ifelse(omics.model.sc$has_ligands == 0, "black", "red")
plot(pr.alldataset.nn_,pch=19, col=colors)
plot(density(pr.alldataset.nn_[omics.model$has_ligands==F]))
lines(density(pr.alldataset.nn_[omics.model$has_ligands==T]), col="red")

}


#' prepare_data_modeling
#' @description fix some variable names (from gene ontology), filter variables, and normalize data
#'
#' @param df data frame containing data to be modeled
#' @param useVariables a vector containing variable names to be used. When NULL, all variables are used.
#'
#' @export
#'
prepare_data_modeling <- function(df, useVariables=NULL){

    omics.model <- df %>% distinct(protein_entry, .keep_all=TRUE)
    omics.rownames <- omics.model$protein_entry

    if(!is.null(useVariables)) omics.model <- omics.model[, useVariables]

    #rename variables for the model
    renaming.vars <- str_extract(names(omics.model), "\\[GO:.*\\]")
    renaming.vars <- gsub("\\[", "", renaming.vars)
    renaming.vars <- gsub("\\]", "", renaming.vars)
    renaming.vars <- gsub(":", "", renaming.vars)

    renaming.vars.index <- which(!is.na(renaming.vars))
    renaming.vars <-renaming.vars[renaming.vars.index]
    names(omics.model)[renaming.vars.index] <- renaming.vars


    # Normalize variables
    omics.model[is.na(omics.model)] <- 0.0

    maxs <- apply(omics.model, 2, max)
    mins <- apply(omics.model, 2, min)

    num.vars <- ncol(omics.model) - 1

    omics.model.sc <- as.data.frame(scale(omics.model, center = mins, scale = maxs - mins))
    rownames(omics.model.sc) <- omics.rownames

    return(omics.model.sc)
}

