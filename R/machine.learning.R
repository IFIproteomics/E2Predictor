
#' prelim.ML
#' @description makes some preliminary steps useful for Machine Learning (A PCA and a correlation analysis)
#'
#' @param df data frame for analysis
#' @param ncp number of dimensions kept in the results (by default 5) of PCA analysis
#' @param dim.reduction.factor Correlation threshold over which correlated variables will be dropped
#' @param save.to when not NULL, folder path where result files will be written
#' @param suffix when save.to is not NULL, it adds a suffix to file names
#'
#' @return list with max four objects: PCA and correlation of variables, and if variable reduction is possible, PCA and correlation of reduced variables.
#'
#' @export
#'
prelim_ML <- function(df, ncp = 5, dim.reduction.factor = 0.7, save.to = NULL, suffix = ""){

    prelim.ML.obj <- list()

    if(!is.null(save.to)) mkdir(save.to)

    pca.summary.file <- ifelse(is.null(save.to),
                               "",
                               file.path(save.to, paste("PCA_summary", suffix , "txt", sep = ".") ) )

    pca.red.summary.file <- ifelse(is.null(save.to),
                                   "",
                                   file.path(save.to, paste("PCA_summary_reduced", suffix , "txt", sep = ".") ) )


    # Only numeric columns can be used in PCA and correlation
    df.num <- df[,  sapply(df, is.numeric)]

    # Do a Principal Component Analysis
    pca.df <- FactoMineR::PCA(df.num, scale.unit=TRUE, ncp=ncp, graph = F)
    prelim.ML.obj$pca <- pca.df
    summary(pca.df, nbelements = Inf , nbind = Inf, file = pca.summary.file)

    #dimdesc(pca.df)

    if(!is.null(save.to)){
        pdf(file.path(save.to, paste("PCA_graphs", suffix, "pdf", sep=".")), width = 10, height = 10)
        devNum = dev.cur()
    }
    plot(pca.df, cex=0.6, choix = "ind", unselect = 0.7, label="ind.sup" )
    plot(pca.df, cex=0.6,  choix = "var")
    if(!is.null(save.to)){
        dev.off(devNum)
    }

    #Calculate variable correlations
    corr.df <- cor(df.num)
    prelim.ML.obj$correlation <- corr.df
    tl.cex <- ifelse(ncol(df.num) > 30, 0.6, 1.0)
    if(!is.null(save.to)){
        pdf(file.path(save.to, paste("Correlation_graphs", suffix, "pdf", sep=".")), width = 10, height = 10)
        devNum = dev.cur()
    }
    corrplot::corrplot(corr.df, order="hclust", tl.cex = tl.cex, method = "ellipse", diag = F, tl.col = "darkblue")
    if(!is.null(save.to)){
        dev.off(devNum)
    }

    # Perform a similar analysis with reduced number of variables
    highlyCor <- caret::findCorrelation(corr.df, dim.reduction.factor)

    if(length(highlyCor) > 0){
        prelim.ML.obj$highly.correlated.vars <- names(df.num[, highlyCor])
        df.red <- df.num[,-highlyCor]
        corr.df.red <- cor(df.red)
        prelim.ML.obj$correlation.reduced <- corr.df.red

        if(!is.null(save.to)){
            pdf(file.path(save.to, paste("Correlation_Reduced_graphs", suffix, "pdf", sep=".")), width = 10, height = 10)
            devNum = dev.cur()
        }
        tl.cex.red <- ifelse(ncol(df.red) > 30, 0.6, 1.0)
        corrplot::corrplot(corr.df.red, order = "hclust", tl.cex = tl.cex.red, method = "ellipse", diag = F, tl.col = "darkblue")
        if(!is.null(save.to)){
            dev.off(devNum)
        }

        pca.df.red <- FactoMineR::PCA(df.red, scale.unit = T, ncp=5, graph = F)
        prelim.ML.obj$pca.reduced <- pca.df.red
        summary(pca.df.red, nbelements = Inf , nbind = Inf, file = pca.red.summary.file)

        if(!is.null(save.to)){
            pdf(file.path(save.to, paste("PCA_graphs_red", suffix, "pdf", sep=".")), width = 10, height = 10)
            devNum = dev.cur()
        }

        plot(pca.df.red, cex=0.6, choix = "ind", unselect = 0.7, label="ind.sup" )
        plot(pca.df.red, cex=0.6,  choix = "var")

        if(!is.null(save.to)){
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
#' @param factor.training percentage of data used for training. It takes the portion factor.training from the smallest class and a equivalent number of elements from the biggest class.
#' @param grep_rowname string containing a regular expression that limits the elements used for the training by choosing matching elements at rowname.
#'
#' @return list with two data frames: trainset and testset
#'
#' @export
#'
get_training_sample <- function(df, class_category, factor.training, grep_rowname=NULL){

    # Sampling for training
    df[, class_category] = as.logical(df[, class_category])

    df_discardedForTraining <- NULL
    if(!is.null(grep_rowname)){
        df_discardedForTraining <- df[grep(grep_rowname, rownames(df), invert = T) , ]
        df                      <- df[grep(grep_rowname, rownames(df)) , ]
    }

    df.class.true <- df[df[,class_category] == T,] # %>% filter_(class_category == T)
    df.class.false<- df[df[,class_category] == F,] # %>% filter_(class_category == F)

    #set the size of the training depending on the smallest dataset
    trainset.size <- round(min(nrow(df.class.true), nrow(df.class.false)) * factor.training, 0)

    trainset.index.true  <- sample(1:nrow(df.class.true),  size = trainset.size, replace = F)
    trainset.index.false <- sample(1:nrow(df.class.false), size = trainset.size, replace = F)

    trainset <- rbind(df.class.true[ trainset.index.true, ], df.class.false[ trainset.index.false, ])
    testset <-  rbind(df.class.true[-trainset.index.true, ], df.class.false[-trainset.index.false,])

    if(!is.null(df_discardedForTraining)) testset <- rbind(testset, df_discardedForTraining)

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
retrieve_data_modeling <- function(cell.line, is_perturbation, ML_method){

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


#' prepare_data_modeling
#' @description fix some variable names (from gene ontology), filter variables, and normalize data
#'
#' @param df data frame containing data to be modeled
#' @param useVariables a vector containing variable names to be used. When NULL, all variables are used.
#'
#' @export
#'
prepare_data_modeling <- function(df, useVariables=NULL){

    omics.model <- df %>% distinct(key, protein_entry, .keep_all=TRUE)

    # Combine in rownames: key and protein name.
    omics.key_protein <- paste(omics.model$key, omics.model$protein_entry, sep="#")

    if(!is.null(useVariables)) omics.model <- omics.model[, useVariables]

    #rename variables for the model
    renaming.vars <- stringr::str_extract(names(omics.model), "\\[GO:.*\\]")
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

    rownames(omics.model.sc) <- omics.key_protein

    return(omics.model.sc)
}

#' evaluate_ML
#' @description evaluate a Machine Learning prediction
#'
#' @param predictor ML predict output object
#' @param testset data set for testing
#' @param class_category classification category
#' @param save.to path to save results. When NULL, results will be displayed online
#'
#' @export
#'
evaluate_ML <- function(predictor, testset, class_category, save.to=NULL){

    # predictor=annFit.2.red
    # testset=omics.combined_tr$testset
    # class_category = class_category
    # save.to = path.ANN.Comb.red

    if(!is.null(save.to)) mkdir(save.to)

    pred.summary <- summary(predictor)

    testset$classification <- testset[, class_category]


    # If this is a neuralnetwork object it works different to other predictors
    if(class(predictor) == "nn" ){

        vals <- as.matrix(testset[, predictor$model.list$variables])
        testset$prediction <- as.numeric(neuralnet::compute(predictor,  vals)$net.result)

    }

    if(!class(predictor) == "nn" ){
        testset$prediction <- as.numeric(predict(predictor, testset))
    }

    f <- as.formula(paste(class_category, "prediction", sep = "~"))

    my.roc <- pROC::roc(classification ~ prediction, data=testset, plot=F, smooth=F)


    if(!is.null(save.to)) pdf(file.path(save.to, "roc_plot.pdf"))
    plot(my.roc)
    if(!is.null(save.to)) dev.off()

    if(!is.null(save.to)) sink(file.path(save.to, "prediction_summary.txt"))
    print(pred.summary)
    if(!is.null(save.to)) sink()

    if(!is.null(save.to)) pdf(file.path(save.to, "separation_plots.pdf"))
    colors <- ifelse(my.roc$response == 0, "black", "red")
    plot(my.roc$predictor,pch=19,cex=0.3, col=colors)
    y = as.numeric(my.roc$predictor)
    x = as.numeric(my.roc$response)
    plot(density(y[x==0]))
    lines(density(y[x>0]), col="red")
    if(!is.null(save.to)) dev.off()

    if(!is.null(save.to)) saveRDS(predictor, file=file.path(save.to, "predictor.Rds"  ))


    if(!is.null(save.to)) sink(file.path(save.to, "prediction_ROC_summary.txt"))

        print(my.roc$call)
        print(my.roc$auc)
        print(cat("Classification levels:", my.roc$levels) )
        print(cat("Direction:", my.roc$direction))
        print(cat("Calculated in percent:", my.roc$percent))
        print(cat("Sensitivities:", my.roc$sensitivities))
        print(cat("Specifities:", my.roc$specificities))
        print(cat("Thresholds:", my.roc$thresholds))

    if(!is.null(save.to)) sink()




    return(my.roc)
}
