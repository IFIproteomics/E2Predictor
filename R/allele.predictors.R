
#' predictNeoEpitopes predicts the most likely alleles for ligands
#'
#' @param ligands a data frame, which should contain at least two variables: Sequence, cell.line
#' @param predictor.value indicates which value is used to estimate the best prediction between different alleles
#'
#' @return netMHCpredictions data.frame containing allele predictions
#'
#' @export
#'
predictNeoEpitopes <- function(ligands, predictor.value = c("rank", "ic50"), seq.min.length = 8, seq.max.length = 14){

    # do for each cell line present in ligands data frame (in most cases there is just one, but...)
    cell.lines = unique(ligands$cell.line)

    predict.cellline <- function(cell.line){
        ligs.cellline <- ligands %>%
                            filter(cell.line == cell.line &
                                       nchar(Sequence) >= seq.min.length &
                                       nchar(Sequence) <= seq.max.length)

        message(paste("Predicting with netMHC 4.0 the alleles of ligands from ", cell.line, sep=" "))
        predictions.cellline <- netMHCpan.alleles(ligs.cellline$Sequence, cell.line)
        predictions.cellline$cell.line <- cell.line
        return(predictions.cellline)
    }

    predictions <- do.call("rbind", lapply(cell.lines, predict.cellline))

    predictions <- predictions %>% group_by(cell.line, Peptide)
    predictions$predictor <- NA

    predictions_min.rank <- predictions %>% top_n(1, desc(Rank))
    predictions_min.rank$predictor <- "rank"

    predictions_min.ic50 <- predictions %>% top_n(1, desc(nM))
    predictions_min.ic50$predictor <- "ic50"

    predictions_min <- predictions %>% filter(cell.line == "IonlyWantAnEmptyDataFrame")

    if("rank" %in% predictor.value){
        predictions_min <- rbind(predictions_min, predictions_min.rank)
    }
    if("ic50" %in% predictor.value){
        predictions_min <- rbind(predictions_min, predictions_min.ic50)
    }

    return(predictions_min)
}


#' exec.netMHCpan executes the command netMHCpan for a set of peptides of SAME length and ONE allele (a single prediction)
#'
#' @param allele allele name (defined as of E2Predictor.Config$hla_alleles)
#' @param peptide peptide sequence string
#' @param outputPath
#'
exec.netMHCpan <- function(allele, peptides, outputPath = file.path(E2Predictor.Config$working.path, "temp")){

    len_pep <- nchar(peptides[1])

    #check peptide length is the same for all sequences
    check_length <- unique(sapply(peptides, nchar))
    if(length(check_length) > 1){
        stop("All peptides supplied to E2Predictor::exec.netMHCpan function must have the same length!")
    }

    mkdir(outputPath)
    allele_txt <- gsub("\\*", "", allele)
    allele_txt <- gsub(":", "", allele_txt)
    output.file <- file.path(outputPath, "netmhcInput.fasta")

    zz <- file(output.file, open = "wt")
    sink(zz, type = "message")
    for(i in 1:length(peptides))
    {
        message(paste0(">", i))
        message(peptides[i])
    }
    sink(type = "message")
    close(zz)

    netMHCpan.command <- paste0(E2Predictor.Config$netMHCpath,
                            " -xls -xlsfile ", file.path(outputPath, "test.xls"),
                            " -l " , len_pep,
                            " -a HLA-", allele_txt,
                            " -f ", output.file
                            )

    message( paste0("Predicting ", len_pep, "-mers for allele ", allele))
    dd <- system(netMHCpan.command, wait = TRUE, ignore.stdout = T)

    netMHC.prediction <- readr::read_tsv(file.path(outputPath, "test.xls"), skip = 1)
    netMHC.prediction$allele <- allele

    return(netMHC.prediction)
}


netMHCpan.alleles <- function(peptides, cell.line){

    alleles <- E2Predictor.Config$hla_alleles[cell.line][[1]]

    #do for each different peptide length in vector peptides
    peptides.df <- data.frame( peptides = peptides, pep.length = nchar(peptides))
    pep.lengths <- sort(unique(peptides.df$pep.length))

    netMHC.length <- function(sel.pep.length){
        peps.l <- as.vector((peptides.df %>% filter(pep.length == sel.pep.length))$peptides)
        peps.l.preds <- do.call("rbind", lapply(alleles, exec.netMHCpan, peps.l))
        return(peps.l.preds)
    }

    netMHC.predictions <- do.call("rbind", lapply(pep.lengths, netMHC.length))

    return(netMHC.predictions)
}

