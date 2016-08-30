#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#library(Biostrings)

##### Read specific file formats and filter out innecessary columns

#' read_file_RPKM reads a file with RPKM readings of a SINGLE measurement (2 columns!)
#'
#' @param thefiles vector of file name(s) containing RPKM readings
#' @param useAverage if TRUE, averaged values of all file names provided are returned
#'
#' @export
#'
read_file_RPKM <- function(thefiles, useAverage = FALSE, minValueRPKM = 0.0, ...){

    df <- do.call("rbind", lapply(thefiles, read_file_add_filename, keepFolderName = F, UseGenericHeaders = T))
    names(df) <- c("gene_entry", "RPKM", "FileName")

    if(useAverage){
        df <- df %>% group_by(gene_entry) %>% summarise(RPKM.avg = mean(RPKM))
    }

    df[df < minValueRPKM] <- NA

    return(df)
}

#' read_file_protein reads a file with protein readings
#'
#' @param thefile file with protein readings
#'
#' @export
#'
read_file_protein <- function(thefile, ...){
    # thefile = "/Users/napedro/CloudStation/protein_expression_data/perturbations/JY/2015-078 Jenny Perturbations JY_user designed 20160306-084126_quantification_report_Values.xlsx"
    skip = 0
    is_excel_file = F
    if(tools::file_ext(thefile) == "xls" | tools::file_ext(thefile) == "xlsx" ){
        skip = 1
        is_excel_file = T
    }

    df <- read_file_add_filename(thefile = thefile, keepFolderName = T, skip = skip)
    names(df) <- make.names(toupper(names(df)))
    df <- df %>% rename(protein_entry = ACCESSION)

    return(df)
}

#' proteomics.select.condition selects the quantitative data from a proteomics (ISOQuant) data frame
#'
#' @param proteindf data frame containing Proteomics data (ISOQuant format)
#' @param condition condition or name you want to use to filter data
#' @param useAverage if TRUE, averaged values of columns found under 'condition' will be returned
#'
#' @export
#'
proteomics.select.condition <- function(proteindf, condition, useAverage = T){

    #proteindf <- proteomics

    proteindf.common <- proteindf[,1:9]

    proteindf.condition <- proteindf %>% select(grep(condition, names(proteindf), ignore.case = T))

    if(useAverage){
        proteindf.condition <- proteindf.condition %>% select(grep("AVERAGE", names(proteindf.condition)))
        proteindf.condition <- rowMeans(proteindf.condition, na.rm = T)
    }



    proteindf.condition <- cbind(proteindf.common, PROTEIN_QUANTITY = proteindf.condition)

    return(proteindf.condition)
}



#' read_ligands_PEAKS reads and combines PEAKS identifications of PEAKS results folders.
#'
#' @param PEAKSresultsFolders folders containing PEAKS results
#' @param minReplicates minimum number of replicates (in different PEAKS results folders) desired for each identified peptide
#' @param ligands_min_length minimum length of ligand sequences
#' @param ligands_max_length maximum length of ligand sequences
#'
#' @export
#'
read_ligands_PEAKS <- function(PEAKSresultsFolders,
                                    minReplicates = 1,
                                    ligands_min_length = 8,
                                    ligands_max_length = 14,
                                    uniprotDB = NULL,
                                    alleles.prediction = NULL){


    ligands <- do.call("rbind", lapply(PEAKSresultsFolders, read_batch_in_folder,
                                                            filepattern = "protein-peptides.csv$",
                                                            keepFolderName = T,
                                                            UseGenericHeaders = F))

    ligands$Sequence <- gsub("\\s*\\([^\\)]+\\)","",as.character(ligands$Peptide))
    ligands$Sequence <- unlist(stringr::str_extract_all(ligands$Sequence, "[A-Z]{2,}"))
    #ligands$Sequence_length <- str_length(ligands$Sequence)

    ligands <- ligands %>% filter(Length >= ligands_min_length & Length <= ligands_max_length)

    ligands.confident <- ligands %>%
                        group_by(Sequence) %>%
                        mutate(num_replicates = n_distinct(FileName)) %>%
                        filter(num_replicates >= minReplicates) %>%
                        arrange(FileName) %>%
                        filter(row_number() == 1) %>%
                        select(`Protein Accession`:Area, Sequence:num_replicates ) %>% ungroup()

    if(!is.null(uniprotDB)){

        # Relate sequences with proteins
        message("Mapping peptides to proteins. This may take a while; if you have a single core system grab a coffee...")
        no_cores <- parallel::detectCores() #- 1

        #estimate time (roughly)
        fraction_estimate <- 0.01
        ds_st <- sample_frac(ligands.confident, size = fraction_estimate)
        cl <- parallel::makeCluster(no_cores, type = "PSOCK") # type = "FORK" is more efficient, but it does not work under Windows
        est.time.fraction <- system.time(parallel::parLapply(cl, ds_st$Sequence, grep, uniprotDB$Sequence, fixed=T))[3]
        td <- lubridate::seconds_to_period(est.time.fraction / fraction_estimate)
        parallel::stopCluster(cl)
        message("Estimated time for mapping:")
        message(sprintf('%02g hour(s) %02g minute(s) %02g second(s)', td@hour, lubridate::minute(td), lubridate::second(td)))


        # Relate sequences with proteins
        cl <- parallel::makeCluster(no_cores, type = "PSOCK") # type = "FORK" is more efficient, but it does not work under Windows
        seqs_prots <- parallel::parLapply(cl, ligands.confident$Sequence, grep, uniprotDB$Sequence, fixed=T)
        parallel::stopCluster(cl)
        # system.time(seqs_prots <- pblapply(ligands.confident$Sequence, grep, uniprotDB$Sequence, fixed=T))

        # link entries to ligands table
        getEntries <- function(indices_list){
            myl <- unlist(indices_list)
            return( paste(uniprotDB$protein_entry[myl], collapse = ";") )
        }

        getNumEntries <- function(indices_list){
            myl <- unlist(indices_list)
            return(length(myl))
        }

        ligands.confident$protein_entries <- sapply(seqs_prots, getEntries)
        ligands.confident$num_protein_entries <- sapply(seqs_prots, getNumEntries)

        message("Mapping done.")
    }

    return(ligands.confident)
}

AAStringSet.to.dataframe <- function(aastringsetobj){
    aastringsetobj <- fasta[1]
    proteinchar <- as.character(aastringsetobj)

    # parse names
    protein_split <- strsplit(names(proteinchar), "\\|")[[1]]
    protein_entry <- protein_split[2]
    protein_accession <- strsplit(protein_split[3], "\\s")[[1]][1]
    protein_description <- gsub(protein_accession, "", protein_split[3])
    protein_sequence <- unname(proteinchar)

    protein_df <- data.frame( protein_entry = c(protein_entry)
                ,protein_accession = c(protein_accession)
                ,protein_description = c(protein_description)
                ,protein_sequence = c(protein_sequence)
                )
    return(protein_df)
}

read_file_database <- function(thefile, ...){
    # thefile = "/Users/napedro/CloudStation/databases/canonical/UPSPhuman130605referenceproteome20251entries.fasta"

    if(tools::file_ext(thefile) == "fasta"){
        fasta <- readAAStringSet(thefile, format="fasta",
                        nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=T)
    }

    my_database <- apply(fasta, 1, AAStringSet.to.dataframe)
}


#' read_file_crawler reads a file produced by crawler script and
#' prepares the table to be used for modeling
#'
#' @param thefile crawler results file
#'
#' @export
#'
read_file_crawler <- function(thefile, ...){

    crawler <- readxl::read_excel(thefile, sheet = 1)

    crawler$accession <- gsub("\\|", "", stringr::str_extract(string = crawler$FastaTitle, pattern = "\\|.*\\|") )
    crawler <- crawler %>% distinct(accession, .keep_all =T)


    crawler <- crawler %>%
        tidyr::spread(CELLO.Location, CELLO.Score, fill=0.0, sep=".") %>% select(-CELLO.Location.NA)

    crawler <- crawler %>%
        tidyr::spread(SUBLOC.Location, SUBLOC.Reliability, fill=0.0, sep=".") %>% select(-SUBLOC.Location.NA)

    crawler <- crawler %>%
        tidyr::spread(TARGETP.Location, TARGETP.Reliability, fill=0.0, sep=".") %>% select(-TARGETP.Location.NA)

    crawler <- crawler %>%
        tidyr::spread(WOLFPSORT.Location, WOLFPSORT.Score, fill=0.0, sep=".") %>% select(-WOLFPSORT.Location.NA)

    crawler2 <- crawler %>% select(-c(accession, id, AccessionNumber, FastaTitle))

    crawler2[is.na(crawler2)] <- 0.0

    maxs <- apply(crawler2[, 2:ncol(crawler2)], 2, max)
    mins <- apply(crawler2[, 2:ncol(crawler2)], 2, min, na.rm=T)

    num.vars <- ncol(crawler2) - 1

    crawler.sc <- as.data.frame(scale(crawler2[, 2:ncol(crawler2)], center = mins, scale = maxs - mins))


    crawler.sc.id <- cbind(protein_entry = crawler$accession, crawler.sc)

    return(crawler.sc.id)
}


#' load_data_modeling
#' @description loads a data modeling table
#'
#' @param cell.line cell line to be loaded
#' @param is_perturbation data from perturbations?
#' @param condition if data from perturbations, condition to be loaded
#' @param addgrepcond an extra parameter to locate the right condition. Example: to distinguish conditions like "Control I" and "Control II"
#' @param requireAllele proteins with ligands require at least one ligand associated to an allele.
#'
#' @return a data.frame containing data for modeling
#'
#' @export
#'
load_data_modeling <- function(cell.line, is_perturbation, condition="", addgrepcond = "", requireAllele = F){

    experim.folder <- ifelse(is_perturbation, "perturbations", "cell.lines")
    modeling.folder <- file.path(E2Predictor.Config$working.path, "tables.for.modeling", experim.folder, cell.line)

    if(is_perturbation) modeling.folder <- file.path(modeling.folder, paste(condition, addgrepcond, sep = "_") )
    modeling.folder <- gsub(" ", "", modeling.folder)
    modeling.file <- file.path(modeling.folder, "omics_table_for_modeling.RData")

    # load data for modeling
    omics <- readRDS(modeling.file)

    omics <- omics %>% mutate(has_ligands = ifelse(is.na(ligands), FALSE, TRUE))  # variable to be predicted

    if(requireAllele) omics$has_ligands <- sapply(omics$predicted.allele, function(x) !all(is.na(x)) )

    omics$key <- paste(cell.line, "perturbation", is_perturbation, condition, addgrepcond, "requireAllele", requireAllele, sep = "_" )

    # Just for the moment, so that files are homogene.
    if("FILENAME" %in% names(omics)) omics$FILENAME <- NULL

    return(omics)
}

