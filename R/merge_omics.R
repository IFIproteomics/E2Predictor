#' merge.omics merges proteomics, genomics, and ligandomics data of an experiment
#'
#' @param selected.cell.line cell line to integrate data.If it is a 3-unit vector, it defines (Proteomics, Genomics, Ligandomics)
#' @param dbfile Uniprot database filename (NOT in fasta format! Use the xls format)
#' @param is_perturbation indicates if the experiment is a perturbation experiment. If it is a 3-unit vector, it defines (Proteomics, Genomics, Ligandomics)
#' @param condition in case of a perturbation, this indicates which condition you want to use. If it is a 3-unit vector, it defines (Proteomics, Genomics, Ligandomics)
#' @param additionalGrepCondition regular expressions(s) defined to find which conditions from perturbations should be used. If it is a 3-unit vector, it defines (Proteomics, Genomics, Ligandomics)
#' @param minReplicates for ligandomics data, the minimum number of biological replicates required to validate a ligand
#' @param allele.predictions data frame containing allele predictions for the ligands. If NULL, they will be estimated
#' @param allele.predictor c("rank", "ic50") predictor to be used.
#' @param nM.threshold minimum nanomolar threshold for validating the predictions
#' @param saveAllResults save allele predictions and merged omics table at the working directory
#' @param GOcategories a string vector containing gene ontology categories you want to include as separate variables
#'
#' @export
#'
merge_omics <- function(selected.cell.line, dbfile,
                        is_perturbation = F, condition = NULL, additionalGrepCondition = "",
                        minReplicates = 1, allele.predictions = NULL, allele.predictor = "ic50",
                        nM.threshold = 1000, saveAllResults = F, GOcategories = NULL){


    if(!length(selected.cell.line) %in% c(1,3)){
        message("Error: length of vector selected.cell.line must be 1 or 3")
        stop()
    }
    if(!length(is_perturbation) %in% c(1,3)){
        message("Error: length of vector is_perturbation must be 1 or 3")
        stop()
    }
    if(!length(condition) %in% c(1,3)){
        message("Error: length of vector condition must be 1 or 3")
        stop()
    }
    if(!length(additionalGrepCondition) %in% c(1,3)){
        message("Error: length of vector additionalGrepCondition must be 1 or 3")
        stop()
    }


    if(length(selected.cell.line)==1) selected.cell.line = rep(selected.cell.line, 3)
    if(length(is_perturbation)==1) is_perturbation = rep(is_perturbation, 3)
    if(length(condition)==1) condition = rep(condition, 3)
    if(length(additionalGrepCondition)==1) additionalGrepCondition = rep(additionalGrepCondition, 3)

    ## Function config #################
    experiment_folder <- ifelse(is_perturbation, "perturbations", "cell_lines")

    proteomics.cellline.folder <- file.path(E2Predictor.Config$working.path, "protein_expression_data", experiment_folder[1])
    genomics.cellline.folder <- file.path(E2Predictor.Config$working.path, "mRNA_expression_data", experiment_folder[2])
    ligandomics.cellline.folder <- file.path(E2Predictor.Config$working.path, "ligandomes", experiment_folder[3])

    allele.predictions.dir <- file.path(E2Predictor.Config$working.path, "allele.predictions", experiment_folder, selected.cell.line[3])
    allele.predictions.filename <- "allele.predictions.Rds"

    omics.table.folder.name.pr <- paste("Prot", selected.cell.line[1], ifelse(is_perturbation[1], paste("Perturb", condition[1], additionalGrepCondition[1], sep="_" ), ""),  sep="_")
    omics.table.folder.name.gn <- paste("Gen", selected.cell.line[2], ifelse(is_perturbation[2], paste("Perturb", condition[2], additionalGrepCondition[2], sep="_" ), ""),  sep="_")
    omics.table.folder.name.lg <- paste("Lig", selected.cell.line[3], ifelse(is_perturbation[3], paste("Perturb", condition[3], additionalGrepCondition[3], sep="_" ), ""),  sep="_")
    omics.table.folder.name <- paste(omics.table.folder.name.pr, omics.table.folder.name.gn, omics.table.folder.name.lg, sep="-")
    omics.table.folder.name <- gsub(" ", "", omics.table.folder.name)

    omics.tables.folder <- file.path(E2Predictor.Config$working.path,"tables.for.modeling", omics.table.folder.name)


    if(is_perturbation[1]){

        proteomics.cellline.folder <- file.path(proteomics.cellline.folder, selected.cell.line[1])
        proteomics.file <- list.files(proteomics.cellline.folder, all.files = F, full.names = T)

    }else{
        proteomics.file <- list.files(proteomics.cellline.folder, all.files = F, full.names = T)
        proteomics.file <- proteomics.file[grep(selected.cell.line[1], proteomics.file)]
    }

    if(is_perturbation[2]){

        genomics.file <- list.files( file.path(genomics.cellline.folder, selected.cell.line[2]) , pattern = condition[2], full.names = T)
    }else{

        genomics.file <- list.files(genomics.cellline.folder, pattern = selected.cell.line[2], full.names = T)
    }

    if(is_perturbation[3]){

        ligandomics.files <- list.dirs(file.path(ligandomics.cellline.folder, selected.cell.line[3]), recursive = F)
        ligandomics.files <- ligandomics.files[grep(condition[3], ligandomics.files)]
        if(nchar(additionalGrepCondition[3]) > 0){
            ligandomics.files <- ligandomics.files[grep(additionalGrepCondition[3], ligandomics.files)]
        }

        allele.predictions.dir <- file.path(allele.predictions.dir,
                                            paste(condition[3], gsub(" ", "", additionalGrepCondition[3]), sep = "_"))

    }else{

        ligandomics.files <- file.path(ligandomics.cellline.folder, selected.cell.line[3])
    }

    ####################################


    ## Read uniprot database ###########
    message("Reading Uniprot database...")
    db <- readr::read_tsv(file = dbfile)
    db <- db %>%
        filter(Status == "reviewed") %>%
        rename(protein_entry = Entry)
    message(paste0("Number of protein entries in Uniprot database: ", nrow(db)))

    ####################################



    ## Read Proteomics data ############
    message("Reading Proteomics data...")
    proteomics <- read_file_protein(proteomics.file)
    if(is_perturbation[1]){
        proteomics <- proteomics.select.condition(proteomics, condition[1])
    }else{
        proteomics <- proteomics %>% rename_( "PROTEIN_QUANTITY" = selected.cell.line[1] )
    }

    # filter protein names from the proteomics data not present at the Uniprot DB
    # (from other organisms, like Yeast or ...)
    proteomics <- proteomics %>% filter(protein_entry %in% db$protein_entry)
    message(paste0("Number of proteins in Proteomics data: ", nrow(proteomics)))
    ####################################

    ## Read Ligandomics data ###########
    message("Reading Ligandomics data...")
    ligandomics <- read_ligands_PEAKS(ligandomics.files, minReplicates = minReplicates, ligands_min_length = 8, ligands_max_length = 14, uniprotDB = db)
    ligandomics$cell.line <- selected.cell.line[3]
    message(paste0("Number of ligands: ", nrow(ligandomics)))
    ####################################

    ## Allele predictions ##############

    if(is.null(allele.predictions)){
        allele.predictions <- predictNeoEpitopes(ligandomics)
    }

    allele.predictions.1pred <- allele.predictions %>% ungroup() %>%
        filter(predictor == allele.predictor, cell.line == selected.cell.line[3]) %>%
        mutate(allele = ifelse(nM < nM.threshold, allele, NA)) %>%
        select(Peptide, allele) %>%
        rename(Sequence = Peptide, predicted.allele = allele)

    ligandomics <- merge(ligandomics, allele.predictions.1pred, by="Sequence", all.x = T)

    ####################################

    ## Read Genomics data ##############
    message("Reading Genomics data...")
    genomics <- read_file_RPKM(genomics.file, useAverage = T, minValueRPKM = 0.0)

    db.genesprots <- db %>%
        select( protein_entry, `Entry name`, `Gene names`, `Protein names`, `Gene ontology (GO)` ) %>%
        mutate(gene_entry = strsplit(`Gene names`, " ")) %>%
        tidyr::unnest(gene_entry)

    # Filter db.genesprots dictionary by known gene entries (Allumina readings)
    db.genesprots <- db.genesprots %>%
        filter(gene_entry %in% genomics$gene_entry)
    message(paste0("Number of genes: ", nrow(genomics)))
    ####################################

    ## Merge different -omics ##########
    ## we use the dictionary db.genesprots as a dictionary, so that we can merge proteomics and genomics data.
    ## there are still some multiple matching proteins to genes: we use sum RPKM values for these cases.
    message("Merging all tables...")

    # merge proteomics
    omics <- merge(db.genesprots, proteomics, by="protein_entry", all = T)

    # merge genomics
    omics <- merge(omics, genomics, by="gene_entry", all.x = T)

    prot_has_ligands <- ligandomics %>%
        mutate(protein_entry = strsplit(protein_entries, ";")) %>%
        tidyr::unnest(protein_entry)

    prot_has_ligands_nest.seqs <- prot_has_ligands %>%
        group_by(protein_entry) %>%
        tidyr::nest(ligands = Sequence) %>%
        rename(ligands = data)


    prot_has_ligands_nest.alleles <- prot_has_ligands %>%
        group_by(protein_entry) %>%
        tidyr::nest(predicted.allele) %>%
        rename(predicted.allele = data)

    prot_has_ligands_nest.numentries <- prot_has_ligands %>%
        group_by(protein_entry) %>%
        tidyr::nest(ligand.num_protein_entries = num_protein_entries) %>%
        rename(ligand.num_protein_entries = data)


    omics <- merge(omics, prot_has_ligands_nest.seqs, by = "protein_entry", all.x = T)
    omics <- merge(omics, prot_has_ligands_nest.alleles, by = "protein_entry", all.x = T)
    omics <- merge(omics, prot_has_ligands_nest.numentries, by = "protein_entry", all.x = T)

    omics <- merge(omics, (db %>% select(protein_entry, `Gene ontology (cellular component)`)))
    message("done.")

    ####################################

    ## add GO variables ################
    addGOvariable <- function(GOcategory, df){
        df.go <- gsub(" ", "", df$`Gene ontology (GO)` )
        df.go.cat <- grepl(GOcategory, df.go, fixed = T)
        df[, GOcategory] <- df.go.cat
        return(df)
    }

    if(!is.null(GOcategories)){
        message("Adding separate GO categories...")
        for(GOcat in 1:length(GOcategories)){
            omics <-  addGOvariable(GOcategories[GOcat], omics)
        }
        message("done.")
    }

    ####################################


    ## Save results ####################
    if(saveAllResults){
        mkdir(allele.predictions.dir[3])
        saveRDS(allele.predictions,file.path(allele.predictions.dir[3], allele.predictions.filename))

        mkdir(omics.tables.folder)
        saveRDS(omics, file = file.path(omics.tables.folder, "omics_table_for_modeling.RData"))
        saveRDS(db, file = file.path(omics.tables.folder, "uniprot.db.RData"))
    }
    ####################################

    return(omics)
}
