library(readr)
library(dplyr)
library(E2Predictor)

joinMods <- function(mods, removeSILAC = FALSE){

    mods <- strsplit(mods, ",")[[1]]
    mods_SILAC <- grepl("SILAC", mods)

    mods_return <- mods
    if(removeSILAC){
        mods_return <- mods[!mods_SILAC]
    }
    if(length(mods_return) == 0L){
        mods_return_txt <- ""
    }else{
        mods_return_txt <- paste(mods_return, sep =", ")
    }

    return(mods_return_txt)
}

silac_file <- "/Volumes/users/00 Samples 2014/Sample Set Information 2014/2014-092/ISOQuant_reports/Waters_SILAC_WF/2014-092 JY Kinetik_ SILAC time course template 20160720-122210_peptide_quantification_report.csv"

silac <- read_csv(silac_file)

silac_heavy <- silac %>% filter( grepl("SILAC", modifier) == TRUE )
silac_heavy$silac_num_mods <- countCharOccurrences(silac_heavy$modifier, "SILAC")
silac_heavy$silac_mods <- sapply(silac_heavy$modifier, joinMods, T) #joinMods(silac_heavy$modifier, T)
silac_heavy$sequence_mods <- paste(silac_heavy$sequence, silac_heavy$silac_mods, sep="_")

silac_light <- silac %>% filter( grepl("SILAC", modifier) == FALSE )
silac_light$sequence_mods <- paste(silac_light$sequence, silac_light$modifier, sep="_")

silac_linked <- merge(silac_light, silac_heavy, by = "sequence_mods", all = T, suffixes = c(".light", ".heavy") )
silac_linked.lim <- merge(silac_light, silac_heavy, by = "sequence_mods", all = F, suffixes = c(".light", ".heavy") )

silac_linked.lim$ratio.1 <- silac_linked.lim$`intensity in Default 1.light` / silac_linked.lim$`intensity in Default 1.heavy`
silac_linked.lim$ratio.2 <- silac_linked.lim$`intensity in Default 2.light` / silac_linked.lim$`intensity in Default 2.heavy`
silac_linked.lim$ratio.3 <- silac_linked.lim$`intensity in Default 3.light` / silac_linked.lim$`intensity in Default 3.heavy`
silac_linked.lim$ratio.4 <- silac_linked.lim$`intensity in Default 4.light` / silac_linked.lim$`intensity in Default 4.heavy`

silac_linked_seq <- merge(silac_light, silac_heavy, by="sequence", all = F)


silac_linked_proteins <- silac_linked.lim %>% group_by(accession.light) %>% summarise()
