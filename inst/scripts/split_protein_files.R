library(E2Predictor)

splitCellLines = F
splitPerturbations = T

datafolder <- "/Users/napedro/CloudStation/protein_expression_data/perturbations/JY"
rest.file <- list.files(datafolder, "xlsx$")
setwd(datafolder)

df <- read_the_file(rest.file)
names(df) <- make.names(names(df))

cell.lines <- names(df)[11:length(names(df))]

dict.cell.lines <- make.names(c("HEK293",
                                 "LCLC-103H",
                                 "NIH-OVCAR-3",
                                 "COLO-699-N",
                                 "JAR",
                                 "MEL-526",
                                 "JY",
                                 "MEL-HO",
                                 "K562_A0201",
                                 "SK-MEL-37"))

# "HEK293", "LCLC.103H", "NIH.OVCAR.3", "COLO.699.N", "JAR", "MEL.526", "JY", "MEL.HO", "K562A2", "SK.MEL.37"

dict.perturbations <- c("Control", "Bortezomib", "IFNg", "Rapamycin", "DMOG")


split_cell_line <- function(cell_line){
     common_columns <- names(df)[1:9]
     df_split_common <- df %>%
                         select(description:entry)
     df_split <- df %>%
         select_(cell_line)

     df_split <- cbind(df_split_common, df_split)

     write_tsv(df_split, paste("protein_TOP3", cell_line, "tsv", sep = "."))

}

if(splitCellLines){
     nix <- sapply(cell.lines, split_cell_line)
}


split_perturbation <- function(perturbation){
     common_columns <- names(df)[1:9]
     df_split_common <- df %>%
         select(description:entry)
     df_split <- df %>%
         select_(matches(perturbation))

     df_split <- cbind(df_split_common, df_split)

     write_tsv(df_split, paste("protein_TOP3", cell_line, "tsv", sep = "."))

}

if(splitPerturbations){
     nix <- sapply(perturbations, split_perturbation)
}

