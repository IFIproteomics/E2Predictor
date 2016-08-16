library(E2Predictor)

datafolder <- "/Users/napedro/CloudStation/mRNA_expression_data/cell_lines/rest"
rest.file <- "all_rpkm_genes"
setwd(datafolder)

df <- read_the_file(rest.file)
names(df) <- make.names(names(df))

cell.lines <- names(df)[2:length(names(df))]

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

split_cell_line <- function(cell_line){
     df_split <- df %>%
                 select_("ID", cell_line)
     write_tsv(df_split, paste("gene_RPKM", cell_line, "tsv", sep = "."))

}

nix <- sapply(cell.lines, split_cell_line)
