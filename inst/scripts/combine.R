library(readr)
library(dplyr)
library(tools)
library(tidyr)
library(E2Predictor)

cell.line <- "JY"

is_perturbation_exp <- F

dbfile <- "/Users/napedro/CloudStation/databases/canonical/uniprot-human.tab.gz"

proteomics.cellline.folder <- "/Users/napedro/CloudStation/protein_expression_data/cell_lines"
proteomics.perturbations.folder <- "/Users/napedro/CloudStation/protein_expression_data/perturbations"

genomics.cellline.folder <- "/Users/napedro/CloudStation/mRNA_expression_data/cell_lines"
genomics.perturbations.folder <- "/Users/napedro/CloudStation/mRNA_expression_data/perturbations"

ligandomics.cellline.folder <- ""
ligandomics.perturbations.folder <- ""

if(is_perturbation_exp){
    proteomics.file <- list.files(file.path(proteomics.perturbations.folder, cell.line), pattern = "xlsx$")[1]
}else{
    proteomics.file <- list.files(proteomics.cellline.folder, full.names = T)
    proteomics.file <- proteomics.file[grep(cell.line, proteomics.file, ignore.case = T)]

    genomics.file <- list.files(genomics.cellline.folder, recursive = F, include.dirs = F, full.names = T)
    genomics.file <- genomics.file[grep(cell.line, genomics.file, ignore.case = T)]

}


db <- read_tsv(file = dbfile)
db <- db %>% filter(Status == "reviewed")

proteomics <- read_file_protein(proteomics.file)
proteomics <- proteomics %>% rename(protein_entry = accession)

genomics <- read_file_RPKM(genomics.file)
genomics <- genomics %>% rename(gene_entry = X1, gene_RPKM = X2)

gene_has_prots <- gen2proteins(db)

protein_sum <- merge(proteomics, gene_has_prots, by = "protein_entry")
protein_sum <- merge(protein_sum, genomics, by = "gene_entry", all = T)

pr_sum_tmp <- protein_sum %>% group_by(protein_entry) %>% summarise(gene_RPKM_sum = sum(gene_RPKM, na.rm=T))
protein_sum <- merge(protein_sum, pr_sum_tmp, by= "protein_entry", all = T)
rm(pr_sum_tmp)
