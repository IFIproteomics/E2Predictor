library(NeuralNetTools)
library(caret)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(tidyr)

library(FactoMineR)
library(corrplot)

omics <- readRDS("/Users/napedro/CloudStation/tables.for.modeling/cell.lines/JY/omics_table_for_modeling.RData")

dbfile <- "/Users/napedro/CloudStation/databases/canonical/uniprot-human.tab.gz"
db <- readr::read_tsv(file = dbfile)
db <- db %>%
    filter(Status == "reviewed") %>%
    rename(protein_entry = Entry)

nn <- readRDS("/Users/napedro/CloudStation/data_analyses/machine_learning/ANN/cell.lines/JY/ANN.trained.GOandprotein.has_ligands_PROTEIN_QUANTITY_RPKM.avg_GO.positive_GO.negative_HL_2_has_ligands_PROTEIN_QUANTITY_RPKM.avg_GO.positive_GO.negative_HL_2.Rds")

plotnet(mod_in = nn)

crawler <- read_file_crawler("/Users/napedro/CloudStation/predictions.crowler/Uniprot_Predictions_Localization_Crawler.xlsx")


crawler.preML <- prelim.ML(df = crawler, dim.reduction.factor = 0.4, saveInFolder = "/Users/napedro/Desktop/loquesea")

corrplot(preML$correlation)
plot(preML$pca, cex=0.6,  choix = "var")




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
GO.index <- grep("^GO.", names(omics.model))
GO.index.pos <- GO.index[1:15]
GO.index.neg <- GO.index[16:30]
names(omics.model)[GO.index.pos] <- paste(names(omics.model)[GO.index.pos], "pos", sep=".")
names(omics.model)[GO.index.neg] <- paste(names(omics.model)[GO.index.neg], "neg", sep=".")

omics.model <- merge(omics.model, crawler.sc.id , by = "protein_entry", all.x = T)
#omics.model$GO.positive <- rowSums(omics.GO.pos)
#omics.model$GO.negative <- rowSums(omics.GO.neg)

omics.model.variables <- c("has_ligands", "PROTEIN_QUANTITY", "RPKM.avg", grep("^GO.*", names(omics.model), value = T)) #
#omics.model.variables <- c("has_ligands", "PROTEIN_QUANTITY", "RPKM.avg", "GO.positive", "GO.negative") #

#omics.model <- omics.model[, omics.model.variables]
omics.model <- omics.model %>% select(-c(protein_entry:DESCRIPTION, FDR.LEVEL:ENTRY,FILENAME, ligands:`Gene ontology (cellular component)`))
# Normalize variables
omics.model[is.na(omics.model)] <- 0.0

maxs <- apply(omics.model, 2, max)
mins <- apply(omics.model, 2, min)

num.vars <- ncol(omics.model) - 1

omics.model.sc <- as.data.frame(scale(omics.model, center = mins, scale = maxs - mins))

corr.omics.model.sc <- cor(omics.model.sc)
corrplot(corr.omics.model.sc, order="hclust", tl.cex = 0.6)

highlyCor <- findCorrelation(corr.omics.model.sc, 0.60)
#Apply correlation filter at 0.70,
#then we remove all the variable correlated with more 0.7.
omics.model.sc.red <- omics.model.sc[,-highlyCor]
corr.omics.model.sc.red <- cor(omics.model.sc.red)
corrplot(corr.omics.model.sc.red, order = "hclust", tl.cex = 0.6)


# PCA with function PCA
lig_idx <- grep("has_ligands", names(omics.model.sc) )
pca <- PCA(omics.model.sc, scale.unit=TRUE, ncp=5, quali.sup = c(lig_idx), graph = F)
summary(pca)
plot(pca, cex=0.6, habillage = lig_idx, choix = "ind", unselect = 0.7, label="ind.sup" )
plot(pca, cex=0.6, habillage = lig_idx, choix = "var")
#scale all the features,  ncp: number of dimensions kept in the results (by default 5)

dimdesc(pca)
#This line of code will sort the variables the most linked to each PC. It is very useful when you have many variables.

lig_idx.red <- grep("has_ligands", names(omics.model.sc.red) )
pca.red <- PCA(omics.model.sc.red, scale.unit = TRUE, ncp=5, quali.sup = c(lig_idx.red), graph = F)
summary(pca.red)
plot(pca.red, cex=0.6, habillage = 1, choix = "var")

