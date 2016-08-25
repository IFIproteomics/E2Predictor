library(E2Predictor)

initConfiguration(netMHCpath = "/Users/napedro/external_sources/netMHC-4.0/netMHC",
                  working.path = "/Users/napedro/CloudStation")

cell.lines <- names(E2Predictor.Config$hla_alleles)
conditions <- c("Control", "Bortezomib", "DMOG", "IFNg", "Rapamycin")
cell.line <- cell.lines[4]

dbfile <- file.path(E2Predictor.Config$working.path, "databases", "canonical", "uniprot-human.tab.gz")
is_perturbation <- T
condition <- conditions[1]
addgrepcond <- "_II_"
minReplicates <- 1
# You don't need to have the allele predictions, it makes the predictions when allele.predictions = NULL
#allele.predictions <- readRDS(file.path(E2Predictor.Config$working.path, "allele.predictions", ifelse(is_perturbation, "perturbations", "cell_lines"), cell.line, "allele.predictions.Rds"))
allele.predictions <- NULL
allele.predictor <- "ic50"

nM.threshold <- 1000
saveAllResults <- T

GOenrichment.pos <- readRDS(file.path(E2Predictor.Config$working.path, "GO.enrichment", "go.enrichment.positive.Rds"))
GOenrichment.neg <- readRDS(file.path(E2Predictor.Config$working.path, "GO.enrichment", "go.enrichment.negative.Rds"))

GOcategories <- c(GOenrichment.pos$GOid, GOenrichment.neg$GOid)

omics <- merge_omics(selected.cell.line = cell.line, dbfile = dbfile,
                       is_perturbation = is_perturbation, condition = condition, additionalGrepCondition = addgrepcond,
                       minReplicates = minReplicates, allele.predictions = allele.predictions,
                       allele.predictor = allele.predictor, nM.threshold = nM.threshold,
                       saveAllResults = saveAllResults, GOcategories = GOcategories)

