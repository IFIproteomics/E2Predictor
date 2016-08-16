library(E2Predictor)
library(dplyr)
library(stringr)


E2Predictor::initConfiguration()
E2Predictor::changeConfigutation(working.path = "/Users/napedro/CloudStation")


cell.lines <- c("JY", "LCLC.103H", "COLO.699.N", "MEL.526", "HEK293", "MEL.HO", "NIH.OVCAR.3", "SK.MEL.37", "K562A2")


uniprot <- readRDS(file.path(E2Predictor.Config$working.path, "tables.for.modeling", "cell.lines", cell.lines[1], "uniprot.db.RData"))

uniprot.go <- uniprot %>%
                select(protein_entry, `Gene ontology (GO)`) %>%
                mutate(GOid = strsplit(`Gene ontology (GO)`, split = ";"))  %>%
                unnest(GOid) %>%
                mutate(GOid = gsub(" ", "", GOid)) %>%
                group_by(GOid) %>%
                summarise(n_norm = n())

omics.tables.folder <- file.path(E2Predictor.Config$working.path,"tables.for.modeling", "cell.lines")

omics.JY          = readRDS( file.path(omics.tables.folder, cell.lines[1], "omics_table_for_modeling.RData")) %>% filter(!is.na(ligands))
omics.LCLC.103H   = readRDS( file.path(omics.tables.folder, cell.lines[2], "omics_table_for_modeling.RData")) %>% filter(!is.na(ligands))
omics.COLO.699.N  = readRDS( file.path(omics.tables.folder, cell.lines[3], "omics_table_for_modeling.RData")) %>% filter(!is.na(ligands))
omics.MEL.526     = readRDS( file.path(omics.tables.folder, cell.lines[4], "omics_table_for_modeling.RData")) %>% filter(!is.na(ligands))
omics.HEK293      = readRDS( file.path(omics.tables.folder, cell.lines[5], "omics_table_for_modeling.RData")) %>% filter(!is.na(ligands))
omics.MEL.HO      = readRDS( file.path(omics.tables.folder, cell.lines[6], "omics_table_for_modeling.RData")) %>% filter(!is.na(ligands))
omics.NIH.OVCAR.3 = readRDS( file.path(omics.tables.folder, cell.lines[7], "omics_table_for_modeling.RData")) %>% filter(!is.na(ligands))
omics.SK.MEL.37   = readRDS( file.path(omics.tables.folder, cell.lines[8], "omics_table_for_modeling.RData")) %>% filter(!is.na(ligands))
omics.K562A2      = readRDS( file.path(omics.tables.folder, cell.lines[9], "omics_table_for_modeling.RData")) %>% filter(!is.na(ligands))

omics.ligands <- rbind(omics.JY,
                       omics.LCLC.103H,
                       omics.COLO.699.N,
                       omics.MEL.526,
                       omics.HEK293,
                       omics.MEL.HO,
                       omics.NIH.OVCAR.3,
                       omics.SK.MEL.37,
                       omics.K562A2) %>%
                    select(protein_entry, `Gene ontology (GO)`) %>%
                    group_by(protein_entry) %>%
                    distinct(protein_entry, .keep_all = TRUE)

go.enrichment <- omics.ligands %>%
                        mutate(GOid = strsplit(`Gene ontology (GO)`, split = ";"))  %>%
                        unnest(GOid) %>%
                        mutate(GOid = gsub(" ", "", GOid)) %>%
                        group_by(GOid) %>%
                        summarise(numIDs = n()) %>%
                        arrange(desc(numIDs))

go.enrichment <- merge(go.enrichment, uniprot.go, by = "GOid") %>%
                    mutate(numIDs_norm = numIDs / n_norm) %>%
                    arrange(desc(numIDs_norm)) %>%
                    filter(n_norm >= 100 & !is.na(GOid)) %>%
                    arrange(numIDs_norm)

go.enrichment.negative <- go.enrichment[1:15,]
go.enrichment.positive <- (go.enrichment %>% arrange(desc(numIDs_norm)))[1:15,]

go.enrichment.folder <- file.path(E2Predictor.Config$working.path, "GO.enrichment")

E2Predictor::mkdir(go.enrichment.folder)

saveRDS(go.enrichment.positive, file.path(go.enrichment.folder, "go.enrichment.positive.Rds"))
saveRDS(go.enrichment.negative, file.path(go.enrichment.folder, "go.enrichment.negative.Rds"))


go.enrich.neg.str <- stringr::str_extract(go.enrichment.negative$GOid, "\\[GO:.*\\]")
go.enrich.pos.str <- stringr::str_extract(go.enrichment.positive$GOid, "\\[GO:.*\\]")

grep(go.enrich.neg.str[1], uniprot$`Gene ontology (GO)`, fixed = T)

