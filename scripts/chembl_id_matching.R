library(tidyverse)
library(taigr)
library(useful)

# A small utility function to collapse equivalent rows 
update.pk <- function(df){
  ixs <- df %>% 
    dplyr::mutate(value = 1) %>% 
    dplyr::filter(is.finite(secondary.key)) %>% 
    reshape2::acast(secondary.key ~ primary.key, fill = 0)
  
  xx <- crossprod(ixs)
  ixs <- colnames(xx)
  ixs1 <- ixs[colSums(xx) > diag(xx)]
  if(length(ixs1) > 1){
    xx <- (xx[ixs1, ixs1] > 0) + 0 
    A <- solve((diag(rep(1, dim(xx)[1])) - xx))
    colnames(A) <- ixs1; rownames(A) <- ixs1
    A %<>% 
      reshape2::melt() %>% 
      dplyr::filter(value > 0) %>% 
      dplyr::group_by(Var1) %>% 
      dplyr::summarise(Var2 = min(Var2)) %>% 
      dplyr::filter(Var1 != Var2) %>% 
      dplyr::rename(primary.key = Var1, primary.key.new = Var2)
    
    B <- setdiff(df$primary.key, A$primary.key)
    B <- tibble(primary.key = B, primary.key.new = B)
    C <- dplyr::bind_rows(A, B)
  } else{
    C <- df %>%
      dplyr::distinct(primary.key) %>% 
      dplyr::mutate(primary.key.new = primary.key)
  }
  
  return(C)
  
}


# Load the existing compound meta-data
compound.metadata <- load.from.taiga(data.name='compound-metadata-de37', data.version=14, data.file='compound_metadata') %>%
  dplyr::select(-InChIKey) %>% dplyr::distinct()

# Iron-out the confusing ID's, Smiles, Drug Name's and Synonyms
# Number of rows drops from 7692 to 7314
# Still far from perfect, example: VENETOCLAX vs VENOTOCLAX
compound.metadata <- compound.metadata %>% 
  dplyr::mutate(cm_ix = 1:n()) %>% 
  tidyr::separate_rows(referenced_datasets, sep = ";") %>% 
  tidyr::separate_rows(IDs, sep = ";") %>% 
  tidyr::separate_rows(Synonyms, sep = ";") %>% 
  tidyr::pivot_longer(cols = c("Drug.Name", "Synonyms"),
                      names_to = "dummy", values_to = "Drug.Name") %>%
  dplyr::select(-dummy) %>% 
  dplyr::filter(Drug.Name != "", IDs != "", referenced_datasets != "") %>% 
  dplyr::distinct() %>%
  dplyr::mutate(ix1 = dplyr::group_indices(., Drug.Name)) %>%
  dplyr::mutate(ix2 = dplyr::group_indices(., IDs)) %>%
  dplyr::mutate(ix3 = dplyr::group_indices(., SMILES)) %>%
  dplyr::mutate(ix3 = ifelse(SMILES == "", NA, ix3)) 

compound.metadata <- compound.metadata %>%
  dplyr::mutate(primary.key = cm_ix,
                secondary.key = ix1) %>% 
  dplyr::distinct(primary.key, secondary.key) %>%
  update.pk() %>% 
  dplyr::rename(cm_ix = primary.key, cm_ix.new = primary.key.new) %>%
  dplyr::right_join(compound.metadata) %>%
  dplyr::mutate(cm_ix = cm_ix.new) %>% 
  dplyr::select(-ix1, -cm_ix.new) %>% 
  dplyr::distinct()

compound.metadata <- compound.metadata %>%
  dplyr::mutate(primary.key = cm_ix,
                secondary.key = ix2) %>% 
  dplyr::distinct(primary.key, secondary.key) %>%
  update.pk() %>% 
  dplyr::rename(cm_ix = primary.key, cm_ix.new = primary.key.new) %>%
  dplyr::right_join(compound.metadata) %>%
  dplyr::mutate(cm_ix = cm_ix.new) %>% 
  dplyr::select(-ix2, -cm_ix.new) %>% 
  dplyr::distinct()

compound.metadata <- compound.metadata %>%
  dplyr::mutate(primary.key = cm_ix,
                secondary.key = ix3) %>% 
  dplyr::distinct(primary.key, secondary.key) %>%
  update.pk() %>% 
  dplyr::rename(cm_ix = primary.key, cm_ix.new = primary.key.new) %>%
  dplyr::right_join(compound.metadata) %>%
  dplyr::mutate(cm_ix = cm_ix.new) %>% 
  dplyr::select(-ix3, -cm_ix.new) %>% 
  dplyr::distinct()

compound.metadata <- compound.metadata %>%
  tidyr::separate_rows(MOA, sep = ";") %>% 
  tidyr::separate_rows(repurposing_target, sep = ";") %>% 
  tidyr::separate_rows(Indication, sep = ";") %>% 
  tidyr::separate_rows(`Disease Area`, sep = ";") %>% 
  tidyr::separate_rows(Phase, sep = ";") %>% 
  dplyr::group_by(cm_ix) %>%
  dplyr::mutate(Synonyms = paste0(sort(unique(Drug.Name)), collapse= ";"),
                IDs = paste0(sort(unique(IDs)), collapse= ";"),
                referenced_datasets = paste0(sort(unique(referenced_datasets)), collapse= ";"),
                SMILES = paste0(sort(unique(SMILES[SMILES != ""])), collapse= ";"),
                MOA = paste0(sort(unique(MOA[MOA != ""])), collapse= ";"),
                repurposing_target = paste0(sort(unique(repurposing_target[repurposing_target != ""])), collapse= ";"),
                Indication = paste0(sort(unique(Indication[Indication != ""])), collapse= ";"),
                `Disease Area` = paste0(sort(unique(`Disease Area`[`Disease Area` != ""])), collapse= ";"),
                Phase = paste0(sort(unique(Phase[Phase != ""])), collapse= ";")) %>%
  dplyr::top_n(1, Drug.Name) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() 



# Load the ChemBL export
chembl <- load.from.taiga(data.name='compound-metadata-de37', data.version=15, data.file='CHEMBL') %>%
  dplyr::distinct(`ChEMBL ID`, Name, Synonyms, Type, `Max Phase`, Smiles, `Inchi Key`)


# list of compounds in the meta-data
compound.list <- compound.metadata  %>% 
  tidyr::separate_rows(Synonyms, sep = ";") %>% 
  tidyr::pivot_longer(cols = c("Drug.Name", "Synonyms"),
                      names_to = "dummy", values_to = "Drug.Name") %>%
  dplyr::select(-dummy) %>% dplyr::filter(Drug.Name != "") %>% 
  dplyr::distinct(cm_ix, Drug.Name, SMILES) %>%
  dplyr::distinct()

# list of compounds in ChemBL
chembl.list <- chembl %>% 
  dplyr::distinct(`ChEMBL ID`, Name, Synonyms, Smiles) %>% 
  tidyr::separate_rows(Synonyms, sep = fixed("\\|")) %>%  
  tidyr::separate_rows(Synonyms, sep = fixed(", ")) %>%  
  tidyr::separate_rows(Smiles, sep = fixed(", ")) %>%  
  tidyr::pivot_longer(cols = c("Name", "Synonyms"),
                      names_to = "dummy", values_to = "Drug.Name") %>%
  dplyr::select(-dummy) %>% dplyr::filter(Drug.Name != "") %>% 
  dplyr::mutate(Drug.Name = toupper(Drug.Name)) %>% 
  dplyr::rename(SMILES = Smiles) %>%
  dplyr::distinct()

# mapping two lists by name
name.mapped <- compound.list %>% 
  dplyr::inner_join(chembl.list, by = "Drug.Name")

# mapping two lists by SMILES
smiles.mapped <- compound.list %>% 
  dplyr::filter(SMILES != "") %>% 
  dplyr::inner_join(chembl.list %>% 
                      dplyr::filter(SMILES != ""),
                    by = "SMILES")

# combine the maps - note some compounds has multiple ChEMBL IDs
map <- name.mapped %>% 
  dplyr::distinct(cm_ix, `ChEMBL ID`) %>% 
  dplyr::bind_rows(smiles.mapped %>% 
                     dplyr::distinct(cm_ix, `ChEMBL ID`)) %>% 
  dplyr::distinct()


# put everything together
# note the each row of the original dataset is maintained as the unique key
compound.metadata.expanded <- compound.metadata %>% 
  dplyr::left_join(map) %>% 
  dplyr::left_join(chembl %>%
                     dplyr::select(-Name, -Synonyms),
                   by = "ChEMBL ID") %>%
  tidyr::separate_rows(Smiles, sep = ",\ ") %>% 
  tidyr::separate_rows(SMILES, sep = ",\ ") %>% 
  dplyr::mutate(Smiles = toupper(Smiles),
                SMILES = toupper(SMILES)) %>% 
  dplyr::group_by(cm_ix) %>% 
  dplyr::mutate(`ChEMBL ID` = paste0(sort(unique(`ChEMBL ID`)), collapse = ";")) %>% 
  tidyr::pivot_longer(cols = c("SMILES", "Smiles"),
                      names_to = "dummy", values_to = "SMILES") %>%
  dplyr::select(-dummy) %>% 
  dplyr::mutate(SMILES = paste0(sort(unique(SMILES[SMILES != ""])), collapse= ";"),
                `Inchi Key` = paste0(sort(unique(`Inchi Key`)), collapse= ";"),
                Type = paste0(sort(unique(Type)), collapse= ";"),
                `Max Phase` = max(`Max Phase`, na.rm = T)) %>% 
  dplyr::mutate(`Max Phase` = ifelse(is.finite(`Max Phase`), `Max Phase` , NA)) %>%
  dplyr::ungroup() %>% 
  dplyr::rename(`Phase(RepHub)` = Phase,
                `Phase(ChEMBL)` = `Max Phase`) %>%
  dplyr::distinct() %>% 
  dplyr::select(cm_ix, Drug.Name, Synonyms, `ChEMBL ID`, IDs, referenced_datasets, 
                everything()) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(`Phase(ChEMBL)` = ifelse(`Phase(ChEMBL)` == 4, "Approved",
                                         ifelse(`Phase(ChEMBL)` == 3, "Phase 3",
                                                ifelse(`Phase(ChEMBL)` == 2, "Phase 2",
                                                       ifelse(`Phase(ChEMBL)` == 1, "Phase 1",
                                                              ifelse(`Phase(ChEMBL)` == 0.5, "Early Phase",
                                                                     ifelse(`Phase(ChEMBL)` == -1, "Unknown",
                                                                            "Precilinical")))))))

# sanity check
compound.metadata.expanded %>% 
  dplyr::group_by(cm_ix) %>% 
  dplyr::filter(n() > 1)

compound.metadata.expanded %>% 
  dplyr::select(-cm_ix) %>% dplyr::distinct() %>% 
  write_csv("data/compound_metadata_expanded.csv")


# Question for Phil: Should we combine rows with the same ChEMBL ids or InChI keys? What about the names?

compound.metadata %>% 
  dplyr::distinct(SMILES, cm_ix) %>%
  dplyr::summarise(n = sum(SMILES != ""))

compound.metadata.expanded %>% 
  dplyr::distinct(SMILES, cm_ix) %>%
  dplyr::summarise(n = sum(SMILES != ""))
  
# Seems like SMILES are slightly different half of the times.
compound.metadata.expanded %>% 
  dplyr::distinct(cm_ix, Drug.Name, SMILES) %>% 
  tidyr::separate_rows(SMILES, sep = ";") %>% 
  dplyr::filter(SMILES != "") %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(cm_ix) %>% 
  dplyr::filter(n() == 1) %>%
  .$cm_ix %>% unique() %>% length
