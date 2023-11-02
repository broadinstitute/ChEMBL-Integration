library(tidyverse)
library(taigr)
library(useful)


# Load the existing compound meta-data
compound.metadata <- load.from.taiga(data.name='compound-metadata-de37', data.version=14, data.file='compound_metadata')

# Iron-out the confusing ID's, Smiles, Drug Name's and Synonyms
# Number of rows drops from 7692 to 7314
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
  dplyr::mutate(ix3 = ifelse(SMILES == "", NA, ix3)) %>% 
  dplyr::group_by(cm_ix) %>% 
  dplyr::mutate(ix1 = min(ix1, na.rm = T), ix2 = min(ix2, na.rm = T), ix3 = min(ix3, na.rm = T)) %>% 
  dplyr::group_by(ix1) %>% 
  dplyr::mutate(cm_ix = min(cm_ix, na.rm = T), ix2 = min(ix2, na.rm = T), ix3 = min(ix3, na.rm = T)) %>% 
  dplyr::group_by(ix2) %>% 
  dplyr::mutate(ix1 = min(ix1, na.rm = T), cm_ix = min(cm_ix, na.rm = T), ix3 = min(ix3, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(ix3 = dplyr::group_indices(., cm_ix, ix3)) %>% 
  dplyr::group_by(ix3) %>% 
  dplyr::mutate(ix1 = min(ix1, na.rm = T), ix2 = min(ix2, na.rm = T), cm_ix = min(cm_ix, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-ix1, -ix2, -ix3) %>% 
  dplyr::distinct() %>% 
  tidyr::separate_rows(MOA, sep = ";") %>% 
  tidyr::separate_rows(repurposing_target, sep = ";") %>% 
  tidyr::separate_rows(Indication, sep = ";") %>% 
  tidyr::separate_rows(`Disease Area`, sep = ";") %>% 
  tidyr::separate_rows(InChIKey, sep = ";") %>% 
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
                InChIKey = paste0(sort(unique(InChIKey[InChIKey != ""])), collapse= ";"),
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
  tidyr::pivot_longer(cols = c("Name", "Synonyms"),
                      names_to = "dummy", values_to = "Drug.Name") %>%
  dplyr::select(-dummy) %>% dplyr::filter(Drug.Name != "") %>% 
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
compound.metadata.expanded <- compound.metadata %>% 
  dplyr::left_join(map) %>% 
  dplyr::left_join(chembl, by = "ChEMBL ID") %>% 
  tidyr::separate_rows(Drug.Name, sep = ";") %>% 
  tidyr::separate_rows(Name, sep = ";") %>% 
  tidyr::separate_rows(Synonyms.x, sep = ";") %>% 
  tidyr::separate_rows(Synonyms.y, sep = fixed("\\|")) %>%  
  tidyr::separate_rows(Synonyms.y, sep = fixed(", ")) %>%  
  tidyr::pivot_longer(cols = c("Name", "Synonyms.x", "Synonyms.y", "Drug.Name"),
                      names_to = "dummy", values_to = "Drug.Name") %>%
  dplyr::select(-dummy) %>% 
  dplyr::filter(Drug.Name != "") %>% 
  dplyr::mutate(Drug.Name = toupper(Drug.Name)) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(cm_ix) %>% 
  dplyr::mutate(Synonyms = paste0(sort(unique(Drug.Name)), collapse = ";"),
                `ChEMBL ID` = paste0(sort(unique(`ChEMBL ID`)), collapse = ";")) %>% 
  dplyr::top_n(1, Drug.Name) %>%  
  dplyr::ungroup() %>% dplyr::distinct() %>% 
  tidyr::pivot_longer(cols = c("SMILES", "Smiles"),
                      names_to = "dummy", values_to = "SMILES") %>%
  dplyr::select(-dummy) %>% 
  dplyr::group_by(cm_ix) %>% 
  dplyr::mutate(SMILES = paste0(sort(unique(SMILES[SMILES != ""])), collapse= ";")) %>%
  dplyr::ungroup() %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(cm_ix) %>% 
  tidyr::pivot_longer(cols = c("InChIKey", "Inchi Key"),
                      names_to = "dummy", values_to = "InChIKey") %>%
  dplyr::select(-dummy) %>% 
  dplyr::group_by(cm_ix) %>% 
  dplyr::mutate(InChIKey = paste0(sort(unique(InChIKey[InChIKey != ""])), collapse= ";"),
                Type = paste0(sort(unique(Type[Type != ""])), collapse= ";")) %>%
  dplyr::mutate(`Max Phase` = max(`Max Phase`, na.rm = T)) %>% 
  dplyr::mutate(`Max Phase` = ifelse(is.finite(`Max Phase`), `Max Phase` , NA)) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(`Phase(RepHub)` = Phase,
                `Phase(ChEMBL)` = `Max Phase`) %>%
  dplyr::select(-cm_ix) %>% 
  dplyr::select(Drug.Name, Synonyms, `ChEMBL ID`, IDs, referenced_datasets, 
                everything()) %>% 
  dplyr::distinct() 



compound.metadata.expanded %>% 
  write_csv("data/compound_metadata_expanded.csv")

