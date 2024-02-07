#!/usr/bin/env Rscript

# Curate everything with Regex

library(ontologyIndex)
library(igraph)
library(tidyverse)
library(stringr)
library(readr)

source("functions.R")

args <- commandArgs(trailingOnly=TRUE)

go_obo_path <- args[[1]]
curated_regs <- args[[2]]
good_terms <- args[[3]]
bad_terms <- args[[4]]

print(args)
ontology <- ontologyIndex::get_OBO(go_obo_path,
                                  extract_tags = "everything",
                                  propagate_relationships = c("is_a", "part_of"))


# Add CTAffiliations column
new_ontology <- processOntology(ontology)
new_ontology_df <- as.data.frame(new_ontology) %>% 
  select(c('id', 'name', 'def', 'CTAffiliations'))
colnames(new_ontology_df) <- c("ID", "Name", "Description", "CTAffiliation")
# Assign the index numbers as row names
rownames(new_ontology_df) <- 1:nrow(new_ontology_df)
# Replace "NA" with NA
new_ontology_df$CTAffiliation[new_ontology_df$CTAffiliation == "NA"] <- NA

new_ontology_df <- new_ontology_df %>%
  filter(!is.na(Name))

stopifnot(new_ontology$CTAffiliations[3] == 'false')


MF_descendents <- get_descendants(ontology, roots = c("GO:0003674"), exclude_roots = TRUE)
BP_descendents <- get_descendants(ontology, roots = c("GO:0008150"), exclude_roots = TRUE)


getNameDescription <- function(descendent, ontology) {
  # Get the name and description for 1 descendent 
  
  #descendent : go_id (str)
  # ontology: ontology index
  
  # Returns, character vector of descendent, name, description
  
  name <- get_term_property(ontology, 'name', term = descendent, as_names = FALSE)
  description <- get_term_property(ontology, 'def', term = descendent, as_names = FALSE)
  CTAffiliation <- get_term_property(ontology, 'CTAffiliations', term = descendent, as_names = FALSE)
  
  
  return (c(descendent, name, description, CTAffiliation))
}

MFWithInfo <-lapply(MF_descendents , FUN = getNameDescription, ontology = new_ontology)
MFWithInfo_df <- data.frame(do.call(rbind, MFWithInfo))
colnames(MFWithInfo_df) <- c("ID", "Name", "Description", "CTAffiliation")

BPWithInfo <-lapply(BP_descendents, FUN = getNameDescription, ontology= new_ontology)
BPWithInfo_df <- data.frame(do.call(rbind, BPWithInfo))
colnames(BPWithInfo_df) <- c("ID", "Name", "Description", "CTAffiliation")


detectTermToDescription <- function(ontologyDF, term){
  ### 
  # Returns a bool vector corresponding to if the term is found in the ontology DF Name
  ###
  loDetected <-  stringr::str_detect(str_to_lower(ontologyDF$Name), str_to_lower(term))
  return (loDetected)
}

addTrueCTAffiliation <- function(ontologyDF, loDetected) {
  ###
  # for the ontology dataframe, make the CTAffiliation column equal to "true" 
  # if the loDetected has TRUE in that position
  ###
  
  for (i in seq_along(loDetected)) {

    if (loDetected[i] == TRUE) {
      ontologyDF[[i, "CTAffiliation"]] <- "true"
    }
  }
  return (ontologyDF)
}

addFalseCTAffiliation <- function(ontologyDF, loDetected) {
  ###
  # for the ontology dataframe, make the CTAffiliation column equal to "true" 
  # if the loDetected has false in that position
  ###
  for (i in seq_along(loDetected)) {

    if (loDetected[i] == TRUE && is.na(ontologyDF[i, "CTAffiliation"])) {
      ontologyDF[i, "CTAffiliation"] <- "false"
    }
  }
  return (ontologyDF)
}


recurseDetection <- function(ontologyDF, terms, n = 0, bool){
  ###
  #update an ontology dataframe with true or false
  # Searches the ontology df's name for a regex string. If that string is found, then make that go term's CTAffiliation equal to a bool
  # if bool == TRUE, then we are using this function for a good terms list, make the CTaffilaition = "true"
  # if bool == FALSE, then we are using this function for a bad terms list, make the CTAffiliation = "false"
  ###
  n <- n + 1
  if (n > length(terms)) {
    return(ontologyDF)
  } else {

    if (bool == TRUE) {
      
      loDetected <- detectTermToDescription(ontologyDF, term = terms[n])

      ontologyDF_new <- addTrueCTAffiliation(ontologyDF = ontologyDF, loDetected = loDetected)
      print(paste(sum(is.na(ontologyDF$CTAffiliation)) - sum(is.na(ontologyDF_new$CTAffiliation)), "Is the number of na terms altered by", terms[n]))
      
      return(recurseDetection(ontologyDF = ontologyDF_new, terms = terms, n = n, bool = bool))
      
    } else if (bool == FALSE) {
      
      loDetected <- detectTermToDescription(ontologyDF, term = terms[n])
      ontologyDF_new <- addFalseCTAffiliation(ontologyDF = ontologyDF, loDetected = loDetected)
      print(paste(sum(is.na(ontologyDF$CTAffiliation)) - sum(is.na(ontologyDF_new$CTAffiliation)), "Is the number of na terms altered by", terms[n]))
      
      return(recurseDetection(ontologyDF = ontologyDF_new, terms = terms, n = n, bool = bool))
    }
  }
}

updateOntology <- function(descendendsWithInfo_df, good_terms, bad_terms) {
  ontologyUpdated<- recurseDetection(descendendsWithInfo_df, terms = good_terms, bool = TRUE)
  ontologyUpdated<- recurseDetection(ontologyUpdated, terms = bad_terms, bool = FALSE)
  
  return (ontologyUpdated)
}

checkRemainingNAs <- function(ontologyDF){
  nas <- filter(ontologyDF, is.na(CTAffiliation))
  print(nas$ID)
  print(nas$Name)
  print(nrow(nas))
}

filterNamesOfOneREGEX <- function(updatedOntology, term) {
  ### 
  #  filter dataframe for only those where a stringr regex match for a term in the rows description matches
  # For purposes of exploration/ensuring your terms are correctly identifying false/true cases
  ###
  
  filtered <- updatedOntology %>%
    filter(stringr::str_detect(str_to_lower(Name), str_to_lower(term)))
  
  return (filtered)
}


good_terms <- readr::read_lines(good_terms)
bad_terms <- readr::read_lines(bad_terms)

#updatedOntology <- updateOntology(new_ontology_df, good_terms = good_terms, bad_terms= bad_terms)
MFupdatedOntology <- updateOntology(MFWithInfo_df, good_terms = good_terms, bad_terms= bad_terms)
sum(MFupdatedOntology$CTAffiliation=="true", na.rm = TRUE)
sum(MFupdatedOntology$CTAffiliation=="false", na.rm = TRUE)
checkRemainingNAs(MFupdatedOntology)



BPupdatedOntology <- updateOntology(BPWithInfo_df, good_terms = good_terms, bad_terms= bad_terms)

sum(BPupdatedOntology$CTAffiliation=="true", na.rm = TRUE)
sum(BPupdatedOntology$CTAffiliation=="false", na.rm = TRUE)

checkRemainingNAs(BPupdatedOntology)

write.csv(MFupdatedOntology, file = "MF_ctRelatedness.csv",  row.names =FALSE)
write.csv(BPupdatedOntology, file = "BP_ctRelatedness.csv",  row.names =FALSE)




#checkNames <- filterNamesOfOneREGEX(updatedOntology = BPupdatedOntology, term = "biosynthetic process")
# print(paste(checkNames$Name, checkNames$CTAffiliation))
# print(checkNames$CTAffiliation)
# 
# false <-  checkNames[checkNames$CTAffiliation == "false", ]
# print(paste(false$Name, false$CTAffiliation))
# 
# true<-  checkNames[checkNames$CTAffiliation == "true", ]
# print(paste(true$Name, true$CTAffiliation))

#write.csv(updatedOntology, file = "bin/ontologyPropagation/curatedAllRegex.csv")
