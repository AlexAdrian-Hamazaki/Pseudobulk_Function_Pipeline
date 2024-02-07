extractCTAffiliation <- function(loOntology) {
  ###
  # For each Ontology term in a list of ontology terms in an ontoloy indexo bject
  # Scan each term for a CTAffiliation
  # Returns: List of bools/None based on a terms CTAffiliation
  ###
  loCTAffiliations <- c()
  loGOIDs = loOntology$id
  loPropValues = loOntology$property_value
  for (index in 1:length(loGOIDs)) {
    #print(index)
    
    CTAffiliation = returnCTAffiliation(loPropValues[index])
    loCTAffiliations <- c(loCTAffiliations, CTAffiliation)
    
  }
  
  names(loCTAffiliations) <- loGOIDs
  
  return (loCTAffiliations)
}

returnCTAffiliation <- function(loPropertyValues) {
  ###
  # Searches the list of property values for 1 GO ontolgoy term
  # If it finds a CTA affiliation return that
  # Otherwise return None
  ###
  for (PropertyValue in loPropertyValues[[1]]) {
    
    if (identical(PropertyValue, character(0))) {
      #print("No Property Value")
      return("None")
    }else if (startsWith(PropertyValue, "CTAffiliation")) {
      #print(paste("CTAffiliation Found: ", PropertyValue ))
      return (PropertyValue)
    }
  }
  #print("CTAffiliation not Found")
  return ("None")
}



giveCTAVAffiliationsToOntology <- function(ontology, loCTAffiliations) {
  # Give the ontology object a new "CTAffiliation" parameter.
  # Set the CTAffiliation according to loCTAffiliations
  
  # True if Go term is CT related
  # False if GO term is not CT related
  # None if GO term has not been curated, or a value has not been propagated to it yet
  ontology[["CTAffiliations"]] <- loCTAffiliations
  return (ontology)
}

processOntology <- function(ontologyIndex){
  # Process a freshly loaded ontology so it has the CTAffiliation property
  # Returns an ontology index with CTAffiliation property
  
  # Extract all CTAffiliations
  loCTAffiliations <- extractCTAffiliation(ontologyIndex)
  
  # Fix CTAffiliations
  loCTAffiliations <- lapply(X = loCTAffiliations, fixCTAValue)
  
  # add CTAffiliations to eaech GO term in the ontologyindex object
  ontology <- giveCTAVAffiliationsToOntology(ontology = ontology, loCTAffiliations = loCTAffiliations)
  
  return (ontology)
}

fixCTAValue <- function(CTAValues) {
  # Parse the output so that the value is either True or False
  # For use in an lapply
  
  result <- str_match(CTAValues, '".*"')
  
  extracted_value <- substr(result[1, 1], 2, nchar(result[1, 1]) - 1)
  
  return (extracted_value)
}



giveCTAtoGOID <- function(GOID, CTAValue, ontology){
  # for 1 goID, CTAvalue pair. Find that goID in the ontologyIndex, and then change that goIDs CTAffiliation
  
  # GOID: str
  # CTAValue: str -> "true" or "false"
  # ontology -> ontologyindex object
  
  
  ontology$CTAffiliations[GOID] <- CTAValue
  
  return (ontology)
}

addCuratedData <- function(curatedDF, ontology) {
  ###
  # For a curated dataframe that has columns "ID", and "CTAffiliation"
  # Search the ontology index for each ID (go id) in the ontology index, and replace that term's CTAffiliation with
  # the CTAffiliation in the dataframe. Only replace if CTAffiliation is == NA
  ###
  
  # Iterate over each row
  for (i in 1:nrow(curatedDF)) {
    row <- curatedDF[i, ]  # Access row i
    
    goID <- row[,"ID"]
    CTAffiliation <- row[,"CTAffiliation"]
    
    ontology <- giveCTAtoGOID(GOID = goID, CTAValue = CTAffiliation, ontology = ontology)
  }
  
  return (ontology)
  
}


getCTAValues <- function(loCTAffiliations) {
  ###
  # Returns the GO id Names and values of the GO ids with CTAffiliations
  # loCTAffiliations: List of CTAffiliations (ontology$CTAffiliations)
  ###
  return (loCTAffiliations[!is.na(loCTAffiliations)])
}






giveTermToDescendents <- function(lodescendents, CTAValue, ontology){
  # for each descendent in a list of descendents. Give the descendent false or True
  # save that result to the ontologyindex object
  
  # loDescendent -> lostrings. strings are GO Ids
  # term -> string
  # ontology -> ontologyindex object
  CTAbool <- CTAValue[[1]]
  
  for (descendent in lodescendents) { # For each descendent. change that descendent's CTAffiliation to the CTAvalue (parents CTAffiliation)
    ontology$CTAffiliations[descendent] <- CTAbool
  }
  
  return (ontology)
}

propagateCTAForGOTerm <- function(loCTAValues, ontology) {
  # Take a list of CTA Values. For each GO id in the list, get the descendents
  # Then for each descendent, give it the CTA bool value from the list of CTA values
  
  # Returns an ontologyindex object where daughter terms have the same CTA values as parent terms
  
  # CTAValues: List of GOIds:CTAValues (bool). 
  # ontology: ontologyindex object
  
  
  # Get descendents for a CTAValue
  
  for (i in seq_along(loCTAValues)) {
    goID <- names(loCTAValues)[i]
    CTAValue <- loCTAValues[[i]]
    
    # Get Descendents of this CTAValue
    lodescendents <- ontologyIndex::get_descendants(ontology = ontology, 
                                                    roots = goID,
                                                    exclude_roots = TRUE)
    
    #Give each of the descendents the same CTA value in the ontology object
    ontology <- giveTermToDescendents(lodescendents, CTAValue, ontology) 
    
    #loCTAffiliations <- extractCTAffiliation(ontology)
    
    # Prints how many ontology terms have CTAffiliations
    # print(paste(length(getCTAValues(ontology$CTAffiliations))), "terms have CTAffiliations")
    
    
  }
  # return new instance of ontology object
  return (ontology)
}




countCTAffiliations <- function(ontology){
  # Count the number of True, False and None CTAffiliations
  
  #ontology -> ontologyindex object
  
  return (table(unlist(newOntology$CTAffiliations)))
}


splitOntologyNameSpaces <- function(ontology) {
  # Split ontology into two ontologies, BP and MF
  
  ontologyDF <- as.data.frame(ontology)
  
  MF <- ontologyDF %>%
    filter(namespace == "molecular_function")
  
  BP <- ontologyDF %>%
    filter(namespace == "biological_process")
  
  return (list(MF, BP))
}
