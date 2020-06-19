# NEXT FUNCTION CREATED BY LEE F STANISH. TESTING USE ONLY FOR NOW #
#### FUNCTION downloadSequenceMetadataRev ####
downloadSequenceMetadataRev <- function(sites='all', startYrMo, endYrMo, 
                                        targetGene= "", dpID = "DP1.10108.001", dir="") {
  # author: Lee Stanish
  # date: 2020-06-16
  # function loads soil marker gene sequencing metadata for target gene, site(s) and date(s)
  # option to download output by providing a valid output directory
  # sites: character vector of valid site ID's, or 'all' for all sites
  # targetGene: '16S' or 'ITS'
  # startYrMo: start date, format YYYY-MM
  # endYrMo: end date, format YYYY-MM
  # dpID: NEON data product of interest. Default is soil marker gene sequences, and currently code only works for this dpID
  # dir (optional): If a local copy of the filtered metadata is desired, provide path to output dir
  
  library(neonUtilities)
  library(plyr)
  library(dplyr)
  
  # check valid data values entered
  ## validate dpID ##
  if(!grepl("DP1", dpID) | !grepl('\\.001', dpID) | !grepl('10108|20280|20282', dpID)) {
    print("Invalid Data Product ID: must follow convention 'DP1.[5-digit value].001' and must be a marker genes data product ID")
    return(NULL)
  } else {
    dpID <- dpID
  }
  
  # validate target gene
  if(!grepl("16S|ITS", targetGene)) {
    print("Invalid targetGene: must be either '16S' or 'ITS'")
    return(NULL)
  } else {
    targetGene <- targetGene
  }
  
  # validate site(s)
  terrSiteList <- c("all","HARV","SCBI","OSBS","GUAN","UNDE","KONZ","ORNL","TALL","WOOD","CPER","CLBJ","YELL","NIWO",
                    "SRER","ONAQ","WREF","SJER","TOOL","BONA","PUUM","BART","BLAN","SERC","SCBI","DSNY","JERC","LAJA",
                    "TREE","STEI","KONA","UKFS","MLBS","GRSM","LENO","DELA","NOGP","DCFS","STER","RMNP","OAES","MOAB",
                    "JORN","ABBY","TEAK","SOAP","BARR","DEJU","HEAL")
  if(any(sites %in% terrSiteList)==FALSE){
    print("Invalid site(s): must be a valid NEON site or 'all'")
    return(NULL)
  } else {
    sites <- sites
  }
  
  print("loading metadata...")
  mmgL1 <- loadByProduct(dpID, sites, package = 'expanded', check.size = F, startdate = startYrMo, enddate = endYrMo) # output is a list of each metadata file
  
  
  # for target data product and targetGene: extract lists into data.frames
  if(grepl("10108", dpID)) {
    if(targetGene=="16S") {
      print("filtering to 16S data")
      seq <- mmgL1$mmg_soilMarkerGeneSequencing_16S
      raw <- mmgL1$mmg_soilRawDataFiles
      
    } else {
      print("filtering to ITS data")
      seq <- mmgL1$mmg_soilMarkerGeneSequencing_ITS
      raw <- mmgL1$mmg_soilRawDataFiles
    }
  }
  
  if(grepl("20280", dpID)) {
    if(targetGene=="16S") {
      print("filtering to 16S data")
      seq <- mmgL1$mmg_benthicMarkerGeneSequencing_16S
      raw <- mmgL1$mmg_benthicRawDataFiles
    } else {
      print("filtering to ITS data")
      seq <- mmgL1$mmg_benthicMarkerGeneSequencing_ITS
      raw <- mmgL1$mmg_benthicRawDataFiles
    }
  }
  
  if(grepl("20282", dpID)) {
    if(targetGene=="16S") {
      print("filtering to 16S data")
      seq <- mmgL1$mmg_swMarkerGeneSequencing_16S
      raw <- mmgL1$mmg_swRawDataFiles
    } else {
      print("filtering to ITS data")
      seq <- mmgL1$mmg_swMarkerGeneSequencing_ITS
      raw <- mmgL1$mmg_swRawDataFiles
    }
  } 
  
  # convert factors to characters (bug in output of loadByProduct)
  i <- sapply(seq, is.factor)
  seq[i] <- lapply(seq[i], as.character)
  j <- sapply(raw, is.factor)
  raw[j] <- lapply(raw[j], as.character)
  
  
  # Join sequencing metadata with raw data files metadata
  if(targetGene=="16S") {
    if(any(grepl("ITS", raw$rawDataFileName))) {
      rawCleaned <- raw[-grep("ITS", raw$rawDataFileName), ]  
    } else {
      rawCleaned <- raw
    }
    joinedTarget <- left_join(rawCleaned, seq, by=c('dnaSampleID', 'sequencerRunID', 'internalLabID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
  } else if(targetGene=="ITS") {
    if(any(grepl("16S", raw$rawDataFileName))) {
      rawCleaned <- raw[-grep("16S", raw$rawDataFileName), ]  
    } else {
      rawCleaned <- raw
    }
    joinedTarget <- left_join(rawCleaned, seq, by=c('dnaSampleID', 'sequencerRunID', 'internalLabID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
  }
  
  # download local copy if user provided output dir path
  if(dir != "") {
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    write.csv(out, paste0(dir, "/mmg_soilMetadata_", targetGene, "_", Sys.Date(), ".csv"),
              row.names=F)
    print(paste0("metadata downloaded to: ", dir, "/mmg_soilMetadata_", targetGene, "_", Sys.Date(), ".csv") )
  }
  return(out)
  
  ### END FUNCTION ###
}
