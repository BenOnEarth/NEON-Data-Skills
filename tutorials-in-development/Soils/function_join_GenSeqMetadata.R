#### FUNCTION downloadSequenceMetadataRev ####
downloadSequenceMetadataRev <- function(sites='all', startYrMo, endYrMo, 
                                        targetGene= "", dpID = "DP1.10108.001",
                                        dir = "", rawDownload = FALSE) {
  # authors: Lee Stanish, Ben Shetterly
  # date: 2020-06-29
  # function loads marker gene sequencing metadata for target data product, gene, site(s) and date(s)
  # option to download output by providing a valid output directory
  # option to download raw sequence archive files
  # sites: character vector of valid site ID's, or 'all' for all sites
  # targetGene: '16S' or 'ITS'
  # startYrMo: start date, format YYYY-MM
  # endYrMo: end date, format YYYY-MM
  # dpID: NEON data product of interest. Default is soil marker gene sequences, but should work for surface water and benthic microbes
  # dir (optional): If a local copy of the filtered metadata and raw sequences is desired, provide path to output dir
  # rawDownload: set to TRUE to download raw microbe marker sequence data. NOTE: THESE ARE VERY LARGE FILES
  library(neonUtilities)
  library(dplyr)
  library(httr)
  
  # check valid argument values entered
  # validate dpID
  if(!grepl("^DP1", dpID) | !grepl('\\.001$', dpID) | !grepl('10108|20280|20282', dpID)) {
    print("Invalid Data Product ID: must follow format 'DP1.#####.001' and must be a marker gene sequences product ID")
    return(NULL)
  } else {
    dpID <- dpID
    dpInfo<-neonUtilities::getProductInfo(dpID = dpID)
  }
  # validate target gene
  if(!grepl("^16S$|^ITS$", toupper(targetGene))) {
    print("Invalid targetGene: must be either '16S' or 'ITS'")
    return(NULL)
  } else targetGene <- toupper(targetGene)
  # validate site data availability
  # Use product info to get DP-specific sites availability
  if("all" %in% tolower(sites)){
    sites <- "all"
  } else {
    availSites <- dpInfo$siteCodes$siteCode
    if(!all(sites %in% availSites)){
      print("Invalid site(s): must be NEON site code(s) associated with the data product, or 'all'")
      return(NULL)
    } else sites <- toupper(sites)
  }
  # Load expanded NEON data product containing metadata into a list using loadByProduct
  print(paste("Downloading metadata for",dpID,dpInfo$productName,"at site(s)", paste(sites,collapse=" "),"start month",startYrMo,"and end",endYrMo,sep=" "))
  mmgL1 <- neonUtilities::loadByProduct(dpID = dpID, site = sites, package = 'expanded', check.size = F, startdate = startYrMo, enddate = endYrMo)
  # use table_types to get the relevant table names corresponding to the DPID, 
  # which have similar but not identical names across DPs.
  tablenames <- neonUtilities::table_types%>%
    filter(grepl(dpID,productID),
           grepl(paste("DnaExtraction$|RawDataFiles$|GeneSequencing_",targetGene,"$",sep = ""),tableName))
  # extract the variables table
  vartable <- mmgL1[[grep("^variables_",names(mmgL1))]]
  
  # assign column classes using neonUtilities readTableNEON and variables table, applied to selected tables
  md<-sapply(mmgL1[tablenames$tableName], neonUtilities::readTableNEON, varFile = vartable, simplify = FALSE)
  
  # extract lists into data.frames
  raw <- md[[grep("RawDataFiles$",names(md))]]
  seq <- md[[grep(paste0("GeneSequencing_",targetGene,"$"),names(md))]]
  ext <- md[[grep("DnaExtraction$",names(md))]]
  
  # character vector of join columns used for metadata on all Microbe marker gene sequence DPs
  join_common <- c("domainID","siteID","namedLocation","collectDate","dnaSampleID","internalLabID")
  # filter soil extraction table to remove non-relevant rows (suggested per DP user guide). 
  # Include plotID column for joining soils extraction sequencing data
  if(dpID=="DP1.10108.001"){
    extCleaned <- ext%>%
      filter(grepl("marker gene|marker gene and metagenomics",sequenceAnalysisType))
    join1_cols<-c(join_common,"plotID")
  } else {
    extCleaned <- ext
    join1_cols <- join_common
  }
  
  # Join extraction with sequencing
  j1 <- left_join(extCleaned, seq, by=join1_cols, suffix=c(".ext",".seq"))
  
  # clean raw data table to remove uneeded files with nontarget gene data
  if(toupper(targetGene)=="16S") nontarget<-"ITS" else if(toupper(targetGene)=="ITS") nontarget<-"16S"
  rawCleaned <- raw%>%
    filter(!grepl(paste0("_",nontarget),rawDataFileName))
  # join extraction_sequencing table with raw data file table
  out <- left_join(j1, rawCleaned, by=c(join_common, "sequencerRunID"), suffix=c(".ext_seq",".raw"))
  
  # download local copy if user provided output dir path
  if(dir != "") {
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    # give informative name to metadata file and display to user
    md_filename<-paste0("mmg_metadata_", substr(dpID,5,9),"_", targetGene, "_", format(Sys.Date(),"%Y%m%d"), ".csv")
    write.csv(out, file = paste0(dir, "/", md_filename), row.names=FALSE)
    print(paste0("metadata downloaded to: ", dir, "/", md_filename))
    # prepare and download raw sequence data files
    if(rawDownload){
      # Some raw data files are not available (404 response) even using
      # the URLs given in _RawDataFiles table. Initialize a table with all the unique URLs/filenames, 
      print("Checking for raw sequence data availability")
      rawfile_status<-out%>%
        select(rawDataFilePath, rawDataFileName)%>%
        mutate(rawdata_available=NA)%>%
        distinct
      # loop to retrieve each header status code and fill table with data availability
      for(u in rawfile_status$rawDataFilePath){
        r<-httr::HEAD(u)
        if(httr::status_code(r) == 200){ rawfile_status[rawfile_status$rawDataFilePath==u,]$rawdata_available <- TRUE
        } else rawfile_status[rawfile_status$rawDataFilePath==u,]$rawdata_available <- FALSE
      }
      # Alert user if some raw sequence data is not available
      if(!all(rawfile_status$rawdata_available)){
        print("NOTE: The following raw sequence files are not available for download:")
        print(rawfile_status[!rawfile_status$rawdata_available,]$rawDataFileName)
      }
      # for zipsByURI, variables file must be 
      # in the same directory as uri data and named "variables.csv"
      if(!file.exists(paste0(dir,"/variables.csv"))){
        write.csv(vartable, file = paste0(dir, "/variables.csv"), row.names=FALSE)
      }
      # get raw data table name, since this depends on data product
      # and zipsByURI needs it to match the variables table.
      raw_tablename<-vartable%>%
        filter(tolower(dataType) == "uri", grepl("RawDataFiles$", table))%>%.[[1,"table"]]
      # use filter to remove records with unavailable raw data files
      # then prepare it for zipsByURI by writing to a file
      write.csv(filter(rawfile_status, rawdata_available), file = paste0(dir, "/", raw_tablename,".csv"), row.names = FALSE)
      # download raw sequence data to the given dir/ECS_zipFiles
      zipsByURI(filepath = dir, check.size = TRUE, unzip = FALSE, saveZippedFiles = TRUE)
    }
  } else if(rawDownload) print("No save directory provided. To download raw sequence data, specify a save path with 'dir = '") # Warns user if no dir provided but rawDownload = TRUE
  return(out)
  
  ### END FUNCTION ###
}