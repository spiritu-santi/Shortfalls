library(tidyverse)
library(magrittr)
library(here)
library(sp)
library(raster,exclude = "select")
library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)
library(showtext)
library(patchwork)
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(ape)

###### PRE-PROCESSING: SETTING UP DATA FILES RETRIEVED FROM GBIF #####
###### USE BASH TO REDUCE FILE SIZE BY SELECTING COLUMNS
### bash commands are commented out!
# unzip 0188599-210914110416597.zip
# rename unzipped file to Traqueos_NeoTropics_raw.csv
# awk -F"\t" '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13"\t"$22"\t"$23"\t"$33}' Traqueos_NeoTropics_raw.csv > Traqueos_NeoTropics_COR.csv
# wc -l Traqueos_NeoTropics_raw.csv
# wc -l Traqueos_NeoTropics_COR.csv

###### Keep records identified to species level or below ########


## Utilities 
total <- function(x) { 
  x %>% summarise(n_tot = n()) %>% pull(n_tot) -> n_tot
  return(n_tot)
}
total_rbcL <- function(x) { 
  x %>% count(rbcL) %>% filter(!is.na(rbcL)) %>% pull(n) -> n_rbcl
  if(length(n_rbcl) == 0) n_rbcl <- 0
  return(n_rbcl)
}
get_mostrecent_year <- function(x){
  x %>% pull(year) %>% .[which.max(.)] -> result
  return(result)
}
roundUp <- function(x,to=10) {to*(x%/%to + as.logical(x%%to))}
font_add_google(name="EB Garamond")
showtext_auto()

#### FUNCTIONS ####
gbif_pow_exactjoin <- function(gbif_data = "Traqueos_NeoTropics_COR.csv",
                               checklist = "wcvp_names.txt",
                               taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),
                               output_file="Joined_base.Rdata") { 
aa <- data.table::fread(here("data",gbif_data)) %>% as_tibble()
aa %<>% mutate(ID=1:nrow(.)) %>% relocate(ID,.before=1) %>% separate(scientificName,into=c("G","S","extra"),extra = "merge",sep=" ") %>% unite("binomial",c(G,S),sep="_",remove=F) %>% unite("scientificName",c(G,S,extra),sep=" ",remove=T) %>% filter(species!="")
list_names <- data.table::fread(here("wcvp_2022",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ",remove = F) %>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble() %>%  filter(taxon_status %in% all_of(taxon_status_cats)) %>% unite("binomial",c(genus,species),sep="_",remove = F) %>% filter(accepted_plant_name_id!="")

## EXACT JOIN 
joined <- aa %>% #rename("taxon_name"=scientificName) %>% 
  left_join(.,list_names,by = "scientificName")
cat("Searching for duplicate IDs due to multiple matches!","\n")
joined %>% pull(ID) %>% duplicated() -> dups
joined$ID[which(dups)] -> id_dups
cat("Found",length(id_dups),"duplicated IDs","\n")
joined %>% filter(ID %in% all_of(id_dups)) -> to_correct
joined <- joined[which(!joined$ID %in% id_dups),]
pow_distributions="wcvp_distribution.txt"
pow_dist <- data.table::fread(here("wcvp_2022",pow_distributions),sep="|",quote="") %>% as_tibble() %>% #left_join(.,list_names,by="accepted_db_id") %>% 
  filter(introduced==0)
joined %<>% mutate(Dist_correction=T)
to_correct %<>% mutate(Dist_correction=T)
cat("First step: resolving by distribution","\n")
for(i in 1:nrow(to_correct)){
  cat("processing",i,"\r")
pow_dist %>% filter(plant_name_id==to_correct$accepted_plant_name_id[i]) %>% pull(continent_code_l1) -> tar_dist
tar_dist %in% c(1:6,9) -> res
if(sum(res)!=0) {to_correct$Dist_correction[i] <- FALSE}
}
to_correct %<>% filter(Dist_correction)
rbind(to_correct,joined) -> joined
joined %>% pull(ID) %>% duplicated() -> dups
joined$ID[which(dups)] -> id_dups
cat("Remaining",length(id_dups),"duplicated IDs","\n")
joined %>% filter(ID %in% all_of(id_dups)) -> to_correct
joined <- joined[which(!joined$ID %in% id_dups),]
to_correct %>% filter(homotypic_synonym=="") %>% group_split(ID) -> to_correct
cat("Second step: resolve manually","\n")
for(i in 1:length(to_correct)){
  cat("processing",i,"\n")
  if(nrow(to_correct[[i]])==1) next
to_correct[[i]] %>% select(scientificName,Accepted_Name,homotypic_synonym,accepted_plant_name_id,taxon_status,geographic_area) %>% print()
to_select <- as.numeric(readline("Which one to keep?"))
if(is.na(to_select)) next
to_correct[[i]] %<>% slice(-all_of(to_select))
}
data.table::rbindlist(to_correct) %>% as_tibble() -> to_correct
rbind(to_correct,joined) -> joined
joined %>% pull(ID) %>% duplicated() -> dups
joined$ID[which(dups)] -> id_dups
cat("Still",length(id_dups),"some unsolved records!","\n")
joined %>% filter(is.na(plant_name_id)) %>% distinct(scientificName) %>% nrow() -> spp
cat(spp,"unresolved names left!","\n")
save(joined,file=here("interim/",output_file))
}

# NEEDS EDIT STILL
gbif_pow_first_fuzzy <- function(data="Joined_base.Rdata",
                                 checklist = "wcvp_names.txt",
                                 taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),
                                 output_file="Joined_fuzzy.Rdata"){ 
load(here("interim",data))
joined %>% filter(is.na(plant_name_id)) %>% distinct(scientificName) -> spp
  spp %>% filter(grepl("×",scientificName)) -> hybrids
  spp %>% filter(grepl("\\?",scientificName)) -> noidea
  spp %>% filter(!grepl("×",scientificName),!grepl("\\?",scientificName)) -> species

list_names <- data.table::fread(here("wcvp_2022",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
  list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ",remove = F)
  list_names %<>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble() %>% filter(taxon_status %in% all_of(taxon_status_cats)) %>% unite("binomial",c(genus,species),sep="_",remove = F) %>% filter(accepted_plant_name_id!="")

joined_nana <- joined %>% filter(is.na(plant_name_id))
joined %<>% filter(!is.na(plant_name_id))

resolved_matches <- list()
non_matches <- list()
multiple_matches <- list()
max.dist=0.15
for (i in 1:nrow(species)){
  species %>% slice(i) %>% pull(scientificName) -> quien
  cat(i,"Attempting match on:", quien,"\n")
  joined_nana %>% filter(scientificName %in% all_of(quien)) -> target
target %>% slice(1) %>% select(binomial.x,scientificName,accepted_plant_name_id,Accepted_Name) %>% rename(binomial=binomial.x) -> sole_target

list_names %>% filter(binomial== sole_target$binomial) -> names_tomatch
fuzzyjoin::stringdist_left_join(species %>% slice(i),names_tomatch,by="scientificName",method="jw",distance_col="distance",max_dist=max.dist) -> result
if(nrow(result)>1){ 
  if(length(unique(result$accepted_plant_name_id))==1){
    cat("        Single match!","\n")
    result[1,] -> transfer
    target$accepted_plant_name_id <- transfer$accepted_plant_name_id
    target$Accepted_Name <- transfer$Accepted_Name
    resolved_matches[[i]] <- target
    } 
  if(length(unique(result$accepted_plant_name_id))>1) {
     if(min(result$distance)<0.01){
       cat("        Single match!","\n")
       result[which.min(result$distance),] -> transfer
       target$accepted_plant_name_id <- transfer$accepted_plant_name_id
       target$Accepted_Name <- transfer$Accepted_Name
       resolved_matches[[i]] <- target
     }
     if(min(result$distance)>0.01){
   cat("        Multiple matches!","\n")
   target$accepted_plant_name_id <- NA
   target$Accepted_Name <- NA
   multiple_matches[[i]] <- target
    }
  }
  next
}
if(is.na(result$plant_name_id)){ cat("       No match!","\n")
  target$accepted_plant_name_id <- NA
  target$Accepted_Name <- NA
  non_matches[[i]] <- target
  next}
if(!is.na(result$plant_name_id)){
  cat("        Single match!","\n")
result %>% select(scientificName.x,accepted_plant_name_id,Accepted_Name) -> transfer
target$accepted_plant_name_id <- transfer$accepted_plant_name_id
target$Accepted_Name <- transfer$Accepted_Name
resolved_matches[[i]] <- target
}
}

resolved_matches %>% data.table::rbindlist(.) %>% as_tibble() -> joined_fuzzy
joined_fuzzy %<>% bind_rows(joined,.)
joined_fuzzy %>% distinct(Accepted_Name)

save(joined_fuzzy,file=here("interim",output_file))
save(non_matches,file=here("interim","non_matches_fuzzy.Rdata"))
save(multiple_matches,file=here("interim","multiples_matches_fuzzy.Rdata"))

}

gbif_pow_second_fuzzy <- function(output_file="Joined_final.v1.Rdata"){ 
load(here("interim/non_matches_fuzzy.Rdata"))
  load("interim/Joined_fuzzy.Rdata")
  non_matches %<>% data.table::rbindlist() %>% as_tibble() 
  non_matches %>% distinct(scientificName) -> species
resolved_matches <- list()
non_matches_final <- list()
multiple_matches_final <- list()
max.dist=0.15
for (i in 1:nrow(species)){
  species %>% slice(i) %>% pull(scientificName) -> quien
  cat(i,"Attempting match on:", quien,"\n")
  non_matches %>% filter(scientificName %in% all_of(quien)) -> target
  target %>% slice(1) %>% select(binomial.x,scientificName,accepted_plant_name_id,Accepted_Name,genus.x) %>% rename(binomial=binomial.x,genus=genus.x) -> sole_target
  
  list_names %>% filter(genus == sole_target$genus) -> names_tomatch
  fuzzyjoin::stringdist_left_join(species %>% slice(i),names_tomatch,by="scientificName",method="jw",distance_col="distance",max_dist=max.dist) -> result
  if(nrow(result)>1){ 
    if(length(unique(result$accepted_plant_name_id))==1){
      cat("        Single match!","\n")
      result[1,] -> transfer
      target$accepted_plant_name_id <- transfer$accepted_plant_name_id
      target$Accepted_Name <- transfer$Accepted_Name
      resolved_matches[[i]] <- target
    } 
    if(length(unique(result$accepted_plant_name_id))>1) {
      if(min(result$distance)<0.01){
        cat("        Single match!","\n")
        result[which.min(result$distance),] -> transfer
        target$accepted_plant_name_id <- transfer$accepted_plant_name_id
        target$Accepted_Name <- transfer$Accepted_Name
        resolved_matches[[i]] <- target
      }
      if(min(result$distance)>0.01){
        cat("        Multiple matches!","\n")
        target$accepted_plant_name_id <- NA
        target$Accepted_Name <- NA
        multiple_matches_final[[i]] <- target
      }
    }
    next
  }
  if(is.na(result$plant_name_id)){ cat("       No match!","\n")
    target$accepted_plant_name_id <- NA
    target$Accepted_Name <- NA
    non_matches_final[[i]] <- target
    next}
  if(!is.na(result$plant_name_id)){
    cat("        Single match!","\n")
    result %>% select(scientificName.x,accepted_plant_name_id,Accepted_Name) -> transfer
    target$accepted_plant_name_id <- transfer$accepted_plant_name_id
    target$Accepted_Name <- transfer$Accepted_Name
    resolved_matches[[i]] <- target
  }
}
non_matches_final %>% data.table::rbindlist(.) %>% as_tibble() 
multiple_matches %>% data.table::rbindlist(.) %>% as_tibble() %>% distinct(scientificName)

resolved_matches %>% data.table::rbindlist(.) %>% as_tibble() -> joined_fuzzy_f
joined_fuzzy %<>% bind_rows(.,joined_fuzzy_f)
joined_fuzzy %>% distinct(Accepted_Name)
save(joined_fuzzy,file=here("output",output_file))
save(non_matches_final,file="output/no_matches_final.v1.Rdata")
save(multiple_matches_final,file="output/multiples_matches_final.v1.Rdata")

}

## This is the same for all data bases
geographic_filter <- function(data="Joined_final.v1.Rdata",output_file = "Joined_finalPOW.v1.Rdata",
                              perform_tests=c("centroids","institutions", "equal", "gbif","capitals", "zeros","seas")) {
  load(here("output",data))
  ### APPLY FILTERS. NOT USING THE OUTLIER TEST, BECAUSE WE BASE THIS ON KEW'S DISTRIBUTIONS.
  cat("Filtering using CoordinateCleaner","\r")
  to_filter <- joined_fuzzy %>% select(decimalLongitude,decimalLatitude) %>% as.data.frame() 
  f <- CoordinateCleaner::clean_coordinates(to_filter, lon = "decimalLongitude",lat = "decimalLatitude",centroids_rad = 1000, centroids_detail = "both",tests=perform_tests,species=NULL,value="flagged",seas_scale = 110)
  joined_fuzzy %<>% mutate(Outlier_Test = f)
  save(joined_fuzzy,file=here("output",output_file))
}
  
powdist_filter <- function(data="Joined_finalPOW.v1.Rdata",output_file = "Joined_finalPOWdist.v2.Rdata",
                              pow_distributions="wcvp_distribution.txt",
                              wgsrpd = "level3/level3.shp") { 
  load(here("output",data))
  pow_dist <- data.table::fread(here("wcvp_2022",pow_distributions),sep="|",quote="") %>% as_tibble() %>% #left_join(.,list_names,by="accepted_db_id") %>% 
    filter(introduced==0)

  cat("....... reading polygons from WGSRPD","\r")
  poly <- rgdal::readOGR(here("wgsrpd-master",wgsrpd))
  pointos <- SpatialPoints(as.data.frame(joined_fuzzy[,11:10]))
  proj4string(pointos)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  points_powo <- sp::over(pointos,poly)
  joined_fuzzy %<>% bind_cols(.,points_powo)
  cat("Filtering using KEW's distributions","\n")
  joined_fuzzy %>% distinct(Accepted_Name)
  joined_fuzzy %>% group_split(Accepted_Name) -> aver
  some <- list()
  for (i in 1:length(aver)){
      cat(i,"out of",scales::comma(length(aver)),"--","Getting POW data for", unique(aver[[i]]$Accepted_Name),"               ","\r")
    temp <- aver[[i]] %>% mutate(POW_distribution = FALSE)
    if(is.na(unique(aver[[i]]$Accepted_Name))) {
      some[[i]] <- temp %>% mutate(POW_distribution = TRUE)
    }
    ok <- pow_dist %>% filter(plant_name_id==unique(temp$accepted_plant_name_id))
    matching <- !is.na(match(temp$LEVEL3_COD,ok %>% pull(area_code_l3),nomatch = NA))
    if(length(which(matching)) == 0) { 
        if(ok %>% nrow() == 0){
          some[[i]] <- temp %>% mutate(POW_distribution = TRUE)
          next
        } else {some[[i]] <- temp %>% mutate(POW_distribution = FALSE)
        next}
      some[[i]] <- temp %>% mutate(POW_distribution = TRUE)
      next
    }
    if(length(which(matching)) != 0) some[[i]] <- temp %>% mutate(POW_distribution = matching)
  }
  rm(aver,joined_fuzzy,pointos,poly)
  some %<>% data.table::rbindlist(.) %>% as_tibble()
  cat("-------------- DONE --------------","\r")
  save(some,file=here("output",output_file))
 }



###### THIS ONE IS MISSING ARGUMENTS TO CONTROL THE PLOT #####

fetch_accessions <- function(target="rbcL OR atpB OR matK OR ndhF OR trnG or trnL",group="Tracheophyta",seq.length="100:7000",db="nuccore",ret.max=850000,output_file="NCBI_v1.txt"){ 
  query <- glue::glue('{target} AND {group}[ORGN] AND {seq.length}[SLEN] AND plants[filter] AND biomol_genomic[PROP] AND is_nuccore[filter]')
  accesions <- reutils::esearch(query,db=db,retmax=ret.max,rettype="uilist",usehistory=T)
  reutils::efetch(accesions,db=db,rettype="acc",outfile=here("data",output_file))
} # FETCH ACCESSION NUMBERS FROM NCBI

ncbi_databases <- function() { 
  taxonomizr::getAccession2taxid(outDir = "data")
  taxonomizr::read.accession2taxid(taxaFiles = "data/nucl_gb.accession2taxid.gz",sqlFile = "data/acc_taxid",indexTaxa=T)
  # CREATE THE TAXA ID NAMES DATABASE
  taxonomizr::prepareDatabase(sqlFile = "nameNode.sqlite",getAccessions = F)
} # CREATE THE ACCESSION TO TAXAID DATABASE

do_format_completeness <- function(){
  data="output/Joined_finalPOWdist.v2.Rdata"
  load(data)
  wgsrpd = "level2/level2.shp"
  poly <- rgdal::readOGR(here("wgsrpd-master",wgsrpd))
  poly <- poly[poly$LEVEL2_NAM %in% c("Mexico","Central America"),]
  cat("....... buffering polygons from WGSRPD","\r")
  spList = vector("list", length(poly))
  for (i in 1:length(poly)) {
    cat(i,"\r")
    a <- rgeos::gBuffer(poly[i,], width = 0.5)
    a$LEVEL2_COD = poly[i,]$LEVEL2_COD
    a$LEVEL2_NAM = poly[i,]$LEVEL2_NAM
    spList[[i]] <- a
  }
  poly <- do.call("rbind", spList)
  
  projection = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  #projection = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs"
  
  #### ONE DEGREE
  g <- raster::raster(nrows=180*1,ncols=360*1,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% projectRaster(.,crs = CRS(projection)) %>% as(., 'SpatialPixels')
  some %>% select(decimalLongitude,decimalLatitude) %>% 
    sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID_1d"))) %>% select(CellID_1d) %>% bind_cols(.,some) -> gridded
  
g@coords %>% as_tibble() %>% mutate(CellID_1d=1:nrow(.)) -> coords
  pointos <- coords %>% select(1:2) %>% as.data.frame(.) %>% SpatialPoints(.)
  proj4string(pointos) <- projection
  points_powo <- sp::over(pointos,poly) %>% as_tibble()
  coords %<>% bind_cols(.,points_powo)
  
  gridded %>%  filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution) %>% select(CellID_1d,Accepted_Name) %>% group_by(CellID_1d,Accepted_Name) %>% summarize(Counts=n()) %>% ungroup() %>% left_join(.,coords,by="CellID_1d") -> one_degree
  
  one_degree %<>% mutate(binomial_accepted=Accepted_Name) %>% mutate(binomial_accepted=sub("× ","",binomial_accepted))
  one_degree %<>% mutate(binomial_accepted=sub(" ","_",binomial_accepted))
  one_degree %<>% separate(binomial_accepted,into=c("binomial"),extra = "warn",sep=" ")
  one_degree
  one_degree %<>% ungroup() %>%  select(binomial,x,y,Counts,LEVEL2_NAM) %>% 
    rename("Species"=binomial,"Longitude"=x,"Latitude"=y)
  
  ### HALF DEGREE
  g <- raster::raster(nrows=180*2,ncols=360*2,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% projectRaster(.,crs = CRS(projection)) %>% as(., 'SpatialPixels')
  gridded %>% select(decimalLongitude,decimalLatitude) %>% 
    sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID_0.5d"))) %>% select(CellID_0.5d) %>% bind_cols(.,gridded) -> gridded
  
  g@coords %>% as_tibble() %>% mutate(CellID_0.5d=1:nrow(.)) -> coords
  pointos <- coords %>% select(1:2) %>% as.data.frame(.) %>% SpatialPoints(.)
  proj4string(pointos) <- projection
  points_powo <- sp::over(pointos,poly) %>% as_tibble()
  coords %<>% bind_cols(.,points_powo)
  
  gridded %>% filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution) %>%select(CellID_0.5d,Accepted_Name) %>% group_by(CellID_0.5d,Accepted_Name) %>% summarize(Counts=n()) %>% ungroup() %>% right_join(.,coords,by="CellID_0.5d") -> half_degree
  
  half_degree %<>% mutate(binomial_accepted=Accepted_Name) %>% mutate(binomial_accepted=sub("× ","",binomial_accepted))
  half_degree %<>% mutate(binomial_accepted=sub(" ","_",binomial_accepted))
  half_degree %<>% separate(binomial_accepted,into=c("binomial"),extra = "warn",sep=" ")
  half_degree
  
  half_degree %<>% ungroup() %>%  select(binomial,x,y,Counts,LEVEL2_NAM) %>% 
    rename("Species"=binomial,"Longitude"=x,"Latitude"=y)

  
  ### FIFTH DEGREE
  g <- raster::raster(nrows=180*5,ncols=360*5,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% projectRaster(.,crs = CRS(projection)) %>% as(., 'SpatialPixels')
  gridded %>% select(decimalLongitude,decimalLatitude) %>% 
    sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID_0.2d"))) %>% select(CellID_0.2d) %>% bind_cols(.,gridded) -> gridded
  
  g@coords %>% as_tibble() %>% mutate(CellID_0.2d=1:nrow(.)) -> coords
  pointos <- coords %>% select(1:2) %>% as.data.frame(.) %>% SpatialPoints(.)
  proj4string(pointos) <- projection
  points_powo <- sp::over(pointos,poly) %>% as_tibble()
  coords %<>% bind_cols(.,points_powo)
  
  gridded %>% filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution) %>%select(CellID_0.2d,Accepted_Name) %>% group_by(CellID_0.2d,Accepted_Name) %>% summarize(Counts=n()) %>% ungroup() %>% right_join(.,coords,by="CellID_0.2d") -> fifth_degree
  
  fifth_degree %<>% mutate(binomial_accepted=Accepted_Name) %>% mutate(binomial_accepted=sub("× ","",binomial_accepted))
  fifth_degree %<>% mutate(binomial_accepted=sub(" ","_",binomial_accepted))
  fifth_degree %<>% separate(binomial_accepted,into=c("binomial"),extra = "warn",sep=" ")
  fifth_degree
  fifth_degree %<>% ungroup() %>%  select(binomial,x,y,Counts,LEVEL2_NAM) %>% 
    rename("Species"=binomial,"Longitude"=x,"Latitude"=y)
  
  ### TENTH DEGREE
  g <- raster::raster(nrows=180*10,ncols=360*10,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% projectRaster(.,crs = CRS(projection)) %>% as(., 'SpatialPixels')
  gridded %>% select(decimalLongitude,decimalLatitude) %>% 
    sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID_0.1d"))) %>% select(CellID_0.1d) %>% bind_cols(.,gridded) -> gridded
  
  g@coords %>% as_tibble() %>% mutate(CellID_0.1d=1:nrow(.)) -> coords
  pointos <- coords %>% select(1:2) %>% as.data.frame(.) %>% SpatialPoints(.)
  proj4string(pointos) <- projection
  points_powo <- sp::over(pointos,poly) %>% as_tibble()
  coords %<>% bind_cols(.,points_powo)
  
  gridded %>% filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution) %>% select(CellID_0.1d,Accepted_Name) %>% group_by(CellID_0.1d,Accepted_Name) %>% summarize(Counts=n()) %>% ungroup() %>% right_join(.,coords,by="CellID_0.1d") -> tenth_degree
  
  tenth_degree %<>% mutate(binomial_accepted=Accepted_Name) %>% mutate(binomial_accepted=sub("× ","",binomial_accepted))
  tenth_degree %<>% mutate(binomial_accepted=sub(" ","_",binomial_accepted))
  tenth_degree %<>% separate(binomial_accepted,into=c("binomial"),extra = "warn",sep=" ")
  tenth_degree
  tenth_degree %<>% ungroup() %>%  select(binomial,x,y,Counts,LEVEL2_NAM) %>% 
    rename("Species"=binomial,"Longitude"=x,"Latitude"=y)
  
  
  for_comp <- list(one_degree=one_degree,half_degree=half_degree,
                   fifth_degree=fifth_degree,tenth_degree=tenth_degree)
  
  save(for_comp,file = here("output/pams_for_completeness_v2.Rdata"))
  

  }
#####

process_completeness <- function(){
  data="output/completeness_1degree/Estimators.RData"
  load(data)
  roads <- rnaturalearth::ne_countries(scale = 110,returnclass = "sf")
  datos <- values %>% as_tibble() %>%
    mutate(CompSlope = cut(Slope,breaks=c(0,0.02,0.1,0.3,1),labels=c("High","Fair","Insufficient","Poor"))) %>% 
    mutate(CompComp = cut(Completeness,breaks=c(0,50,75,90,100),labels=c("Poor","Insufficient","Fair","High"))) %>% 
    mutate(CompRatio = cut(Ratio,breaks=c(0,3,10,15,100),labels=c("Poor","Insufficient","Fair","High"))) %>% 
    mutate(Category = "Fair") %>% 
    mutate(Category = ifelse(CompSlope%in%c("High","Fair") & CompComp%in%c("High","Fair") & CompRatio%in%c("High","Fair"),"Fair",Category)) %>% 
    mutate(Category = ifelse(CompSlope%in%c("Insufficient","Poor") & CompComp%in%c("Insufficient","Poor") & CompRatio%in%c("Insufficient","Poor"),"Insufficient",Category)) %>% 
    mutate(Category = ifelse(CompSlope=="High" & CompComp=="High" & CompRatio=="High","High",Category)) %>% 
    mutate(Category = ifelse(CompSlope=="Poor" & CompComp=="Poor" & CompRatio=="Poor","Poor",Category)) %>% 
    mutate(Category = ifelse(is.na(Category),"Poor",Category))
  
  datos %>% summarise(n=sum(Records)) %>% pull(n) -> final
  
  datos <- datos %>% mutate(Category= factor(Category,levels=c("High","Fair","Insufficient","Poor")))
  
  
  
  
  ggplot() +
    geom_sf(data=roads,colour=NA,fill=adjustcolor(MetBrewer::met.brewer("Demuth",direction=-1)[1],alpha.f = 0.2),size=0.5) + xlim(c(-120,-70)) + ylim(5,35) + 
    theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),legend.position = c(-.1,0.3),legend.background = element_rect(fill = NA),legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"),legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
    geom_tile(data=datos,aes(x=Longitude,y=Latitude,fill=Category)) +
    scale_fill_manual(values=MetBrewer:::met.brewer("Hiroshige",n=12)[c(1,4,8,12)],name="") +
    labs(x="",y="",title="Digital botanical knowledge in Mesoamerica",subtitle=paste0(formatC(final,format="f",digits = 0, big.mark=",")," occurrence records from GBIF")) +
    standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
    NULL
}

estimate_metrics <- function(
  data="Joined_finalPOWdist.v2.Rdata",
  resolution = 60,
  data_comp="output/completeness_60secs/Estimators.RData"){ 
  load(here("output",data))

  some %<>% mutate(binomial_accepted=Accepted_Name) %>% mutate(binomial_accepted=sub("× ","",binomial_accepted))
  some %<>% mutate(binomial_accepted=sub(" ","_",binomial_accepted))
  some %<>% separate(binomial_accepted,into=c("binomial"),extra = "warn",sep=" ")
  
  wgsrpd = "level2/level2.shp"
  poly <- rgdal::readOGR(here("wgsrpd-master",wgsrpd))
  poly <- poly[poly$LEVEL2_NAM %in% c("Mexico","Central America"),]
  cat("....... buffering polygons from WGSRPD","\r")
  spList = vector("list", length(poly))
  for (i in 1:length(poly)) {
    cat(i,"\r")
    a <- rgeos::gBuffer(poly[i,], width = 0.5)
    a$LEVEL2_COD = poly[i,]$LEVEL2_COD
    a$LEVEL2_NAM = poly[i,]$LEVEL2_NAM
    spList[[i]] <- a
  }
  poly <- do.call("rbind", spList)
  
  projection = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  res <- 60/resolution
  g <- raster::raster(nrows=180*res,ncols=360*res,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% projectRaster(.,crs = CRS(projection)) %>% as(., 'SpatialPixels')
  some %>% select(decimalLongitude,decimalLatitude) %>% 
    sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID"))) %>% select(CellID) %>% bind_cols(.,some) -> gridded
  
  g@coords %>% as_tibble() %>% mutate(CellID=1:nrow(.)) -> coords
  pointos <- coords %>% select(1:2) %>% as.data.frame(.) %>% SpatialPoints(.)
  proj4string(pointos) <- projection
  points_powo <- sp::over(pointos,poly) %>% as_tibble()
  coords %<>% bind_cols(.,points_powo)

  # NEED TO PROCESS THE DATA AS IN DO_FORMAT_COMPLETENESS TO GET THE IDS AND THEN DO ALL THE ANALYSES:
  #     - SPECIES RICHNESS
  gridded %>% select(binomial,CellID,year,Outlier_Test,POW_distribution) %>% 
    filter(!is.na(binomial),Outlier_Test,POW_distribution) %>% 
    distinct(binomial,CellID) %>% count(CellID) %>% 
    left_join(.,gridded %>% distinct(CellID,.keep_all = T),by="CellID") %>% 
    select(CellID,n,LEVEL3_NAM) %>% rename(SR=n) -> SR
  SR %<>% mutate(SR_std = (SR - min(SR,na.rm=TRUE)) / (max(SR,na.rm=TRUE) - min(SR,na.rm=TRUE)))
  #     - NUMBER OF RECORDS
  gridded %>% filter(!is.na(binomial),Outlier_Test,POW_distribution) %>% count(CellID) %>% 
    left_join(.,gridded %>% distinct(CellID,.keep_all = T),by="CellID") %>% select(CellID,n) %>% inner_join(.,SR,by="CellID") %>% mutate(n_std = (n - min(n,na.rm=TRUE)) / (max(n,na.rm=TRUE) - min(n,na.rm=TRUE))) -> SR
  nombre <-  paste0("RichnessRecords_",resolution,".Rdata")
  richness <- SR %>% left_join(.,coords,by="CellID")
  
 #  ctyshp = "paises.shp"
 #  poly_cty <- rgdal::readOGR(here("data/countries",ctyshp))
 #  #poly_cty <- poly_cty[poly_cty$NAME_0 %in% c("Mexico","United States","Guatemala","El Salvador","Belize","Nicaragua","Honduras","Costa Rica","Panama","Colombia","Cuba","Jamaica","Puerto Rico"),]
 #  g_ras <- raster::raster(nrows=180*res,ncols=360*res,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% projectRaster(.,crs = CRS(projection))
 #  
 # #  if(resolution==15) {
 # # points_ctys_1 <- sp::over(pointos[1:150000],poly_cty) %>% as_tibble()
 # # points_ctys_2 <- sp::over(pointos[150001:250000],poly_cty) %>% as_tibble()
 # # points_ctys_3 <- sp::over(pointos[250001:350000],poly_cty) %>% as_tibble()
 # # points_ctys_4 <- sp::over(pointos[350001:450000],poly_cty) %>% as_tibble()
 # # points_ctys_5 <- sp::over(pointos[450001:550000],poly_cty) %>% as_tibble()
 # # points_ctys_6 <- sp::over(pointos[550001:650000],poly_cty) %>% as_tibble()
 # # points_ctys_7 <- sp::over(pointos[650001:750000],poly_cty) %>% as_tibble()
 # # points_ctys_8 <- sp::over(pointos[750001:850000],poly_cty) %>% as_tibble()
 # # points_ctys_9 <- sp::over(pointos[850001:950000],poly_cty) %>% as_tibble()
 # # points_ctys_10 <- sp::over(pointos[950001:1036800],poly_cty) %>% as_tibble()
 # # points_ctys <- list(points_ctys_1,points_ctys_2,points_ctys_3,points_ctys_4,points_ctys_5,points_ctys_6,points_ctys_7,points_ctys_8,points_ctys_9,points_ctys_10) %>% do.call("rbind", .)
 # # } else points_ctys <- sp::over(pointos,poly_cty) %>% as_tibble()
 # #  coords %<>% bind_cols(.,points_ctys)
 #  ras_cty <- rasterize(poly_cty,g_ras)
 #  ras_cty[is.na(ras_cty)] <- 9999
 # 
 #  ras_cty[24562]
 #  
 #  richness %>% left_join(.,cty_tib,by="CellID") %>% count(NAME_0) %>% view
 #  
  save(richness,file=here("interim",nombre))
  
#     - AGE OF COLLECTIONS
gridded %>% select(CellID,binomial,year,Outlier_Test,POW_distribution) %>% filter(Outlier_Test,POW_distribution,year >=1922) %>% group_by(CellID,binomial) %>% nest() -> test
 test <- test %>% ungroup() %>% rowwise() %>% mutate(MostRecent = (get_mostrecent_year(data)))
 nombre <-  paste0("Ages_",resolution,".Rdata")
 families <- test %>% 
   left_join(.,gridded %>% select(binomial,family.x,CellID,year,Outlier_Test,POW_distribution) %>% 
               filter(!is.na(binomial),Outlier_Test,POW_distribution) %>% 
               distinct(family.x,binomial),by="binomial") 
  families %<>% left_join(.,coords,by="CellID")
  save(families,file=here("interim",nombre))

  test %>% group_by(CellID) %>% summarise(age= 2022 - median(MostRecent,na.rm=T)) -> ages
  SR <- ages %>% inner_join(.,richness,by="CellID")

  #     - NCBI
  ## JOIN RECORDS WITH RBCL DATA
  data = "NCBI_v1.txt"
  data.table::fread(here("data",data),header = F) %>% as_tibble() -> accs
  accs %>% pull() -> query_accs
  taxonomizr::accessionToTaxa(accessions = query_accs, sqlFile = "data/acc_taxid") %>% as_tibble() %>% bind_cols(.,query_accs) -> aver
  aver %>% rename_with(~c("TaxaID","accession")) %>% filter(!is.na(TaxaID)) -> aver
  taxonomizr::getTaxonomy(aver %>%  pull(TaxaID),sqlFile = "data/nameNode.sqlite", desiredTaxa = c("family","genus","species")) -> taxonomy
  aver %>% bind_cols(.,taxonomy %>% as_tibble()) %>% mutate(region="rbcL",.after=TaxaID) -> aver
  data.table::fread("data/names.dmp") -> names
  names %>% as_tibble() %>% filter(grepl("authority",V4)) %>% rename(TaxaID=V1,authority=V2) %>% left_join(aver,.,by="TaxaID") %>% select(-starts_with("V")) %>% mutate(authority=gsub("\t","",authority)) -> names
  
  checklist = "wcvp_names.txt"
  taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced")
  list_names <- data.table::fread(here("wcvp_2022",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
  list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ") %>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble()
  list_names %<>% filter(infraspecific_rank=="",taxon_status %in% all_of(taxon_status_cats))
  list_names %<>% unite("binomial",c(genus,species),sep=" ",remove = F) 
  names %<>% filter(!grepl("x",species)) %>% filter(!grepl("]",species)) %>% filter(!grepl(")",species))  %>%
    #distinct(TaxaID,.keep_all = T) %>% 
    rename(binomial=species) %>% left_join(.,list_names,by = "binomial") %>% select(TaxaID,region,binomial,accepted_plant_name_id,Accepted_Name)
  names
  
 genebanks <-  gridded %>% select(CellID,year,Outlier_Test,POW_distribution,accepted_plant_name_id) %>% 
    mutate(TaxaID=names$TaxaID[match(gridded$accepted_plant_name_id,names$accepted_plant_name_id)]) %>%
    mutate(ncbi=names$region[match(gridded$accepted_plant_name_id,names$accepted_plant_name_id)]) %>% filter(Outlier_Test,POW_distribution)  %>% distinct(CellID,accepted_plant_name_id,.keep_all = T) %>% group_by(CellID) %>% 
    mutate(value=ifelse(is.na(ncbi),0,1)) %>% summarise(seqs = sum(value)/n())
 
 SR %<>% inner_join(.,genebanks,by="CellID")
 
 # - COMPLETENESS
 load(data_comp)
 datos <- values %>% as_tibble()
 datos %>% select(Longitude,Latitude) %>% 
   sf::st_as_sf(x = ., coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID"))) %>% select(CellID) %>% bind_cols(.,datos) -> datos
 SR %<>% rename(Longitude = x, Latitude = y) %>% right_join(.,datos,by="CellID")
 datos <- SR
 nombre <-  paste0("Allmetrics_",resolution,".Rdata")
 save(datos,file=here("interim",nombre))
 
}

plot_mappy_things <- function(data="interim/Allmetrics_60.Rdata",
                              resolution=60,
                              rich_bin=200,
                              comp_bin=5,
                              seq_bin=0.05,
                              age_bin=5){
load(here(data))
  if(resolution==60) { datos %<>% rename(x.y=Longitude.y,y.y=Latitude.y)}
roads <- rnaturalearth::ne_countries(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
 counts <- rnaturalearth::ne_coastline(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
 xlimits = c(-125,-75)
 ylimits = c(0,35)
 ###

 ####
 
 ##  richness
 upperlim <-  roundUp(max(datos$SR,na.rm=T),to=1000)
 richness <- ggplot() + xlim(xlimits) + ylim(ylimits) + 
   theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
         legend.position = "none",
         legend.background = element_rect(fill = NA),
         legend.spacing.y = unit(c(0.2),"cm"),
         legend.title = element_text(family="EB Garamond"),
         legend.text = element_text(family="EB Garamond"),
         legend.margin=margin(t=-25),
         legend.key.size = unit(0.1,"cm"),
         plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
   geom_tile(data = datos,aes(x=x.y,y=y.y,fill=SR,color=SR)) +
   scale_fill_stepsn(colors=(MetBrewer::met.brewer("Derain",direction=1)),name="Species richness",
                     breaks=seq(0,upperlim,rich_bin),labels= scales::comma,
                     na.value="NA",limits=c(0,upperlim),
                     guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5)) +
   scale_color_stepsn(colors=(MetBrewer::met.brewer("Derain",direction=1)),name="Species richness",
                      breaks=seq(0,upperlim,rich_bin),labels= scales::comma,
                      na.value="NA",limits=c(0,upperlim),
                      guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5)) +
   labs(x="",y="") +
  # standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
   geom_sf(data=roads,colour="black",fill=NA,size=0.2) +
   NULL

 n_bins <- length(ggplot2:::bin_breaks_width(range(0,upperlim), width = rich_bin)$breaks) - 1L
 
hist_rich <- ggplot() + 
  geom_histogram(data=datos,aes(x=SR),binwidth=rich_bin,fill=MetBrewer::met.brewer("Derain",direction=1,n=n_bins)) +
   theme(panel.background = element_blank(),
         panel.grid = element_blank(),
         axis.line = element_line(),
         legend.position = "none",
         legend.background = element_rect(fill = NA),
         legend.spacing.y = unit(c(0.2),"cm"),
         legend.title = element_text(family="EB Garamond"),
         legend.text = element_text(family="EB Garamond"),
         axis.title.y = element_text(family="EB Garamond",size=10),
         axis.title.x = element_text(family="EB Garamond",size=10),
         axis.text =  element_text(family="EB Garamond",size = 8),
         legend.margin=margin(t=-25),
         legend.key.size = unit(0.1,"cm"),
         plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Species richness",y="Number of cells") +
  scale_x_continuous(breaks=seq(0,upperlim,1000),limits = c(0,upperlim)) +
 NULL

richness <- richness +
  inset_element(hist_rich,left=0.01,right=0.6,bottom=0.01,top=0.5,on_top=FALSE)
 

### completeness
upperlim <-  roundUp(max(datos$Completeness,na.rm=T),to=100)
 completeness  <- ggplot() + 
   geom_sf(data=roads,colour=NA,fill=adjustcolor(MetBrewer::met.brewer("Demuth",direction=-1)[1],alpha.f = 0.1),size=0.5) + xlim(xlimits) + ylim(ylimits) + 
   theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
         legend.position = "none",
         legend.background = element_rect(fill = NA),
         legend.key.size = unit(0.1,"cm"),
         legend.spacing.y = unit(c(0.2),"cm"),
         legend.title = element_text(family="EB Garamond"),
         legend.text = element_text(family="EB Garamond"),
         legend.margin=margin(t=-25),
         plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
   geom_tile(data = datos,aes(x=x.y,y=y.y,fill=Completeness),color="white") +
   scale_fill_stepsn(colors=(MetBrewer::met.brewer("Hokusai3",direction=1)),name="Completeness",
                     breaks=seq(0,upperlim,comp_bin),labels= scales::comma,
                     na.value="NA",limits=c(0,upperlim),
                     guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5,keyheight=0.4)) +
   labs(x="",y="") +
   #standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
   geom_sf(data=roads,colour="black",fill=NA,size=0.2) +
   NULL
 n_bins <- length(ggplot2:::bin_breaks_width(range(0,upperlim), width = comp_bin)$breaks) - 1L
 
 hist_comp <- ggplot() + 
   geom_histogram(data=datos,aes(x=Completeness),fill=MetBrewer::met.brewer("Hokusai3",direction=1,n=n_bins),binwidth=comp_bin,color="white") +
   theme(panel.background = element_blank(),
         panel.grid = element_blank(),
         axis.line = element_line(),
         legend.position = "none",
         legend.background = element_rect(fill = NA),
         legend.spacing.y = unit(c(0.2),"cm"),
         legend.title = element_text(family="EB Garamond"),
         legend.text = element_text(family="EB Garamond"),
         axis.title.y = element_text(family="EB Garamond",size=10),
         axis.title.x = element_text(family="EB Garamond",size=10),
         axis.text =  element_text(family="EB Garamond",size = 8),
         legend.margin=margin(t=-25),
         legend.key.size = unit(0.1,"cm"),
         plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
   labs(x="Completeness",y="Number of cells") +
   scale_x_continuous(breaks=seq(0,upperlim,10),limits = c(0,upperlim)) +
   NULL
 completeness <- completeness +
   inset_element(hist_comp,left=0.01,right=0.6,bottom=-0.01,top=0.5,on_top=FALSE)
 
 # age
 upperlim <-  roundUp(max(datos$age,na.rm=T),to=100)
 age <-  ggplot() + geom_sf(data=roads,colour=NA,fill=adjustcolor(MetBrewer::met.brewer("Demuth",direction=-1)[1],alpha.f = 0.1),size=0.5) + xlim(xlimits) + ylim(ylimits) + 
   theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
         legend.position = "none",
         legend.background = element_rect(fill = NA),
         legend.key.size = unit(0.1,"cm"),
         legend.spacing.y = unit(c(0.2),"cm"),
         legend.title = element_text(family="EB Garamond"),
         legend.text = element_text(family="EB Garamond"),
         legend.margin=margin(t=-25),
         plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
   geom_tile(data = datos,aes(x=x.y,y=y.y,fill=age),color="white") +
   scale_fill_stepsn(colors=(MetBrewer::met.brewer("Demuth",direction=1)),name="Age",breaks=seq(0,upperlim,age_bin),labels= scales::comma,limits=c(0,upperlim),
                     na.value="NA",guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5)) +
   labs(x="",y="") +
   #standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
   geom_sf(data=roads,colour="black",fill=NA,size=0.2) +
   NULL
 n_bins <- length(ggplot2:::bin_breaks_width(range(0,upperlim), width = age_bin)$breaks) - 1L
 
 hist_age <- ggplot() + 
   geom_histogram(data=datos,aes(x=age),fill=MetBrewer::met.brewer("Demuth",direction=1,n=n_bins),binwidth=age_bin,color="white") +
   theme(panel.background = element_blank(),
         panel.grid = element_blank(),
         axis.line = element_line(),
         legend.position = "none",
         legend.background = element_rect(fill = NA),
         legend.spacing.y = unit(c(0.2),"cm"),
         legend.title = element_text(family="EB Garamond"),
         legend.text = element_text(family="EB Garamond"),
         axis.title.y = element_text(family="EB Garamond",size=10),
         axis.title.x = element_text(family="EB Garamond",size=10),
         axis.text =  element_text(family="EB Garamond",size = 8),
         legend.margin=margin(t=-25),
         legend.key.size = unit(0.1,"cm"),
         plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
   labs(x="Age",y="Number of cells") +
   scale_x_continuous(breaks=seq(0,upperlim,10),limits = c(0,upperlim)) +
   NULL
 
 age <- age +
   inset_element(hist_age,left=0.01,right=0.6,bottom=-0.01,top=0.5,on_top=FALSE)
 
 # sequenced
 upperlim <-  roundUp(max(datos$seqs,na.rm=T),to=1)
 #n_bins <- length(ggplot2:::bin_breaks_width(range(0,upperlim), width = seq_bin)$breaks) - 1L
 n_bins = 12
 colas <- rev(MetBrewer::met.brewer("Tam",direction=-1,n=n_bins))
 colas <- c(rep(colas[1],7),colas,rep(colas[length(colas)],2))
 sequenced <- ggplot() + 
   xlim(xlimits) + ylim(ylimits) + 
   theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
         legend.position = "none",
         legend.background = element_rect(fill = NA),
         legend.key.size = unit(0.1,"cm"),
         legend.spacing.y = unit(c(0.2),"cm"),
         legend.title = element_text(family="EB Garamond"),
         legend.text = element_text(family="EB Garamond"),
         legend.margin=margin(t=-25),
         plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),
         plot.margin = unit(c(.0,.0,.0,.0), "cm"),
         legend.direction="horizontal") +
   geom_tile(data = datos,aes(x=x.y,y=y.y,fill=seqs),color="white") +
   scale_fill_stepsn(colors=colas,
                     name="Sequenced",
                     breaks = seq(0,upperlim,seq_bin),
                     #c(0.2,0.5,seq(0.55,upperlim,seq_bin)),
                     labels= scales::comma,limits=c(0.2,upperlim),
                     guide=guide_legend(reverse = TRUE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5)) +
   labs(x="",y="") +
   #standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
   geom_sf(data=roads,colour="black",fill=NA,size=0.2) +
   NULL

 hist_sequenced <- ggplot() + 
   geom_histogram(data=datos,aes(x=seqs),
                  binwidth=seq_bin,
                  fill=colas,color="white") +
   theme(panel.background = element_blank(),
         panel.grid = element_blank(),
         axis.line = element_line(),
         legend.position = "none",
         legend.background = element_rect(fill = NA),
         legend.spacing.y = unit(c(0.2),"cm"),
         legend.title = element_text(family="EB Garamond"),
         legend.text = element_text(family="EB Garamond"),
         axis.title.y = element_text(family="EB Garamond",size=10),
         axis.title.x = element_text(family="EB Garamond",size=10),
         axis.text =  element_text(family="EB Garamond",size = 8),
         legend.margin=margin(t=-25),
         legend.key.size = unit(0.1,"cm"),
         plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
   labs(x="Sequenced",y="Number of cells") +
   scale_x_continuous(breaks=seq(0.2,upperlim,seq_bin*2),limits = c(0,upperlim)) +
   NULL

 sequenced <- sequenced +
   inset_element(hist_sequenced,left=0.01,right=0.6,bottom=-0.01,top=0.5,on_top=FALSE)

 
pp <- wrap_plots(list(richness,completeness,age,sequenced)) +
   plot_layout(byrow=T,nrow = 2,guides="keep",tag_level = 'new') & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"))
resolution
nombre <- paste0("Maps_",resolution,".pdf")
ggsave(filename = here("figures",nombre),plot = pp)

}
 
plot_phylo_things <- function(phylo = "data/Carta_tofamily.tre",
                              age="interim/Ages_60.Rdata",
                              data="interim/Allmetrics_60.Rdata",
                              resolution=60){ 
  load(here(data))
  load(here(age))
  datos %>% pull(CellID) ->  cells
  families %<>% select(-data) %>% group_by(binomial) %>% summarise(Youngest=max(MostRecent),family=first(family.x)) %>% group_by(family) #%>% summarise(Youngest_mean_age = median(Youngest,na.rm=T)) 
  tree <- ape::read.tree(here(phylo))
  tree <- ape::drop.tip(tree,tip=which(is.na(match(tree$tip.label,families %>% distinct(family) %>% pull()))))
  families %<>% mutate(Record_age = 2022 - Youngest)
  families %<>% mutate(median=mean(Record_age))
  p <- ggtree(tree, layout="fan", open.angle=10) + 
    geom_cladelab(node=ape::getMRCA(tree,c("Equisetaceae","Cyatheaceae")),offset=100,label="F",barsize=1, 
                  barcolor=MetBrewer::met.brewer("Kandinsky")[1],
                  textcolor="white",
                  fill=MetBrewer::met.brewer("Kandinsky")[1],offset.text=10,geom='label') +
    geom_cladelab(node=ape::getMRCA(tree,c("Cupressaceae","Zamiaceae")),label="G",offset=100,barsize=1, 
                  barcolor=MetBrewer::met.brewer("Kandinsky")[3],
                  fill=MetBrewer::met.brewer("Kandinsky")[3],
                  textcolor="white",offset.text=10,geom='label') +
    geom_cladelab(node=ape::getMRCA(tree,c("Cabombaceae","Lamiaceae")),label="A",offset=100,barsize=1, 
                  barcolor=MetBrewer::met.brewer("Kandinsky")[2],
                  fill=MetBrewer::met.brewer("Kandinsky")[2],
                  textcolor="white",
                  offset.text=55,geom='label') +
    geom_cladelab(node=ape::getMRCA(tree,c("Selaginellaceae","Isoetaceae")),label="L",offset=100,barsize=1, 
                  barcolor=MetBrewer::met.brewer("Kandinsky")[4],
                  fill=MetBrewer::met.brewer("Kandinsky")[4],
                  textcolor="white",
                  offset.text=10,geom='label') 
  p <- p +
    geom_fruit(
      data=families,
      geom=geom_point,
      mapping = aes(
        y=family,
        x=Record_age,
        color=Record_age,
        fill=Record_age),
      alpha=0.9,size=.1,
      outlier.size=0.1,
      outlier.stroke=0.08,
      outlier.shape=21,
      axis.params=list(
        axis       = "x",
        text.size  = 3,
        hjust      = 0,
        vjust      = 0.5,
        nbreak     = 3,
      ),show.legend=TRUE,
      grid.params=list()
    ) + scale_fill_stepsn(colors=MetBrewer::met.brewer("Demuth"),n.breaks=20,
                          guide=guide_legend(label.position = "bottom", title.hjust = 0, nrow=1,title.position = "top"
                                             ,override.aes = list(size=4,shape=15)
                          ),
                          name="Age in years") +
    scale_color_stepsn(colors=MetBrewer::met.brewer("Demuth"),n.breaks=20,
                       guide=guide_legend(label.position = "bottom", title.hjust = 0,nrow=1,title.position = "top"
                                          ,override.aes = list(size=4,shape=15)
                                          
                       ),
                       name="Age in years") +
    #guides(color = guide_legend(override.aes = list(size = 0)),fill = guide_legend(override.aes = list(size = 0)))+
    theme(legend.position=c(0.5,0.06),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0,0,0), "cm"),
          text = element_text(family="EB Garamond"),
          legend.title=element_text(size=10,family="EB Garamond"), 
          legend.text=element_text(size=7,family="EB Garamond"),
          plot.title=element_text(family="EB Garamond",face="bold",hjust = 0,size=18,vjust = -2),
          plot.caption = element_text(family="EB Garamond",face="italic",hjust = 0,size=8,vjust=5),
          legend.spacing.x = unit(0.2, 'cm')
          # ,legend.key.height = unit(2, 'lines')
    ) + #labs(title="Age of the most recent records for vascular plants \nspecies across Mesoamerica", caption = "Occurrence data: GBIF\nTaxonomy: Kew's Checklist of Vascular Plants\nPhylogeny: Carta et al. (2022)") +
    NULL
  ggsave(filename = here("figures","Age_by_family.pdf"),plot = p)

  }


more_phylo_things <- function( ){
data = "NCBI_v1.txt"
data.table::fread(here("data",data),header = F) %>% as_tibble() -> accs
accs %>% pull() -> query_accs
taxonomizr::accessionToTaxa(accessions = query_accs, sqlFile = "data/acc_taxid") %>% as_tibble() %>% bind_cols(.,query_accs) -> aver
aver %>% rename_with(~c("TaxaID","accession")) %>% filter(!is.na(TaxaID)) -> aver
taxonomizr::getTaxonomy(aver %>%  pull(TaxaID),sqlFile = "data/nameNode.sqlite", desiredTaxa = c("class","family")) -> taxonomy
load(here("interim/Ages_60.Rdata"))
families %<>% ungroup() %>% mutate(family.x = ifelse(family.x == "Turneraceae","Passifloraceae",family.x ))%>%   mutate(family.x = ifelse(family.x == "Cochlospermaceae","Bixaceae",family.x )) %>% 
  mutate(family.x = ifelse(family.x == "Vivianiaceae","Francoaceae",family.x )) %>% 
  mutate(family.x = ifelse(family.x == "Stixaceae","Resedaceae",family.x )) %>% 
  mutate(family.x = ifelse(family.x == "Quiinaceae","Ochnaceae",family.x )) %>% 
  mutate(family.x = ifelse(family.x == "Hypseocharitaceae","Geraniaceae",family.x )) %>% 
  mutate(family.x = ifelse(family.x %in% c("Malesherbiaceae","Parnassiaceae"),"Celastraceae",family.x )) %>% 
  mutate(family.x = ifelse(family.x == "Tetradiclidaceae","Nitrariaceae",family.x )) %>% 
  mutate(family.x = ifelse(family.x == "Hydnoraceae","Hydnoraceae",family.x )) 
taxonomy %<>% as_tibble() 
a <- match(families$family.x,taxonomy$family)
families %<>% ungroup() %>% mutate(class = taxonomy$class[a],.after=1) 
families %<>% 
  mutate(class = replace_na(class,"Magnoliopsida")) %>% 
  mutate(class = ifelse(class %in% c("Cycadopsida","Gnetopsida") ,"Pinopsida",class ))
families %<>% separate(binomial,into=c("Genus","species"),sep="_",remove=FALSE) 
families %<>% 
  mutate(class = ifelse(class =="Pinopsida" ,"Gymnosperms",class))%>% 
  mutate(class = ifelse(class =="Lycopodiopsida","Lycopods",class))%>% 
  mutate(class = ifelse(class =="Magnoliopsida" ,"Angiosperms",class))%>% 
  mutate(class = ifelse(class =="Polypodiopsida" ,"Ferns",class))

families %>% filter(!is.na(LEVEL2_NAM)) %>%
  rowwise() %>% mutate(Records = nrow(data)) %>% select(-data) %>% 
  group_by(binomial) %>% summarise(Class= first(class),Family=first(family.x),Genus=first(Genus),Youngest=max(MostRecent),Records=sum(Records)) %>% group_by(Class) %>% mutate(Record_age = 2022 - Youngest) %>%
  relocate(binomial,.after=Genus) %>% 
  mutate(Class=factor(Class,levels=c("Lycopods","Ferns","Gymnosperms","Angiosperms"))) %>% 
  write.csv(.,file=here("MS/Supplementary_data.csv"))

families %>% filter(!is.na(LEVEL2_NAM)) %>% rowwise() %>% mutate(Records = nrow(data)) %>%  select(-data) %>% 
  group_by(Genus) %>% summarise(Class= first(class),Family=first(family.x),Genus=first(Genus),
                                Youngest=max(MostRecent),Richness=n(),Records=sum(Records)) %>% relocate(Genus,.after=Family) %>% mutate(Record_age = 2022 - Youngest) %>%
  mutate(Class=factor(Class,levels=c("Lycopods","Ferns","Gymnosperms","Angiosperms"))) %>% 
  write.csv(.,file=here("MS/Supplementary_data.csv"))
  

#ungroup() %>% 
  #filter(Record_age >= 22) %>% count() %>% mutate(n/32021)
  #summarize((1-ecdf(Record_age)(21))*100)

families %>% filter(!is.na(LEVEL2_NAM)) %>% 
  rowwise() %>% select(-data) %>% 
  group_by(binomial) %>% summarise(Youngest=max(MostRecent),class= first(class)) %>% group_by(class) %>% mutate(Record_age = 2022 - Youngest + 0.01) %>% mutate(class=factor(class,levels=c("Lycopods","Ferns","Gymnosperms","Angiosperms"))) %>% 
  ggplot(aes(y=class,x=Record_age,fill=class)) +
  geom_vline(xintercept = c(10,20,30,40,50,60,70,80,90,100),linetype="dashed",alpha=0.5)+
  geom_point(aes(color=class),position = position_jitter(height = 0.1),alpha=.8,size=0.7,shape=19) +
  scale_x_continuous(trans="log10",limits=c(1,100)) +
  ggdist::stat_halfeye(color="black",height = 0.8,alpha=0.8,position = position_nudge(y = 0.11)) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Kandinsky",direction=1,n=4)[c(4,1,3,2)]),guide=guide_legend(reverse = FALSE,ncol = 1,title.position="top",label.position="left", title.hjust = 0)) +
  scale_color_manual(values=(MetBrewer::met.brewer("Kandinsky",direction=1,n=4)[c(4,1,3,2)]),guide=guide_legend(reverse = FALSE,ncol = 1,title.position="top",label.position="left", title.hjust = 0)) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.title.y = element_blank(),
        axis.text =  element_text(family="EB Garamond",size = 8),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  NULL

ggsave(filename = here("figures/AgeSpecies_class.pdf"),last_plot())



load("interim/Allmetrics_60.Rdata")
datos
resolution = 60
datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -108,"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -97,"Mexico Southwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 20,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y > -83,"Panama",LEVEL3_NAM))
datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama"))))

families %>% left_join(.,datos %>% select(CellID,LEVEL3_NAM),by="CellID") %>% 
  filter(!is.na(LEVEL2_NAM)) %>% 
  rowwise() %>% select(-data) %>% 
  group_by(binomial,LEVEL3_NAM) %>% summarise(Youngest=max(MostRecent),class= first(class)) %>% group_by(class,LEVEL3_NAM) %>% mutate(Record_age = 2022 - Youngest) %>% mutate(class=factor(class,levels=c("Lycopods","Ferns","Gymnosperms","Angiosperms"))) %>% 
  ggplot(aes(y=LEVEL3_NAM,x=Record_age,fill=class,color=LEVEL3_NAM)) +
  facet_wrap(~class) +
 # geom_point(aes(color=class),position = position_jitter(height = 0.1),alpha=.5,size=0.7,shape=15)+
  xlim(c(0,100)) +
  #ggdist::stat_interval(.width = seq(0,1,0.25)) +
  ggdist::stat_pointinterval(.width = seq(0,1,0.25)) +
scale_fill_manual(values=(MetBrewer::met.brewer("Kandinsky",direction=1,n=4)[c(4,1,3,2)]),guide=guide_legend(reverse = FALSE,ncol = 1,title.position="top",label.position="left", title.hjust = 0)) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.title.y = element_blank(),
        text=element_text(family="EB Garamond"),
        axis.text =  element_text(family="EB Garamond",size = 8),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  NULL

ggsave(here("figures/AgeSpecies_RegionClass.pdf"),last_plot())
}
####


main_map <- function(){ 
  data="Joined_finalPOWdist.v2.Rdata"
  resolution=30
  load(here("output",data))
some %<>% filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution)
wgsrpd = "level2/level2.shp"
poly <- rgdal::readOGR(here("wgsrpd-master",wgsrpd))
poly <- poly[poly$LEVEL2_NAM %in% c("Mexico","Central America"),]
cat("....... buffering polygons from WGSRPD","\r")
spList = vector("list", length(poly))

for (i in 1:length(poly)) {
  cat(i,"\r")
  a <- rgeos::gBuffer(poly[i,], width = 0.5)
  a$LEVEL2_COD = poly[i,]$LEVEL2_COD
  a$LEVEL2_NAM = poly[i,]$LEVEL2_NAM
  spList[[i]] <- a
}
poly <- do.call("rbind", spList)

projection = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
res <- 60/resolution
g <- raster::raster(nrows=180*res,ncols=360*res,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% projectRaster(.,crs = CRS(projection)) %>% as(., 'SpatialPixels')
some %>% select(decimalLongitude,decimalLatitude) %>% 
  sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID"))) %>% select(CellID) %>% bind_cols(.,some) -> gridded

g@coords %>% as_tibble() %>% mutate(CellID=1:nrow(.)) -> coords
pointos <- coords %>% select(1:2) %>% as.data.frame(.) %>% SpatialPoints(.)
proj4string(pointos) <- projection
points_powo <- sp::over(pointos,poly) %>% as_tibble()
coords %<>% bind_cols(.,points_powo)

counts <- rnaturalearth::ne_coastline(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
xlimits = c(-125,-75)
ylimits = c(0,35)
###

####

gridded %>% distinct(CellID,Accepted_Name,.keep_all = TRUE) %>% 
  count(CellID) %>% left_join(.,coords,by="CellID") -> datos_gr
datos_gr %>% select(CellID,n) -> SR
gridded %>% #distinct(CellID,Accepted_Name,.keep_all = TRUE) %>% 
  count(CellID) %>% left_join(.,coords,by="CellID") -> datos_gr
datos_gr %>% select(CellID,x,y,n) %>% left_join(.,SR %>% rename(SR=n)) %>% 
  filter(x <= -75 & x >= -125 & y <= 35 & y>=0) %>% 
  ggplot(aes(x=n,y=SR),alpha=0.6) + 
  geom_point() +
  scale_x_continuous(trans="log10",breaks = c(1,10,100,1000,10000,100000),labels= scales::comma) +
  scale_y_continuous(trans="log10",breaks = c(1,10,100,1000,10000,100000),labels= scales::comma) +
  geom_smooth(color=MetBrewer::met.brewer("Hiroshige")[1])+
  labs(x="Ocurrence records",y="Species richness") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = c(1,0.5),
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.1),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  NULL
ggsave(filename = here("figures/SR_n.pdf"),last_plot())

##  richness
upperlim <-  roundUp(max((datos_gr$n),na.rm=T),to=100)
rich_bin = 1

main <- datos_gr %>% 
  ggplot(aes(x=x,y=y,fill=(n)),color="NA") +
  xlim(c(-200,-30)) + ylim( c(-70,80)) + 
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),#panel.border = element_rect(fill=NA),
        legend.position = c(1,0.5),
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.1),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        #legend.spacing = unit(0.001,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  geom_tile() +
  scale_fill_stepsn(colors=(MetBrewer::met.brewer("Hiroshige",direction=-1)),name="Species\nrichness",trans="log10",
                    #breaks=c(0,1,10,20,50,100,500,1000,2000,5000,10000,upperlim),
                    breaks=c(0,1,10,100,1000,5000,10000,20000,50000,75000,upperlim),
                    labels= scales::comma,
                    na.value="NA",
                    guide=guide_legend(reverse = FALSE,ncol = 1,title.position="top",label.position="left", title.hjust = 0)) + 
    labs(x="",y="") +
  geom_sf(data=counts,colour="black",fill=NA,size=0.1,alpha=0.5,inherit.aes = FALSE) +
  NULL


meso <- datos_gr %>% 
  ggplot(aes(x=x,y=y,fill=(n)),color="NA") +
  coord_cartesian(xlim=c(-125,-75),ylim= c(0,35)) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.1),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        #legend.spacing = unit(0.001,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  geom_tile() +
  scale_fill_stepsn(colors=(MetBrewer::met.brewer("Hiroshige",direction=-1)),name="Species\nrichness",trans="log10",
                   # breaks=c(0,1,10,20,50,100,500,1000,2000,5000,10000,upperlim),
                    breaks=c(0,1,10,100,1000,5000,10000,20000,50000,75000,upperlim),
                    labels= scales::comma,
                    na.value="NA",
                    guide=guide_legend(reverse = FALSE,ncol = 1,title.position="top",label.position="left", title.hjust = 0)) + 
  labs(x="",y="") +
  geom_sf(data=counts,colour="black",fill=NA,size=0.3,alpha=0.5,inherit.aes = FALSE) +
  #xlim(c(-125,-75)) + ylim( c(0,35)) + 
  coord_sf(
    xlim = c(-125,-75),
    ylim = c(0,35),
    expand = FALSE
  )+
  NULL
main  +
  geom_rect(
    xmin = -125,
    ymin = 0,
    xmax = -75,
    ymax = 35,
    fill = NA, 
    colour = "black",
    size = 0.4
  )  +
  inset_element(meso,left= -0.1,right=0.55,bottom= -0.1,top=0.6,on_top=TRUE) +
  NULL
  
ggsave(filename = here("figures/GeneralOccur_Map.pdf"),last_plot())


 }

plot_polygons_maps <- function(){
counts <- rnaturalearth::ne_countries(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
xlimits = c(-125,-75)
ylimits = c(0,35)

load(here("interim/Ages_30.Rdata"))
families %>% ungroup() %>% summarise(median(2022-MostRecent),quantile(2022-MostRecent,probs=c(0.10,0.90)))
load(here("interim/Ages_60.Rdata"))
families %>% ungroup() %>% summarise(median(2022-MostRecent),quantile(2022-MostRecent,probs=c(0.10,0.90)))
load(here("interim/Ages_15.Rdata"))
families %>% ungroup() %>% summarise(median(2022-MostRecent),quantile(2022-MostRecent,probs=c(0.10,0.90)))
}
###### 30 MINUTES ######
by_region_30 <- function(){ 
load(here("interim/Allmetrics_30.Rdata"))
resolution = 30
datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -106,"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -97,"Mexico Southwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -94,"Mexico Gulf",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 17.5,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -90,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y > -84,"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y > -85.5,"Honduras",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y > -88,"Nicaragua",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 15,"Belize",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 15,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
datos %>% filter(is.na(LEVEL3_NAM))
datos %>% write_csv(.,file=here("figures/temp.csv"))
aqui <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = (c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  group_by(LEVEL3_NAM) %>% 
  summarize(across(c(Completeness,age,seqs),
                   list(Median=function(x)median(x,na.rm=T),
                        Q0.1=function(x)quantile(x,na.rm=T,probs=0.1),
                        Q0.9=function(x)quantile(x,na.rm=T,probs=0.9)),.names = "{.fn}_{.col}")) %>% mutate(Resolution=resolution)
  save(aqui,file = "interim/summary_30.Rdata")
 

lm(age ~ 1 + Completeness,data=datos)  %>% summary()
lm(Completeness ~ 1 + SR,data=datos)  %>% summary()
plot(Completeness ~  Slope,data=datos)

age <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=age,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=3,size=0,geom="density_ridges_gradient",rel_min_height=0.01) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_blank(),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Age",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,100)) +
  NULL

completeness <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=Completeness,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=5,size=0,geom="density_ridges_gradient",rel_min_height=0.001) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=10),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_text(family="EB Garamond",size = 10),legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Completeness",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,100)) +
  NULL

seqs <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=seqs,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=0.02,size=0.1,geom="density_ridges_gradient",rel_min_height=0.001) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_blank(),legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Sequenced",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,1)) +
  NULL

xlimits = c(-123,-75)
ylimits = c(0,35)

map <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=Longitude.x,y=Latitude.x,color=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  geom_tile() + 
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        #legend.spacing = unit(0.001,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="",y="") +
  # standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
  geom_sf(data=counts,colour="black",fill=NA,size=0.2,inherit.aes = FALSE) +
  xlim(xlimits) + ylim(ylimits) +
  NULL

seq_comp <- datos %>% #mutate(Completeness = replace_na(Completeness,0)) %>% 
  group_by(LEVEL3_NAM) %>% #summarise(across(c(2,23,16),median,na.rm=T)) %>% 
  ggplot(aes(x=Completeness,y=seqs,fill=LEVEL3_NAM,color=LEVEL3_NAM))+
  geom_point(size=0.2,color="black",fill="black",alpha=0.2) +
  #geom_smooth(method="lm",se = TRUE,alpha=0.2,size=0.1) +
  geom_smooth(method="lm",se = FALSE,alpha=1,size=0.3) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        axis.title = element_text(family="EB Garamond",size=8),
        axis.text = element_text(family="EB Garamond",size=8),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  labs(x="Completeness",y="Sequenced") +
  coord_cartesian(ylim=c(0.2,1)) +
  #coord_flip() +
  NULL

age_comp <- datos %>% #mutate(Completeness = replace_na(Completeness,0)) %>% 
  group_by(LEVEL3_NAM) %>% #summarise(across(c(2,23,16),median,na.rm=T)) %>% 
  ggplot(aes(x=Completeness,y=age,fill=LEVEL3_NAM,color=LEVEL3_NAM))+
  geom_point(size=0.2,color="black",fill="black",alpha=0.2) +
  #geom_smooth(method="lm",se = TRUE,alpha=0.2,size=0.1) +
  geom_smooth(method="lm",se = FALSE,alpha=1,size=0.3) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        axis.title = element_text(family="EB Garamond",size=8),
        axis.text = element_text(family="EB Garamond",size=8),
        #legend.spacing = unit(0.001,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  labs(x="Completeness",y="Age") + 
  #coord_flip() +
  NULL



lm(seqs ~ 1 + Completeness,data=datos)  %>% broom.mixed::tidy(.)

layout <- "
ABC
DEF"
pp <- wrap_plots(design = layout,list(map,age_comp,seq_comp,completeness,age,seqs)) +
  plot_layout(byrow=T,nrow = 2,guides="keep",tag_level = 'new') & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"))

nombre <- paste0("Ridges_",resolution,".pdf")
ggsave(filename = here("figures",nombre),plot = pp)

}

bi_plot <- function(){ 
  roads <- rnaturalearth::ne_countries(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  counts <- rnaturalearth::ne_coastline(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  xlimits = c(-125,-75)
  ylimits = c(0,35)
  load(here("interim/Allmetrics_15.Rdata"))
  resolution = 15
  datos %>% filter(is.na(LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -106,"Mexico Northwest",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -98,"Mexico Southwest",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 25,"Mexico Northeast",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 16 & Longitude.y < -93,"Mexico Gulf",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y > -83,"Panama",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 11,"Costa Rica",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 14 & Longitude.y < -92,"Mexico Southeast",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 14 & Longitude.y < -88,"El Salvador",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 13,"Nicaragua",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 17 & Longitude.y > -88,"El Salvador",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 17 & Longitude.y < -90,"Mexico Southeast",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 18,"Belize",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 18 & Longitude.y > -90,"Mexico Southeast",LEVEL3_NAM))
  datos %>% filter(is.na(LEVEL3_NAM))
  datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) -> esto
  
  esto %>% ggplot(aes(y=(SR),x=seqs,color=LEVEL3_NAM,fill=LEVEL3_NAM)) +
geom_point() +
#ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=200,size=0.1,geom="density_ridges_gradient",rel_min_height=0.001) +
scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
theme(panel.background = element_blank(),panel.grid = element_blank(),
axis.line = element_line(),
legend.position = "none",
legend.background = element_rect(fill = NA),
legend.spacing.y = unit(c(0.2),"cm"),
legend.title = element_text(family="EB Garamond"),
legend.text = element_text(family="EB Garamond"),
axis.title.y = element_blank(),
axis.title.x = element_text(family="EB Garamond",size=10),
axis.text =  element_text(family="EB Garamond",size = 8),
#axis.text.y =  element_blank(),
legend.margin=margin(t=-25),
legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
labs(y="Species Richness",x="Temporal shortfall") +
#geom_hline(yintercept = 1:13,size=0.05)  +
#coord_cartesian(xlim=c(0,8000)) +
#xlim(c(0,8000)) +
scale_y_continuous(trans="log10") +
scale_x_continuous(limits=c(0,1)) +
    coord_polar()+
geom_smooth(method="loess",se = FALSE) +
NULL

  
  
library(biscale)
  colors <- bi_class(esto %>% filter(!is.na(Completeness),!is.na(SR)), x = Completeness, y = SR, style = "fisher", dim = 4)
colors %>% count(bi_class)
  legend <- bi_legend(pal = "PurpleGrn",dim = 4,
                      xlab = "Wallacean shortfall",ylab = "Species richness",
                      size = 10,rotate_pal = FALSE,flip_axes = FALSE,
breaks = bi_class_breaks(esto %>% filter(!is.na(Completeness),!is.na(SR)), x = Completeness, y = SR, style = "fisher", dim = 4))
colors %>% 
    ggplot(aes(x=Longitude.x,y=Latitude.x,fill=bi_class)) +
    #geom_point(shape=15) +
    xlim(xlimits) + ylim(ylimits) + 
    geom_tile(colour="white") +
    #scale_fill_manual(values=MetBrewer::met.brewer("Cassatt2",n=16))
    bi_scale_fill(pal = "PurpleGrn", dim = 4,rotate_pal=FALSE,flip_axes = FALSE) +
    #bi_scale_color(pal = "GrPink", dim = 4) +
    theme(panel.background = element_blank(),panel.grid = element_blank(),
          legend.position = "none",
          legend.background = element_rect(fill = NA),
          legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          axis.title = element_blank(),
          axis.text =  element_blank(),
          axis.line=element_blank(),
          axis.ticks = element_blank(),
          legend.margin=margin(t=-25),
          legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
    geom_sf(data=roads,inherit.aes = FALSE,colour="black",fill=NA,size=0.2) +
    #cowplot::draw_plot(legend,x=-120,y=8,width = 12,height = 12) +
    inset_element(legend,left=0.0,right=0.6,bottom=-0.2,top=0.5,on_top=FALSE) +
    NULL
  
colors <- bi_class(esto %>% filter(!is.na(SR)), x = age, y = SR, style = "jenks", dim = 4)
colors %>% count(bi_class)
legend <- bi_legend(pal = "PurpleGrn",
                    dim = 4,
                    xlab = "Temporal shortfall",
                    ylab = "Species richness",
                    size = 10,rotate_pal = FALSE,flip_axes = FALSE,
                    breaks = bi_class_breaks(esto %>% filter(!is.na(SR)), x = age, y = SR, style = "jenks", dim = 4))

colors %>% 
 ggplot(aes(x=Longitude.x,y=Latitude.x,fill=bi_class)) +
  #geom_point(shape=15) +
  xlim(xlimits) + ylim(ylimits) + 
  geom_tile(colour="white") +
  #scale_fill_manual(values=MetBrewer::met.brewer("Cassatt2",n=16))
  bi_scale_fill(pal = "PurpleGrn", dim = 4,rotate_pal=FALSE,flip_axes = FALSE) +
  #bi_scale_color(pal = "GrPink", dim = 4) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title = element_blank(),
        axis.text =  element_blank(),
        axis.line=element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  geom_sf(data=roads,inherit.aes = FALSE,colour="black",fill=NA,size=0.2) +
  #cowplot::draw_plot(legend,x=-120,y=8,width = 12,height = 12) +
  inset_element(legend,left=0.0,right=0.6,bottom=-0.2,top=0.5,on_top=FALSE) +
  NULL

colors <- bi_class(esto %>% filter(!is.na(SR)), x = seqs, y = SR, style = "jenks", dim = 4)
colors %>% count(bi_class)
legend <- bi_legend(pal = "PurpleGrn",
                    dim = 4,
                    ylab = "Species richness",
                    xlab = "Darwinian shortfall",
                    size = 10,rotate_pal = FALSE,flip_axes = FALSE,
                    breaks = bi_class_breaks(esto %>% filter(!is.na(SR)) , x = seqs, y = SR, style = "jenks", dim = 4))

colors %>% 
  ggplot(aes(x=Longitude.x,y=Latitude.x,fill=bi_class)) +
  #geom_point(shape=15) +
  xlim(xlimits) + ylim(ylimits) + 
  geom_tile(colour="white") +
  #scale_fill_manual(values=MetBrewer::met.brewer("Cassatt2",n=16))
  bi_scale_fill(pal = "PurpleGrn", dim = 4,rotate_pal=FALSE,flip_axes = FALSE) +
  #bi_scale_color(pal = "GrPink", dim = 4) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title = element_blank(),
        axis.text =  element_blank(),
        axis.line=element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  geom_sf(data=roads,inherit.aes = FALSE,colour="black",fill=NA,size=0.2) +
  #cowplot::draw_plot(legend,x=-120,y=8,width = 12,height = 12) +
  inset_element(legend,left=0.0,right=0.6,bottom=-0.2,top=0.5,on_top=FALSE) +
  NULL


}

######

##### 60 MINUTES ####
by_region_60 <- function(){
load(here("interim/Allmetrics_60.Rdata"))
  resolution = 60
  datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -108,"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -97,"Mexico Southwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 20,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y > -83,"Panama",LEVEL3_NAM))
datos %>% filter(is.na(LEVEL3_NAM))
datos %>% write_csv(.,file=here("figures/temp.csv"))

aqui <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = (c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  group_by(LEVEL3_NAM) %>% 
  summarize(across(c(Completeness,age,seqs),
                   list(Median=function(x)median(x,na.rm=T),
                        Q0.1=function(x)quantile(x,na.rm=T,probs=0.1),
                        Q0.9=function(x)quantile(x,na.rm=T,probs=0.9)),.names = "{.fn}_{.col}")) %>% mutate(Resolution=resolution)
save(aqui,file = "interim/summary_60.Rdata")



lm(age ~ 1 + Completeness,data=datos)  %>% summary()

age <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=age,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=3,size=0.2,geom="density_ridges_gradient",rel_min_height=0.01) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_blank(),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Age",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,100)) +
  NULL

completeness <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=Completeness,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=5,size=0.2,geom="density_ridges_gradient",rel_min_height=0.001) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=10),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_text(family="EB Garamond",size = 10),legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Completeness",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,100)) +
  NULL

seqs <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=seqs,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=0.02,size=0.2,geom="density_ridges_gradient",rel_min_height=0.001) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_blank(),legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Sequenced",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,1)) +
  NULL

xlimits = c(-123,-75)
ylimits = c(0,35)

map <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=Longitude.x,y=Latitude.x,color=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  geom_tile() + 
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        #legend.spacing = unit(0.001,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="",y="") +
  # standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
  geom_sf(data=counts,colour="black",fill=NA,size=0.2,inherit.aes = FALSE) +
  xlim(xlimits) + ylim(ylimits) +
  NULL

seq_comp <- datos %>% #mutate(Completeness = replace_na(Completeness,0)) %>% 
  group_by(LEVEL3_NAM) %>% #summarise(across(c(2,23,16),median,na.rm=T)) %>% 
  ggplot(aes(x=Completeness,y=seqs,fill=LEVEL3_NAM,color=LEVEL3_NAM))+
  geom_point(size=0.2,color="black",fill="black",alpha=0.2) +
  #geom_smooth(method="lm",se = TRUE,alpha=0.2,size=0.1) +
  geom_smooth(method="lm",se = FALSE,alpha=1,size=0.3) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        axis.title = element_text(family="EB Garamond",size=8),
        axis.text = element_text(family="EB Garamond",size=8),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  labs(x="Completeness",y="Sequenced") +
  coord_cartesian(ylim=c(0.2,1)) +
  #coord_flip() +
  NULL

age_comp <- datos %>% #mutate(Completeness = replace_na(Completeness,0)) %>% 
  group_by(LEVEL3_NAM) %>% #summarise(across(c(2,23,16),median,na.rm=T)) %>% 
  ggplot(aes(x=Completeness,y=age,fill=LEVEL3_NAM,color=LEVEL3_NAM))+
  geom_point(size=0.2,color="black",fill="black",alpha=0.2) +
  #geom_smooth(method="lm",se = TRUE,alpha=0.2,size=0.1) +
  geom_smooth(method="lm",se = FALSE,alpha=1,size=0.3) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        axis.title = element_text(family="EB Garamond",size=8),
        axis.text = element_text(family="EB Garamond",size=8),
        #legend.spacing = unit(0.001,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  labs(x="Completeness",y="Age") + 
  #coord_flip() +
  NULL



lm(seqs ~ 1 + Completeness,data=datos)  %>% broom.mixed::tidy(.)

layout <- "
ABC
DEF"
pp <- wrap_plots(design = layout,list(map,age_comp,seq_comp,completeness,age,seqs)) +
  plot_layout(byrow=T,nrow = 2,guides="keep",tag_level = 'new') & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"))

nombre <- paste0("Ridges_",resolution,".pdf")
ggsave(filename = here("figures",nombre),plot = pp)

}

#####

#### 15 MINUTES ####
by_region_15 <- function(){
  load(here("interim/Allmetrics_15.Rdata"))
  resolution = 15
  datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -106,"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -98,"Mexico Southwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 25,"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 16 & Longitude.y < -93,"Mexico Gulf",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y > -83,"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 11,"Costa Rica",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 14 & Longitude.y < -92,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 14 & Longitude.y < -88,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 13,"Nicaragua",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 17 & Longitude.y > -88,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 17 & Longitude.y < -90,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 18,"Belize",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 18 & Longitude.y > -90,"Mexico Southeast",LEVEL3_NAM))
datos %>% filter(is.na(LEVEL3_NAM))

datos %>% write_csv(.,file=here("figures/temp.csv"))

aqui <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = (c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  group_by(LEVEL3_NAM) %>% 
  summarize(across(c(Completeness,age,seqs),
                   list(Median=function(x)median(x,na.rm=T),
                        Q0.1=function(x)quantile(x,na.rm=T,probs=0.1),
                        Q0.9=function(x)quantile(x,na.rm=T,probs=0.9)),.names = "{.fn}_{.col}")) %>% mutate(Resolution=resolution)
save(aqui,file = "interim/summary_15.Rdata")



lm(seqs ~ 1 + Completeness, data=datos)  %>% summary()


age <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=age,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=3,size=0.2,geom="density_ridges_gradient",rel_min_height=0.01) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_blank(),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Age",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,100)) +
  NULL

completeness <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=Completeness,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=5,size=0.2,geom="density_ridges_gradient",rel_min_height=0.001) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=10),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_text(family="EB Garamond",size = 10),legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Completeness",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,100)) +
  NULL

seqs <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=seqs,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=0.02,size=0.2,geom="density_ridges_gradient",rel_min_height=0.001) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_blank(),legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Sequenced",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,1)) +
  NULL

xlimits = c(-123,-75)
ylimits = c(0,35)

map <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=Longitude.x,y=Latitude.x,color=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  geom_tile() + 
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        #legend.spacing = unit(0.001,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="",y="") +
  # standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
  geom_sf(data=counts,colour="black",fill=NA,size=0.2,inherit.aes = FALSE) +
  xlim(xlimits) + ylim(ylimits) +
  NULL

seq_comp <- datos %>% #mutate(Completeness = replace_na(Completeness,0)) %>% 
  group_by(LEVEL3_NAM) %>% #summarise(across(c(2,23,16),median,na.rm=T)) %>% 
  ggplot(aes(x=Completeness,y=seqs,fill=LEVEL3_NAM,color=LEVEL3_NAM))+
  geom_point(size=0.2,color="black",fill="black",alpha=0.2) +
  #geom_smooth(method="lm",se = TRUE,alpha=0.2,size=0.1) +
  geom_smooth(method="lm",se = FALSE,alpha=1,size=0.3) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        axis.title = element_text(family="EB Garamond",size=8),
        axis.text = element_text(family="EB Garamond",size=8),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  labs(x="Completeness",y="Sequenced") +
  coord_cartesian(ylim=c(0.2,1)) +
  #coord_flip() +
  NULL

age_comp <- datos %>% #mutate(Completeness = replace_na(Completeness,0)) %>% 
  group_by(LEVEL3_NAM) %>% #summarise(across(c(2,23,16),median,na.rm=T)) %>% 
  ggplot(aes(x=Completeness,y=age,fill=LEVEL3_NAM,color=LEVEL3_NAM))+
  geom_point(size=0.2,color="black",fill="black",alpha=0.2) +
  #geom_smooth(method="lm",se = TRUE,alpha=0.2,size=0.1) +
  geom_smooth(method="lm",se = FALSE,alpha=1,size=0.3) +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = "none",
        legend.background = element_rect(fill = NA),
        #legend.key.width = unit(0.01, "cm"),
        #legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),
        axis.title = element_text(family="EB Garamond",size=8),
        axis.text = element_text(family="EB Garamond",size=8),
        #legend.spacing = unit(0.001,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
  labs(x="Completeness",y="Age") + 
  #coord_flip() +
  NULL


layout <- "
ABC
DEF"
pp <- wrap_plots(design = layout,list(map,age_comp,seq_comp,completeness,age,seqs)) +
  plot_layout(byrow=T,nrow = 2,guides="keep",tag_level = 'new') & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"))

nombre <- paste0("Ridges_",resolution,".pdf")
ggsave(filename = here("figures",nombre),plot = pp)

}
#####

median_plots <- function(){ 
load("interim/summary_15.Rdata")
tabla <- aqui
load("interim/summary_30.Rdata")
tabla <- aqui %>% bind_rows(tabla,.)

load("interim/summary_60.Rdata")
tabla <- aqui %>% bind_rows(tabla,.)
tabla_glob <- tabla

tabla %>% relocate(Resolution,.after = 1) %>% mutate(across(2:ncol(.),round,2)) %>% 
  mutate(Completeness = paste0(Median_Completeness,"\n\U005B",Q0.1_Completeness,"\U2014",Q0.9_Completeness,"\U005D")) %>% 
  mutate(Age = paste0(Median_age,"\n\U005B",Q0.1_age,"\U2014",Q0.9_age,"\U005D")) %>% 
  mutate(Sequenced = paste0(Median_seqs,"\n\U005B",Q0.1_seqs,"\U2014",Q0.9_seqs,"\U005D")) %>% 
  select(LEVEL3_NAM,Resolution,Completeness,Age,Sequenced) #%>% 
  #flextable::flextable() %>% flextable::save_as_docx(.,path=here("figures/Table_2.docx"))

tabla <- tabla_glob %>% select(LEVEL3_NAM,Median_Completeness,Q0.1_Completeness,Q0.9_Completeness,Resolution) %>% pivot_longer(col=2:4,values_to = "val",names_to = "measure") %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels=rev(levels(LEVEL3_NAM)))) %>% group_by(LEVEL3_NAM)

comp <- ggplot() +
  geom_segment(data = tabla %>% filter(measure=="Median_Completeness")%>% group_by(LEVEL3_NAM) %>% slice_min(order_by = val) %>% ungroup() ,aes(x = 0, xend = val, y = LEVEL3_NAM, yend = LEVEL3_NAM),linetype = "dotted", size = 0.5, color = "gray80") +
  #geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==15),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.3),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==30),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.4),alpha=0.5) + 
 # geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.4),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.4),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==60),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.5),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.5),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.5),alpha=0.7) +
  geom_segment(data = tabla %>% filter(measure=="Median_Completeness") %>% group_by(LEVEL3_NAM) %>% summarise(start = range(val)[1], end = range(val)[2]) %>% ungroup(),aes(x = start, xend = end, y = LEVEL3_NAM, yend = LEVEL3_NAM),color = "gray80", size = 1)  +
  geom_point(data = tabla %>% filter(measure=="Median_Completeness"),  aes(val, LEVEL3_NAM, fill = factor(Resolution)),size = 3,shape=21) +
  geom_point(data = tabla %>% filter(measure=="Median_Completeness"),  aes(val, LEVEL3_NAM, fill = factor(Resolution)),size = 3,shape=21) +
  scale_fill_manual(values=MetBrewer::met.brewer("Derain",direction=1)[c(1,4,7)],name="Spatial resolution") +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        #legend.position = c(0.94,0.6),legend.direction = "vertical",
        legend.background = element_blank(),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"),
        legend.text = element_text(family="EB Garamond",size=8),
        legend.title = element_text(family="EB Garamond"),
        legend.key = element_rect(fill="transparent"),
        legend.spacing = unit(0.1,"cm"),
        legend.box.background = element_blank(),
        axis.text = element_text(family="EB Garamond",size=10),
        axis.title = element_text(family="EB Garamond",size=12),
        plot.title = element_text(family="EB Garamond",face="bold",size=20,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) + 
 # guides(color = guide_colorbar(title.position="top",title.hjust = 0, ticks = F),alpha="none",size = "none"
         #guide_legend(reverse=T,ticks=F,override.aes = list(alpha=1,fill="grey50"))) +
  labs(x="Median Completeness",y="",title="") +
  NULL

tabla <- tabla_glob %>% select(LEVEL3_NAM,Median_age,Q0.1_age,Q0.9_age,Resolution) %>% pivot_longer(col=2:4,values_to = "val",names_to = "measure") %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels=rev(levels(LEVEL3_NAM)))) %>% group_by(LEVEL3_NAM)

age <- ggplot() +
  geom_segment(data = tabla %>% filter(measure=="Median_age")%>% group_by(LEVEL3_NAM) %>% slice_min(order_by = val) %>% ungroup() ,aes(x = 0, xend = val, y = LEVEL3_NAM, yend = LEVEL3_NAM),linetype = "dotted", size = 0.5, color = "gray80") +
  #geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==15),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.3),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==30),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.4),alpha=0.5) + 
  # geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.4),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.4),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==60),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.5),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.5),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.5),alpha=0.7) +
  geom_segment(data = tabla %>% filter(measure=="Median_age") %>% group_by(LEVEL3_NAM) %>% summarise(start = range(val)[1], end = range(val)[2]) %>% ungroup(),aes(x = start, xend = end, y = LEVEL3_NAM, yend = LEVEL3_NAM),color = "gray80", size = 1)  +
  geom_point(data = tabla %>% filter(measure=="Median_age"),  aes(val, LEVEL3_NAM, fill = factor(Resolution)),size = 3,shape=21) +
  geom_point(data = tabla %>% filter(measure=="Median_age"),  aes(val, LEVEL3_NAM, fill = factor(Resolution)),size = 3,shape=21) +
  scale_fill_manual(values=MetBrewer::met.brewer("Derain",direction=1)[c(1,4,7)],name="Spatial resolution") +  theme(panel.background = element_blank(),
        axis.line = element_line(),
        #legend.position = c(0.94,0.6),legend.direction = "vertical",
        legend.background = element_blank(),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"),
        legend.text = element_text(family="EB Garamond",size=8),
        legend.title = element_text(family="EB Garamond"),
        legend.key = element_rect(fill="transparent"),
        legend.spacing = unit(0.1,"cm"),
        legend.box.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(family="EB Garamond",size=10),
        axis.title = element_text(family="EB Garamond",size=12),
        plot.title = element_text(family="EB Garamond",face="bold",size=20,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) + 
  # guides(color = guide_colorbar(title.position="top",title.hjust = 0, ticks = F),alpha="none",size = "none"
  #guide_legend(reverse=T,ticks=F,override.aes = list(alpha=1,fill="grey50"))) +
  labs(x="Median age",y="",title="") +
  NULL

tabla <- tabla_glob %>% select(LEVEL3_NAM,Median_seqs,Q0.1_seqs,Q0.9_age,Resolution) %>% pivot_longer(col=2:4,values_to = "val",names_to = "measure") %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels=rev(levels(LEVEL3_NAM)))) %>% group_by(LEVEL3_NAM)


seqs <- ggplot() +
  geom_segment(data = tabla %>% filter(measure=="Median_seqs")%>% group_by(LEVEL3_NAM) %>% slice_min(order_by = val) %>% ungroup() ,aes(x = 0, xend = val, y = LEVEL3_NAM, yend = LEVEL3_NAM),linetype = "dotted", size = 0.5, color = "gray80") +
  #geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==15),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.3),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==30),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.4),alpha=0.5) + 
  # geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.4),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.4),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==60),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.5),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.5),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.5),alpha=0.7) +
  geom_segment(data = tabla %>% filter(measure=="Median_seqs") %>% group_by(LEVEL3_NAM) %>% summarise(start = range(val)[1], end = range(val)[2]) %>% ungroup(),aes(x = start, xend = end, y = LEVEL3_NAM, yend = LEVEL3_NAM),color = "gray80", size = 1)  +
  geom_point(data = tabla %>% filter(measure=="Median_seqs"),  aes(val, LEVEL3_NAM, fill = factor(Resolution)),size = 3,shape=21) +
  geom_point(data = tabla %>% filter(measure=="Median_seqs"),  aes(val, LEVEL3_NAM, fill = factor(Resolution)),size = 3,shape=21) +
  scale_fill_manual(values=MetBrewer::met.brewer("Derain",direction=1)[c(1,4,7)],name="Spatial resolution") +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        #legend.position = c(0.94,0.6),legend.direction = "vertical",
        legend.background = element_blank(),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"),
        legend.text = element_text(family="EB Garamond",size=8),
        legend.title = element_text(family="EB Garamond"),
        legend.key = element_rect(fill="transparent"),
        legend.spacing = unit(0.1,"cm"),
        legend.box.background = element_blank(),
        axis.text = element_text(family="EB Garamond",size=10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(family="EB Garamond",size=12),
        plot.title = element_text(family="EB Garamond",face="bold",size=20,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) + 
  # guides(color = guide_colorbar(title.position="top",title.hjust = 0, ticks = F),alpha="none",size = "none"
  #guide_legend(reverse=T,ticks=F,override.aes = list(alpha=1,fill="grey50"))) +
  labs(x="Median Sequenced",y="",title="") +
  NULL


layout <- "
ABC"
pp <- wrap_plots(design = layout,list(comp,age,seqs)) +
  plot_layout(byrow=T,nrow = 2,guides="collect",tag_level = 'new') & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"),legend.position = "bottom")
pp

nombre <- paste0("Dumbbells.pdf")
ggsave(filename = here("figures",nombre),plot = pp)

}

hotcold_spots <- function(){
base <- raster("output/raster15_CompletenessG.tif")
base[base==999] <- NA
base %>% getValues(.) -> comp; length(comp)
base <- raster("output/raster15_ageG.tif")
base %>% getValues(.) -> age; length(age)
base <- raster("output/raster15_seqsG.tif")
base %>% getValues(.) -> seqs; length(seqs)
hfp <- raster(here("data/LUL/LULC_1920__CMIP6.1.tif"))
hfp <- crop(hfp,extent(base))
getValues(hfp) -> lul_1920; length(lul_1920)

hfp <- raster(here("data/LUL/LULC_1950__CMIP6.1.tif"))
hfp <- crop(hfp,extent(base))
getValues(hfp) -> lul_1950; length(lul_1950)

hfp <- raster(here("data/LUL/LULC_1980__CMIP6.1.tif"))
hfp <- crop(hfp,extent(base))
getValues(hfp) -> lul_1980; length(lul_1980)

hfp <- raster(here("data/LUL/LULC_2010__CMIP6.1.tif"))
hfp <- crop(hfp,extent(base))
getValues(hfp) -> lul_2010; length(lul_2010)

hfp <- raster(here("data/Human_footprint/HFP_gw_MESO_15ok.tif"))
getValues(hfp) -> HFP; length(HFP)

load("interim/Allmetrics_15.Rdata")
base <- raster(here("data/base_rasters/baseRaster_15.tif"))
datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -106,"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -98,"Mexico Southwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 25,"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 16 & Longitude.y < -93,"Mexico Gulf",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y > -83,"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 11,"Costa Rica",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 14 & Longitude.y < -92,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 14 & Longitude.y < -88,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 13,"Nicaragua",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 17 & Longitude.y > -88,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 17 & Longitude.y < -90,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 18,"Belize",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 18 & Longitude.y > -90,"Mexico Southeast",LEVEL3_NAM))
datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(lvl3=as.factor(LEVEL3_NAM),.after=1)
rasterize(datos[,c("Longitude.y","Latitude.y")],base,
 field=as.numeric(datos$lvl3)) %>% crop(.,extent(hfp)) %>% getValues(.) -> regions; length(regions)

aver <- raster(here("output/sumaHotspot.tif"))
aver %>% crop(.,extent(hfp)) %>% getValues(.) -> zonation
aver %>% crop(.,extent(hfp)) %>% coordinates(.) -> coords

tibble(x=coords[,1],y=coords[,2],zone=zonation,comp=comp,age=age,seqs=seqs,lul_2010=lul_2010,lul_1980=lul_1980,lul_1950=lul_1950,lul_1920=lul_1920,HFP=HFP,regions=regions)  %>% 
  mutate(across(starts_with("lul"), ~case_when(.x %in% c(1:6, 9, 12) ~"Human use",.x %in% c(7,8)~"Primary",.x %in% c(10,11)~"Secondary"))) %>% 
  filter(!is.na(HFP)) ->  tiba

tiba <- tiba %>% mutate(regions=factor(regions,levels=1:13))
levels(tiba$regions) <- levels(datos$lvl3)
tiba %>% mutate(spot_comp = ifelse(comp <= -1.96,"Coldspot",ifelse(comp >= 1.96,"Hotspot",NA))) %>% mutate(spot_age = ifelse(age <= -1.96,"Coldspot",ifelse(age >= 1.96,"Hotspot",NA))) %>% mutate(spot_seqs = ifelse(seqs <= -1.96,"Coldspot",ifelse(seqs >= 1.96,"Hotspot",NA))) -> datas

datas %>% write_csv(.,file=here("figures/temp.csv"))

datas %>% filter(!is.na(regions)) %>% select(regions,spot_seqs,spot_comp,spot_age) %>% 
pivot_longer(cols=2:4,names_to = "variable",values_to ="type" ) %>% 
  group_by(regions,variable,type) %>% 
  count() %>% group_by(regions) %>% mutate(prop=n/sum(n)*100) %>% select(-n) %>% 
  filter(!is.na(type)) %>% 
pivot_wider(names_from = c(type,variable),values_from = prop) %>% 
  mutate(across(everything(),replace_na,0)) %>% flextable::flextable() %>% flextable::save_as_docx(path=here("figures/table_3.docx"))

datas %>% select(-HFP) %>%  pivot_longer(cols = starts_with("lul"),names_to = "lul",values_to = "type") %>% filter(lul=="lul_2010") %>% filter(!is.na(type)) %>% 
  group_by(type,spot_seqs) %>% count() %>% rename_with(~c("type","spot","n")) %>% 
  group_by(spot) %>% 
  mutate(n= n / sum(n)) %>% 
  filter(!is.na(spot)) %>% 
  ggplot(aes(y=n,x=spot,fill=factor(type,levels=c("Primary","Secondary","Human use")))) + #facet_wrap(~lul) +
  geom_bar(stat="identity",position="dodge") + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = "bottom",
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=12),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.x =  element_text(family="EB Garamond",size = 12),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.6,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="vertical")  + labs(title="",y="Percent area",x="") +
  scale_y_continuous(limits=c(0,0.8)) +
  scale_fill_manual(values = MetBrewer::met.brewer("Kandinsky")[c(2:4)],name="") +
  NULL


datas %>% summarise(quantile(HFP,probs = c(0.15,0.85)))
human_cat <- datas %>% mutate(CatHFP = ifelse(HFP <= quantile(HFP,probs=c(0.15),na.rm=T),"Lowest impact",ifelse(HFP >= quantile(HFP,probs=c(0.85),na.rm=T),"Highest impact",
 ifelse(HFP >= quantile(HFP,probs=c(0.6),na.rm=T)& HFP < quantile(HFP,probs=c(0.85),na.rm=T),"High impact",
  ifelse(HFP > quantile(HFP,probs=c(0.15),na.rm=T) & HFP <= quantile(HFP,probs=c(0.4),na.rm=T),"Low impact","Intermediate impact"))))) %>% 
  select(-HFP) %>%
  pivot_longer(cols = starts_with("Cat"),names_to = "lul",values_to = "type") %>% 
  filter(!is.na(type)) %>% 
  group_by(lul,type,spot_seqs) %>% count() %>% group_by(spot_seqs,lul) %>% 
  mutate(n= n / sum(n)) %>% 
  filter(!is.na(spot_seqs)) %>% 
  filter(type %in% c("Highest impact","Lowest impact")) %>% 
  ggplot(aes(y=n,x= spot_comp,fill=factor(type,levels=c("Highest impact",
                                                  #"High impact",
                                                  #"Intermediate impact",
                                                  #"Low impact",
                                                  "Lowest impact"
                                                  )))) +
  geom_bar(stat="identity",position="dodge") + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = "bottom",
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=12),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.x =  element_text(family="EB Garamond",size = 12),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.6,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="vertical")  + labs(title="",y="Percent area",x="") +
  scale_fill_manual(values = c(MetBrewer::met.brewer("Kandinsky")[c(3:4)]),name="") +
  NULL
human_cat
#wrap_plots(list(human_cat,landuse))
ggsave(filename = here("figures/Comp_Antro.pdf"),last_plot())

library(UpSetR)
datas %>% select(starts_with("spot")) %>%
  mutate(across(everything(),replace_na,"Neutral")) %>% group_by(spot_comp,spot_age,spot_seqs) %>% count() %>% filter(spot_comp!="Coldspot",spot_age!="Coldspot",spot_seqs!="Coldspot") %>% ungroup() %>% summarise( (sum(n)-1024)/3436 )

datas %>% select(starts_with("spot")) %>% mutate(ID=1:nrow(.),.before=1) %>% pivot_longer(cols = -1,names_to = "type",values_to = "spots") %>% 
  mutate(type=sub("^spot_","",type)) %>% 
  mutate(val=1) %>% pivot_wider(names_from = c(spots,type),values_from = val) %>%  
mutate(across(everything(),replace_na,0)) %>% select(-starts_with("NA_")) -> forup
upset(as.data.frame(forup),nsets = 6,sets = rev(c("Hotspot_comp","Hotspot_age","Hotspot_seqs","Coldspot_comp","Coldspot_age","Coldspot_seqs")),keep.order = TRUE)


ggsave(filename = here("figures/upset_hotspots.pdf"),last_plot())

}

base <- raster("output/raster15_CompletenessG.tif")
base[base==999] <- NA
base %>% getValues(.) -> comp; length(comp)
base <- raster("output/raster15_ageG.tif")
base %>% getValues(.) -> age; length(age)
base <- raster("output/raster15_seqsG.tif")
base %>% getValues(.) -> seqs; length(seqs)
hfp <- raster(here("data/LUL/LULC_1920__CMIP6.1.tif"))
hfp <- crop(hfp,extent(base))
getValues(hfp) -> lul_1920; length(lul_1920)

hfp <- raster(here("data/LUL/LULC_1950__CMIP6.1.tif"))
hfp <- crop(hfp,extent(base))
getValues(hfp) -> lul_1950; length(lul_1950)

hfp <- raster(here("data/LUL/LULC_1980__CMIP6.1.tif"))
hfp <- crop(hfp,extent(base))
getValues(hfp) -> lul_1980; length(lul_1980)

hfp <- raster(here("data/LUL/LULC_2010__CMIP6.1.tif"))
hfp <- crop(hfp,extent(base))
getValues(hfp) -> lul_2010; length(lul_2010)

hfp <- raster(here("data/Human_footprint/HFP_gw_MESO_15ok.tif"))
getValues(hfp) -> HFP; length(HFP)

load("interim/Allmetrics_15.Rdata")
base <- raster(here("data/base_rasters/baseRaster_15.tif"))
datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -106,"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y < -98,"Mexico Southwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 25,"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 16 & Longitude.y < -93,"Mexico Gulf",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude.y > -83,"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 11,"Costa Rica",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 14 & Longitude.y < -92,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 14 & Longitude.y < -88,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 13,"Nicaragua",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 17 & Longitude.y > -88,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 17 & Longitude.y < -90,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y < 18,"Belize",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude.y > 18 & Longitude.y > -90,"Mexico Southeast",LEVEL3_NAM))
datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(lvl3=as.factor(LEVEL3_NAM),.after=1)
rasterize(datos[,c("Longitude.y","Latitude.y")],base,
          field=as.numeric(datos$lvl3)) %>% crop(.,extent(hfp)) %>% getValues(.) -> regions; length(regions)
datos %<>% mutate(across(SR,replace_na,0))
rasterize(datos[,c("Longitude.y","Latitude.y")],base,
          field=datos$SR) %>% crop(.,extent(hfp)) %>% getValues(.) -> SR; length(SR)

aver <- raster(here("output/sumaColdspot.tif"))
aver %>% crop(.,extent(hfp)) %>% getValues(.) -> zonation
aver %>% crop(.,extent(hfp)) %>% coordinates(.) -> coords

tibble(x=coords[,1],y=coords[,2],SR=SR,zone=zonation,comp=comp,age=age,seqs=seqs,lul_2010=lul_2010,lul_1980=lul_1980,lul_1950=lul_1950,lul_1920=lul_1920,HFP=HFP,regions=regions)  %>% 
  mutate(across(starts_with("lul"), ~case_when(.x %in% c(1:6, 9, 12) ~"Human use",.x %in% c(7,8) ~ "Primary",.x %in% c(10,11) ~ "Secondary"))) %>% 
  filter(!is.na(HFP)) ->  tiba
tiba <- tiba %>% mutate(regions=factor(regions,levels=1:13))
levels(tiba$regions) <- levels(datos$lvl3)
tiba %>% mutate(spot_comp = ifelse(comp <= -1.96,"Coldspot",ifelse(comp >= 1.96,"Hotspot",NA))) %>% mutate(spot_age = ifelse(age <= -1.96,"Coldspot",ifelse(age >= 1.96,"Hotspot",NA))) %>% mutate(spot_seqs = ifelse(seqs <= -1.96,"Coldspot",ifelse(seqs >= 1.96,"Hotspot",NA))) -> datas

xlimits = c(-125,-75)
ylimits = c(0,35)
colas <- sample(1:7,7,replace=FALSE)
datas %>% mutate(zone=ifelse(zone==0,NA,
                             ifelse(zone==1,"Comp",
                             ifelse(zone==10,"Seq",
                             ifelse(zone==11,"Comp + Seq",
                             ifelse(zone==20,"Age",
                             ifelse(zone==21,"Comp + Age",
                             ifelse(zone==30,"Seq + Age",
                             ifelse(zone==31,"Comp + Seq + Age",NA))))))))) %>% 
  filter(!is.na(zone)) %>% 
  mutate(zone=factor(zone,levels=c("Comp","Age","Seq","Comp + Age","Comp + Seq","Seq + Age","Comp + Seq + Age"))) %>% 
  ggplot(aes(x=x,y=y,color=zone)) + 
  geom_point(size=1) +
  coord_cartesian(xlim=c(-120,-75),ylim = c(5,35)) +
  scale_fill_manual(values=MetBrewer::met.brewer("Archambault",direction=1,n=7)[colas]) +
  scale_color_manual(values=MetBrewer::met.brewer("Archambault",direction=1,n=7)[colas]) +
  theme(panel.background = element_blank(),
        #panel.grid = element_line(size=0.2,color="grey30"),
        panel.border = element_rect(fill="NA"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        #legend.position = c(0.94,0.6),legend.direction = "vertical",
        legend.background = element_blank(),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"),
        legend.text = element_text(family="EB Garamond",size=8),
        legend.title = element_text(family="EB Garamond"),
        legend.key = element_rect(fill="transparent"),
        legend.spacing = unit(0.1,"cm"),
        legend.box.background = element_blank(),
        legend.position=c(0.1,0.3),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(family="EB Garamond",face="bold",size=20,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
  geom_sf(data=roads,colour="black",fill=NA,size=0.2,inherit.aes = FALSE) +
  xlim(xlimits) + ylim(ylimits) +
  NULL

datas %>% select(SR,spot_comp) %>% filter(SR!=0) %>% rename_with(~c("SR","spot")) %>% 
  filter(!is.na(spot)) %>% 
ggplot(aes(y=SR,side=ifelse(spot=="Hotspot","right","left"),color=spot)) +
  #geom_jitter(alpha=0.3,size=0.5) +
ggdist::stat_dist_dots(shape=15,binwidth=0.033) +
  #geom_boxplot(fill="NA",outlier.stroke = 0,outlier.size = 0) + 
  scale_y_continuous(trans="log10")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = "right",
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=12),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.x =  element_text(family="EB Garamond",size = 12),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.6,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="vertical")  + 
  labs(title="",y="Species richness",x="Wallacean shortfall") +
  
NULL




#ggsave(filename = here("figures/zone_cold.pdf"),last_plot())

datas %>% mutate(zone=ifelse(zone==0,"nonspot","spot")) %>% 
  mutate(zone=ifelse(is.na(zone),"nonspot",zone)) %>% 
mutate(CatHFP = ifelse(HFP <= quantile(HFP,probs=c(0.2),na.rm=T),"Lowest impact",ifelse(HFP >= quantile(HFP,probs=c(0.8),na.rm=T),"Highest impact",ifelse(HFP >= quantile(HFP,probs=c(0.6),na.rm=T)& HFP < quantile(HFP,probs=c(0.8),na.rm=T),"High impact",ifelse(HFP > quantile(HFP,probs=c(0.2),na.rm=T) & HFP <= quantile(HFP,probs=c(0.4),na.rm=T),"Low impact","Intermediate impact"))))) %>% 
  group_by(zone,CatHFP) %>% count() %>% group_by(zone) %>% 
  mutate(n= n / sum(n)) %>% 
  #filter(!is.na(lul_2010)) %>% 
  filter(zone=="spot") %>% 
  filter(CatHFP %in% c("Highest impact","Lowest impact")) %>% 
  ggplot(aes(y=n,x= zone,
             #fill=factor(lul_2010, levels=c("Primary","Secondary","Human use"))
             fill=factor(CatHFP,levels=c("Highest impact","Lowest impact"))
             )) +
  geom_bar(stat="identity",position="dodge") + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = "bottom",
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=12),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.x =  element_text(family="EB Garamond",size = 12),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.6,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="vertical")  + labs(title="",y="Percent area",x="") +
  scale_y_continuous(limits=c(0,0.3)) +
  scale_fill_manual(values = c(MetBrewer::met.brewer("Kandinsky")[c(1,5:4)]),name="") +
  NULL





localG <- function(){
  xlimits = c(-125,-75)
  ylimits = c(0,35)
resolution = 15 
distance = 60
var = "seqs"
hfp <- raster(here("data/Human_footprint/HFP_gw_meso_rcl_15.tif"))
load("interim/Allmetrics_15.Rdata")
base <- raster("data/base_rasters/baseRaster_15.tif")
datos %<>% filter(!is.na(seqs))
nb <- dnearneigh(as.matrix(datos[,c("Longitude.y","Latitude.y")]),0,distance,longlat = TRUE)
aver <- spdep::localG(1-datos[[glue::glue("{var}")]],nb2listw(nb,style = "B",zero.policy=TRUE))
datos$localG <- as.vector(aver)
prefix = paste0("raster",resolution,"_",var,".tif")
rasterize(datos[,c("Longitude.y","Latitude.y")],base,field=datos[[glue::glue("{var}")]]) %>% crop(.,extent(hfp)) %>%   writeRaster(.,filename = here("output",prefix),overwrite=TRUE)
prefix = paste0("raster",resolution,"_",var,"G",".tif")
#datos$localG[which(is.na(datos$localG))] <- 999
rasterize(datos[,c("Longitude.y","Latitude.y")],base,field=datos$localG) %>% crop(.,extent(hfp)) %>% writeRaster(.,filename = here("output",prefix),overwrite=TRUE)
prefix = paste0("raster",resolution,"_","Slope",".tif")
rasterize(datos[,c("Longitude.y","Latitude.y")],base,field=datos$Slope) %>% crop(.,extent(hfp)) %>%   writeRaster(.,filename = here("output",prefix))

range(datos$localG)
#datos %<>% mutate(across(localG,replace_na,999))
datos %>% ggplot(aes(x=Longitude.y,y=Latitude.y,color=seqs,fill=seqs)) +
  geom_point(size=1) +
  coord_cartesian(xlim=c(-120,-75),ylim = c(5,35)) +
  scale_fill_stepsn(colours=MetBrewer::met.brewer("Hiroshige",direction=-1),breaks=c(-3,-2.5,-1.96,-1.65,0,1.65,1.96,2.5,3.0),limits=c(-9,9)) +
  scale_color_stepsn(colours=MetBrewer::met.brewer("Hiroshige",direction=-1),breaks=c(-3,-2.5,-1.96,-1.65,0,1.65,1.96,2.5,3.0),limits=c(-9,9)) +
  theme(panel.background = element_blank(),
        #panel.grid = element_line(size=0.2,color="grey30"),
        panel.border = element_rect(fill="NA"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        #legend.position = c(0.94,0.6),legend.direction = "vertical",
        legend.background = element_blank(),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"),
        legend.text = element_text(family="EB Garamond",size=8),
        legend.title = element_text(family="EB Garamond"),
        legend.key = element_rect(fill="transparent"),
        legend.spacing = unit(0.1,"cm"),
        legend.box.background = element_blank(),
        legend.position=c(0.1,0.3),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(family="EB Garamond",face="bold",size=20,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
  geom_sf(data=counts,colour="black",fill=NA,size=0.2,inherit.aes = FALSE) +
  xlim(xlimits) + ylim(ylimits)



}

map_hotspots <- function(){ 
counts <- rnaturalearth::ne_countries(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
xlimits = c(-125,-75)
ylimits = c(0,35)
resolution = 15 
distance = 60
var = "seqs"
hfp <- raster("output/raster15_ageG.tif")
load("interim/Allmetrics_15.Rdata")
base <- raster("data/base_rasters/baseRaster_15.tif")
datos %<>% filter(!is.na(seqs))
nb <- dnearneigh(as.matrix(datos[,c("Longitude.y","Latitude.y")]),0,distance,longlat = TRUE)
aver <- spdep::localG(1-datos[[glue::glue("{var}")]],nb2listw(nb,style = "B",zero.policy=TRUE))
datos$localG <- as.vector(aver)
datos %>% mutate(spots = ifelse(localG >= 1.65,"HC",ifelse(localG <= -1.65,"HC","N"))) %>% 
  ggplot(aes(x=Longitude.y,y=Latitude.y,color=localG,alpha=spots)) +
  scale_alpha_discrete(range=c(1,0.0)) +
  geom_point(size=0.7) +
  scale_color_stepsn(colors=(MetBrewer::met.brewer("Hiroshige",direction=-1)),n.breaks=20,name="local G",limits=c(-9,9)) +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        legend.position = c(0.2,0.3),
        axis.ticks=element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"),
        legend.text = element_text(family="EB Garamond",size=8),
        legend.title = element_text(family="EB Garamond"),
        legend.key = element_rect(fill="transparent"),
        legend.spacing = unit(0.1,"cm"),
        legend.box.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(family="EB Garamond",face="bold",size=20,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) + 
  geom_sf(data=counts,colour="black",fill=NA,size=0.1,alpha=0.5,inherit.aes = FALSE) +
  xlim(xlimits) + ylim(ylimits) + 
  guides(alpha="none") +
NULL

ggsave(filename = here("figures/Seqs_hotcoldspots.pdf"),last_plot())
}






library(ggalluvial)
as_tibble(Titanic)
datos %>% mutate(across(Completeness, replace_na,0)) %>% 
  mutate(CompCat=cut(Completeness,breaks=c(0,50,75,90,100),include.lowest = TRUE, 
                             labels=c("Incomplete","Fairly\nIncomplete","Fairly\nComplete","Complete")),.before=x) %>% 
          mutate(SeqCat=cut(seqs,breaks=c(0,0.50,0.7,1),include.lowest = TRUE, 
                            labels=c("< 0.5","0.5-0.7",">  0.7")),.before=x) %>% 
          mutate(AgeCat=cut(age,breaks=c(100,40,20,0),include.lowest = TRUE, 
                            labels=c("> 40","20-40","< 20")),.before=x) %>% 
  group_by(CompCat,SeqCat,AgeCat) %>% summarise(n=n()) %>% ungroup() %>% 
  ggplot(aes(y=n,axis1 = CompCat, axis3 = SeqCat, axis2 = AgeCat)) +
   #geom_alluvium(aes(fill = CompCat),width = 0,knot.pos=1,alpha=0.9,reverse = TRUE,curve_type="sine",knot.prop = TRUE) +
  geom_flow(aes(fill = CompCat),curve_type = "sine") +
  scale_x_discrete(limits = c("Completeness", "Age (years)","Sequencing (%)"),expand = c(.2, .05)) +
  geom_stratum(width = 0.4,aes(fill=CompCat),size=0.2) +
  scale_fill_manual(values=MetBrewer::met.brewer("Archambault",n=6)) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum),family="EB Garamond"),color="black")  +
  theme(panel.background = element_blank(),
        #panel.grid = element_line(size=0.2,color="grey30"),
        #panel.border = element_rect(),
        axis.line = element_line(),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"),
        legend.text = element_text(family="EB Garamond",size=8),
        legend.title = element_text(family="EB Garamond"),
        legend.key = element_rect(fill="transparent"),
        legend.spacing = unit(0.1,"cm"),
        legend.box.background = element_blank(),
        axis.text = element_text(family="EB Garamond",size=16),
        axis.title = element_text(family="EB Garamond",size=18),
        plot.title = element_text(family="EB Garamond",face="bold",size=20,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5))

  # ternary_plot %>% 
  # ggplot(aes(fill=completeness,x=age,y=(sr),size=(n),alpha=completeness),color="grey50") +
  #   geom_point(shape=21,stroke=0.1) +
  #   scale_size_binned(name="Size of\ncollection") +
  #   scale_fill_stepsn(colors=(MetBrewer::met.brewer("Hiroshige",direction=-1)),name="Completeness",n.breaks=10,labels=function(x) sprintf("%.2f", (x))) +
  #   scale_y_continuous(trans="log",labels= function(x) sprintf("%.0f", x)) + 
  #   scale_x_continuous(trans="identity",labels= function(x) sprintf("%.0f", x)) +
  #   theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = c(0.94,0.6),legend.direction = "vertical",
  #         legend.background = element_blank(),
  #         legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"),
  #         legend.text = element_text(family="EB Garamond",size=8),
  #         legend.title = element_text(family="EB Garamond"),
  #         legend.key = element_rect(fill="transparent"),
  #         legend.spacing = unit(0.1,"cm"),
  #         legend.box.background = element_blank(),
  #         axis.text = element_text(family="EB Garamond"),
  #         axis.title = element_text(family="EB Garamond",size=12),
  #         plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),
  #         plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)
  #   ) + guides(color = guide_colorbar(title.position="top", title.hjust = 0, ticks = F),alpha="none",size=guide_legend(reverse=T,ticks=F,override.aes = list(alpha=1,fill="grey50"))) +
  #   labs(y="Species Richness (log scale)",x="Median age of collections",title="Digitally accesible botanical knowledge in the Americas",subtitle=paste0(formatC(final, big.mark=",")," occurrence records from GBIF")) +
  #   standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
  #   annotate("text", x = 150, y = 2600, size = 4, color = "gray20", lineheight = .9,
  #   label = "Re-exploration:\ncomplete but old",family="EB Garamond") +
  #   geom_curve(aes(x = 145, y = 3300, xend = 92, yend = 7500),inherit.aes = F,arrow = arrow(length = unit(0.07, "inch")), size = 0.2,color = "gray50", curvature = 0.2) +
  # annotate("text", x = 130, y = 60, size = 4, color = "gray20", lineheight = .9,
  #            label = "Discovery:\nincomplete and old",family="EB Garamond") +
  #   geom_curve(aes(x = 140, y = 70, xend = 161, yend = 300),inherit.aes = F,arrow = arrow(length = unit(0.07, "inch")), size = 0.2,color = "gray50", curvature = 0.2) +
  #   NULL
  
  
  









#####
some2 %>% distinct(Accepted_Name,CellID,.keep_all = T) %>% filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution) %>% group_by(CellID) %>% nest() -> some2

some2 %>% ungroup() %>% rowwise() %>% mutate(n_tot = total(data)) %>% mutate(n_rbcL = total_rbcL(data)) -> SR
g@coords %>% as_tibble() %>% mutate(CellID=1:nrow(.)) %>% right_join(.,SR,by="CellID") -> SR
roads <- rnaturalearth::ne_countries(scale = 110,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs"))
counts <- rnaturalearth::ne_coastline(scale = 110,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs"))
xlimits = c(-12650052,-1774452)
ylimits = c(-6774199,8606756)
SR %>% mutate(Darwinian_Gap = 1 - (n_rbcL / n_tot)) -> SR
ggplot() + geom_sf(data=roads,colour=NA,fill=adjustcolor(MetBrewer::met.brewer("Hiroshige",direction=-1)[1],alpha.f = 0.2),size=0.5) + xlim(xlimits) + ylim(ylimits) + theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),legend.position = "right") +
  #geom_sf(data=counts,size=0.2) +
  geom_tile(data = SR ,aes(x=x,y=y,fill=log(n_tot),color=log(n_tot))) +
  #scale_fill_gradientn(colors=(MetBrewer::met.brewer("Hokusai2",direction=1)[-1]),name="species\nrichness (log)") +
  #scale_color_gradientn(colors=(MetBrewer::met.brewer("Hokusai2",direction=1)[-1]),name="species\nrichness (log)") +
  scale_fill_stepsn(colors=(MetBrewer::met.brewer("Hiroshige",n = 10,direction=-1)),n.breaks=10,name="Gap") +
  scale_color_stepsn(colors=(MetBrewer::met.brewer("Hiroshige",n = 10,direction=-1)),n.breaks=10,name="Gap") +
  labs(x="",y="",title="The Darwinian Gap",
       subtitle="Proportion of species without sequence data") +
  standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
  NULL





roads <- rnaturalearth::ne_countries(scale = 110,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
counts <- rnaturalearth::ne_coastline(scale = 110,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

SR %>% summarise(range(decimalLongitude))
xlimits = c(-115,-83)
ylimits = c(10.7,25.4)
ggplot() + geom_sf(data=roads,colour="grey95",fill="grey95",size=0.5) + xlim(xlimits) + ylim(ylimits) + theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA)) +
  geom_sf(data=counts,size=0.2) +
  geom_tile(data = SR ,aes(x=decimalLongitude,y=decimalLatitude,fill=log(n),color=log(n)),size=0.5) +
  scale_fill_gradientn(colors=(MetBrewer::met.brewer("Hokusai2"))) +
  scale_color_gradientn(colors=(MetBrewer::met.brewer("Hokusai2"))) +
  labs(x="",y="") +
  
  NULL

#### PROCESSED SCRIPT ####
gbif_pow_exactjoin(gbif_data = "Traqueos_NeoTropics_COR.csv",
                             checklist = "wcs_and_atoz_genera_species_infraspecies_2019_apr.txt",
                             taxon_status = c("Accepted","Orthographic","Synonym","Unplaced"),
                             output=here("output","Joined_base.Rdata"))
query_Checklist(bmm_checklist = "Lista_SPP.txt",output = here("output","JoinedVilla_base.Rdata"))
geographic_filter(data="JoinedVilla_base.Rdata",output = "Joined_filtered.Rdata",
                              pow_distributions="wcs_and_atoz_genera_species_infraspecies_dist_2019.txt",
                              perform_tests=c("seas","centroids","institutions", "equal", "gbif","capitals", "zeros"),
                              wgsrpd = "level3/level3.shp")

#### EDITING ####
load(here("output","Joined_POWdist.Rdata"))
some

checklist = "wcvp_names.txt"
taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced")
list_names <- data.table::fread(here("wcvp_2022",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ") %>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble() %>%  filter(infraspecific_rank=="",taxon_status %in% all_of(taxon_status_cats))
  

## Fuzzy join
# PERHAPS DO THE EXACT JOIN WITH THE BINOMIAL, AND IDENTIFY RECORDS WITH MORE THAN ONE HIT..... THEN DO JOIN WITH AUTHORS....

some %>% filter(is.na(Accepted_Name)) %>% distinct(scientificName,.keep_all = T) -> target 
target %>% distinct(genus.x) %>% pull(genus.x) -> genus_list
list_names %>% filter(genus %in% all_of(genus_list)) %>% unite(binomial,c(genus,species),sep="_") -> list_names_tofuzzy

target %>% slice(1:5) %>% fuzzyjoin::stringdist_left_join(.,list_names_tofuzzy,by="binomial",method="lv",distance_col="distance") -> aver
aver %>% select(ID:Accepted_Name.x,Accepted_Name.y,distance) %>% view







load("~/Documents/1.PROYECTOS/10.BMM/CF_Project/input_data/Historical_data.R")
reduct %>% #distinct(Resolved_ACCEPTED,CellID,year,.keep_all = T) %>% 
  count(year) %>% 
  filter(year >= 1921) %>% 
  mutate(n=n) -> original

MetBrewer::met.brewer("Hiroshige") %>% as.vector() -> colores
some4 %>% #distinct(AcceptedName,CellID,year,.keep_all = T) %>% 
  count(year) %>% 
  filter(year >= 1921 & year < 2010) %>% 
  mutate(n=n) %>% 
  ggplot() + 
  geom_bar(aes(x=year,y= n),stat="identity",fill=colores[1]) + 
  #geom_bar(data=original,aes(x=year,y= n), stat="identity",fill=colores[length(colores)/2+1]) +
  theme(panel.grid = element_blank(),panel.background = element_blank(),axis.line = element_line()) +
  geom_vline(xintercept = c(1976,1986))  + 
  NULL


some4 %>% count(year,species) %>% filter(year >= 1950 & year < 2011) %>% 
  mutate(n = ifelse(n <= median(n), n, median(n)))-> N

some4 %>% group_by(year,species) %>% filter(year >= 1950 & year < 2011) %>% nest() %>% left_join(.,N,by=c("year","species")) %>% mutate(data=map2(data,n,replace=F,sample_n)) %>% unnest(data) -> rarified

some4 %>% #distinct(AcceptedName,CellID,year,.keep_all = T) %>% 
  count(year) %>% 
  filter(year >= 1950 & year < 2011) %>% 
  mutate(n=n) %>% 
  ggplot() + 
  geom_vline(xintercept = seq(1950,2010,5),alpha=0.2) +
  geom_bar(aes(x=year,y= -n),stat="identity",fill=colores[1]) + 
  geom_bar(data= rarified  %>% count(year),aes(x=year,y= n), stat="identity",fill=colores[length(colores)]) +
  theme(panel.grid = element_blank(),panel.background = element_blank(),axis.line = element_line()) + scale_y_continuous(breaks = c(-15000,-10000,-5000,0,5000),labels=c(c(15000,10000,5000,0,5000))) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1950,2010,5))  +
  labs(x="Year of collection",y="Number of collection records",title="Dowsampling of yearly collection effort by species",subtitle="3,198 species") +
  geom_label(aes(x=1958, y=-10000, label="Complete dataset"), color=colores[1]) + 
  geom_label(aes(x=1958, y=5000, label="Downsampled dataset"), color=colores[length(colores)]) +
  ylim(c(-18000,5500)) +
  NULL




some4 %>% filter(year >= 1950 & year < 2011) %>% count(AcceptedName) %>% filter(n >= 50)
rarified %>% ungroup() %>% count(AcceptedName) %>% filter(n >= 50)





reduct %>% ggplot(aes(x=decimallongitude,y=decimallatitude)) + geom_point()
some2 %>% count(AcceptedName) %>% filter(n >= 100) %>% summarise(sum(n))
some3 %>% count(AcceptedName) %>% filter(n >= 100) %>% summarise(sum(n))

#### These are three cut-offs based on lat/lon: it takes out species (not only occurrences) that extend beyond these limits
spp <- unique(p_test$Resolved_ACCEPTED[which(p_test$decimalLongitude < -135 | p_test$decimalLongitude > -30)])
cat("Number of species beyound America:", spp %>% length(.) ,"\n") #### This is the number of species in Villaseñor's list
p_test_2 <- p_test[-which(p_test$Resolved_ACCEPTED %in% spp),]
length(unique(p_test$Resolved_ACCEPTED)) ######## This is the number of species in the database
length(unique(p_test_2$Resolved_ACCEPTED)) ######## This is the number of species without the ones extending past the limits
b$name_ACCEPT %>% .[which(!is.na(.))] %>% length(.) #### This is the number of species in Villaseñor's list





setwd("~/Documents/1.PROYECTOS/10.BMM/")
load("data/BMM_POW_clean1.Rdata")
p_test_clean
p_test_clean %>% distinct(Resolved_ACCEPTED,decimallatitude,decimallongitude) %>% group_by(Resolved_ACCEPTED) %>% summarise(Count=n()) %>% filter(Count >= 50) %>% select(Count) %>% sum()

### The species list was processed this way:
### 1) Vetted against Villaseñor
### 2) Corrected with POW
### 3) Records in America only
### 4) Species with less than 2.5% of records in 'more-arid' biomes
### 5) Species with more than 50 records

g <- raster::raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% as(., 'SpatialPixels')
r1 <- SpatialPoints(p_test_clean[,8:7],proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% 
  sp::over(.,g) %>% enframe(.,name="name") %>% mutate(.,name = p_test_clean$Resolved_ACCEPTED) %>% rename_with(.,~all_of(c("Species","CellID")))
p_test_clean <- p_test_clean %>% bind_cols(CellID=r1$CellID,.)
p_test_clean %>% distinct(Resolved_ACCEPTED,CellID,.keep_all = T) %>% 
  dplyr::select(Resolved_ACCEPTED )%>% table() -> tabla
which(tabla >= 50) %>% length()
p_test_clean %>% distinct(Resolved_ACCEPTED,CellID,.keep_all = T) %>% 
  dplyr::select(family_ACCEPTED_POW,class,order) -> tabla
tabla %>% distinct(family_ACCEPTED_POW) %>% dim()
tabla %>% distinct(order) %>% dim()
tabla %>% distinct(class) %>% dim()








## PULL UNMATCHED RECORDS TO ATTEMPT A FUZZY JOIN: SPLIT INTO GROUPS DUE TO MEMORY ISSUES
joined %>% filter(is.na(db_id)) %>%  # get unmatched records
  filter(!grepl("×",scientificName) ) %>%   # get rid of hybrids, vars, and subsp.
  filter(!grepl("var.",scientificName) ) %>% 
  filter(!grepl("subsp.",scientificName) ) %>% 
  select(1:12) -> to_fuzzy

to_fuzzy %>% group_by(scientificName) %>% group_split(.keep=T) -> split_unmatched
to_fuzzy %>%  distinct(genus) %>% pull(genus) -> list_genus
list_names %>% separate(binomial,into=c("genus","species"),sep="_",remove=F) %>% 
  filter(genus %in% all_of(list_genus)) %>% select(-genus,-species) -> Sub_list_names
split_fuzzymatched <- list()
n = 0
for(i in 1:length(split_unmatched)) { 
  ## REDUCE THE CHECKLIST TO THE TARGET GENERA: exact match!
  cat("Attempting fuzzy join----",length(split_unmatched) - n,"species left","\n")
  n = n + 1
 ## Fuzzy join using lv distance
  split_unmatched[[i]] %>% slice(1) %>% 
    fuzzyjoin::stringdist_left_join(.,Sub_list_names,by="scientificName",method="lv") -> aver
  aver %>% filter(!is.na(db_id)) %>% nrow() -> x
  if(x > 0) {
    cat("---- SUCCESS! -----","\n")
  split_fuzzymatched[[i]] <-  split_unmatched[[i]] %>% bind_cols(.,aver %>% select(13:25))
  split_unmatched[[i]] <- NULL
  }
  if(x == 0) cat("---- FAIL! -----","\n")
  
  }






