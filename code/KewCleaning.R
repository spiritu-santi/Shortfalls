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

######################################################################################################
#####         Last update: January 20, 2023
##### R script to perform all necessary processing for the geographic data as obtained
##### from GBIF and compared against WCVP. There is a pre-processing step in bash to reduce
##### file size, but row ID can be traced back to the original data as downloaded from GBIF with 
##### extra information.
##### The WCVP data can be directly downloaded from Kew's Plant of the World webpage.
######################################################################################################

###### PRE-PROCESSING: SETTING UP DATA FILES RETRIEVED FROM GBIF #####
###### USE BASH TO REDUCE FILE SIZE BY SELECTING COLUMNS
### bash commands are commented out!
# unzip 0188599-210914110416597.zip
# rename unzipped file to Traqueos_NeoTropics_raw.csv
# awk -F"\t" '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13"\t"$22"\t"$23"\t"$33}' Traqueos_NeoTropics_raw.csv > Traqueos_NeoTropics_COR.csv
# wc -l Traqueos_NeoTropics_raw.csv
# wc -l Traqueos_NeoTropics_COR.csv
# awk -F"\t" '{print $13"\t"$22"\t"$23"\t"$37}' Traqueos_NeoTropics_raw.csv > Traqueos_NeoTropics_CODES.csv

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

### DEPRECATED FUNCTIONS ####
#### This function is now integrated into estimate_metrics()
phylo_shortfall_deprecated <-function(data="Joined_finalPOWdist.v2.Rdata",checklist = "wcvp_names.txt",taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced")) { 
  
  load(here("output",data))
  some %<>% mutate(binomial_accepted=Accepted_Name) %>% mutate(binomial_accepted=sub("× ","",binomial_accepted)) %>%  mutate(binomial=sub(" ","_",binomial_accepted)) %>% separate(binomial,into=c("binomial"),extra = "warn",sep=" ") 
  some %<>% filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution)
  
  list.files("data/Trees_originals/") -> arboles
  cat("reading",arboles[2],"\n")
  ALLMB <- ape::read.tree(paste0("data/Trees_originals/",arboles[2]))
  ferns <- ape::extract.clade(ALLMB,node=ape::getMRCA(ALLMB,tip=ALLMB$tip.label[c(74461,74359)]))
  ALLMB <- ape::drop.tip(ALLMB, tip=ferns$tip.label)
  cat("Tree size (number of species):",length(ALLMB$tip.label),"\n")
  ## Matching the original names to the POW names 'Resolved_ACCEPTED'
  cat("Matching species without correction to tree","\n")
  qui <- tibble(binomial=ALLMB$tip.label,in_tree="SB") %>% left_join(.,list_names %>% filter(taxon_rank=="Species"),by="binomial") %>% filter(!is.na(binomial)) %>% select(binomial,in_tree,Accepted_Name)
  qui <- ftolr::accessions_long %>% select(species,sci_name) %>% mutate(in_tree = "NITTA") %>% distinct(sci_name,.keep_all = TRUE) %>% rename(binomial=species) %>% left_join(.,list_names %>% filter(taxon_rank=="Species"),by="binomial") %>% filter(!is.na(Accepted_Name)) %>% select(binomial,in_tree,Accepted_Name) %>% bind_rows(qui,.)
  
  some <- qui %>% distinct(Accepted_Name,.keep_all = TRUE) %>% right_join(.,some,by="Accepted_Name")
  some %>% select(Accepted_Name,in_tree,LEVEL2_COD) %>% #filter(LEVEL2_COD %in% c("79","80")) %>% 
    distinct(Accepted_Name,.keep_all = TRUE) %>% filter(!is.na(in_tree))
  
  
  # This are the species in the tree and in GBIF
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
  
  some %>% select(decimalLongitude,decimalLatitude) %>% sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS( "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) %>% sp::over(.,poly) %>% as_tibble() -> gridded
  gridded %<>% bind_cols(some %>% select(-LEVEL2_COD),.)
  gridded %>% select(Accepted_Name,in_tree,LEVEL2_COD,year) %>% 
    filter(!is.na(year)) %>% 
    filter(LEVEL2_COD %in% c("79","80")) #%>% 
  #distinct(Accepted_Name,.keep_all = TRUE) %>% 
  filter(!is.na(in_tree))
  
  
  
  
} 
#### This function is now integrated into map_hotspots()
localG_deprecated <- function(resolution = 15,distance = 75,var = "Completeness"){
  require(spdep)
  xlimits = c(-125,-75)
  ylimits = c(0,35)
  
  load("interim/Allmetrics_15.Rdata")
  base <- raster("data/base_rasters/baseRaster_15.tif")
  datos %<>% mutate(Completeness = 100 - Completeness, seqs= 1- seqs, phylo = 1 - phylo)
  datos <- datos[!is.na(datos[[glue::glue('{var}')]]),]
  datos <- datos[!is.na(datos[["x"]]),]
  nb <- dnearneigh(as.matrix(datos[,c("x","y")]),0,distance,longlat = TRUE)
  aver <- spdep::localG(datos[[glue::glue('{var}')]],nb2listw(nb,style = "B",zero.policy=TRUE))
  datos$localG <- as.vector(aver)
  prefix = paste0("raster",resolution,"_",var,".tif")
  rasterize(datos[,c("x","y")],base,field=datos[[glue::glue("{var}")]]) %>% crop(.,extent(hfp)) %>%   writeRaster(.,filename = here("output",prefix),overwrite=TRUE)
  prefix = paste0("raster",resolution,"_",var,"G",".tif")
  #datos$localG[which(is.na(datos$localG))] <- 999
  rasterize(datos[,c("x","y")],base,field=datos$localG) %>% crop(.,extent(hfp)) %>% writeRaster(.,filename = here("output",prefix),overwrite=TRUE)
  prefix = paste0("raster",resolution,"_","Slope",".tif")
  rasterize(datos[,c("x","y")],base,field=datos$Slope) %>% crop(.,extent(hfp)) %>%   writeRaster(.,filename = here("output",prefix),overwrite=TRUE)
  
}
### --- ####

###### WORKING FUNCTIONS: NCBI ####
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
} # CREATE THE ACCESSION TO TAXAID DATABASE (THIS PRODUCES A BIG FILE!!)


#### WORKING FUNCTIONS: FILTERING ####
gbif_pow_exactjoin <- function(gbif_data = "Traqueos_NeoTropics_COR.csv",checklist = "wcvp_names.txt",taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),output_file="Joined_base.Rdata") { 
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

gbif_pow_first_fuzzy <- function(data="Joined_base.Rdata",checklist = "wcvp_names.txt",taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),output_file="Joined_fuzzy.Rdata"){ 
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

geographic_filter <- function(data="Joined_final.v1.Rdata",output_file = "Joined_finalPOW.v1.Rdata",perform_tests=c("centroids","institutions", "equal", "gbif","capitals", "zeros","seas")) {
  load(here("output",data))
  ### APPLY FILTERS. NOT USING THE OUTLIER TEST, BECAUSE WE BASE THIS ON KEW'S DISTRIBUTIONS.
  cat("Filtering using CoordinateCleaner","\r")
  to_filter <- joined_fuzzy %>% select(decimalLongitude,decimalLatitude) %>% as.data.frame() 
  f <- CoordinateCleaner::clean_coordinates(to_filter, lon = "decimalLongitude",lat = "decimalLatitude",centroids_rad = 1000, centroids_detail = "both",tests=perform_tests,species=NULL,value="flagged",seas_scale = 110)
  joined_fuzzy %<>% mutate(Outlier_Test = f)
  save(joined_fuzzy,file=here("output",output_file))
}
  
powdist_filter <- function(data="Joined_finalPOW.v1.Rdata",output_file = "Joined_finalPOWdist.v2.Rdata",pow_distributions="wcvp_distribution.txt",wgsrpd = "level3/level3.shp") { 
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

#### WORKING FUNCTIONS: ANALYSES ####
do_format_completeness <- function(data="output/Joined_finalPOWdist.v2.Rdata",wgsrpd = "level2/level2.shp"){
  load(data)
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

#### Get institution codes
get_codes_GBIF <- function(gbif_codes = "Traqueos_NeoTropics_CODES.csv",data="Joined_finalPOWdist.v2.Rdata"){
  aa <- data.table::fread(here("data",gbif_data)) %>% as_tibble()
  load(here("output",data))
  some %>% filter(LEVEL2_COD %in% c(79,80),!is.na(Accepted_Name),Outlier_Test,POW_distribution) %>% pull(ID) -> idds
  aa %>% slice(idds) %>% #filter(institutionCode!="") %>% 
    group_by(institutionCode) %>% 
    summarise(n=n(),issue=first(issue)) %>%
    mutate(prop=n/sum(n)*100) -> codes
  codes %>% arrange(desc(prop)) %>% mutate(cum = cumsum(prop)) %>% mutate(Note = ifelse(grepl("INSTITUTION_MATCH_NONE",issue),"INSTITUTION_MATCH_NONE","NA")) %>% select(-issue) %>% write.table(.,file=here("MS_REV_23/tables","table_Codes.csv"),sep=",",row.names=FALSE,col.names=TRUE,quote = FALSE) 
  
  codes %>% arrange(desc(prop)) %>% mutate(cum = cumsum(prop)) %>% slice(1:20) %>% 
    ## This has been set after manual check on insitution codes.... have to get this automated!
    mutate(Country = c("US","MEX","MEX","CR","US","MEX","US","CR","US","US","US","US","US","MEX","MEX","US","MEX","US","US","US")) %>% group_by(Country) %>% summarize(sum(prop))
}

# Estimate metrics per grid-cell
estimate_metrics <- function(data="Joined_finalPOWdist.v2.Rdata",checklist = "wcvp_names.txt",taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),resolution = 15,data_comp="output/completeness_15secs/Estimators.RData",  projection = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"){ 
  
list_names <- data.table::fread(here("wcvp_2022",checklist),sep="|",quote="") %>% as_tibble()
list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ",remove = F) %>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble() %>%  filter(taxon_status %in% all_of(taxon_status_cats)) %>% unite("binomial",c(genus,species),sep="_",remove = F) %>% filter(accepted_plant_name_id!="")
  
load(here("output",data))
  some %<>% mutate(binomial_accepted=Accepted_Name) %>% mutate(binomial_accepted=sub("× ","",binomial_accepted)) %>%
      mutate(binomial=sub(" ","_",binomial_accepted)) %>% separate(binomial,into=c("binomial"),extra = "warn",sep=" ")
  some %<>% filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution)
some %<>% mutate(family.y = ifelse(is.na(family.y),family.x,family.y))
  list.files("data/Trees_originals/") -> arboles
  cat("reading",arboles[3],"\n")
  ALLMB <- ape::read.tree(paste0("data/Trees_originals/",arboles[3]))
  one <- grep("Equisetum",ALLMB$tip.label)[1]
  two <- grep("Cyathea",ALLMB$tip.label)[1]
  ferns <- ape::extract.clade(ALLMB,node=ape::getMRCA(ALLMB,tip=ALLMB$tip.label[c(one,two)]))
  ALLMB <- ape::drop.tip(ALLMB, tip=ferns$tip.label)
  cat("Tree size (number of species):",length(ALLMB$tip.label),"\n")
  ## Matching the original names to the POW names 'Resolved_ACCEPTED'
  cat("Matching species without correction to tree","\n")
  qui <- tibble(binomial=ALLMB$tip.label,in_tree="SB") %>% left_join(.,list_names %>% filter(taxon_rank=="Species"),by="binomial") %>% filter(!is.na(binomial)) %>% select(binomial,in_tree,Accepted_Name)
  qui <- ftolr::accessions_long %>% select(species,sci_name) %>% mutate(in_tree = "NITTA") %>% distinct(sci_name,.keep_all = TRUE) %>% rename(binomial=species) %>% left_join(.,list_names %>% filter(taxon_rank=="Species"),by="binomial") %>% filter(!is.na(Accepted_Name)) %>% select(binomial,in_tree,Accepted_Name) %>% bind_rows(qui,.)
  some <- qui %>% distinct(Accepted_Name,.keep_all = TRUE) %>% right_join(.,some,by="Accepted_Name")

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
  
  some %>% select(decimalLongitude,decimalLatitude) %>% sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,poly) %>% as_tibble() -> gridded
  some <- gridded %>% bind_cols(some %>% select(-LEVEL2_COD),.)

  res <- 60/resolution
  g <- raster::raster(nrows=180*res,ncols=360*res,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% projectRaster(.,crs = CRS(projection)) %>% as(., 'SpatialPixels')
  some %>% select(decimalLongitude,decimalLatitude) %>% 
    sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID"))) %>% select(CellID) %>% bind_cols(.,some) -> gridded
  
  g@coords %>% as_tibble() %>% mutate(CellID=1:nrow(.)) -> coords
  pointos <- coords %>% select(1:2) %>% as.data.frame(.) %>% SpatialPoints(.)
  proj4string(pointos) <- projection
  points_powo <- sp::over(pointos,poly) %>% as_tibble()
  coords %<>% bind_cols(.,points_powo)

  # PROCESS THE DATA TO DO ALL THE ANALYSES:
  #     - SPECIES RICHNESS
  gridded %>% select(Accepted_Name,CellID,year,Outlier_Test,POW_distribution) %>% 
    filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution) %>% 
    distinct(Accepted_Name,CellID) %>% count(CellID) %>% 
    left_join(.,gridded %>% distinct(CellID,.keep_all = T),by="CellID") %>% 
    select(CellID,n,LEVEL3_NAM) %>% rename(SR=n) -> SR
  SR %<>% mutate(SR_std = (SR - min(SR,na.rm=TRUE)) / (max(SR,na.rm=TRUE) - min(SR,na.rm=TRUE)))
  #     - NUMBER OF RECORDS
  gridded %>% filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution) %>% count(CellID) %>% 
    left_join(.,gridded %>% distinct(CellID,.keep_all = T),by="CellID") %>% select(CellID,n) %>% inner_join(.,SR,by="CellID") %>% mutate(n_std = (n - min(n,na.rm=TRUE)) / (max(n,na.rm=TRUE) - min(n,na.rm=TRUE))) -> SR
  nombre <-  paste0("RichnessRecords_",resolution,".Rdata")
  richness <- SR %>% left_join(.,coords,by="CellID")
  save(richness,file=here("interim",nombre))
  
#     - AGE OF COLLECTIONS
  gridded %>% select(CellID,Accepted_Name,year,Outlier_Test,POW_distribution) %>% filter(Outlier_Test,POW_distribution,year >= 1922) %>% group_by(CellID,Accepted_Name) %>% nest() -> test
test <- test %>% ungroup() %>% rowwise() %>% mutate(MostRecent = (get_mostrecent_year(data)))
 nombre <-  paste0("Ages_",resolution,".Rdata")
 families <- test %>% 
   left_join(.,gridded %>% select(Accepted_Name,family.y,CellID,year,Outlier_Test,POW_distribution) %>% 
               filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution) %>% 
               distinct(family.y,Accepted_Name),by="Accepted_Name") 
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
names %<>% filter(!grepl(" x ",species)) %>% filter(!grepl("]",species)) %>% filter(!grepl(")",species))  %>%
    #distinct(TaxaID,.keep_all = T) %>% 
    rename(binomial=species) %>% left_join(.,list_names,by = "binomial") %>% select(TaxaID,region,binomial,accepted_plant_name_id,Accepted_Name)
  names
  
 genebanks <-  gridded %>% select(CellID,year,Outlier_Test,POW_distribution,accepted_plant_name_id,in_tree) %>% 
    mutate(TaxaID=names$TaxaID[match(gridded$accepted_plant_name_id,names$accepted_plant_name_id)]) %>%
    mutate(ncbi=names$region[match(gridded$accepted_plant_name_id,names$accepted_plant_name_id)]) %>% 
    filter(Outlier_Test,POW_distribution)  %>% distinct(CellID,accepted_plant_name_id,.keep_all = T) %>% group_by(CellID) %>% 
    mutate(value=ifelse(is.na(ncbi),0,1)) %>% summarise(seqs = sum(value)/n())
 
 gridded %>% select(CellID,year,Outlier_Test,POW_distribution,accepted_plant_name_id,in_tree) %>% 
   mutate(TaxaID=names$TaxaID[match(gridded$accepted_plant_name_id,names$accepted_plant_name_id)]) %>%
   mutate(ncbi=names$region[match(gridded$accepted_plant_name_id,names$accepted_plant_name_id)]) %>% 
   filter(Outlier_Test,POW_distribution) 
 
 SR %<>% inner_join(.,genebanks,by="CellID")
 
phylos <- gridded %>% select(CellID,year,Outlier_Test,POW_distribution,accepted_plant_name_id,in_tree) %>% 
   filter(Outlier_Test,POW_distribution)  %>% distinct(CellID,accepted_plant_name_id,.keep_all = T) %>% group_by(CellID) %>% 
   mutate(value=ifelse(is.na(in_tree),0,1)) %>% summarise(phylo = sum(value)/n())

SR %<>% inner_join(.,phylos,by="CellID")


 # - COMPLETENESS
 load(data_comp)
 datos <- values %>% as_tibble()
 datos %>% select(Longitude,Latitude) %>% 
   sf::st_as_sf(x = ., coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID"))) %>% select(CellID) %>% bind_cols(.,datos) -> datos
 SR %<>% right_join(.,datos,by="CellID")
 datos <- SR
 nombre <-  paste0("Allmetrics_",resolution,".Rdata")
 save(datos,file=here("interim",nombre))
 
 
 
 
  }

plot_mappy_things <- function(data="interim/Allmetrics_15.Rdata",resolution=60,rich_bin=200,comp_bin=5,seq_bin=0.05,age_bin=5){
load(here(data))
roads <- rnaturalearth::ne_countries(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
 counts <- rnaturalearth::ne_coastline(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
 xlimits = c(-125,-75)
 ylimits = c(0,35)
 
 ## Species richness
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
   geom_tile(data = datos,aes(x=x,y=y,fill=SR),color="white") +
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
 if(resolution==30) {n_bins=27}
 if(resolution==60) {n_bins=38}
 n_tot <- datos %>% filter(!is.na(SR)) %>%  nrow(.)
hist_rich <- datos %>% ggplot(aes(x=SR)) + 
  geom_histogram(binwidth=rich_bin,fill=MetBrewer::met.brewer("Derain",direction=1,n=n_bins)) +
  stat_bin(binwidth=rich_bin,bins=n_bins, geom="text", aes(label=round(..count../n_tot*100,2)),fill="black",vjust=-1.5) +
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
  scale_x_continuous(breaks=seq(0,upperlim,1000)) +
  #scale_y_continuous(trans="log10") +
 NULL

richness <- richness + hist_rich

### Completeness: spatial shortfall
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
   geom_tile(data = datos,aes(x=x,y=y,fill=100-Completeness),color="white") +
   scale_fill_stepsn(colors=(MetBrewer::met.brewer("Hokusai3",direction=1)),name="Spatial shortfall (incompleteness)",
                     breaks=seq(0,upperlim,comp_bin),labels= scales::comma,
                     na.value="NA",limits=c(0,upperlim),
                     guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5,keyheight=0.4)) +
   labs(x="",y="") +
   #standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
   geom_sf(data=roads,colour="black",fill=NA,size=0.2) +
   NULL
 n_bins <- length(ggplot2:::bin_breaks_width(range(0,upperlim), width = comp_bin)$breaks) - 1L
 
 n_tot <- datos %>% filter(!is.na(Completeness)) %>%  nrow(.)
 hist_comp <- datos %>% ggplot(aes(x=100-Completeness)) + 
   geom_histogram(fill=MetBrewer::met.brewer("Hokusai3",direction=1,n=n_bins),binwidth=comp_bin,color="white") +
   stat_bin(binwidth=comp_bin,geom="text", aes(label=round(..count../n_tot*100,1)),fill="black",vjust=-1.5) +
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
   labs(x="Spatial shortfall (incompleteness)",y="Number of cells") +
   scale_x_continuous(breaks=seq(0,upperlim,10),limits = c(0,upperlim)) +
   NULL
 completeness <- completeness + hist_comp

 # Occurrence age: temporal shortfall
 upperlim <-  roundUp(max(datos$age,na.rm=T),to=100)
 n_bins <- length(ggplot2:::bin_breaks_width(range(0,upperlim), width = age_bin)$breaks) - 1L
 n_bins = 11
 colas <- rev(MetBrewer::met.brewer("OKeeffe2",direction=-1,n=n_bins))
 if(resolution==30) { colas <- (c(rep(colas[1],2),colas,rep(colas[length(colas)],7)))}
 if(resolution==60) { colas <- (c(rep(colas[1],1),colas))} else colas <- (c(rep(colas[1],2),colas,rep(colas[length(colas)],8)))
 age <-  ggplot() +   
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
  geom_tile(data = datos,aes(x=x,y=y,fill=age),color="white") +
#scale_fill_stepsn(colors=(MetBrewer::met.brewer("OKeeffe2",direction=1)),name="Age",breaks=seq(0,upperlim,age_bin),labels= scales::comma,limits=c(0,upperlim),na.value="NA",guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5)) +
  scale_fill_stepsn(colors=(colas),
                    name="Temporal shortfall (occurrence age)",
                    breaks = seq(0,upperlim,seq_bin),
                    #c(0.2,0.5,seq(0.55,upperlim,seq_bin)),
                    labels= scales::comma,limits=c(0.2,upperlim),
                    guide=guide_legend(reverse = TRUE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5)) +
   labs(x="",y="") +
   #standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
   geom_sf(data=roads,colour="black",fill=NA,size=0.2) +
   NULL
 
 n_tot <- datos %>% filter(!is.na(age)) %>%  nrow(.)
 hist_age <- datos %>% ggplot(aes(x=age)) + 
   geom_histogram(fill=colas,binwidth=age_bin,color="white") +
   #geom_histogram(fill=MetBrewer::met.brewer("OKeeffe2",direction=1,n=n_bins),binwidth=comp_bin,color="white") +
   stat_bin(binwidth=comp_bin,geom="text", aes(label=round(..count../n_tot*100,1)),fill="black",vjust=-1.5) +
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
   labs(x="Temporal shortfall (occurrence age)",y="Number of cells") +
   scale_x_continuous(breaks=seq(0,upperlim,10)) +
   NULL
 
  age <- age +hist_age
  
 # Proportion sequenced: phylogenetic shortfall (first component)
 upperlim <-  roundUp(max(datos$seqs,na.rm=T),to=1)
 #n_bins <- length(ggplot2:::bin_breaks_width(range(0,upperlim), width = seq_bin)$breaks) - 1L
 n_bins = 12
 colas <- rev(MetBrewer::met.brewer("Tam",direction=-1,n=n_bins))
 colas <- (c(rep(colas[1],7),colas,rep(colas[length(colas)],2)))
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
   geom_tile(data = datos,aes(x=x,y=y,fill=1-seqs),color="white") +
   scale_fill_stepsn(colors=(colas),
                     name="Phylogenetic shortfall (proportion not sequenced)",
                     breaks = seq(0,upperlim,seq_bin),
                     #c(0.2,0.5,seq(0.55,upperlim,seq_bin)),
                     labels= scales::comma,limits=c(0.2,upperlim),
                     guide=guide_legend(reverse = TRUE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5)) +
   labs(x="",y="") +
   #standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
   geom_sf(data=roads,colour="black",fill=NA,size=0.2) +
   NULL

 n_tot <- datos %>% filter(!is.na(seqs)) %>%  nrow(.)
 hist_sequenced <- datos %>% ggplot(aes(x=1-seqs)) + 
   geom_histogram(fill=colas,binwidth=seq_bin,color="white") +
   stat_bin(binwidth=seq_bin,geom="text", aes(label=round(..count../n_tot*100,1)),fill="black",vjust=-1.5) +
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
   labs(x="Phylogenetic shortfall (proportion not sequenced)",y="Number of cells") +
   scale_x_continuous(breaks=seq(0.2,upperlim,seq_bin*2),limits = c(0,upperlim)) +
   NULL

 sequenced <- sequenced +hist_sequenced

# Proporion in phylogeny: phylogenetic shortfall (second component)
 upperlim <-  roundUp(max(datos$phylo,na.rm=T),to=1)
 #n_bins <- length(ggplot2:::bin_breaks_width(range(0,upperlim), width = seq_bin)$breaks) - 1L
 n_bins = 13
 colas <- rev(MetBrewer::met.brewer("Tam",direction=-1,n=n_bins))
 if(resolution==60) { colas <- (c(rep(colas[1],1),colas))} else colas <- c(rep(colas[1],3),colas,rep(colas[length(colas)],5))
 phylo <- ggplot() + 
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
   geom_tile(data = datos,aes(x=x,y=y,fill=1-phylo),color="white") +
   scale_fill_stepsn(colors=(colas),
                     name="Phylogenetic shortfall (proportion not in phylogeny)",
                     breaks = seq(0,upperlim,seq_bin),
                     #c(0.2,0.5,seq(0.55,upperlim,seq_bin)),
                     labels= scales::comma,limits=c(0,upperlim),
                     guide=guide_legend(reverse = TRUE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0.5)) +
   labs(x="",y="") +
   #standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
   geom_sf(data=roads,colour="black",fill=NA,size=0.2) +
   NULL
 
n_tot <- datos %>% filter(!is.na(phylo)) %>%  nrow(.)
hist_phylo <-  datos %>% ggplot(aes(x=1-phylo)) + 
   geom_histogram(fill=colas,binwidth=seq_bin,color="white") +
   stat_bin(binwidth=seq_bin,geom="text", aes(label=round(..count../n_tot*100,1)),fill="black",vjust=-1.5) +
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
   labs(x="Phylogenetic shortfall (proportion not in phylogeny)",y="Number of cells") +
   scale_x_continuous(breaks=seq(0,upperlim,seq_bin*2)) +
   NULL
 
 phylo <- phylo + hist_phylo

pp <- wrap_plots(list(richness,completeness,age,phylo)) +
   plot_layout(byrow=T,nrow = 4,guides="keep",tag_level = 'new') & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"))
resolution
nombre <- paste0("Maps_",resolution,".pdf")
ggsave(filename = here("figures",nombre),plot = pp)

}

plot_phylo_things <- function(phylo = "data/Carta_tofamily.tre",age="interim/Ages_15.Rdata",data="interim/Allmetrics_15.Rdata",resolution=15){ 
  load(here(data))
  load(here(age))
  datos %>% pull(CellID) ->  cells
  families %<>% select(-data) %>% filter(!is.na(LEVEL2_COD)) %>%  group_by(Accepted_Name) %>% summarise(Youngest=max(MostRecent),family=first(family.y)) %>% group_by(family) #%>% summarise(Youngest_mean_age = median(Youngest,na.rm=T)) 
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
load(here("interim/Ages_15.Rdata"))
#Some manual edits to family names
families %>% filter(Accepted_Name=="Stylotrichium sucrei R.M.King & H.Rob.")

families %<>% ungroup() %>% mutate(family = ifelse(family.y == "Turneraceae","Passifloraceae",family.y )) %>% 
  mutate(family = ifelse(family == "Cochlospermaceae","Bixaceae",family )) %>% 
  mutate(family = ifelse(family == "Vivianiaceae","Francoaceae",family )) %>% 
  mutate(family = ifelse(family == "Stixaceae","Resedaceae",family )) %>% 
  mutate(family = ifelse(family == "Quiinaceae","Ochnaceae",family )) %>% 
  mutate(family = ifelse(family == "Hypseocharitaceae","Geraniaceae",family )) %>% 
  mutate(family = ifelse(family %in% c("Malesherbiaceae","Parnassiaceae"),"Celastraceae",family )) %>% 
  mutate(family = ifelse(family == "Tetradiclidaceae","Nitrariaceae",family )) %>% 
  mutate(family = ifelse(family == "Hydnoraceae","Hydnoraceae",family )) 
taxonomy %<>% as_tibble() %>% distinct(family,.keep_all = TRUE)

a <- match(families$family,taxonomy$family)
families %<>% ungroup() %>% mutate(class = taxonomy$class[a],.after=1) 
# Additional manual edits
families %<>% mutate(class = ifelse(family.y %in% c("Viburnaceae","Peltantheraceae","Hydnoraceae",
                                                   "Ruppiaceae","Apodanthaceae","Balanophoraceae",
                                                   "Mitrastemonaceae","Cytinaceae"),"Magnoliopsida",class))

families %<>% 
  mutate(class = ifelse(class %in% c("Cycadopsida","Gnetopsida") ,"Pinopsida",class ))
families %<>% separate(Accepted_Name,into=c("Genus","species","author"),sep=" ",remove=FALSE) 
families %<>% 
  mutate(class = ifelse(class =="Pinopsida" ,"Gymnosperms",class))%>% 
  mutate(class = ifelse(class =="Lycopodiopsida","Lycopods",class))%>% 
  mutate(class = ifelse(class =="Magnoliopsida" ,"Angiosperms",class))%>% 
  mutate(class = ifelse(class =="Polypodiopsida" ,"Ferns",class))

families %>% filter(!is.na(LEVEL2_NAM)) %>%
  rowwise() %>% mutate(Records = nrow(data)) %>% select(-data) %>% 
  group_by(Accepted_Name) %>% summarise(Class= first(class),Family=first(family),Genus=first(Genus),Youngest=max(MostRecent),Records=sum(Records)) %>% group_by(Class) %>% mutate(Record_age = 2022 - Youngest) %>%
  mutate(Class=factor(Class,levels=c("Lycopods","Ferns","Gymnosperms","Angiosperms"))) %>% 
  write.csv(.,file=here("MS_REV_23/Supplementary_data_Species_Ages.csv"))

families %>% filter(!is.na(LEVEL2_NAM)) %>% 
  rowwise() %>% mutate(Records = nrow(data)) %>%  select(-data) %>% 
  group_by(Genus) %>% summarise(Class= first(class),Family=first(family),Genus=first(Genus),Youngest=max(MostRecent),Richness=n(),Records=sum(Records)) %>% mutate(Record_age = 2022 - Youngest) %>%
  mutate(Class=factor(Class,levels=c("Lycopods","Ferns","Gymnosperms","Angiosperms"))) %>% 
  write.csv(.,file=here("MS_REV_23/Supplementary_data_Genus_Ages.csv"))
  
families %>% filter(!is.na(LEVEL2_NAM)) %>% filter(!is.na(family)) %>% 
  rowwise() %>% select(-data) %>% 
  group_by(family.y) %>% summarise(Youngest=max(MostRecent),class= first(class)) %>% group_by(class) %>% mutate(Record_age = 2022 - Youngest + 0.01) %>% mutate(class=factor(class,levels=c("Lycopods","Ferns","Gymnosperms","Angiosperms"))) %>% 
  ggplot(aes(y=class,x=Record_age,fill=class)) +
  geom_vline(xintercept = c(1,5,10,20,30,40,50,60,70,80,90,100),linetype="dashed",alpha=0.5) +
  geom_point(aes(color=class),position = position_jitter(height = 0.1),alpha=.8,size=0.7,shape=19) +
  scale_x_continuous(trans="log10",limits=c(1,100)) +
  ggdist::stat_slabinterval(color="black",height = 0.8,alpha=0.8,position = position_nudge(y = 0.11),
                            .width = c(0.67, 0.89),normalize="panels",slab_type="pdf",adjust=0.5,
                            interval_size_range = c(1, 2),breaks=20) + 
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

ggsave(filename = here("figures/AgeFamily_class.pdf"),last_plot())


}

plot_main_map <- function(data="Joined_finalPOWdist.v2.Rdata", resolution=15){ 
  
 
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

bi_plot <- function(resolution = 15){ 
  require(biscale)
  load(here("interim/Allmetrics_15.Rdata"))
  roads <- rnaturalearth::ne_countries(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  counts <- rnaturalearth::ne_coastline(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  xlimits = c(-125,-75)
  ylimits = c(0,35)
  datos %>% filter(is.na(LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -106,"Mexico Northwest",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -98,"Mexico Southwest",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 25,"Mexico Northeast",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 16 & Longitude < -93,"Mexico Gulf",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude > -83,"Panama",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 11,"Costa Rica",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 14 & Longitude < -92,"Mexico Southeast",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 14 & Longitude < -88,"El Salvador",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 13,"Nicaragua",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 17 & Longitude > -88,"El Salvador",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 17 & Longitude < -90,"Mexico Southeast",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 18,"Belize",LEVEL3_NAM))
  datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 18 & Longitude > -90,"Mexico Southeast",LEVEL3_NAM))
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
  
  esto %<>% mutate(Completeness = 100 - Completeness, phylo=1-phylo, seqs=1-seqs)
  
colors <- bi_class(esto %>% filter(!is.na(Completeness),!is.na(SR)), x = Completeness, y = SR, style = "fisher", dim = 4)
  
colors %>% count(bi_class)
legend <- bi_legend(pal = "PurpleGrn",dim = 4,
                      xlab = "Spatial shortfall",ylab = "Documented species richness",
                      size = 10, rotate_pal = FALSE,flip_axes = FALSE, breaks = bi_class_breaks(esto %>% filter(!is.na(Completeness),!is.na(SR)), x = Completeness, y = SR, style = "fisher", dim = 4))
  comp_bi <- colors %>% 
    ggplot(aes(x=x,y=y,fill=bi_class)) +
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
  
comp_bi
  
colors <- bi_class(esto %>% filter(!is.na(SR)), x = age, y = SR, style = "jenks", dim = 4)
  colors %>% count(bi_class)
  legend <- bi_legend(pal = "PurpleGrn",
                      dim = 4,
                      xlab = "Temporal shortfall",
                      ylab = "Documented species richness",
                      size = 10,rotate_pal = FALSE,flip_axes = FALSE,
                      breaks = bi_class_breaks(esto %>% filter(!is.na(SR)), x = age, y = SR, style = "jenks", dim = 4))
  
  age_bi <- colors %>% 
    ggplot(aes(x=x,y=y,fill=bi_class)) +
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
  
  age_bi
  
  colors <- bi_class(esto %>% filter(!is.na(SR)), x = phylo, y = SR, style = "jenks", dim = 4)
  colors %>% count(bi_class)
  legend <- bi_legend(pal = "PurpleGrn",
                      dim = 4,
                      ylab = "Documented species richness",
                      xlab = "Phylogenetic shortfall",
                      size = 10,rotate_pal = FALSE,flip_axes = FALSE,
                      breaks = bi_class_breaks(esto %>% filter(!is.na(SR)) , x = phylo, y = SR, style = "jenks", dim = 4))
  
  phylo_bi <- colors %>% 
    ggplot(aes(x=x,y=y,fill=bi_class)) +
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
  phylo_bi
  
  colors <- bi_class(esto %>% filter(!is.na(SR)), x = seqs, y = SR, style = "jenks", dim = 4)
  colors %>% count(bi_class)
  legend <- bi_legend(pal = "PurpleGrn",
                      dim = 4,
                      ylab = "Documented species richness",
                      xlab = "Phylogenetic shortfall (NCBI)",
                      size = 10,rotate_pal = FALSE,flip_axes = FALSE,
                      breaks = bi_class_breaks(esto %>% filter(!is.na(SR)) , x = seqs, y = SR, style = "jenks", dim = 4))
  
  seqs_bi <- colors %>% 
    ggplot(aes(x=x,y=y,fill=bi_class)) +
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
  
  seqs_bi
  
  layout <- "
AB
CD"
  pp <- wrap_plots(design = layout,list(comp_bi,age_bi,phylo_bi,seqs_bi)) +
    plot_layout(byrow=T,nrow = 2,guides="keep",tag_level = 'new') & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"))
  
  
  nombre <- paste0("Biplots_",resolution,".pdf")
  ggsave(filename = here("figures",nombre),plot = pp)
  
  
  
}

plot_ridges_by_region <- function(data="Allmetrics_15.Rdata",resolution = 15){ 
load(here("interim",data))
if(resolution==30) {datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -106,"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -97,"Mexico Southwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -94,"Mexico Gulf",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 17.5,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -90,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude > -84,"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude > -85.5,"Honduras",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude > -88,"Nicaragua",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 15,"Belize",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 15,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))}
if(resolution==60) {datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -108,"Mexico Northwest",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -97,"Mexico Southwest",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 20,"Mexico Southeast",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude > -83,"Panama",LEVEL3_NAM))}
if(resolution==15) {datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -106,"Mexico Northwest",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -98,"Mexico Southwest",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 25,"Mexico Northeast",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 16 & Longitude < -93,"Mexico Gulf",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude > -83,"Panama",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 11,"Costa Rica",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 14 & Longitude < -92,"Mexico Southeast",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 14 & Longitude < -88,"El Salvador",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 13,"Nicaragua",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 17 & Longitude > -88,"El Salvador",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 17 & Longitude < -90,"Mexico Southeast",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 18,"Belize",LEVEL3_NAM))
    datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 18 & Longitude > -90,"Mexico Southeast",LEVEL3_NAM))}
cat("Just checking....")
print(datos %>% filter(is.na(LEVEL3_NAM)))
datos %<>% mutate(Completeness = 100 - Completeness, seqs = 1 - seqs, phylo = 1 - phylo) 
#datos %>% write_csv(.,file=here("figures/temp.csv"))
aqui <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = (c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  group_by(LEVEL3_NAM) %>% 
  summarize(across(c(Completeness,age,seqs,phylo),
                   list(Median=function(x)median(x,na.rm=T),
                        Q0.1=function(x)quantile(x,na.rm=T,probs=0.1),
                        Q0.9=function(x)quantile(x,na.rm=T,probs=0.9)),.names = "{.fn}_{.col}")) %>% mutate(Resolution=resolution) #%>% pivot_longer(cols = 1:ncol(.),names_to = "Variable",values_to = "Values")
save(aqui,file = here("interim",paste0("summary_",resolution,".Rdata")))

age <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=age,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  #ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=3,size=0.2,geom="density_ridges_gradient",rel_min_height=0.01) +
  ggridges::geom_density_ridges(quantile_lines = TRUE, quantiles = 2,binwidth=2.5,
      size=0.2,rel_min_height=0.0,stat="binline",scale=0.99) +
  ggstance::stat_summaryh(fun.x=median, geom="text", aes(label=sprintf("%1.1f", ..x..)),
                          position=position_nudge(y=0.5,x=20), colour="black", size=3.5) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,
                                binwidth=0.025,size=0.2,scale=0.99,rel_min_height=0.0,fill=NA) +
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
  labs(x="Temporal shortfall",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,100)) +
  NULL

completeness <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=Completeness,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  #ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=5,size=0.2,geom="density_ridges_gradient",rel_min_height=0.001) +
  ggridges::geom_density_ridges(quantile_lines = TRUE, quantiles = 2,
                                binwidth=2.5,size=0.2,scale=0.99,rel_min_height=0.0,stat="binline") +
  ggstance::stat_summaryh(fun.x=median, geom="text", aes(label=sprintf("%1.1f", ..x..)),
                          position=position_nudge(y=0.5,x=30), colour="black", size=3.5) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,
                                binwidth=0.025,size=0.2,scale=0.99,rel_min_height=0.0,fill=NA) +
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
  labs(x="Spatial shortfall",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,100)) +
  NULL

phylo <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=phylo,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  #ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=0.02,size=0.2,geom="density_ridges_gradient",rel_min_height=0.001) +
  ggridges::geom_density_ridges(quantile_lines = TRUE, quantiles = 2,
                                binwidth=0.025,size=0.2,scale=0.99,rel_min_height=0.0,stat="binline") +
  ggstance::stat_summaryh(fun.x=median, geom="text", aes(label=sprintf("%1.2f", ..x..)),
                          position=position_nudge(y=0.5,x=0.2), colour="black", size=3.5) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,
                                binwidth=0.025,size=0.2,scale=0.99,rel_min_height=0.0,fill=NA) +
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
  labs(x="Phylogenetic shortfall",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,1)) +
  NULL

seqs <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
  ggplot(aes(x=seqs,y=LEVEL3_NAM,fill=LEVEL3_NAM)) +
  #ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,bandwidth=0.02,size=0.2,geom="density_ridges_gradient",rel_min_height=0.001) +
  ggridges::geom_density_ridges(quantile_lines = TRUE, quantiles = 2,
                                binwidth=0.025,size=0.2,scale=0.99,rel_min_height=0.0,stat="binline") +
  scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
  ggstance::stat_summaryh(fun.x=median, geom="text", aes(label=sprintf("%1.2f", ..x..)),
                position=position_nudge(y=0.5,x=0.1), colour="red", size=3.5) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2,
                                binwidth=0.025,size=0.2,scale=0.99,rel_min_height=0.0,fill=NA) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(),legend.position = "none",legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(family="EB Garamond",size=10),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.y =  element_blank(),legend.margin=margin(t=-25),
        legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  labs(x="Phylogenetic shortfall",y="") +
  geom_hline(yintercept = 1:13,size=0.05)  +
  coord_cartesian(xlim=c(0,1)) +
  NULL



xlimits = c(-123,-75)
ylimits = c(0,35)

# map <- datos %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels = rev(c("Mexico Northwest","Mexico Northeast","Mexico Southwest","Mexico Central","Mexico Gulf","Mexico Southeast","Guatemala","Belize","El Salvador","Honduras","Nicaragua","Costa Rica","Panama")))) %>% 
#   ggplot(aes(x=x,y=y,color=LEVEL3_NAM,fill=LEVEL3_NAM)) +
#   geom_tile() + 
#   scale_fill_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
#   scale_color_manual(values=(MetBrewer::met.brewer("Hiroshige",direction=1,n=13))) +
#   theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
#         legend.position = "none",
#         legend.background = element_rect(fill = NA),
#         #legend.key.width = unit(0.01, "cm"),
#         #legend.key.height = unit(0.5,"cm"),
#         legend.spacing.y = unit(c(0.2),"cm"),
  #       legend.title = element_text(family="EB Garamond"),
  #       legend.text = element_text(family="EB Garamond"),
  #       legend.margin=margin(t=-25),
  #       legend.key.size = unit(0.1,"cm"),
  #       #legend.spacing = unit(0.001,"cm"),
  #       plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal")  +
  # labs(x="",y="") +
  # # standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
  # geom_sf(data=counts,colour="black",fill=NA,size=0.2,inherit.aes = FALSE) +
  # xlim(xlimits) + ylim(ylimits) +
  # NULL

phylo_comp <- datos %>% #mutate(Completeness = replace_na(Completeness,0)) %>% 
  group_by(LEVEL3_NAM) %>% #summarise(across(c(2,23,16),median,na.rm=T)) %>% 
  ggplot(aes(x=Completeness,y=phylo,fill=LEVEL3_NAM,color=LEVEL3_NAM)) +
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
  labs(x="Spatial shortfall",y="Phylogenetic shortfall") +
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
  labs(x="Spatial shortfall",y="Temporal shortfall") + 
  #coord_flip() +
  NULL


layout <- "
B
C"
pp <- completeness | age | phylo | (age_comp/phylo_comp + plot_layout(design=layout)) +
  plot_layout(byrow=T,nrow = 2,guides="keep",tag_level = 'new') & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"))
nombre <- paste0("Ridges_",resolution,".pdf")
ggsave(filename = here("figures",nombre),plot = pp)
pp
sink(here("interim",paste0("linearfits_",resolution,".txt")))
lm(phylo ~ 1 + Completeness, data = datos)  %>% summary()
lm(age ~ 1 + Completeness, data = datos)  %>% summary()
lm(phylo ~ 1 + age, data = datos)  %>% summary()
sink()
}

median_plots <- function(){
load("interim/summary_15.Rdata")
tabla <- aqui
load("interim/summary_30.Rdata")
tabla <- aqui %>% bind_rows(tabla,.)
load("interim/summary_60.Rdata")
tabla <- aqui %>% bind_rows(tabla,.)
tabla_glob <- tabla
tabla %>% relocate(Resolution,.after = 1) %>% mutate(across(3:ncol(.),round,2)) %>% 
  mutate(Completeness = paste0(Median_Completeness,"\n\U005B",Q0.1_Completeness,"\U2014",Q0.9_Completeness,"\U005D")) %>% 
  mutate(Age = paste0(Median_age,"\n\U005B",Q0.1_age,"\U2014",Q0.9_age,"\U005D")) %>% 
  mutate(Sequenced = paste0(Median_seqs,"\n\U005B",Q0.1_seqs,"\U2014",Q0.9_seqs,"\U005D")) %>%
  mutate(Phylogenetic = paste0(Median_phylo,"\n\U005B",Q0.1_phylo,"\U2014",Q0.9_phylo,"\U005D")) %>%
  select(LEVEL3_NAM,Resolution,Completeness,Age,Phylogenetic,Sequenced) %>% 
  flextable::flextable() %>% flextable::save_as_docx(.,path=here("interim/table_medians_region.docx"))

tabla <- tabla_glob %>% select(LEVEL3_NAM,Median_Completeness,Q0.1_Completeness,Q0.9_Completeness,Resolution) %>% pivot_longer(col=2:4,values_to = "val",names_to = "measure") %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels=rev(levels(LEVEL3_NAM)))) %>% group_by(LEVEL3_NAM)

comp <- ggplot() +
  geom_segment(data = tabla %>% filter(measure=="Median_Completeness")%>% group_by(LEVEL3_NAM) %>% slice_min(order_by = val) %>% ungroup() ,aes(x = min(val), xend = val, y = LEVEL3_NAM, yend = LEVEL3_NAM),linetype = "dotted", size = 0.5, color = "gray80") +
  #geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==15),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.3),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==30),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.4),alpha=0.5) + 
 # geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.4),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.4),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==60),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.5),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.5),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.5),alpha=0.7) +
  geom_segment(data = tabla %>% filter(measure=="Median_Completeness") %>% group_by(LEVEL3_NAM) %>% summarise(start = range(val)[1], end = range(val)[2]) %>% ungroup(),aes(x = start, xend = end, y = LEVEL3_NAM, yend = LEVEL3_NAM),color = "gray80", size = 1)  +
  geom_point(data = tabla %>% filter(measure=="Median_Completeness"),  aes(val, LEVEL3_NAM, fill = factor(Resolution,levels = c(15,30,60))),size = 3,shape=21) +
  geom_point(data = tabla %>% filter(measure=="Median_Completeness"),  aes(val, LEVEL3_NAM, fill = factor(Resolution,levels = c(15,30,60))),size = 3,shape=21) +
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
  labs(x="Median spatial shortfall",y="",title="") +
  NULL

tabla <- tabla_glob %>% select(LEVEL3_NAM,Median_age,Q0.1_age,Q0.9_age,Resolution) %>% pivot_longer(col=2:4,values_to = "val",names_to = "measure") %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels=rev(levels(LEVEL3_NAM)))) %>% group_by(LEVEL3_NAM)
age <- ggplot() +
  geom_segment(data = tabla %>% filter(measure=="Median_age")%>% group_by(LEVEL3_NAM) %>% slice_min(order_by = val) %>% ungroup() ,aes(x = min(val), xend = val, y = LEVEL3_NAM, yend = LEVEL3_NAM),linetype = "dotted", size = 0.5, color = "gray80") +
  #geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==15),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.3),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==30),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.4),alpha=0.5) + 
  # geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.4),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.4),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==60),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.5),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.5),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.5),alpha=0.7) +
  geom_segment(data = tabla %>% filter(measure=="Median_age") %>% group_by(LEVEL3_NAM) %>% summarise(start = range(val)[1], end = range(val)[2]) %>% ungroup(),aes(x = start, xend = end, y = LEVEL3_NAM, yend = LEVEL3_NAM),color = "gray80", size = 1)  +
  geom_point(data = tabla %>% filter(measure=="Median_age"),  aes(val, LEVEL3_NAM, fill = factor(Resolution,levels = c(15,30,60))),size = 3,shape=21) +
  geom_point(data = tabla %>% filter(measure=="Median_age"),  aes(val, LEVEL3_NAM, fill = factor(Resolution,levels = c(15,30,60))),size = 3,shape=21) +
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
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(family="EB Garamond",size=10),
        axis.title = element_text(family="EB Garamond",size=12),
        plot.title = element_text(family="EB Garamond",face="bold",size=20,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) + 
  # guides(color = guide_colorbar(title.position="top",title.hjust = 0, ticks = F),alpha="none",size = "none"
  #guide_legend(reverse=T,ticks=F,override.aes = list(alpha=1,fill="grey50"))) +
  labs(x="Median temporal shortfall",y="",title="") +
  NULL

tabla <- tabla_glob %>% select(LEVEL3_NAM,Median_seqs,Q0.1_seqs,Q0.9_seqs,Resolution) %>% pivot_longer(col=2:4,values_to = "val",names_to = "measure") %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels=rev(levels(LEVEL3_NAM)))) %>% group_by(LEVEL3_NAM)

seqs <- ggplot() +
  geom_segment(data = tabla %>% filter(measure=="Median_seqs")%>% group_by(LEVEL3_NAM) %>% slice_min(order_by = val) %>% ungroup() ,aes(x = min(val), xend = val, y = LEVEL3_NAM, yend = LEVEL3_NAM),linetype = "dotted", size = 0.5, color = "gray80") +
  #geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==15),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.3),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==30),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.4),alpha=0.5) + 
  # geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.4),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.4),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==60),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.5),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.5),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.5),alpha=0.7) +
  geom_segment(data = tabla %>% filter(measure=="Median_seqs") %>% group_by(LEVEL3_NAM) %>% summarise(start = range(val)[1], end = range(val)[2]) %>% ungroup(),aes(x = start, xend = end, y = LEVEL3_NAM, yend = LEVEL3_NAM),color = "gray80", size = 1)  +
  geom_point(data = tabla %>% filter(measure=="Median_seqs"),  aes(val, LEVEL3_NAM, fill = factor(Resolution,levels = c(15,30,60))),size = 3,shape=21) +
  geom_point(data = tabla %>% filter(measure=="Median_seqs"),  aes(val, LEVEL3_NAM, fill = factor(Resolution,levels = c(15,30,60))),size = 3,shape=21) +
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
  labs(x="Median phylogenetic shortfall (NCBI)",y="",title="") +
  NULL

tabla <- tabla_glob %>% select(LEVEL3_NAM,Median_phylo,Q0.1_phylo,Q0.9_phylo,Resolution) %>% pivot_longer(col=2:4,values_to = "val",names_to = "measure") %>% mutate(LEVEL3_NAM=factor(LEVEL3_NAM,levels=rev(levels(LEVEL3_NAM)))) %>% group_by(LEVEL3_NAM)


phylo <- ggplot() +
  geom_segment(data = tabla %>% filter(measure=="Median_phylo")%>% group_by(LEVEL3_NAM) %>% 
  slice_min(order_by = val) %>% ungroup() ,aes(x = min(val), xend = val, y = LEVEL3_NAM, yend = LEVEL3_NAM),linetype = "dotted", size = 0.5, color = "gray80") +
  #geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==15),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.3),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==15),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.3),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==30),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.4),alpha=0.5) + 
  # geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.4),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==30),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.4),alpha=0.7) +geom_segment(data =tabla %>% pivot_wider(names_from = measure,values_from = val) %>% filter(Resolution==60),aes(x = Q0.1_Completeness, xend = Q0.9_Completeness, y = LEVEL3_NAM, yend = LEVEL3_NAM,color=factor(Resolution)), size = 1,position=position_nudge(y=0.5),alpha=0.5) + 
  #geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.1_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=15,position=position_nudge(y=0.5),alpha=0.7) +geom_point(data = tabla %>% pivot_wider(names_from = measure,values_from = val)%>% filter(Resolution==60),aes(Q0.9_Completeness, LEVEL3_NAM, color = factor(Resolution)),size = 2,shape=30,position=position_nudge(y=0.5),alpha=0.7) +
  geom_segment(data = tabla %>% filter(measure=="Median_phylo") %>% group_by(LEVEL3_NAM) %>% summarise(start = range(val)[1], end = range(val)[2]) %>% ungroup(),aes(x = start, xend = end, y = LEVEL3_NAM, yend = LEVEL3_NAM),color = "gray80", size = 1)  +
  geom_point(data = tabla %>% filter(measure=="Median_phylo"),  aes(val, LEVEL3_NAM, fill = factor(Resolution,levels=c(15,30,60))),size = 3,shape=21) +
  geom_point(data = tabla %>% filter(measure=="Median_phylo"),  aes(val, LEVEL3_NAM, fill = factor(Resolution,levels=c(15,30,60))),size = 3,shape=21) +
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
  labs(x="Median phylogenetic shortfall",y="",title="") +
  NULL



layout <- "
AB
CD"
pp <- wrap_plots(design = layout,list(comp,age,phylo,seqs)) +
  plot_layout(byrow=T,nrow = 2,guides="collect",tag_level = 'new') & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"),legend.position = "bottom")
pp

nombre <- paste0("Dumbbells.pdf")
ggsave(filename = here("figures",nombre),plot = pp)

}

map_hotspots <- function(data = "Allmetrics_15.Rdata",var = "Completeness",resolution = 15,distance = 75,threshold = 1.96){
  require(spdep)
  counts <- rnaturalearth::ne_countries(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  xlimits = c(-125,-75)
  ylimits = c(0,35)
  nombre = paste0("Spots",var,"_",resolution,".pdf")
  load(here("interim",data))
  datos %<>% mutate(Completeness = 100 - Completeness, seqs= 1- seqs, phylo = 1 - phylo)
  datos <- datos[!is.na(datos[[glue::glue('{var}')]]),]
  datos <- datos[!is.na(datos[["x"]]),]
  nb <- dnearneigh(as.matrix(datos[,c("x","y")]),0,distance,longlat = TRUE)
  aver <- localG(datos[[glue::glue('{var}')]],nb2listw(nb,style = "B",zero.policy=TRUE))
  datos$localG <- as.vector(aver)
  datos %>% mutate(spots = ifelse(localG >= threshold,"HC",ifelse(localG <= -threshold,"HC","N"))) %>% 
    ggplot(aes(x=x,y=y,color=localG,alpha=spots)) +
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
  ggsave(filename = here("figures",nombre),last_plot())
  
base <- raster("data/base_rasters/baseRaster_15.tif")
load("data/base_extent.r")
  prefix = paste0("raster",resolution,"_",var,".tif")
  rasterize(datos[,c("x","y")],base,field=datos[[glue::glue("{var}")]]) %>% crop(.,extent(x)) %>%   writeRaster(.,filename = here("output",prefix),overwrite=TRUE)
  prefix = paste0("raster",resolution,"_",var,"G",".tif")
  if(length(which(is.na(datos$localG))!=0)) {datos$localG[which(is.na(datos$localG))] <- -999}
  rasterize(datos[,c("x","y")],base,field=datos$localG) %>% crop(.,extent(x)) %>% writeRaster(.,filename = here("output",prefix),overwrite=TRUE)
  #prefix = paste0("raster",resolution,"_","Slope",".tif")
  #rasterize(datos[,c("x","y")],base,field=datos$Slope) %>% crop(.,extent(x)) %>%   writeRaster(.,filename = here("output",prefix),overwrite=TRUE)
  
}

hotcold_spots <- function(data="interim/Allmetrics_15.Rdata",extent="data/base_extent.r"){
load(extent)
base <- raster("output/raster15_CompletenessG.tif") %>% crop(.,x)
base[base %in% c(999,-999)] <- NA
base %>% getValues(.) -> compHC; length(compHC)
base <- raster("output/raster15_ageG.tif") %>% crop(.,x)
base %>% getValues(.) -> ageHC; length(ageHC)
base <- raster("output/raster15_seqsG.tif") %>% crop(.,x)
base %>% getValues(.) -> seqsHC; length(seqsHC)
base <- raster("output/raster15_phyloG.tif") %>% crop(.,x)
base %>% getValues(.) -> phyloHC; length(phyloHC)

base <- raster("dbVoCC/tmax/Velocidad/d50.tif") %>% crop(.,x)
base %>% getValues(.) -> vocc_tmax; length(vocc_tmax)
base <- raster("dbVoCC/ppt/Velocidad/d50.tif") %>% crop(.,x)
base %>% getValues(.) -> vocc_ppt; length(vocc_ppt)
base <- raster("dbVoCC/tmean/Velocidad/dv50.tif") %>% crop(.,x)
base %>% getValues(.) -> vocc_tmean; length(vocc_tmean)
base <- raster("dbVoCC/tmin/Velocidad/dv50.tif") %>% crop(.,x)
base %>% getValues(.) -> vocc_tmin; length(vocc_tmin)

hfp <- raster(here("data/LUL/LULC_1940__CMIP6.1.tif")) %>% crop(.,x)
getValues(hfp) -> lul_1940; length(lul_1940)
hfp <- raster(here("data/LUL/LULC_1950__CMIP6.1.tif")) %>% crop(.,x)
getValues(hfp) -> lul_1950; length(lul_1950)
hfp <- raster(here("data/LUL/LULC_1960__CMIP6.1.tif")) %>% crop(.,x)
getValues(hfp) -> lul_1960; length(lul_1960)
hfp <- raster(here("data/LUL/LULC_1970_CMIP6.1.tif")) %>% crop(.,x)
getValues(hfp) -> lul_1970; length(lul_1970)
hfp <- raster(here("data/LUL/LULC_1980__CMIP6.1.tif")) %>% crop(.,x)
getValues(hfp) -> lul_1980; length(lul_1980)
hfp <- raster(here("data/LUL/LULC_1990__CMIP6.1.tif")) %>% crop(.,x)
getValues(hfp) -> lul_1990; length(lul_1990)
hfp <- raster(here("data/LUL/LULC_2000__CMIP6.1.tif")) %>% crop(.,x)
getValues(hfp) -> lul_2000; length(lul_2000)
hfp <- raster(here("data/LUL/LULC_2015__CMIP6.1.tif")) %>% crop(.,x)
getValues(hfp) -> lul_2015; length(lul_2015)

load(data)
base <- raster(here("data/base_rasters/baseRaster_15.tif")) %>% crop(.,x)
datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("California","Arizona"),"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("New Mexico","Texas"),"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(LEVEL3_NAM%in%c("Colombia"),"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -106,"Mexico Northwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude < -98,"Mexico Southwest",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 25,"Mexico Northeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 16 & Longitude < -93,"Mexico Gulf",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Longitude > -83,"Panama",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 11,"Costa Rica",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 14 & Longitude < -92,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 14 & Longitude < -88,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 13,"Nicaragua",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 17 & Longitude > -88,"El Salvador",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 17 & Longitude < -90,"Mexico Southeast",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude < 18,"Belize",LEVEL3_NAM))
datos %<>% mutate(LEVEL3_NAM = ifelse(is.na(LEVEL3_NAM) & Latitude > 18 & Longitude > -90,"Mexico Southeast",LEVEL3_NAM))
datos %>% filter(is.na(LEVEL3_NAM))
datos %<>% mutate(lvl3=as.factor(LEVEL3_NAM),.after=1)
datos %<>% mutate(Completeness = ifelse(is.na(Completeness),-999,Completeness))
rasterize(datos[,c("x","y")],base,
 field=as.numeric(datos$lvl3))  %>% getValues(.) -> regions; length(regions)
rasterize(datos[,c("x","y")],base,
          field=as.numeric(datos$age))  %>% getValues(.) -> age; length(age)
rasterize(datos[,c("x","y")],base,
          field=as.numeric(datos$Completeness))  %>% getValues(.) -> comp; length(comp)
rasterize(datos[,c("x","y")],base,
          field=as.numeric(datos$phylo)) %>% getValues(.) -> phylo; length(phylo)
rasterize(datos[,c("x","y")],base,
          field=as.numeric(datos$seqs)) %>% getValues(.) -> seqs; length(seqs)

rasterize(datos[,c("x","y")],base,
          field=as.numeric(datos$SR))  %>% getValues(.) -> SR; length(SR)

aver <- raster(here("output/sumaHotspot.tif")) %>% crop(.,x)
aver %>% getValues(.) -> zonation; length(zonation)
base  %>% coordinates(.) -> coords; length(coords)/2

tibble(x=coords[,1],y=coords[,2],zone=zonation,vocc_tmax=vocc_tmax,
       vocc_ppt=vocc_ppt,vocc_tmean=vocc_tmean,vocc_tmin=vocc_tmin,
       compHC=compHC,ageHC=ageHC,seqsHC=seqsHC,phyloHC=phyloHC,age=age,comp=comp,phylo=phylo,seqs=seqs,
       lul_1940=lul_1940,lul_1950=lul_1950,lul_1960=lul_1960,lul_1970=lul_1970,
       lul_1980=lul_1980,lul_1990=lul_1990,lul_2000=lul_2000,lul_2015=lul_2015,
       regions=regions)  %>% 
  mutate(across(starts_with("lul"), ~case_when(.x %in% c(1:6, 9, 12) ~"Altered",.x %in% c(7,8) ~ "Vegetation",.x %in% c(10,11) ~ "Altered"))) ->  tiba
tiba %<>% mutate(comp = 100 - comp, phylo = 1 - phylo,seqs = 1 - seqs)
tiba <- tiba %>% mutate(regions=factor(regions,levels=1:13))
levels(tiba$regions) <- levels(datos$lvl3)
tiba %>% mutate(spot_comp = ifelse(compHC <= -1.96,"Coldspot",ifelse(compHC >= 1.96,"Hotspot",NA)),.after=compHC) %>% 
  mutate(spot_age = ifelse(ageHC <= -1.96,"Coldspot",ifelse(ageHC >= 1.96,"Hotspot",NA))) %>% 
  mutate(spot_seqs = ifelse(seqsHC <= -1.96,"Coldspot",ifelse(seqsHC >= 1.96,"Hotspot",NA))) %>% 
  mutate(spot_phylo = ifelse(phyloHC <= -1.96,"Coldspot",ifelse(phyloHC >= 1.96,"Hotspot",NA))) -> datas


library(UpSetR)
datas %>% select(regions,starts_with("spot")) %>% select(-spot_seqs) %>% 
  
  mutate(ID=1:nrow(.),.before=1) %>% pivot_longer(cols = -c(1:2),names_to = "type",values_to = "spots") %>% 
  mutate(type=sub("^spot_","",type)) %>% 
  mutate(val=1) %>% pivot_wider(names_from = c(spots,type),values_from = val) %>%  
mutate(across(-2,replace_na,0)) -> forup

forup %>% group_by(regions) %>% summarise(across(-ID,sum)) %>% 
  pivot_longer(-1,names_to = "Variable",values_to = "Values") %>% 
  filter(!is.na(regions)) %>% separate(Variable, into=c("Type","Dimension"),sep="_",remove = FALSE) %>% 
  
  group_by(regions,Dimension) %>% mutate(Values = Values/sum(Values)*100) %>% 
  filter(!grepl("NA_",Variable)) %>% 
  mutate(Variable = factor(Variable,levels=c("Coldspot_comp","Coldspot_age","Coldspot_phylo",
                                             "Hotspot_comp","Hotspot_age","Hotspot_phylo"))) %>% 
  ggplot(aes(x=regions,y=Values,fill=Variable)) + #facet_wrap(~Variable) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip() +
  scale_fill_manual(values=MetBrewer::met.brewer("Hiroshige",n=6,direction=-1),
                    guide = guide_legend(reverse = TRUE) ) +
  labs(x="",y="Spatial units") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.line.y = element_blank(),
        legend.position = c(0.9,0.5),
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
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="vertical") +
  geom_hline(yintercept = 0) +
  NULL
  
ggsave(filename = here("figures/HotCold_region.pdf"),last_plot())


upset(as.data.frame(forup %>% select(ID,starts_with("Coldspot"))),nsets = 4,sets = rev(c("Coldspot_comp","Coldspot_age","Coldspot_phylo")),keep.order = TRUE)

upset(as.data.frame(forup %>% select(ID,starts_with("Hotspot"))), nsets = 3, sets = rev(c("Hotspot_comp","Hotspot_age","Hotspot_phylo")),keep.order = TRUE)

datas %>% select(x,y,starts_with("spot_"),starts_with("lul")) %>% 
  pivot_longer(cols = starts_with("lul"),names_to = "lul",values_to = "type") %>% 
  filter(!is.na(type)) %>% 
  pivot_longer(cols = starts_with("spot_"),names_to = "Dimension",values_to = "spot") %>% 
  group_by(type,spot,lul) %>% 
  count() %>%
  group_by(spot,lul) %>% 
  mutate(n= n / sum(n)) %>%
  filter(!is.na(spot)) %>%
  #filter(Dimension!="spot_seqs") %>% 
  mutate(lul = sub("lul_","",lul)) %>% 
  ggplot(aes(y=n,x=lul,fill=factor(type,levels=c("Vegetation","Altered")))) + 
  facet_wrap(~spot,scales = "fixed") +
  geom_bar(stat="identity",position="stack") +
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
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="vertical")  + labs(title="",y="Percent area",x="") +
  #scale_y_continuous(limits=c(0,0.8)) +
  scale_fill_manual(values = MetBrewer::met.brewer("Kandinsky")[c(2:4)],name="") +
  NULL


ggsave(filename = here("figures/LULtrans_spots.pdf"),last_plot())


datas %>%  filter(!is.na(age)) %>% mutate(LU_change = paste0(lul_1960,"_",lul_2015),.after=1) %>% 
  mutate(LU_change = ifelse(LU_change=="NA_NA",NA,LU_change)) %>% 
  mutate(LU_change = ifelse(LU_change=="Altered_Vegetation", "Altered_Altered",LU_change)) %>% 
  select(x,y,LU_change,starts_with("spot_")) %>% 
  pivot_longer(cols = starts_with("spot_"),names_to = "dimension",values_to = "spot") %>% 
  filter(!is.na(LU_change)) %>% 
  group_by(spot,LU_change,dimension) %>% count() %>% 
  group_by(spot, dimension) %>% mutate(n= n / sum(n)) %>% 
  filter(dimension!="spot_seqs") %>% 
  ggplot(aes(y=n,x=dimension,fill=factor(LU_change))) +
  facet_wrap(~spot) +
  geom_bar(stat="identity",position="stack") +
  geom_text(aes(label=round(n,2)),position=position_stack(vjust = 0.5),family="EB Garamond",fontface="bold",color="white",) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = "right",
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=12),
        axis.title.x = element_text(family="EB Garamond",size=12),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.x =  element_text(family="EB Garamond",size = 10),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.6,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="vertical") +
  scale_fill_manual(values=MetBrewer::met.brewer("Kandinsky")[c(2,3,1)],name="Land-use changes\n(1960-2015)") + 
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_hline(yintercept = 0) +
  ylim(c(0,1)) +
  labs(y="Proportion of grid cells",x="shortfall") +
  NULL

ggsave(filename = here("figures/LUL_spots.pdf"),last_plot())

datas %>% filter(!is.na(age)) %>% select(x,y,vocc_ppt,starts_with("spot_"),-spot_seqs) %>% 
  pivot_longer(cols = starts_with("spot_"),names_to = "dimension",values_to = "spot") %>% 
  mutate(vcc_prec = ifelse(is.na(vocc_ppt),"non_analogue",ifelse(vocc_ppt==0,"analogue","within_dist"))) %>% 
  group_by(spot,dimension,vcc_prec) %>% count() %>% 
  group_by(spot, dimension) %>% mutate(n= n / sum(n)) %>% 
  ggplot(aes(y=n,x=dimension,fill=factor(vcc_prec,levels=c("non_analogue","within_dist","analogue")))) +
  facet_wrap(~spot) +
geom_bar(stat="identity",position="stack") +
  geom_text(aes(label=round(n,2)),position=position_stack(vjust = 0.5),family="EB Garamond",fontface="bold",color="white",) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = "right",
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=12),
        axis.title.x = element_text(family="EB Garamond",size=12),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.x =  element_text(family="EB Garamond",size = 10),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.6,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="vertical") +
  scale_fill_manual(values=MetBrewer::met.brewer("Hiroshige")[c(6,7,9)],name="Climate change\n(1910-2015)") + 
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_hline(yintercept = 0) +
  ylim(c(0,1)) +
  labs(y="Proportion of grid cells",x="shortfall") +
  NULL

ggsave(filename = here("figures/vocc_ppt_spots.pdf"),last_plot())

datas %>% filter(!is.na(age)) %>% select(x,y,vocc_tmean,starts_with("spot_"),-spot_seqs) %>% 
  pivot_longer(cols = starts_with("spot_"),names_to = "dimension",values_to = "spot") %>% 
  mutate(vcc_prec = ifelse(is.na(vocc_tmean),"non_analogue",ifelse(vocc_tmean==0,"analogue","within_dist"))) %>% 
  group_by(spot,dimension,vcc_prec) %>% count() %>% 
  group_by(spot, dimension) %>% mutate(n= n / sum(n)) %>% 
  ggplot(aes(y=n,x=dimension,fill=factor(vcc_prec,levels=c("non_analogue","within_dist","analogue")))) +
  facet_wrap(~spot) +
  geom_bar(stat="identity",position="stack") +
  geom_text(aes(label=round(n,2)),position=position_stack(vjust = 0.5),family="EB Garamond",fontface="bold",color="white",) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = "right",
        legend.background = element_rect(fill = NA),
        legend.spacing.y = unit(c(0.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond"),
        axis.title.y = element_text(family="EB Garamond",size=12),
        axis.title.x = element_text(family="EB Garamond",size=12),
        axis.text =  element_text(family="EB Garamond",size = 8),
        axis.text.x =  element_text(family="EB Garamond",size = 10),
        legend.margin=margin(t=-25),
        legend.key.size = unit(0.6,"cm"),
        plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),
        plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="vertical") +
  scale_fill_manual(values=MetBrewer::met.brewer("Hiroshige")[c(1,3,4)],name="Climate change\n(1910-2015)") + 
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_hline(yintercept = 0) +
  ylim(c(0,1)) +
  labs(y="Proportion of grid cells",x="shortfall") +
  NULL

ggsave(filename = here("figures/vocc_tmean_spots.pdf"),last_plot())

circlize::chordDiagram(x = sub %>% select(-N), 
                       grid.col = my_col,transparency = alphas,directional = 1,
                       direction.type = c("arrows"), 
                       link.arr.type = "big.arrow",
                       link.arr.length = 0.1,link.sort = T,
                       link.largest.ontop = TRUE,grid.border="black",
                       link.border=adjustcolor("black",0.5),
                       link.lwd=0.2#,order=names(my_col)
)

}

#### EXPERIMENTAL FUNCTIONS !!!! ##### NOT PROPERLY TESTED (JAN 2023)
age_description <- function(){ 
  checklist = "wcvp_names.txt"
  data="Joined_finalPOWdist.v2.Rdata"
  load(here("output",data))
  some %<>% filter(!is.na(Accepted_Name),Outlier_Test,POW_distribution)
  list_names <- data.table::fread(here("wcvp_2022",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
  list_names %<>% filter(taxon_status=="Accepted") %>% mutate(first_published = str_replace_all(first_published, "\\*|\\(|\\)", "")) %>% 
    mutate(first_published =  str_trim(first_published)) %>% 
    mutate(first_published_num= as.numeric(substr(first_published,1, 4)))
  list_names %<>% select(accepted_plant_name_id,first_published_num) %>% 
    group_by(accepted_plant_name_id) %>% mutate(Age_spp = min(first_published_num)) %>% 
    distinct(accepted_plant_name_id,.keep_all = TRUE)
  some %>% left_join(.,list_names,by="accepted_plant_name_id") -> alguito
  
  res <- 60/60
  projection = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  g <- raster::raster(nrows=180*res,ncols=360*res,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% projectRaster(.,crs = CRS(projection)) %>% as(., 'SpatialPixels')
  alguito %>% select(decimalLongitude,decimalLatitude) %>% 
    sf::st_as_sf(x = ., coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs = proj4string(g)) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID"))) %>% 
    select(CellID) %>% bind_cols(.,alguito) -> alguito
  
  g@coords %>% as_tibble() %>% mutate(CellID=1:nrow(.)) -> coords
  datos <- alguito %>% group_by(CellID) %>% mutate(Age_spp= median(2022 - Age_spp,na.rm=T),n=n())
  datos %<>% distinct(CellID,.keep_all = TRUE) %>% left_join(.,coords,by="CellID")
  range(datos$Age_spp,na.rm=TRUE)
  roads <- rnaturalearth::ne_countries(scale = 10,returnclass = "sf") %>% sf::st_transform(., crs = CRS(projection))
  xlimits = c(-12500000,-2500000)
  ylimits = c(-7000000,9000000)
  ggplot() + 
    xlim(xlimits) + ylim(ylimits) + 
    theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
          legend.position = c(0.15,0.3),
          legend.background = element_rect(fill = NA),
          #legend.spacing.y = unit(c(0.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          #legend.margin=margin(t=-25),
          legend.key.size = unit(0.5,"cm"),
          plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),
          plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),
          plot.caption = element_text(family="EB Garamond",size=8,hjust = 1),
          legend.direction="horizontal")  +
    geom_tile(data = datos %>% filter(!is.na(Age_spp),n>=10),aes(x=x,y=y,fill=Age_spp,color=Age_spp),color="white") +
    scale_fill_stepsn(colors=(MetBrewer::met.brewer("Hiroshige",direction=-1)),name="Age (years)",
                      breaks=c(10,20,40,60,80,100,120,140,160,180,200,240,269),
                      labels= scales::comma,na.value="NA",
                      guide=guide_legend(reverse = FALSE,ncol = 1,title.position="top",label.position="left", title.hjust = 0.5)) +
    scale_color_stepsn(colors=(MetBrewer::met.brewer("Hiroshige",direction=-1)),name="Age (years)",
                       breaks=c(10,20,40,60,80,100,120,140,160,180,200,240,269),
                       labels= scales::comma,na.value="NA",
                       guide=guide_legend(reverse = FALSE,ncol = 1,title.position="top",label.position="left", title.hjust = 0.5)) +
    labs(x="",y="",title="The geography of plant species ages in America",subtitle="Based on the year of species' publication",caption="* Age was estimated as the mean age of publication (relative to 2022)\nof species ocurring in a given spatial pixel") +
    # standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
    geom_sf(data=roads,colour="black",fill=NA,size=0.01) +
    NULL
} #### This function is still experimental!
