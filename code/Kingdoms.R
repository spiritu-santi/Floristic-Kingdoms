library(data.table)
library(tidyverse)
library(sp)
library(furrr)
library(CoordinateCleaner)
library(ape)

######################################################################################################
#####         ALL THE FOLLOWING IS NOT AUTOMATED (last update: NOVEMBER 10, 2020)                #####
##### The following script does all the necessary processing for the geographic data as obtained #####
##### from GBIF. There is a pre-processing step in bash to reduce file size, but row ID can be   #####
##### traced back to the original data as downloaded from GBIF with extra information.           #####
##### The script includes the necessary code to perform all the analyses and figures presented   #####
##### in Ramírez-Barahona et al. (2020)                                                          #####
######################################################################################################
######################################################################################################
###### PRE-PROCESSING: SETTING UP DATA FILES RETRIEVED FROM GBIF #####
###### USE BASH TO REDUCE FILE SIZE BY SELECTING COLUMNS
### bash commands are commented out!
# awk -F"\t" '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$13"\t"$22"\t"$23"\t"$33}' TRAQUEOS_raw.csv > traqueos_COR.txt
# wc -l TRAQUEOS_raw.csv
# wc -l traqueos_COR.txt
##### PREPARE FILE FOR PYTHON SCRIPT
##### the GBIF file has been renamed to: TRAQUEOS_raw
pre.process <- function(GBIF.data){ 
  py <- fread(GBIF.data);dim(py);colnames(py)
  py <- cbind(py[,7:9],c(1:dim(py)[1])) #### ADD A NUMERIC ID FOR REFEENCE TO THE ORIGINAL DATA
  py[which(py$species!=""),]->py
  cat(dim(py),"\n")
  file.out <- toupper(paste(sub("_.*","",GBIF.data)))
  write.table(py,paste("data/",file.out,"_python.csv",sep=""),col.names = F,row.names = F,quote = F,sep=",")
}
###### FIRST STEP: PERFORM CLEANING USING PYTHON SCRIPT OF Edwards et al. #######
# allHerbaria_ADM1_badCoords.txt is provided! 
# python code/cleanGbifCoords.1.0.py data/TRAQUEOS_python.csv data/allHerbaria_ADM1_badCoords.txt data/TRAQUEOS_good.csv data/TRAQUEOS_bad.csv
###### SECOND STEP: ELIMANTE UNIDENTIFIED RECORDS ########
a <- fread("TRAQUEOS_good.csv"); dim(a) ####### LIST OF GOOD OCURRENCES (ID, LAT, LON, EXTRA)
#bad<-fread("TRAQUEOS_bad.csv"); dim(bad) ####### LIST OF BAD OCURRENCES (ID, LAT, LON, EXTRA)
aa <- fread("traqueos_COR.txt"); dim(aa) ###### ORIGINAL DATA (COMPLETE COLUMNS)
iden<-as.vector(a$V4); b <- aa[iden,] ##### EXTRA COLUMN HAS THE ROWNAME FROM THE ORIGINAL
dim(b); dim(a) ##### CHECK IF LENGTH IS OK!

## Prepare names
aa<-b; aa$species<-sub(" ","_",aa$species)
aa$scientificName <- sub(" ","_",aa$scientificName); aa$scientificName <- sub(" ","_",aa$scientificName)
aa$scientificName_binomial <- paste(unlist(lapply(strsplit(aa$scientificName,"_"),"[[",1)),unlist(lapply(strsplit(aa$scientificName,"_"),"[[",2)),sep="_")
dim(b); dim(aa)
rm(list = c("a","b")); gc()
###### THIRD STEP: NAME FILTERING AND HOMOGENEIZATION:  POW ########
aa<-aa[-which(aa$class==""),]
fams<-read.table("data/APW_synonyms.txt",sep=",") ###### LIST OF FAMILY NAMES
colnames(fams)<-c("Name","Synonym","Order")
fams$ACCEPTED<-rep(NA,dim(fams)[1])
for(i in 1:dim(fams)[1]) {if((fams$Synonym[i]=="")==TRUE){fams$ACCEPTED[i]<-as.character(fams$Name[i])}
  else (fams$ACCEPTED[i]<-as.character(fams$Synonym[i]))}
match(aa$family,fams$Name)->iden
aa$family_ACCEPTED<-fams$ACCEPTED[iden]
list_names<-fread("data/wcs_and_atoz_genera_species_infraspecies_2019_apr.txt") ## data provided directly from Kew's IT (we do not provide it here)
# Name matching
list_names <- list_names[-which(list_names$species==""),]
list_names$Complete_NAME<-paste(list_names$genus,list_names$species,sep="_")
list_names$Complete_NAME_author<-paste(list_names$genus,list_names$species,list_names$b_authors,sep="_")
list_names$Complete_NAME_author <- gsub('._$', '', list_names$Complete_NAME_author)
list_names$New_NAME_ACCEPTED<-list_names$Complete_NAME[match(list_names$accepted_db_id,list_names$db_id,nomatch = NA)]
list_names$New_NAME_ACCEPTED_author<-list_names$Complete_NAME_author[match(list_names$accepted_db_id,list_names$db_id,nomatch = NA)]
list_names[is.na(list_names$New_NAME_ACCEPTED)]$New_NAME_ACCEPTED<-list_names[is.na(list_names$New_NAME_ACCEPTED)]$Complete_NAME
list_names[is.na(list_names$New_NAME_ACCEPTED_author)]$New_NAME_ACCEPTED_author<-list_names[is.na(list_names$New_NAME_ACCEPTED_author)]$Complete_NAME_author

match(aa$scientificName_binomial, list_names$Complete_NAME) -> iden
length(which(is.na(iden)))
aa$species_ACCEPTED_POW <- list_names$New_NAME_ACCEPTED[iden]
match(aa$scientificName,list_names$Complete_NAME_author)->iden
length(which(is.na(iden)))
aa$species_ACCEPTED_POW_author <- list_names$New_NAME_ACCEPTED[iden]
save(aa,file="data/TRAQUEOS_corNAMES.r") ##### DATA BASE

###### FOURTH STEP: MINOR EDITS AND SAVE ########
load("data/TRAQUEOS_corNAMES.r")
length(unique(aa$family_ACCEPTED))
aa[which(aa$species=="Xyris_araracuare"),"species_ACCEPTED_TPL"]<-"Xyris_araracuarae"
aa[which(aa$species=="Xyris_araracuare"),"species_ACCEPTED_POW"]<-"Xyris_araracuarae"
aa[which(aa$species=="Xyris_araracuare"),"family_ACCEPTED"]<-"Xyridaceae"
aa[which(aa$species=="Xyris_araracuare"),"family"]<-"Xyridaceae"
aa[which(aa$genus=="Afrothismia"),"family_ACCEPTED"]<-"Afrothismieae"
##### REMOVE SOME FOSSILS
aa<-aa[-c(which(aa$genus=="Liliacidites")),]
aa<-aa[-c(which(aa$species=="Lindernia_confusa")),]
aa<-aa[-c(which(aa$genus=="Thalassotaenia")),]
aa<-aa[-c(which(aa$genus=="Polypodiisporites")),]
aa<-aa[-c(which(aa$genus=="Otozamites")),]
aa<-aa[-c(which(aa$genus=="Nilssoniopteris")),]
aa<-aa[-c(which(aa$genus=="Pterophyllum")),]
aa<-aa[-c(which(aa$genus=="Ctenophyllum")),]
aa<-aa[-c(which(aa$genus=="Dictyozamites")),]
aa<-aa[-c(which(aa$genus=="Elatocladus")),]
aa<-aa[-c(which(aa$genus=="Cunninghamites")),]
aa[which(aa$family==""),]
save(aa,file="data/TRAQUEOS_corNAMES.r") ##### DATA BASE
write.table(sort(unique(aa$family_ACCEPTED)),file="data/Fams_corNAMES.txt",row.names=F,quote=F)

###### FIFTH STEP: FAMILY RE-CIRCUMSCRIPTION TO MATCH CLASSIFICATION USED IN ANGIOS TIME-TREE 2.0 ########
load("data/TRAQUEOS_corNAMES.r")
aa$family_ACCEPTED_2 <- aa$family_ACCEPTED
tipas<- read.table("data/TIPAS_v.7.1.csv",sep=",",header=T); fa_ti<-unique(tipas$Family)
fams <- sort(unique(aa$family_ACCEPTED))
fa_ti[which(is.na(match(fa_ti,fams)))] ######## IDENTIFY FAMILIES IN PHYLOGENY THAT ARE NOT IN THE DATA BASE
####### THIS IS A BIT OF A NASTY CODE..... NEED TO HAVE A FILE WITH THE NAMES
######## TACCACEA
aa[grep("Tacca_",aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Taccaceae"
aa[grep("Tacca_",aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Taccaceae"
######## CODONACEAE + COLDENIACEAE + HOPLESTIGMATACEAE
aa[grep("Codon_",aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Codonaceae"
aa[grep("Coldenia_",aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Coldeniaceae"
aa[grep("Hoplestigma_",aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Hoplestigmataceae"
aa[grep("Codon_",aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Codonaceae"
aa[grep("Coldenia_",aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Coldeniaceae"
aa[grep("Hoplestigma_",aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Hoplestigmataceae"
####### CORDIACEAE
names.in <- c("Cerdana_","Cienkowskya_","Cordia_","Catonia_","Carpiphea_","Calyptracordia_",
              "Bourgia_","Auxemma_","Ascania_","Acnadena_","Varroniopsis_","Varronia_","Ulmarronia_",
              "Toquera_","Topiaris_","Sebestena_","Salimori_","Saccellium","Quarena_","Rhabdocalyx_",
              "Plethostephia_","Piloisia_","Pilicordia_","Physoclada_","Patagonula_","Paradigmia_",
              "Novella_","Myxa_","Montjolya_","Macria_","Maciela_","Lithocardium_","Hymenesthes_",
              "Hemigymnia_","Gynaion_","Gerascanthus_","Firensia_","Ectemis_","Diacoria_","Cordiopsis_",
              "Cordiada_","Coilanthera_","Collococcus_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Cordiaceae"
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Cordiaceae"
###### EHRETIACEAE
names.in <- c("Ammobroma_","Antrophora_","Beurreria_","Bourreria_","Carmona_","Corallophyllum_","Cortesia_",
              "Crematomia_","Desmophyla_","Diplostylus_","Eddya_","Ehretia_","Galapagoa_","Gaza_","Halgania_",
              "Hilsenbergia_","Lennoa_","Lepidocordia_","Lithothamnus_","Lutrostylis_","Menais_","Monomesia_",
              "Morelosia_","Ptilocalyx_","Pholisma_","Rhabdia_","Rhabdia_","Rotula_","Stegnocarpus_","Subrisia_",
              "Tetracoccus_","Tiquilia_","Tiquiliopsis_","Traxilum_","Zombiana_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Ehretiaceae"
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Ehretiaceae"
####### WELLSTEDIACEAE
aa[grep("Wellstedia_",aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Wellstediaceae"
aa[grep("Wellstedia_",aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Wellstediaceae"
###### NAMACEAE
names.in <- c("Andropus_","Conanthus_","Eriodictyon_","Ernstamra_","Lemmonia_","Marilaunidium_","Nama_",
              "Turricula_","Wigandia_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Namaceae"
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Namaceae"
###### HELIOTROPIACEAE
names.in <- c("Argusia_","Beruniella_","Bourjotia_","Bucanion_","Ceballosia_","Cochranea_",
              "Dialion_","Eliopia_","Euploca_","Heliophytum_","Heliotropium_","Hieranthemum_","Hilgeria_",
              "Ixorhea_","Lithococca_","Mallotonia_","Meladendron_","Messerschmidia_","Messersmidia_","Myriopus_",
              "Nogalia_","Oskampia_","Parabouchetia_","Peristema_","Piptoclaina_","Pittonia_",
              "Schleidenia_","Schobera_","Scorpianthes_","Scorpiurus_","Spilocarpus_","Synzistachium_","Tetrandra_",
              "Tiaridium_","Tournefortia_","Valentina_","Valentiniella_","Verrucaria_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Heliotropiaceae"
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Heliotropiaceae"
######## MAUNDIACEAE
aa[grep("Maundia_",aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Maundiaceae"
aa[grep("Maundia_",aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Maundiaceae"
####### MYODOCARPACEAE
names.in <- c("Delarbrea_","Myodocarpus_","Porospermum_","Pseudosciadium_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Myodocarpaceae"
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Myodocarpaceae"
####### MYSTROPETALACEAE
names.in <- c("Dactylanthus_","Hachettea_","Mystropetalon_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Mystropetalaceae"
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Mystropetalaceae"
###### PELTANTHERACEAE
aa[grep("Peltanthera_",aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Peltantheraceae"
aa[grep("Peltanthera_",aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Peltantheraceae"
###### PETIVERIACEA
names.in <- c("Flueckigera_","Gallesia_","Hilleria_","Ledenbergia_","Mohlana_","Monococcus_",
              "Petiveria_","Rivina_","Schindleria_","Seguieria_","Trichostigma_","Villamilla_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Petiveriaceae"
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Petiveriaceae"
###### THISMIACEAE + AFROTHISMIA
names.in <- c("Bagnisia_","Cymbocarpa_","Geomitra_","Glaziocharis_","Haplothismia_","Mamorea_",
              "Myostoma_","Ophiomeris_","Oxygyne_","Saionia_","Sarcosiphon_","Scaphiophora_","Thismia_",
              "Tiputinia_","Tribrachys_","Triscyphus_","Triurocodon_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Thismiaceae"
aa[grep("Afrothismia_",aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Afrothismieae"
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Thismiaceae"
aa[grep("Afrothismia_",aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Afrothismieae"
######## NYSSACEAE
names.in <- c("Camptotheca_","Davidia_","Diplopanax_","Mastixia_","Nyssa_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_TPL),"family_ACCEPTED"]<-"Nyssaceae"
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED_POW),"family_ACCEPTED_2"]<-"Nyssaceae"

####### LAST EDITS
aa[grep("Stixaceae",aa$family),"family_ACCEPTED"] <- "Resedaceae"
aa[grep("Vivianiaceae",aa$family),"family_ACCEPTED"]<-"Francoaceae"
aa[grep("Stixaceae",aa$family),"family_ACCEPTED_2"]<-"Resedaceae"
aa[grep("Vivianiaceae",aa$family),"family_ACCEPTED_2"]<-"Francoaceae"
familias <- c("Plagiogyriaceae", "Osmundaceae", "Desmophlebiaceae", "Didymochlaenaceae", 
              "Rhachidosoridaceae", "Tempskyaceae", "Braithwaiteaceae", "Lepidodendraceae", 
              "Calamitaceae", "Sigillariaceae")
replace <- c("Plagiogyriaceae", "Osmundaceae", "Aspleniaceae", "NA", 
             "NA", "NA", "NA", "NA","NA", "NA")
for (i in 1:length(familias)){
  aa[which(aa$family==familias[i]),"family_ACCEPTED"] <- replace[i]
}
save(aa,file="data/TRAQUEOS_Corrected.r") ########## SAVE FILE

##### POW FLAGGING #####
load("data/TRAQUEOS_Corrected.r")
p_clean <- as_tibble(aa)
names(p_clean)[7:8] <- c("decimallatitude","decimallongitude")
p_clean <- p_clean[-c(which(is.na(p_clean$decimallatitude)),which(is.na(p_clean$decimallongitude))),]
p_clean$IDs <- 1:dim(p_clean)[1]
names(p_clean)[15] <- "family_ACCEPTED_POW"
no_match <- which((p_clean$species_ACCEPTED_POW == p_clean$species_ACCEPTED_POW_author)==FALSE)
p_clean$Resolved_ACCEPTED <- p_clean$species_ACCEPTED_POW
p_clean[no_match,"Resolved_ACCEPTED"] <- p_clean$species_ACCEPTED_POW_author[no_match]
length(unique(p_clean$Resolved_ACCEPTED))
length(which(is.na(p_clean$Resolved_ACCEPTED)))

poly <- rgdal::readOGR("data/wgsrpd-master/level3/level3.shp")
pointos <- SpatialPoints(as.data.frame(p_clean[,8:7]))
proj4string(pointos)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
points_powo <- sp::over(pointos,poly)
p_clean <- bind_cols(p_clean,points_powo)
names(p_clean)

pow_dist <- list_names<-fread("wcs_and_atoz_genera_species_infraspecies_dist_2019.txt")
list_names <- fread("wcs_and_atoz_genera_species_infraspecies_2019_apr.txt")
list_names <- list_names[-which(list_names$species==""),]
list_names <- list_names[-which(list_names$infraspecific_rank!=""),]
list_names$Complete_NAME <- paste(list_names$genus,list_names$species,sep="_")
list_names$Complete_NAME_author<-paste(list_names$genus,list_names$species,list_names$b_authors,sep="_")
list_names$Complete_NAME_author <- gsub('._$', '', list_names$Complete_NAME_author)
list_names$New_NAME_ACCEPTED <- list_names$Complete_NAME[match(list_names$accepted_db_id,list_names$db_id,nomatch = NA)]
list_names$New_NAME_ACCEPTED_author <- list_names$Complete_NAME_author[match(list_names$accepted_db_id,list_names$db_id,nomatch = NA)]
list_names[is.na(list_names$New_NAME_ACCEPTED)]$New_NAME_ACCEPTED <- list_names[is.na(list_names$New_NAME_ACCEPTED)]$Complete_NAME
list_names[is.na(list_names$New_NAME_ACCEPTED_author)]$New_NAME_ACCEPTED_author <- list_names[is.na(list_names$New_NAME_ACCEPTED_author)]$Complete_NAME_author
pow_dist$Accepted_NAME <- list_names$New_NAME_ACCEPTED[match(pow_dist$accepted_db_id,list_names$accepted_db_id)]
pow_dist$Accepted_NAME_author <- list_names$New_NAME_ACCEPTED_author[match(pow_dist$accepted_db_id,list_names$accepted_db_id)]
pow_dist <- pow_dist[which(pow_dist$introduced==0),]
rm("points_powo","poly","aa","pointos","no_match")

list.spp <- unique(p_clean$Resolved_ACCEPTED)
p_clean$Correct_POW_dist <- FALSE
p_clean %>% dplyr::arrange(Resolved_ACCEPTED) -> p_clean

f <- rep(seq_len(ceiling(dim(p_clean)[1] / 100000)),each = 100000,length.out = ceiling(dim(p_clean)[1]))
splita <- split(p_clean,f)
counter <- length(splita)

master_clean <- lapply(splita,function(lis) { 
  list.spp <- unique(lis$Resolved_ACCEPTED)
  counter <<- counter - 1
  cat("  ", "\n")
  cat("Chunks to go:", counter, "\n")
  cat("  ", "\n")
for (i in 1:length(list.spp)) { 
  cat(i,"--","Getting POW data for", list.spp[i],"               ","\r")
  iden <- which(lis$Resolved_ACCEPTED==list.spp[i])
  if(length(which(!is.na(match(as_vector(lis$LEVEL3_COD[iden]),
                               as_vector(pow_dist[which(pow_dist$Accepted_NAME==list.spp[i]),"area_code_l3"])))))==0) { 
    lis$Correct_POW_dist[iden] <- TRUE;
    next}
  lis$Correct_POW_dist[iden][which(!is.na(match(as_vector(lis$LEVEL3_COD[iden]),as_vector(pow_dist[which(pow_dist$Accepted_NAME==list.spp[i]),"area_code_l3"]))))] <- TRUE
}
  return(lis)
}
)
save(master_clean,file="data/master_clean.Rdata")

###### SIXTH STEP: GEOGRAPHICAL OUTLIER DETECTION BY FAMILY ########
load("data/master_clean.Rdata")
p_clean <- data.table::rbindlist(master_clean)
rm("master_clean")
p_clean <- as_tibble(p_clean);p_clean
aa_na <- p_clean[which(is.na(p_clean$Resolved_ACCEPTED)),]
p_clean$Correct_geographic <- NA
fami <- (unique(p_clean$family_ACCEPTED))
data("buffland")
#### FIRST THE BASIC TESTS
#### THIS CAN BE QUICKER.....
for (i in 1:length(fami)){ 
  fama <- p_clean[which(p_clean$family_ACCEPTED == fami[i]),]
  cat(toupper(fami[i]),"--",i,"--","\n")
  outas_1 <- clean_coordinates(x=fama,lon="decimallongitude",lat="decimallatitude",
                               test=c("duplicates","centroids","capitals","gbif","seas","institutions","equal","zeros"),
                               outliers_method = "distance",inst_rad=0.05,seas_ref=buffland,
                               capitals_rad=0.05,centroids_rad = 0.05,zeros_rad = 0.2,
                               species="Resolved_ACCEPTED",value="flagged",verbose=T)
  p_clean[which(p_clean$family_ACCEPTED == fami[i]),"Correct_geographic"] <- outas_1 
}
table(p_clean$Correct_geographic)

###### EXTRA OUTLIER CORRECTION ON TWO FAMILIES
p_clean[which(p_clean$family_ACCEPTED=="Tetracarpaeaceae"),"Resolved_ACCEPTED"]<- "Tetracarpaea_tasmannica"
p_clean[which(p_clean$family_ACCEPTED=="Kewaceae"),"Resolved_ACCEPTED"] <- p_clean[which(p_clean$family_ACCEPTED=="Kewaceae"),"species"]
fami <- c("Kewaceae","Tetracarpaeaceae")
for (i in 1:length(fami)){ 
  fama <- p_clean[which(p_clean$family_ACCEPTED == fami[i]),]
  cat(toupper(fami[i]),"--",i,"--","\n")
  outas_1 <- clean_coordinates(x=fama,lon="decimallongitude",lat="decimallatitude",
                               test=c("duplicates","centroids","capitals","gbif","seas","institutions","equal","zeros"),
                               outliers_method = "distance",inst_rad=0.05,seas_ref=buffland,
                               capitals_rad=0.05,centroids_rad = 0.05,zeros_rad = 0.2,
                               species="Resolved_ACCEPTED",value="flagged",verbose=T)
  p_clean[which(p_clean$family_ACCEPTED == fami[i]),"Correct_geographic"] <- outas_1 
}
dim(p_clean);table(p_clean$Correct_geographic)
p_clean_2 <- p_clean
table(p_clean_2$Correct_POW_dist)
table(p_clean_2$Correct_geographic)
fami <- (unique(p_clean_2$family_ACCEPTED))
length(unique(p_clean_2$Resolved_ACCEPTED))
p_clean_2$Second_test <- NA
#### THIS CAN BE QUICKER.....
for (i in 1:length(fami)){ 
  fama <- p_clean_2[which(p_clean_2$family_ACCEPTED == fami[i]),]
  cat(toupper(fami[i]),"--",i,"--","\n")
  outas_1 <- clean_coordinates(x=fama,lon="decimallongitude",lat="decimallatitude",
                               test=c("duplicates","centroids","capitals","gbif","seas","institutions","equal","zeros"),
                               outliers_method = "distance",outliers_td=500,inst_rad=0.05,seas_ref=buffland,
                               capitals_rad=0.05,centroids_rad = 0.05,zeros_rad = 0.2,
                               species="Resolved_ACCEPTED",value="flagged",verbose=T)
  p_clean_2[which(p_clean_2$family_ACCEPTED == fami[i]),"Second_test"] <- outas_1
}

p_clean_2[which(p_clean_2$decimallatitude <= -60),"Correct_geographic"] <- FALSE
f <- which(!p_clean_2$Correct_geographic | !p_clean_2$Correct_POW_dist)
save(p_clean_2,file="data/TRAQUEOS_Corrected_POWO_v1.r")




load("data/TRAQUEOS_Corrected_POWO_v1.r")
#### read database
p_clean_2 %>% filter(Correct_POW_dist,Correct_geographic)
p_clean_2 %>% filter(Correct_POW_dist,Correct_geographic) %>% dplyr::select(class) %>% unique()
ape::read.tree("data/Trees_Corrected/GBOTB_TS.tre")->a
p_clean_2 %>% filter(Correct_POW_dist,Correct_geographic) %>% dplyr::select(Resolved_ACCEPTED) %>% unique() -> b
p_clean_2 %>% filter(Correct_POW_dist,Correct_geographic) -> p_clean
p_clean %>% select(species,Resolved_ACCEPTED,class) %>% distinct() -> target_names

pow_dist <- list_names<-fread("wcs_and_atoz_genera_species_infraspecies_dist_2019.txt")
list_names <- fread("wcs_and_atoz_genera_species_infraspecies_2019_apr.txt")
list_names <- list_names[-which(list_names$species==""),]
list_names <- list_names[-which(list_names$infraspecific_rank!=""),]
list_names$Complete_NAME <- paste(list_names$genus,list_names$species,sep="_")
list_names$Complete_NAME_author<-paste(list_names$genus,list_names$species,list_names$b_authors,sep="_")
list_names$Complete_NAME_author <- gsub('._$', '', list_names$Complete_NAME_author)
list_names$New_NAME_ACCEPTED <- list_names$Complete_NAME[match(list_names$accepted_db_id,list_names$db_id,nomatch = NA)]
list_names$New_NAME_ACCEPTED_author <- list_names$Complete_NAME_author[match(list_names$accepted_db_id,list_names$db_id,nomatch = NA)]
list_names[is.na(list_names$New_NAME_ACCEPTED)]$New_NAME_ACCEPTED <- list_names[is.na(list_names$New_NAME_ACCEPTED)]$Complete_NAME
list_names[is.na(list_names$New_NAME_ACCEPTED_author)]$New_NAME_ACCEPTED_author <- list_names[is.na(list_names$New_NAME_ACCEPTED_author)]$Complete_NAME_author
pow_dist$Accepted_NAME <- list_names$New_NAME_ACCEPTED[match(pow_dist$accepted_db_id,list_names$accepted_db_id)]
pow_dist$Accepted_NAME_author <- list_names$New_NAME_ACCEPTED_author[match(pow_dist$accepted_db_id,list_names$accepted_db_id)]
pow_dist <- pow_dist[which(pow_dist$introduced==0),]

#### read trees and rename with POWs database
list.files("data/Trees_originals/") -> arboles
arboles <- arboles[-c(6,7)]
total <- c()
p_clean$class %>% unique() %>% length()
grupos <- matrix(NA,nrow=11,ncol = length(arboles))
rownames(grupos) <- p_clean$class %>% unique()
colnames(grupos) <- arboles
grupos_after <- grupos
grupos_after_thought <- grupos
for (i in 1:length(arboles)){ 
  cat("reading",arboles[i],"\n")
ALLMB <- ape::read.tree(paste0("Trees/",arboles[i]))
c(total,length(ALLMB$tip.label)) -> total
cat("Tree size (number of species):",length(ALLMB$tip.label),"\n")
## Matching the original names to the POW names 'Resolved_ACCEPTED'
cat("Matching species without correction:",length(which(!is.na(match(ALLMB$tip.label,target_names$Resolved_ACCEPTED)))),"\n")
match(ALLMB$tip.label,target_names$Resolved_ACCEPTED) %>% as_tibble() %>% filter(!is.na(value)) %>% as_vector() -> id
target_names$class[id] %>% table() -> before
grupos[match(names(before),rownames(grupos)),i] <- before

## Changing names in tree directly from POW database, and then matching.
list_names$New_NAME_ACCEPTED[match(ALLMB$tip.label,list_names$Complete_NAME)] -> jaja
jaja <- jaja[!is.na(jaja)]
ALLMB$tip.label[which(!is.na(match(ALLMB$tip.label,list_names$Complete_NAME)))] <- jaja
cat("Matching species after second correction:",length(which(!is.na(match(ALLMB$tip.label,target_names$Resolved_ACCEPTED)))),"\n")
target_names$Resolved_ACCEPTED[match(ALLMB$tip.label,target_names$Resolved_ACCEPTED)] -> resolvas
resolvas <- resolvas[!is.na(resolvas)]
ALLMB$tip.label[which(!is.na(match(ALLMB$tip.label,target_names$Resolved_ACCEPTED)))] <- resolvas
match(ALLMB$tip.label,target_names$Resolved_ACCEPTED) %>% as_tibble() %>% filter(!is.na(value)) %>% as_vector() -> id
target_names$class[id] %>% table() -> after
grupos_after_thought[match(names(after),rownames(grupos_after_thought)),i] <- after
drop.tip(ALLMB,tip = which(duplicated(ALLMB$tip.label))) -> ALLMB
write.tree(ALLMB,file=paste0("data/Trees_Corrected/",arboles[i],"_namesCorrected.tre"))

}

## RANDOM CODE??
cbind(Grupo=rownames(grupos),grupos) %>% as_tibble() ->  grupos
#cbind(Grupo=rownames(grupos_after),grupos_after) %>% as_tibble() ->  grupos_after
cbind(Grupo=rownames(grupos_after_thought),grupos_after_thought) %>% as_tibble() ->  grupos_after_thought
clase <- c("Magnoliophyta","Magnoliophyta","Pinophyta","Monilophyta","Monilophyta","Monilophyta", "Pinophyta","Lycophyta",
           "Pinophyta","Monilophyta", "Pinophyta" )
grupos %>% bind_cols(Class=clase,.) %>% mutate(across(.cols = 3:8,.fns = as.numeric)) %>% group_by(Class) %>% 
  summarise(across(.cols = 2:7,sum)) -> sum_before
#grupos_after %>% bind_cols(Class=clase,.) %>% mutate(across(.cols = 3:8,.fns = as.numeric)) %>% group_by(Class) %>% 
#  summarise(across(.cols = 2:7,sum)) -> sum_after
grupos_after_thought %>% bind_cols(Class=clase,.) %>% mutate(across(.cols = 3:8,.fns = as.numeric)) %>% group_by(Class) %>% 
  summarise(across(.cols = 2:7,sum)) -> sum_after_thought
names(total) <- arboles
total %>% as_tibble() %>% bind_cols(Tree=arboles) %>% rename("Total_spp"=value)
bind_rows(sum_before,sum_after_thought) %>% mutate(Process=c(rep("Original",4),rep("Correction_POW",4))) %>% 
  select(Class,Process,ALLMB.tre:TS_namesCorrected.tre) %>% arrange(Class) %>% write.table(.,file="Numbers_by_tree.csv",sep=",",
                                                                                           quote = F,row.names = F,col.names = T)
#### search in database
library(ape)
one <- ape::read.tree("data/Trees_Corrected/GBOTB.extended.tree_namesCorrected.tre")
target_node <- c(grep("Ophioglossum",one$tip.label)[1],grep("Cyathea",one$tip.label)[1]) %>% 
getMRCA(one,.)
two <- ape::read.tree("Trees_Corrected/TS_namesCorrected.tre_namesCorrected.tre")
two$tip.label[which(two$tip.label=="Dendrobium_moniliforme")] <- "Onychium_japonicum"
to_merge <- c(grep("Equisetum",two$tip.label)[1],grep("Cyathea",two$tip.label)[1]) %>% 
  getMRCA(two, .) %>%  extract.clade(two,.)
max(branching.times(to_merge))
extract.clade(one,target_node) %>% branching.times() %>% max() -> to_scale
rescaleTree<-function(tree,scale){
  tree$edge.length<-
    tree$edge.length/max(nodeHeights(tree)[,2])*scale
  return(tree)
}
to_merge$root.edge <- 0
to_merge <- rescaleTree(to_merge,extract.clade(one,target_node) %>% branching.times() %>% max())
max(branching.times(to_merge))
extract.clade(one,target_node) %>% branching.times() %>% max()

sub_tre <- extract.clade(one,target_node)
which(one$tip.label %in% sub_tre$tip.label) -> id
one$tip.label[id] <- paste(one$tip.label[id],"XXXXX",sep="X")
to_merge$root.edge <- 0
bind.tree(one,to_merge,where = target_node) -> binded

drop.tip(binded,grep("XXXXX",binded$tip.label)) -> binded

drop.tip(binded,which(!is.na(as.numeric(binded$tip.label)))) %>%  write.tree(.,file="data/GBOTB_TS.tre")



### RANDOM CODE???
binded$tip.label[which(duplicated(binded$tip.label))]
random <- LETTERS
binded_sub <- extract.clade(binded,node=getMRCA(binded,tip=c("Cyathea_mutica","Alsophila_celebica")))
plot(binded_sub,cex=0.5)
for (i in 1:length(LETTERS)){
hh <- runif(1,0.01,0.09)
r_node <- sample(1:length(binded_sub$node.label))
tip <- random[i]
all_tips <- binded_sub$tip.label[grep("Cyathea_",binded_sub$tip.label)]
sister <- all_tips[sample(length(all_tips),1)]
binded_sub <- bind.tip(binded_sub,tip,where=which(binded_sub$tip.label==sister),
                       edge.length = hh*binded_sub$edge.length[which(binded_sub$edge[,2]==
                                                                        which(binded_sub$tip.label==sister))],
               position=hh*binded_sub$edge.length[which(binded_sub$edge[,2]==
                                                     which(binded_sub$tip.label==sister))])
}
plot(binded_sub,cex=0.5)
is.ultrametric(binded_sub)
