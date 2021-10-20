### Brain morphology analysis 

### load appropriate packages 
library(tidyverse)
library(geomorph)
library(ggplot2)
library(cowplot)
library(ape)
library(phytools)
library(geiger)
library(tidytree)
library(picante)
library(pca3d)
library(Rphylopars)
library(dplyr)


### Table of Contents
### 1) read data
### 2) remove snakes from dataframe, leaving only lizards
### 3) convert 2D array to 3D array
### 4) create individual dataframes for each major brain region 


### Read in tree
s1.tree <- read.nexus("BrainGeomorph.nex")

s1.tree$tip.label <- tolower(s1.tree$tip.label) # make species Id's lowercase

plot(s1.tree)

### Read in trait data
s1.traits <- read.csv("./trait_data/Lizard_Trait_Data.csv", stringsAsFactors = F)

s1.traits$Id <- tolower(s1.traits$Id) #make species Id's lowercase
rownames(s1.traits) <- s1.traits$Id

head(s1.traits)

### impute continuous trait data 

s1.traits.cont <- s1.traits[,c("Id", "maxSVL", "fSVL", "hatchlingSVL", "clutch_size")]

names(s1.traits.cont)[1] <- "species" # rename column from Id to species b/c phylpar requires it


s1.traits.imputed <- phylopars(trait_data = s1.traits.cont, tree = s1.tree)

str(s1.traits.imputed)
s1.traits.imputed$anc_recon

s1.traits.imputed.values <- s1.traits.imputed$anc_recon
s1.traits.imputed.values <- s1.traits.imputed.values[c(1:29),]

# s1.traits.imputed.values isn't a data frame, so needed to manually add imputed values to s1.traits
s1.traits["amphisbaena_scutigerum","fSVL"] = 325.891 
s1.traits["amphisbaena_scutigerum", "hatchlingSVL"] <- 98.704
s1.traits["amphisbaena_scutigerum","clutch_size"] <- 5.507

s1.traits["plestiodon_marginatus","hatchlingSVL"] <- 33.62
s1.traits["plestiodon_marginatus","clutch_size"] <- 4.78

s1.traits["rieppeleon_brevicaudatus","hatchlingSVL"] <- 25.52


### Read in brain data 
s1 <- read.csv("./brain_data/whole_brain_coords_S1.csv", stringsAsFactors = F)

head(s1)

### Remove snakes from brain dataframe 
snakes <- c("Xerotyphlops_vermicularis",
            "Python_regius",
            "Epicrates_cenchria",
            "Eryx_jaculus",
            "Cerastes_cerastes",
            "Hydrophis_platurus",
            "Boaedon_fuliginosus",
            "Chrysopelea_ornata",
            "Dendrelaphis_pictus",
            "Dasypeltis_gansi",
            "Pantherophis_guttatus")

s1.trim <- s1[!((s1$Id) %in% snakes),]

s1.trim$Id <- tolower(s1.trim$Id) # make species Id's lowercase

### order species in alphabetical order

s1.sorted <- s1.trim[order(s1.trim$Id),]

### reassign row numbers based on alphabetical  order
rownames(s1.sorted) <- seq(length=nrow(s1.sorted))
### rename row numbers with the species name
rownames(s1.sorted) <- s1.sorted$Id

### Convert 2D array to 3D array

s1.trim.3D <- as.matrix(s1.sorted[,-(1)])
s1.trim.3D <- arrayspecs(s1.trim.3D, 61, 3)

head(s1.trim.3D)

      ### the 3D array is [a, b, c]
      ### where 
      ### a = landmark id (e.g., landmark 1, landmark 2, landmark 3, etc.) with a max of 61 values
      ### b = the 3D coordinates (e.g., x coord, y coord, z coord) with a max of 3 values
      ### c = species id coded by position in dataframe (e.g., Actonias_meleagris is 1)
      ###  species id       #    , , 1 
      
      ###  3D coordinates   #    [,1]     [,2]     [,3]
      ###  landmark 1       #    [1,] 9104.361 14396.24 7085.975
      ###  landmark 2       #    [2,] 9740.752 13144.82 6818.917


### Create individual dataframes for each major brain region 
### source for subsetting 3D arrays: http://adv-r.had.co.nz/Subsetting.html

### Which landmarks ids (1 through 61) correspond to what major brain region
### *3 and 45 are shared between Diencephalon, Mesencephalon, and Medulla Oblongata
### (even though 3 isn't duplicated between them in their data - weird) 

      # Telencephalon
      # 1:5, 10:14, 19:26
      # Anterior-most extent of the olfactory bulb 1-2
      # Lateral-most extent of the olfactory bulb 3-4
      
      # Diencephalon
      # 6, 15, 43:45, 54:55, 61
      # Optic chiasm - mid-sagittal plane 20
      
      # Mesencephalon
      # 7:8, 16:17, 27:28, 45, 48:49, 52:53, 58:59
      
      # Cerebellum
      # 30:42, 50:51
      
      # Medulla oblongota
      # 9, 18, 29, 45:47, 56:57, 60 

tel <- s1.trim.3D[c(1:5, 10:14, 19:26), 1:3, 1:29]
dien <- s1.trim.3D[c(6, 15, 43:45, 54:55, 61), 1:3, 1:29]
mes <- s1.trim.3D[c(7:8, 16:17, 27:28, 45, 48:49, 52:53, 58:59), 1:3, 1:29]
cere <-s1.trim.3D[c(30:42, 50:51), 1:3, 1:29]
medob <- s1.trim.3D[c(9, 18, 29, 45:47, 56:57, 60), 1:3, 1:29]
# Note: s1.trim.3D is whole brain 



# Have example of all the analysis with the whole brian data
# Then create a function to test same things but on the separate brain regions 
# we shouldn't need to run modularity or integration on the individual brain regions, I don't think 
### Geomorphometric analyses using package 'geomorph'

### Generalized Procrustes Analysis
### how we get shape variables from landmark data (refer to gpagen details)
S1.GPA <- gpagen(s1.trim.3D, ProcD = T, verbose = T)

S1.GPA$procD # Procrustes distance matrix for all specimens 
S1.GPA$points.VCV # variance-covariance matrix among Procrustes shape variables 
plot(S1.GPA) # ? does this plot all the species on top of each other?

head(S1.GPA)

S1.gdf <- geomorph.data.frame(S1.GPA, phy = s1.tree) 

### Calculate phylogenetic signal 
physignal(S1.GPA$coords, s1.tree)
physignal(S1.GPA$Csize, s1.tree)

### Test the integration of the brain structures ('quantifying the degree of morphological integration between modular partitions of shape data')
### landmarks of the brain assigned to brain region partitions
    # tel  = A
    # dien = B (the shared point 45 is classified under B)
    # mes  = C
    # cere = D
    # medob= E

brain.regions <- c('A', 'A', 'A', 'A', 'A',
                   'B',
                   'C', 'C',
                   'E',
                   'A', 'A', 'A', 'A', 'A',
                   'B',
                   'C', 'C', 
                   'E',
                   'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
                   'C', 'C',
                   'E',
                   'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
                   'B', 'B', 'B',
                   'E', 'E', 
                   'C', 'C',
                   'D', 'D',
                   'C', 'C',
                   'B', 'B',
                   'E', 'E',
                   'C', 'C',
                   'E',
                   'B')

Brain.IT <- integration.test(S1.GPA$coords, partition.gp = brain.regions, iter = 999) #two-block partial least squares analysis; correlation between each individual component 
summary(Brain.IT)
plot(Brain.IT)

Brain.IT

Brain.Mod <- modularity.test(S1.GPA$coords, partition.gp = brain.regions, iter = 999)
summary(Brain.Mod)
plot(Brain.Mod)

Brain.Mod$CR.mat # modularity scores 

# refer to adam dean papers on how to interpret the results
    # Adams, D.C. and M.L. Collyer. 2019. 
    # Comparing the strength of modular signal, 
    # and evaluating alternative modular hypotheses, 
    # using covariance ratio effect sizes with morphometric data. 
    # Evolution. 73:2352-2367.


# Phylo-integration and Phylo-modularity 



Brain.phylo.IT <- phylo.integration(S1.gdf$coords, phy = s1.tree, partition.gp = brain.regions, iter = 999)
summary(Brain.phylo.IT)
plot(Brain.phylo.IT)

str(Brain.phylo.IT)

Brain.phylo.Mod <- phylo.modularity(S1.gdf$coords, phy = s1.tree, partition.gp = brain.regions, iter = 999)
summary(Brain.phylo.Mod)
plot(Brain.phylo.Mod)

###### PCA

S1.PCA <- gm.prcomp(S1.GPA$coords, phy = s1.tree, GLS = TRUE)
summary(S1.PCA)
plot(S1.PCA)
plot(S1.PCA, phylo=TRUE)

str(S1.PCA)
pca3d(S1.PCA$x)



Brain.pgls <- procD.pgls(coords ~ Csize, phy=s1.tree, data=S1.gdf)
summary(Brain.pgls)


##### Is brain shape different between foraging modes, reproductive mode, latitude

trait.gdf <- geomorph.data.frame(S1.GPA,
                                 bio_realm = s1.traits$main_biogeographic_realm,
                                 lat = s1.traits$Latitude,
                                 lon = s1.traits$longitude,
                                 ie = s1.traits$insular_endemic,
                                 maxSVL = s1.traits$maxSVL,
                                 fSVL = s1.traits$fSVL,
                                 hSVL = s1.traits$hatchlingSVL,
                                 hm = s1.traits$habitat_modes,
                                 at = s1.traits$activity_time,
                                 substrate = s1.traits$substrate,
                                 microhab = s1.traits$microhabitat,
                                 fm = s1.traits$foraging_mode, 
                                 rm = s1.traits$reproductive_mode,
                                 clutch = s1.traits$clutch_size,
                                 IUCNra = s1.traits$IUCN_redlist_assessment,
                                 IUCNpt = s1.traits$IUCN_population_trend)

#foraging.gdf      <- geomorph.data.frame(S1.GPA, fm = s1.traits$foraging_mode)
#repro.mode.gdf    <- geomorph.data.frame(S1.GPA, rm = s1.traits$reproductive_mode)
#latitude.gdf      <- geomorph.data.frame(S1.GPA, latitude = s1.traits$Latitude)
#biogeography.gdf  <- geomorph.data.frame(S1.GPA, bio_realm = s1.traits$main_biogeographic_realm)
#substrate.gdf    <- geomorph.data.frame(S1.GPA, substrate = s1.traits$substrate)
#activitytime.gdf <- geomorph.data.frame(S1.GPA, activity = s1.traits$activity_time)

fm.pgls <- procD.pgls(S1.GPA$coord ~ S1.GPA$Csize + fm, phy=s1.tree, SS.type="I", data=foraging.gdf)
summary(fm.pgls)


rm.pgls <- procD.pgls(S1.GPA$coord ~ S1.GPA$Csize + rm, phy=s1.tree, SS.type="I", data=repro.mode.gdf)
summary(rm.pgls)

latitude.pgls <- procD.pgls(S1.GPA$coord ~ S1.GPA$Csize + latitude, phy=s1.tree, SS.type="III", data=latitude.gdf, iter=9999)
summary(latitude.pgls)

bio_realm.pgls <- procD.pgls(S1.GPA$coord ~ S1.GPA$Csize + bio_realm, phy=s1.tree, SS.type="III", 
                             data=biogeography.gdf, iter=9999)
summary(bio_realm.pgls)

substrate.pgls <- procD.pgls(S1.GPA$coord ~ S1.GPA$Csize + substrate, phy=s1.tree, SS.type="III", data=substrate.gdf, iter=9999)
summary(substrate.pgls)

activitytime.pgls <- procD.pgls(S1.GPA$coord ~ S1.GPA$Csize + activity, phy=s1.tree, SS.type="III", data=activitytime.gdf, iter=9999)
summary(activitytime.pgls)

#############
# analysis above but with the individual brain regions

cere.GPA <- gpagen(cere, ProcD = T, verbose = T)

cere.GPA$procD # Procrustes distance matrix for all specimens 
cere.GPA$points.VCV # variance-covariance matrix among Procrustes shape variables 
plot(cere.GPA)

cere.trait.gdf <- geomorph.data.frame(cere.GPA,
                                 bio_realm = s1.traits$main_biogeographic_realm,
                                 lat = s1.traits$Latitude,
                                 lon = s1.traits$longitude,
                                 ie = s1.traits$insular_endemic,
                                 maxSVL = s1.traits$maxSVL,
                                 fSVL = s1.traits$fSVL,
                                 hSVL = s1.traits$hatchlingSVL,
                                 hm = s1.traits$habitat_modes,
                                 at = s1.traits$activity_time,
                                 substrate = s1.traits$substrate,
                                 microhab = s1.traits$microhabitat,
                                 fm = s1.traits$foraging_mode, 
                                 rm = s1.traits$reproductive_mode,
                                 clutch = s1.traits$clutch_size,
                                 IUCNra = s1.traits$IUCN_redlist_assessment,
                                 IUCNpt = s1.traits$IUCN_population_trend)

fm.cere.pgls <- procD.pgls(cere.GPA$coord ~ cere.GPA$Csize + fm, phy=s1.tree, SS.type="I", data=cere.trait.gdf)
summary(fm.cere.pgls)

rm.cere.pgls <- procD.pgls(cere.GPA$coord ~ cere.GPA$Csize + rm, phy=s1.tree, SS.type="I", data=cere.trait.gdf)
summary(rm.cere.pgls) # trend 

clutch.cere.pgls <- procD.pgls(cere.GPA$coord ~ cere.GPA$Csize + clutch, phy=s1.tree, SS.type="I", data=cere.trait.gdf)
summary(clutch.cere.pgls)

##################

tel.GPA <- gpagen(tel, ProcD = T, verbose = T)

tel.GPA$procD # Procrustes distance matrix for all specimens 
tel.GPA$points.VCV # variance-covariance matrix among Procrustes shape variables 
plot(tel.GPA)

tel.trait.gdf <- geomorph.data.frame(tel.GPA,
                                      bio_realm = s1.traits$main_biogeographic_realm,
                                      lat = s1.traits$Latitude,
                                      lon = s1.traits$longitude,
                                      ie = s1.traits$insular_endemic,
                                      maxSVL = s1.traits$maxSVL,
                                      fSVL = s1.traits$fSVL,
                                      hSVL = s1.traits$hatchlingSVL,
                                      hm = s1.traits$habitat_modes,
                                      at = s1.traits$activity_time,
                                      substrate = s1.traits$substrate,
                                      microhab = s1.traits$microhabitat,
                                      fm = s1.traits$foraging_mode, 
                                      rm = s1.traits$reproductive_mode,
                                      clutch = s1.traits$clutch_size,
                                      IUCNra = s1.traits$IUCN_redlist_assessment,
                                      IUCNpt = s1.traits$IUCN_population_trend)

fm.tel.pgls <- procD.pgls(tel.GPA$coord ~ tel.GPA$Csize + fm, phy=s1.tree, SS.type="I", data=tel.trait.gdf)
summary(fm.tel.pgls)

rm.tel.pgls <- procD.pgls(tel.GPA$coord ~ tel.GPA$Csize + rm, phy=s1.tree, SS.type="I", data=tel.trait.gdf)
summary(rm.tel.pgls) # trend

clutch.tel.pgls <- procD.pgls(tel.GPA$coord ~ tel.GPA$Csize + clutch, phy=s1.tree, SS.type="I", data=tel.trait.gdf)
summary(clutch.tel.pgls)

# Csize has an effect on teh tel but not the trait itself



#################
dien.GPA <- gpagen(dien, ProcD = T, verbose = T)

dien.GPA$procD # Procrustes distance matrix for all specimens 
dien.GPA$points.VCV # variance-covariance matrix among Procrustes shape variables 
plot(dien.GPA)

dien.trait.gdf <- geomorph.data.frame(dien.GPA,
                                      bio_realm = s1.traits$main_biogeographic_realm,
                                      lat = s1.traits$Latitude,
                                      lon = s1.traits$longitude,
                                      ie = s1.traits$insular_endemic,
                                      maxSVL = s1.traits$maxSVL,
                                      fSVL = s1.traits$fSVL,
                                      hSVL = s1.traits$hatchlingSVL,
                                      hm = s1.traits$habitat_modes,
                                      at = s1.traits$activity_time,
                                      substrate = s1.traits$substrate,
                                      microhab = s1.traits$microhabitat,
                                      fm = s1.traits$foraging_mode, 
                                      rm = s1.traits$reproductive_mode,
                                      clutch = s1.traits$clutch_size,
                                      IUCNra = s1.traits$IUCN_redlist_assessment,
                                      IUCNpt = s1.traits$IUCN_population_trend)

fm.dien.pgls <- procD.pgls(dien.GPA$coord ~ dien.GPA$Csize + fm, phy=s1.tree, SS.type="I", data=dien.trait.gdf)
summary(fm.dien.pgls) # trend

rm.dien.pgls <- procD.pgls(dien.GPA$coord ~ dien.GPA$Csize + rm, phy=s1.tree, SS.type="I", data=dien.trait.gdf)
summary(rm.dien.pgls) # sign effect of rm 

clutch.dien.pgls <- procD.pgls(dien.GPA$coord ~ dien.GPA$Csize + clutch, phy=s1.tree, SS.type="I", data=dien.trait.gdf)
summary(clutch.dien.pgls)


