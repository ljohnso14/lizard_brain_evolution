### Read in Packages
###### Preparation


pacman::p_load(tidyverse, here, geomorph, ggplot2, cowplot, ape, phytools, geiger, tidytree, picante, pca3d)

# Tidy the whole brain coords data 
# 1) Remove all lizards 
# 2) Create separate data frames for each major brain region of interest
# 2a) consider also separating on smaller region scales, like the olfactory, occiptal lobe etc. 
# 3) Convert into a 3D array 

### Read in the data
s1 <- read.csv(here::here("Dropbox", "LaurenBrainAnalysis", "Data", "whole_brain_coords_no_snakes_S1.csv"))
head(s1)

row.names(s1) <- s1$Id


### Read in the tree

s1.tree <- read.nexus(here::here("Dropbox", "LaurenBrainAnalysis", "Data", "BrainGeomorph.nex"))
plot(s1.tree)


#### Read in the Trait data

s1.traits <- read.csv(here::here("Dropbox", "LaurenBrainAnalysis", "Data", "Lizard_Trait_data.csv"))




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

lizards <- (c("Acontias_meleagris",
              "Anguis_fragilis",
              "Basiliscus_vittatus",
              "Blanus_cinereus", 
              "Bradypodion_pumilum",
              "Chalcides_chalcides",
              "Bachia_flavescens",
              "Trioceros_jacksonii",
              "Dasia_olivacea",
              "Draco_volans",
              "Eublepharis_macularius",
              "Plestiodon_marginatus",
              "Gekko_gecko",
              "Amphisbaena_scutigerum",
              "Lygodactylus_picturatus",
              "Melanoseps_loveridgei",
              "Ophiodes_fragilis",
              "Pseudopus_apodus",
              "Lepidothyris_fernandi",
              "Takydromus_sexlineatus",
              "Teratoscincus_scincus",
              "Tropidurus_torquatus",
              "Ablepharus_kitaibelii",
              "Agama_agama",
              "Chalcides_sepsoides",
              "Hemiergis_quadrilineata",
              "Phelsuma_grandis",
              "Pogona_vitticeps",
              "Rieppeleon_brevicaudatus"))

# removing snake species
s1.trim <- s1[!((s1$Id) %in% snakes),]

s1.sorted <- s1.trim[order(s1.trim$Id),]

row.names(s1.sorted) <- seq(length=nrow(s1.sorted)) 
View(s1.trim)
# convert to 3D array
s1.trim.3D <- as.matrix(s1[,-(1)])
s1.trim.3D <- arrayspecs(s1.trim.3D, 61, 3)
head(s1.trim.3D)
# now that i have the 3D array with the 61 landmarks, let's try to then separate those by region!
# hope this works better and reduces the chances of data getting deleted
# http://adv-r.had.co.nz/Subsetting.html

##########################################################
# TESTING SUBSETTING THE ARRAY #
s1.trim.3D[1:3, 1:3, 1:29]
s1.trim.3D[1:2, 1:3, 1:2]
# 3d array [landmark, x/y/z coords, species]
# species          #    , , 1 

# x/y/z coord      #    [,1]     [,2]     [,3]
# landmark 1       #    [1,] 9104.361 14396.24 7085.975
# landmark 2       #    [2,] 9740.752 13144.82 6818.917
##########################################################


# Telencephalon
# 1:5, 10:14, 19:26

# Diencephalon
# 6, 15, 43:45, 54:55, 61

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

##### Generalized Procrustes Analysis

S1.GPA <- gpagen(s1.trim.3D, ProcD=TRUE, verbose = TRUE)

S1.GPA$procD
S1.GPA$points.VCV
plot(S1.GPA)

head(S1.GPA)

S1.gdf <- geomorph.data.frame(S1.GPA, phy=s1.tree)

##### Calculate Phylogenetic Signal

physignal(S1.GPA$coord, s1.tree)
physignal(S1.GPA$Csize, s1.tree)




#### Test of Integration of Brain Structures

brain.regions <- c("A","A","A","A","A",
                   "B","C","C","E","A","A","A","A","A",
                   "B","C","C","E","A","A","A","A","A",
                   "A","A","A","C","C","E","D","D","D","D","D","D","D","D","D",
                   "D","D","D","D","B","B",
                   "B","E","E","C","C","D","D","C","C","B","B","E","E","C",
                   "C","E","B")
Brain.IT <- integration.test(S1.GPA$coords, partition.gp = brain.regions, iter = 999)
summary(Brain.IT)
plot(Brain.IT)

Brain.IT

Brain.Mod <- modularity.test(S1.GPA$coords, partition.gp = brain.regions, iter = 999)
summary(Brain.Mod)
plot(Brain.Mod)

Brain.Mod$CR.mat

#### Phylo-integration and Phylo-modularity


Brain.phylo.IT <- phylo.integration(S1.gdf$coords, phy=s1.tree, partition.gp = brain.regions, iter = 999)
summary(Brain.phylo.IT)
plot(Brain.phylo.IT)

str(Brain.phylo.IT)

Brain.phylo.Mod <- phylo.modularity(S1.GPA$coords, phy=s1.tree, partition.gp = brain.regions, iter = 999)
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

foraging.gdf      <- geomorph.data.frame(S1.GPA, fm = s1.traits$foraging_mode)
repro.mode.gdf    <- geomorph.data.frame(S1.GPA, rm = s1.traits$repro_mode)
latitude.gdf      <- geomorph.data.frame(S1.GPA, latitude = s1.traits$latitude)
biogeography.gdf  <- geomorph.data.frame(S1.GPA, bio_realm = s1.traits$main_biogeographic_realm)
substrate.gdf    <- geomorph.data.frame(S1.GPA, substrate = s1.traits$substrate)


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

str(substrate.pgls)
