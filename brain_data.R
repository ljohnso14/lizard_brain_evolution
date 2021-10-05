# brain data test 




library(tidyverse)
library(geomorph)
library(ggplot2)



getwd()
setwd("C:/Users/laure/OneDrive/Documents/1.WUSTL/LososLab/Lizard_Brains_LEJ_DBM/data/processedforR")

# Tidy the whole brain coords data 
# 1) Remove all lizards 
# 2) Create separate data frames for each major brain region of interest
# 2a) consider also separating on smaller region scales, like the olfactory, occiptal lobe etc. 
# 3) Convert into a 3D array 

# Read in brain data 
s1 <- read.csv("whole_brain_coords_S1.csv", stringsAsFactors = F)

# Read in tree data 
# read.nexus format 
# estimate of branch lengths 
# brain.geomorph.trees.nex

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


# convert to 3D array

s1.trim.3D <- as.matrix(s1.trim[,-(1)])
s1.trim.3D <- arrayspecs(s1.trim.3D, 61, 3)

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

# Generalized Procrustes Analysis 
S1.GPA <- gpagen(s1.trim.3D)
plot(S1.GPA)

# Define modules
brain.modules <- define.modules(s1.trim.3D, 5) # this isn't working

# Test of Integration of Brain structures 
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

# tel = a
# dien = b
# mes = c
# cere = d
# medob = e

Brain.IT <- integration.test(S1.GPA$coords, partition.gp = brain.regions, iter = 999)
summary(Brain.IT)
picknplot.shape(plot(Brain.IT))

Brain.Mod <- modularity.test(S1.GPA$coords, partition.gp = brain.regions, iter = 999)
summary(Brain.Mod)
plot(Brain.Mod)

# refer to adam dean papers on how to interpret the results 


####### 10/5/2021 #######
traits <- read.csv("Lizard_Trait_Data.csv", stringsAsFactors = F)

