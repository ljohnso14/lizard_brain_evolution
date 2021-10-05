### Brain morphology analysis 

### load appropriate packages 
library(tidyverse)
library(geomorph)
library(ggplot2)


### Table of Contents
### 1) read data
### 2) remove snakes from dataframe, leaving only lizards
### 3) convert 2D array to 3D array
### 4) create individual dataframes for each major brain region 


### Read in brain data 
s1 <- read.csv("./brain_data/whole_brain_coords_S1.csv", stringsAsFactors = F)

### Remove snakes from dataframe 
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


### Convert 2D array to 3D array

s1.trim.3D <- as.matrix(s1.trim[,-(1)])
s1.trim.3D <- arrayspecs(s1.trim.3D, 61, 3)

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



### Geomorphometric analyses using package 'geomorph'

### Generalized Procrustes Analysis - how we get shape variables from landmark data (refer to gpagen details)
S1.GPA <- gpagen(s1.trim.3D)
plot(S1.GPA) # ? does this plot all the species on top of each other?

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

Brain.IT <- integration.test(S1.GPA$coords, partition.gp = brain.regions, iter = 999)
summary(Brain.IT)
picknplot.shape(plot(Brain.IT))

Brain.Mod <- modularity.test(S1.GPA$coords, partition.gp = brain.regions, iter = 999)
summary(Brain.Mod)
plot(Brain.Mod)

# refer to adam dean papers on how to interpret the results
    # Adams, D.C. and M.L. Collyer. 2019. 
    # Comparing the strength of modular signal, 
    # and evaluating alternative modular hypotheses, 
    # using covariance ratio effect sizes with morphometric data. 
    # Evolution. 73:2352-2367.


