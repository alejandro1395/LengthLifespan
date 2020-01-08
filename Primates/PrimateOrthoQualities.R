##########################################
########DESCRIPTIVE PLOTS FROM ORTHOLOG###
############QUALITIES#####################
##########################################

library(phytools)
packageVersion("phytools")

#READ THE SPECIES TREE
PrimatesTree <- read.nexus("../../../Downloads/consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues_26.nex")

#READ THE TRAIT VARIABLE
trait_matrix <- as.matrix(read.csv("qualitiesless0.5.tsv",header=F,sep=" ",row.names=1))[,1]
names(trait_matrix) <- c("Pongo_abelii", "Hapalemur_simus", "Cebus_capucinus", "Macaca_fascicularis",
                            "Propithecus_coquereli", "Cercocebus_torquatus_atys", "Pan_troglodytes_troglodytes",
                            "Otolemur_garnettii", "Aotus_nancymaae", "Callithrix_jacchus", "Nomascus_leucogenys",
                            "Saimiri_boliviensis", "Tarsius_syrichta", "Gorilla_gorilla_gorilla", "Colobus_angolensis",
                            "Chlorocebus_pygerythrus", "Microcebus_murinus", "Rhinopithecus_roxellana",
                            "Piliocolobus_tephrosceles", "Macaca_nemestrina", "Papio_anubis", "Mandrillus_leucophaeus",
                            "Rhinopithecus_bieti", "Pan_paniscus", "Macaca_mulatta", "Theropithecus_gelada")
trait_matrix <- log(trait_matrix)

#PLOTS
obj<-contMap(PrimatesTree,-trait_matrix, direction="rightwards", legend = FALSE)
obj<-setMap(obj,colors=c("darkred","yellow","darkgreen"))
plot(obj)

obj<-contMap(PrimatesTree,-trait_matrix, direction="rightwards")
plot(obj,lwd=7,xlim=c(-0.2,3.6))
errorbar.contMap(obj)


#Now let's study the overall length for each one of the species

