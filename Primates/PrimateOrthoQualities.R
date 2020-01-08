##########################################
########DESCRIPTIVE PLOTS FROM ORTHOLOG###
############QUALITIES#####################
##########################################

library(phytools)
library(ggplot2)
library(ggnewscale)
library(ggtree)
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

GeneLengths <- read.csv("GeneLengths.tsv",header=TRUE,sep="\t")
GeneLengths[5, ]

Primate_names <- c("Macaca_mulatta", "Tarsius_syrichta", "Otolemur_garnettii",
                   "Saimiri_boliviensis", "Pan_troglodytes_troglodytes", "Nomascus_leucogenys", 
                   "Microcebus_murinus", "Pongo_abelii", "Pan_paniscus", "Callithrix_jacchus",
                   "Gorilla_gorilla_gorilla", "Propithecus_coquereli", "Macaca_fascicularis", 
                   "Rhinopithecus_roxellana", "Mandrillus_leucophaeus", "Hapalemur_simus",
                   "Theropithecus_gelada", "Cebus_capucinus", "Colobus_angolensis",
                   "Rhinopithecus_bieti", "Chlorocebus_pygerythrus", "Papio_anubis",
                   "Aotus_nancymaae", "Cercocebus_torquatus_atys", "Macaca_nemestrina", 
                   "Piliocolobus_tephrosceles")
PrimateSumLengths<- rep(0, 26)
names(PrimateSumLengths) <- Primate_names

for (i in seq(0, nrow(GeneLengths)-1, by=26)){
  SubsetLength <- GeneLengths[(i+1):(i+26),]
  PrimatesList <- SubsetLength$Species
  PrimatesGenesID <- SubsetLength$Species.Gene.ID
  PrimatesLengths <- SubsetLength$Species.Gene.Length
  LengthPrimatesDataframe <- data.frame(PrimatesList, PrimatesGenesID, PrimatesLengths)
  row.names(LengthPrimatesDataframe) <- c("Macaca_mulatta", "Tarsius_syrichta", "Otolemur_garnettii",
                                          "Saimiri_boliviensis", "Pan_troglodytes_troglodytes", "Nomascus_leucogenys", 
                                          "Microcebus_murinus", "Pongo_abelii", "Pan_paniscus", "Callithrix_jacchus",
                                          "Gorilla_gorilla_gorilla", "Propithecus_coquereli", "Macaca_fascicularis", 
                                          "Rhinopithecus_roxellana", "Mandrillus_leucophaeus", "Hapalemur_simus",
                                          "Theropithecus_gelada", "Cebus_capucinus", "Colobus_angolensis",
                                          "Rhinopithecus_bieti", "Chlorocebus_pygerythrus", "Papio_anubis",
                                          "Aotus_nancymaae", "Cercocebus_torquatus_atys", "Macaca_nemestrina", 
                                          "Piliocolobus_tephrosceles")
  PrimateSumLengths <- PrimateSumLengths + as.numeric(as.character(PrimatesLengths))
}

#Now let's print the whole length for each species
SumLengths_matrix <- as.matrix(PrimateSumLengths)[,1]

obj2<-contMap(PrimatesTree,SumLengths_matrix, direction="rightwards", legend = FALSE)
obj2<-setMap(obj2,colors=c("darkred","yellow","darkblue"))
plot(obj2)




















####################################################################
######################################################################
######################################################################
#Now I want to look whether there is correlation between LQ and and 
#all the average gene length of orthologs correcting by phylogeny


###############
###############
#PHENOTIPIC DB#
###############

#Import Phenotypes from AnAge
Anagedb <- read.csv("anage_data.txt", header = TRUE, sep = "\t")


#Vector of Species wanted
species <- c("Chimpanzee", "Gorilla", "Orangutan", "Greater bamboo lemur", "Golden snub-nosed monkey", "Drill",
             "Bolivian squirrel monkey", "White-tufted-ear marmoset", "Gray mouse lemur", "Angolan colobus",
             "Rhesus monkey", "White-cheeked gibbon", "Human", "Pygmy chimpanzee or bonobo",
             "Philippine tarsier", "Gelada baboon", "Vervet", "Long-tailed macaque",
             "White-faced capuchin", "Small-eared galago", "Pigtail macaque", "Coquerel's sifaka", "Olive baboon",
             "Shooty mangabey")


#Select only Primate Species and include ADW info

AnAge_primates <- Anagedb[which(Anagedb$Order == "Primates"),]
SelectedPrimates <- AnAge_primates[which(AnAge_primates$Common.name %in% species),]
SelectedPrimates$Common.name <- as.character(SelectedPrimates$Common.name)
SelectedPrimates$Species <- as.character(SelectedPrimates$Species)
SelectedPrimates[nrow(SelectedPrimates) + 1,] = c(NA, "Animalia", "Chordata", "Mammalia", "Primates", "Indriidae",
                                                  "Propithecus", "coquereli", "Coquerel's sifaka", NA, NA, NA, NA, NA,
                                                  NA, NA, NA, NA, NA, NA, 30, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
SelectedPrimates[nrow(SelectedPrimates) + 1,] = c(NA, "Animalia", "Chordata", "Mammalia", "Primates", "Cercopithecidae",
                                                  "Papio", "anubis", "Olive baboon", NA, NA, NA, NA, NA,
                                                  NA, NA, NA, NA, NA, NA, 37.5, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
SelectedPrimates[nrow(SelectedPrimates) + 1,] = c(NA, "Animalia", "Chordata", "Mammalia", "Primates", "Cercopithecidae",
                                                  "Cercocebus", "atys", "Sooty Mangabey", NA, NA, NA, NA, NA,
                                                  NA, NA, NA, NA, NA, NA, 18, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
SelectedPrimates[nrow(SelectedPrimates) + 1,] = c(NA, "Animalia", "Chordata", "Mammalia", "Primates", "Cercopithecidae",
                                                  "Rhinopithecus", "bieti", "Black snub-nosed monkey", NA, NA, NA, NA, NA,
                                                  NA, NA, NA, NA, NA, NA, 14.7, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
View(SelectedPrimates)

#hist(as.numeric(SelectedPrimates$Maximum.longevity..yrs.), main="Dsitribution of max lifespan",
col=rainbow(length(SelectedPrimates$Species)), xlab="Ages")

#REPLACE WRONG entries matchin common species names
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Rhesus monkey"] <- "Macaque"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Philippine tarsier"] <- "Tarsier"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Small-eared galago"] <- "Bushbaby"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "White-cheeked gibbon"] <- "Gibbon"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "White-tufted-ear marmoset"] <- "Marmoset"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Long-tailed macaque"] <- "Crab-eating macaque"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Gelada baboon"] <- "Gelada"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "White-faced capuchin"] <- "Capuchin"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Angolan colobus"] <- "Angola colobus"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Sooty Mangabey"] <- "Sooty mangabey"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Pigtail macaque"] <- "Pig-tailed macaque"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Pygmy chimpanzee or bonobo"] <- "Bonobo"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Gray mouse lemur"] <- "Mouse Lemur"
SelectedPrimates$Common.name[SelectedPrimates$Common.name == "Vervet"] <- "Vervet-AGM"


#######################
#######################
#BODY MASS INDEX INFO##
#######################

SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Macaque"] <- 8000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Gibbon"] <- 5700
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Crab-eating macaque"] <- 5000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Gelada"] <- 17000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Capuchin"] <- 3000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Angola colobus"] <- 8900
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Sooty mangabey"] <- 9500
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Pig-tailed macaque"] <- 9600
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Bonobo"] <- 44000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Mouse Lemur"] <- 60
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Vervet-AGM"] <- 4000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Bolivian squirrel monkey"] <- 750
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Drill"] <- 18250
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Golden snub-nosed monkey"] <- 13000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Gorilla"] <- 275000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Chimpanzee"] <- 48000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Orangutan"] <- 60000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Greater bamboo lemur"] <- 2250
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Coquerel's sifaka"] <- 4000
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Olive baboon"] <- 19500
SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == "Black snub-nosed monkey"] <- 10000


#######################
#PGLS OVERALL ANALYSIS#
#######################

Species_names <- c("Callithrix_jacchus", "Cebus_capucinus",
                   "Saimiri_boliviensis", "Chlorocebus_pygerythrus", 
                   "Colobus_angolensis", "Macaca_fascicularis",
                   "Macaca_mulatta", "Macaca_nemestrina",
                   "Mandrillus_leucophaeus", "Rhinopithecus_roxellana",
                   "Theropithecus_gelada", "Microcebus_murinus", "Otolemur_garnettii",
                   "Gorilla_gorilla_gorilla", "Pan_paniscus", 
                   "Pan_troglodytes_troglodytes", "Pongo_abelii",
                   "Nomascus_leucogenys", "Hapalemur_simus",
                   "Tarsius_syrichta", "Propithecus_coquereli",
                   "Papio_anubis", "Cercocebus_torquatus_atys", "Rhinopithecus_bieti")
Species_LQs <- log10(as.numeric(as.character(SelectedPrimates$Maximum.longevity..yrs.))/
                       (4.88*(as.numeric(as.character(SelectedPrimates$Body.mass..g.)))^0.153))
Species_LQs <- Species_LQs[-15]

#Only those ones with lifespan information
PrimatesLifespanTree <- read.nexus("../../../Downloads/consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")
PrimatesLifespanTree <- drop.tip(PrimatesLifespanTree, "Homo_sapiens")

#Keep only primates with lengths summed and lifespan info too
SumLengths_matrix <- SumLengths_matrix[names(SumLengths_matrix) != "Aotus_nancymaae"]
SumLengths_matrix <- SumLengths_matrix[names(SumLengths_matrix) != "Piliocolobus_tephrosceles"]
SumLengths_matrix <- SumLengths_matrix[order(factor(names(SumLengths_matrix), levels = Species_names))]

plotDataFrame <- as.data.frame(cbind(SumLengths_matrix, Species_LQs))
rownames(plotDataFrame) <- Species_names

circ <- ggtree(PrimatesLifespanTree, layout = "circular")

pgp1 <- gheatmap(circ, plotDataFrame[, "SumLengths_matrix", drop=F], width=2, offset = 1,
                 colnames_angle = 90, colnames_offset_y = .25, font.size = 0.00001) + 
  scale_fill_viridis_c(option="E", name = "Lengths", begin = 1, end = 0) 


# Add gheatmap 2 for assay 2
pgp2 <- pgp1 + new_scale_fill()

pgp_all <- gheatmap(pgp2, plotDataFrame[, "Species_LQs", drop=F], offset = 1.5, width = 1.1,
                    colnames_angle = 90, colnames_offset_y = .05, high = "darkgreen",
                    low = "darkred", legend_title = "LQ", font.size = 0.00001) + geom_tiplab(size=2)

