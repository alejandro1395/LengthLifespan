##########################################
##########################################
###RScript to make correlations between###
###Data of Primates from AnAge for their###
###lifespan and the computed ortholog ####
####lengths from ensembl#################
###########################################

#IMPORT libraries

library(ape)
library(geiger)
library(nlme)
library(phytools)
library(ggplot2)
library(qqman)
library(qualityTools)
library(Haplin)

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


#Select only Primate Species

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

#ADD BM index of each one of the organisms

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


View(SelectedPrimates)

#######################
###GENETIC DATABASE####
#######################

#Import Lengths from Primate genes
#Lengthdb <- read.csv("/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/Primates/GeneLengths.tsv", header = TRUE, sep = "\t")
Lengthdb <- read.csv("AllIntronContent.tsv", header = TRUE, sep = "\t")
Lengthdb <- Lengthdb[,-1]
print(Lengthdb[1:26,])

#Let's remove rows of species we do not want
Lengthdb <- Lengthdb[!grepl("Ma's night monkey", Lengthdb$Species),]
Lengthdb <- Lengthdb[!grepl("Ugandan red Colobus", Lengthdb$Species),]
print(Lengthdb[1:24,])

#Get tree of primates
PrimatesTree <- read.nexus("~/Downloads/consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")
PrimatesTree$tip.label
#plot(PrimatesTree)

# R function for ordering dataset for PGLS
order_species_subset = function(x, vec) {
  vec <- append(vec, x)
  return(vec)
} 

order_species_len = function(x, vec) {
  vec <- append(vec, as.character(x))
  return(vec)
} 

pvalues <- c()
gene_ids <- c()
#LOOP FOR EACH GENE
for (i in seq(0, nrow(Lengthdb)-1, by=24)){
  #DATABASE CONSTRUCTION
  SubsetLength <- Lengthdb[(i+1):(i+24),]
  SpeciesList <- c()
  GenesID <- c()
  Lengths <- c()
  SpeciesList[1] <- "Human"
  hum_len <- as.vector(SubsetLength$Human.Total.Intron.Length[1])
  species_len <- as.vector(SubsetLength$Species.Total.Intron.Length)
  GenesID[1] <- as.character(SubsetLength$Human.Gene.ID[1])
  Lengths[1] <- SubsetLength$Human.Total.Intron.Length[1]
  PrimatesList <- SubsetLength$Species
  PrimatesGenesID <- SubsetLength$Species.Gene.ID
  PrimatesLengths <- SubsetLength$Species.Total.Intron.Length
  Total_SpeciesList <- c(as.character(SpeciesList), as.character(PrimatesList))
  Total_GenesID <- c(as.character(GenesID), as.character(PrimatesGenesID))
  Total_Lengths <- log10(as.numeric(as.character(c(hum_len, species_len))))
  LengthPrimatesDataframe <- data.frame(Total_SpeciesList, Total_GenesID, Total_Lengths)
  names(LengthPrimatesDataframe) <- c("Total_SpeciesList", "Total_GenesID", "Total_Lengths")
  
  #Loop for changing common names to match and including longevities
  
  ordered_longevities <- c()
  ordered_bodymass <- c()
  for ( i in 1:length(LengthPrimatesDataframe$Total_SpeciesList)){
    ordered_longevities[i]  <- ifelse(LengthPrimatesDataframe$Total_SpeciesList[i] %in% SelectedPrimates$Common.name, 
                                      as.numeric(SelectedPrimates$Maximum.longevity..yrs[SelectedPrimates$Common.name == LengthPrimatesDataframe$Total_SpeciesList[i]]),
                                      NA )
    ordered_bodymass[i] <- ifelse(LengthPrimatesDataframe$Total_SpeciesList[i] %in% SelectedPrimates$Common.name, 
                                  as.numeric(SelectedPrimates$Body.mass..g.[SelectedPrimates$Common.name == LengthPrimatesDataframe$Total_SpeciesList[i]]),
                                  NA )
  }
  LengthPrimatesDataframe$MaxLifespan <- log10(ordered_longevities)
  LengthPrimatesDataframe$LQ <- log10(ordered_longevities/(4.88*(ordered_bodymass)^0.153))
  
  row.names(LengthPrimatesDataframe) <- c("Homo_sapiens", "Macaca_mulatta", "Tarsius_syrichta", "Otolemur_garnettii",
                                          "Saimiri_boliviensis", "Pan_troglodytes_troglodytes", "Nomascus_leucogenys", 
                                          "Microcebus_murinus", "Pongo_abelii", "Pan_paniscus", "Callithrix_jacchus",
                                          "Gorilla_gorilla_gorilla", "Propithecus_coquereli", "Macaca_fascicularis", 
                                          "Rhinopithecus_roxellana", "Mandrillus_leucophaeus", "Hapalemur_simus",
                                          "Theropithecus_gelada", "Cebus_capucinus", "Colobus_angolensis",
                                          "Rhinopithecus_bieti", "Chlorocebus_pygerythrus", "Papio_anubis",
                                          "Cercocebus_torquatus_atys", "Macaca_nemestrina")
  
  print(as.character(SubsetLength$Human.Gene.ID[1]))
  #PGLS
  tree.coorel <-corBrownian(phy=PrimatesTree)
  if (class(try(gls(as.numeric(LQ) ~ as.numeric(Total_Lengths), 
                    data = LengthPrimatesDataframe, correlation = tree.coorel,
                    method = "ML"))) == "try-error"){
    print("Singular data")
    
  }
  else{
    pglsModel <- gls(as.numeric(LQ) ~ as.numeric(Total_Lengths), 
                     data = LengthPrimatesDataframe, correlation = tree.coorel,
                     method = "ML")
    pv <- anova(pglsModel)[2,3]
    pvalues<- append(pvalues, pv)
    gene_ids <- append(gene_ids, as.character(SubsetLength$Human.Gene.ID[1]))
  }
}

#Correction according to multiple testing of genes
print(length(pvalues))
print(length(gene_ids))

names(pvalues) <- gene_ids
corrected_pvalues <- p.adjust(pvalues, method = "fdr", n=length(pvalues))
sorted_corrected_pvalues <- sort(corrected_pvalues , decreasing = FALSE)

best <- sorted_corrected_pvalues[sorted_corrected_pvalues <= 0.05]
length(best)

gene_length_total_candidates <- c("ENSG00000213889", "ENSG00000108797", "ENSG00000116062", "ENSG00000138829",
                                  "ENSG00000146842", "ENSG00000125434", "ENSG00000119698", "ENSG00000172663",
                                  "ENSG00000269313", "ENSG00000158113", "ENSG00000159216", "ENSG00000139546", 
                                  "ENSG00000127824", "ENSG00000147606",  "ENSG00000143951","ENSG00000162419", 
                                  "ENSG00000155052", "ENSG00000119537", "ENSG00000048405", "ENSG00000164134", 
                                  "ENSG00000109339",  "ENSG00000076706", "ENSG00000039537", "ENSG00000117569", 
                                  "ENSG00000122034", "ENSG00000138964", "ENSG00000212443", "ENSG00000131473",  
                                  "ENSG00000196476", "ENSG00000136048")
intron_length_candidates <- names(best)

intron_length_candidates[intron_length_candidates %not in% gene_length_total_candidates]
setdiff(intron_length_candidates, gene_length_total_candidates)
#qqplot and lambda values
qq(pvalues, main = "Q-Q plot of correlation p-values", 
   pch = 18, col = "blue4", cex = 1.5, las = 1) #NOMINAL pvalues

lambda1 <-median(qchisq(1-pvalues,1))/qchisq(0.5,1)

pvalues_lambda_corrected <- pvalues*lambda1

pQQ(pvalues_lambda_corrected, conf = 0.95)

