###########################################################
###########################################################
#########RScript to make correlations between##############
########Data of Primates from AnAge for their##############
##########lifespan and the computed ortholog###############
#######lengths from ensembl - Weighted by Quality##########
###########################################################

#IMPORT libraries

library(ape)
library(geiger)
library(nlme)
library(phytools)
library(ggplot2)
library(qqman)
library(qualityTools)
library(Haplin)
library(ggpubr)
library(gridExtra)
library(grid)


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



#######################
###GENETIC DATABASE####
#######################

#Import Lengths from Primate genes
#Lengthdb <- read.csv("/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/Primates/GeneLengths.tsv", header = TRUE, sep = "\t")
Lengthdb <- read.csv("GeneLengthsQual.tsv", header = TRUE, sep = "\t")
Lengthdb <- Lengthdb[,-1]
print(Lengthdb[1:26,])

#Let's remove rows of species we do not want
Lengthdb <- Lengthdb[!grepl("Ma's night monkey", Lengthdb$Species),]
Lengthdb <- Lengthdb[!grepl("Ugandan red Colobus", Lengthdb$Species),]
print(Lengthdb[1:24,])

#Get tree of primates
PrimatesTree <- read.nexus("../../../Downloads/consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")
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
  hum_len <- as.vector(SubsetLength$Human.Gene.Length[1])
  hum_qual <- 1
  species_len <- as.vector(SubsetLength$Species.Gene.Length)
  GenesID[1] <- as.character(SubsetLength$Human.Gene.ID[1])
  Lengths[1] <- SubsetLength$Human.Gene.Length[1]
  PrimatesList <- SubsetLength$Species
  PrimatesGenesID <- SubsetLength$Species.Gene.ID
  PrimatesLengths <- SubsetLength$Species.Gene.Length
  PrimatesQualities <- SubsetLength$Species.Quality
  Total_SpeciesList <- c(as.character(SpeciesList), as.character(PrimatesList))
  Total_GenesID <- c(as.character(GenesID), as.character(PrimatesGenesID))
  Total_Lengths <- log10(as.numeric(as.character(c(hum_len, species_len))))
  Total_Qualities <- as.numeric(as.character(c(hum_qual, PrimatesQualities)))
  LengthPrimatesDataframe <- data.frame(Total_SpeciesList, Total_GenesID, Total_Lengths, Total_Qualities)
  names(LengthPrimatesDataframe) <- c("Total_SpeciesList", "Total_GenesID", "Total_Lengths", "Total_Qualities")
  
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
                    method = "ML", weights = ~as.numeric(as.character(Total_Qualities))))) == "try-error"){
    print("Singular data")
    
  }
  else{
    pglsModel <- gls(as.numeric(LQ) ~ as.numeric(Total_Lengths), 
                     data = LengthPrimatesDataframe, correlation = tree.coorel,
                     weights = ~as.numeric(as.character(Total_Qualities)),
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
corrected_pvalues <- p.adjust(pvalues, method = "bonf", n=length(pvalues))
sorted_corrected_pvalues <- sort(corrected_pvalues , decreasing = FALSE)

best <- sorted_corrected_pvalues[sorted_corrected_pvalues <= 0.05]
length(best)

#qqplot and lambda values
qq(pvalues, main = "Q-Q plot of correlation p-values", 
   pch = 18, col = "blue4", cex = 1.5, las = 1) #NOMINAL pvalues

lambda1 <-median(qchisq(1-pvalues,1))/qchisq(0.5,1)

pvalues_lambda_corrected <- pvalues*lambda1

pQQ(pvalues, conf = 0.95)


#qqplot2
lambda1 = median(qchisq(1-pvalues,1))/qchisq(0.5,1)

pvalues_lambda_corrected <- pvalues*lambda1
qqplot(-log10(ppoints(length(pvalues_lambda_corrected))), -log10(pvalues))
abline(0,1)

#PLOT FOR BEST CANDIDATES
p <- list()
num <- 0
pvalues <- c()
gene_ids <- c()
candidates <- names(sorted_corrected_pvalues)[1:30]
par(mfcol=c(5,6), oma=c(1,1,0,0), mar=c(1,1,1,0),
    
    tcl=-0.1, mgp=c(0,0,0))
for (value in candidates){
  num<- num +1
  for (i in seq(0, nrow(Lengthdb)-1, by=24)){
    #DATABASE CONSTRUCTION
    SubsetLength <- Lengthdb[(i+1):(i+24),]
    SpeciesList <- c()
    GenesID <- c()
    Lengths <- c()
    SpeciesList[1] <- "Human"
    hum_len <- as.vector(SubsetLength$Human.Gene.Length[1])
    hum_qual <- 1
    species_len <- as.vector(SubsetLength$Species.Gene.Length)
    GenesID[1] <- as.character(SubsetLength$Human.Gene.ID[1])
    if (GenesID[1] == value){
    Lengths[1] <- SubsetLength$Human.Gene.Length[1]
    PrimatesList <- SubsetLength$Species
    PrimatesGenesID <- SubsetLength$Species.Gene.ID
    PrimatesLengths <- SubsetLength$Species.Gene.Length
    PrimatesQualities <- SubsetLength$Species.Quality
    Total_SpeciesList <- c(as.character(SpeciesList), as.character(PrimatesList))
    Total_GenesID <- c(as.character(GenesID), as.character(PrimatesGenesID))
    Total_Lengths <- log10(as.numeric(as.character(c(hum_len, species_len))))
    Total_Qualities <- as.numeric(as.character(c(hum_qual, PrimatesQualities)))
    LengthPrimatesDataframe <- data.frame(Total_SpeciesList, Total_GenesID, Total_Lengths, Total_Qualities)
    names(LengthPrimatesDataframe) <- c("Total_SpeciesList", "Total_GenesID", "Total_Lengths", "Total_Qualities")
    
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
                      method = "ML", weights = ~as.numeric(as.character(Total_Qualities))))) == "try-error"){
      print("Singular data")
      
    }
    else{
      pglsModel <- gls(as.numeric(LQ) ~ as.numeric(Total_Lengths), 
                       data = LengthPrimatesDataframe, correlation = tree.coorel,
                       weights = ~as.numeric(as.character(Total_Qualities)),
                       method = "ML")
      pv <- anova(pglsModel)[2,3]
      pvalues<- append(pvalues, pv)
      gene_ids <- append(gene_ids, as.character(SubsetLength$Human.Gene.ID[1]))
    }
      candidate_gene <- LengthPrimatesDataframe}}
  
  #Correction according to multiple testing of genes
  length_gen <- candidate_gene[,"Total_Lengths"]
  max_lifesp <- candidate_gene[,"LQ"]
  
  names(max_lifesp) <- names(length_gen) <- rownames(candidate_gene)
  
  hPic <- pic(length_gen, PrimatesTree)
  aPic <- pic(max_lifesp, PrimatesTree)
  
  # Make a model
  picModel <- lm(hPic ~ aPic - 1)
  #plot(hPic ~ aPic, main = value)
  #abline(a = 0, b = coef(picModel))
  
  grob <- grobTree(textGrob(paste("R = ", round(cor(hPic, aPic),2)), x=0.1,  y=0.95, hjust=0,
                            gp=gpar(col="blue", fontsize=13, fontface="italic")))
  
  
  p[[num]] <- ggplot(as.data.frame(cbind(hPic, aPic)), aes(x = hPic, y = aPic)) + 
    geom_point() +  geom_smooth(method='lm',formula=y~x) + annotation_custom(grob) + xlab(value)
  
}

do.call(grid.arrange,p)
