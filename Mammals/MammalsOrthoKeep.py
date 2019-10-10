#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from itertools import combinations
#import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
import pickle
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import copy
from scipy.stats import sem, t
from scipy import mean
import re
import os

"""
CREATING DATABASE WITHT PRIMATE OTHOLOGS
"""

#SYS.VARIABLES

"""
if len(sys.argv) == 2:
else:
    sys.exit("The command line should be the following I am going to explain here")
"""

pd.options.mode.chained_assignment = None  # default='warn'

#VARIABLES
PATH = "/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/data/"
Mammals = ["Chimpanzee", "Bolivian squirrel monkey", "Gibbon", "Gorilla", "Bonobo", "Orangutan",
"Marmoset", "Mouse Lemur", "Tarsier", "Drill", "Gelada", "Macaque", "Angola colobus", "Capuchin",
"Black snub-nosed monkey", "Bushbaby", "Coquerel's sifaka", "Crab-eating macaque", "Greater bamboo lemur",
"Golden snub-nosed monkey", "Vervet-AGM", "Ugandan red Colobus", "Sooty mangabey", "Pig-tailed macaque",
"Olive baboon", "Ma's night monkey", "Algerian mouse", "Alpaca", "Alpine marmot", "American beaver", "American bison",
"American black bear" , "American mink", "Arctic ground squirrel", "Armadillo",
"Brazilian guinea pig", "Cat", "Chinese hamster CHOK1GS", "Chinese hamster CriGri", "Cow", "Damara mole rat",
"Daurian ground squirrel", "Degu", "Dingo", "Dolphin", "Donkey", "Elephant", "Ferret", "Goat", "Golden Hamster",
"Guinea Pig", "Hedgehog", "Horse", "Hyrax", "Kangaroo rat", "Leopard", "Lesser Egyptian jerboa", "Lesser hedgehog tenrec",
"Long-tailed chinchilla", "Megabat", "Microbat", "Mongolian gerbil", "Mouse Lemur", "Naked mole-rat female",
"Northern American deer mouse", "Opossum", "Panda", "Pig", "Pika", "Platypus", "Prairie vole", "Rabbit", "Rat",
"Ryukyu mouse", "Sheep", "Shrew mouse", "Sloth", "Squirrel", "Steppe mouse", "Tasmanian devil", "Tiger", "Tree Shrew",
"Upper Galilee mountains blind mole rat", "Wallaby"]
Human_genes_id = {}
#MAIN

#Create the dataframe


for filename in os.listdir(PATH):
    print(filename)
    if filename.startswith("Mammal"):
        ortholog_dataset = pd.read_csv(PATH + filename, sep='\t', low_memory=False)
        current_mammals = []
        for species in Mammals:
            for col in ortholog_dataset:
                if col.startswith(species + " "):
                    current_mammals.append(species)
        current_mammals = set(current_mammals)
        print(current_mammals)
        if current_mammals:
            for i, row in ortholog_dataset.iterrows():
                if row["Gene stable ID"] not in Human_genes_id:
                    Human_genes_id[row["Gene stable ID"]] = {}
                for organism2 in current_mammals:
                    if row[organism2 + " homology type"] == "ortholog_one2one":
                        Human_genes_id[row["Gene stable ID"]][organism2] = row[organism2 + " gene stable ID"]
        else:
            continue


#PRINT THE DATASET

Mammals_Orthos = pd.DataFrame(columns=["Human Gene ID", "Species", "Species Gene ID"])
for key1 in Human_genes_id:
    if len(Human_genes_id[key1]) >= 0.7*len(Mammals):
        for key2 in Human_genes_id[key1]:
            Mammals_Orthos = Mammals_Orthos.append({"Human Gene ID": key1, "Species": key2, "Species Gene ID": Human_genes_id[key1][key2]}, ignore_index=True)
Mammals_Orthos.to_csv(sys.argv[1], sep="\t")
