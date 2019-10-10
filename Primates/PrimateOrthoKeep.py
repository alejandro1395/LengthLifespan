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
Primates = ["Chimpanzee", "Bolivian squirrel monkey", "Gibbon", "Gorilla", "Bonobo", "Orangutan",
"Marmoset", "Mouse Lemur", "Tarsier", "Drill", "Gelada", "Macaque", "Angola colobus", "Capuchin",
"Black snub-nosed monkey", "Bushbaby", "Coquerel's sifaka", "Crab-eating macaque", "Greater bamboo lemur",
"Golden snub-nosed monkey", "Vervet-AGM", "Ugandan red Colobus", "Sooty mangabey", "Pig-tailed macaque",
"Olive baboon", "Ma's night monkey"]

#
Human_genes_id = {}
#MAIN

#Create the dataframe


for filename in os.listdir(PATH):
    print(filename)
    if filename.startswith("Mammal"):
        ortholog_dataset = pd.read_csv(PATH + filename, sep='\t', low_memory=False)
        current_primates = []
        for primate in Primates:
            for col in ortholog_dataset:
                if col.startswith(primate):
                    current_primates.append(primate)
        current_primates = set(current_primates)
        print(current_primates)
        if current_primates:
            for i, row in ortholog_dataset.iterrows():
                if row["Gene stable ID"] not in Human_genes_id:
                    Human_genes_id[row["Gene stable ID"]] = {}
                for primate2 in current_primates:
                    if row[primate2 + " homology type"] == "ortholog_one2one":
                        Human_genes_id[row["Gene stable ID"]][primate2] = row[primate2 + " gene stable ID"]
        else:
            continue


#PRINT THE DATASET

Primates_Orthos = pd.DataFrame(columns=["Human Gene ID", "Species", "Species Gene ID"])
for key1 in Human_genes_id:
    for key2 in Human_genes_id[key1]:
        Primates_Orthos = Primates_Orthos.append({"Human Gene ID": key1, "Species": key2, "Species Gene ID": Human_genes_id[key1][key2]}, ignore_index=True)
Primates_Orthos.to_csv(sys.argv[1], sep="\t")
