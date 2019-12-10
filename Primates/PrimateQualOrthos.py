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
PATH1 = "/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/data/"
PATH2 = "/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/"

cols = pd.read_csv(PATH2 + "Primates/GeneLengths.tsv", sep='\t', nrows=1).columns
Primate_ortho_lengths = pd.read_csv(PATH2 + "Primates/GeneLengths.tsv", sep='\t', low_memory=False, usecols=cols[1:])
Primates = ["Chimpanzee", "Bolivian squirrel monkey", "Gibbon", "Gorilla", "Bonobo", "Orangutan",
"Marmoset", "Mouse Lemur", "Tarsier", "Drill", "Gelada", "Macaque", "Angola colobus", "Capuchin",
"Black snub-nosed monkey", "Bushbaby", "Coquerel's sifaka", "Crab-eating macaque", "Greater bamboo lemur",
"Golden snub-nosed monkey", "Vervet-AGM", "Ugandan red Colobus", "Sooty mangabey", "Pig-tailed macaque",
"Olive baboon", "Ma's night monkey"]

#
Human_genes_id = {}
Primate_ortho_lengths["Species Quality"] = np.nan
matching_species_ids = Primate_ortho_lengths['Species Gene ID'].tolist()
#MAIN

#Create the dataframe with qualities

for filename in os.listdir(PATH1):
    if filename.startswith("primates"):
        quality_dataset = pd.read_csv(PATH1 + filename, sep='\t', low_memory=False)
        current_primates = []
        for primate in Primates:
            for col in quality_dataset:
                if col.startswith(primate):
                    current_primates.append(primate)
        current_primates = set(current_primates)
        if current_primates:
            for i, row in quality_dataset.iterrows():
                for primate2 in current_primates:
                    qual_parm1 = str(row[primate2 + " Gene-order conservation score"])
                    qual_parm2 = str(row[primate2 + " Whole-genome alignment coverage"])
                    if qual_parm1 != "nan" and qual_parm2 != "nan":
                        score = (float(qual_parm1)+float(qual_parm2))/200
                    elif qual_parm1 != "nan":
                        score = float(qual_parm1)/100
                    elif qual_parm2 != "nan":
                        score = float(qual_parm2)/100
                    else:
                        score = 0
                    print(score)
                    if str(row[primate2 + " gene stable ID"]) in matching_species_ids:
                        Primate_ortho_lengths.loc[Primate_ortho_lengths['Species Gene ID'] == row[primate2 + " gene stable ID"], 'Species Quality'] = score

#PRINT THE DATASET
Primate_ortho_lengths.to_csv(sys.argv[1], sep="\t")
