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
from time import sleep
import sqlite3


#VARIABLES

PATH="/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/"
DATASETS_PATH="/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/data/Ensembl_genotypes/"
OUTPATH= PATH + "Mammals/GeneLengths.tsv"
"""
First we open the ortholog list of the species groups we work with
"""
cols = pd.read_csv(PATH + "Mammals/ortho1to1.tsv", sep='\t', nrows=1).columns
Mammals_orthos = pd.read_csv(PATH + "Mammals/ortho1to1.tsv", sep='\t', low_memory=False, usecols=cols[1:])


"""
API FUNCTION TO RETRIEVE SEQUENCE LENGTH
"""
def retrieve_length_of_sequence(ID_seq, database):
    print(database)
    conn = sqlite3.connect(database)
    ensembl_id = ID_seq
    c = conn.cursor()
    for row in c.execute('SELECT GeneID FROM annotation WHERE GeneID=?', (ensembl_id,)):
        print(row)
    conn.close()


"""
Then, for each gene family we loop and print the lengths
"""

#Put human gene lengths
"""
human_lengths = []
ref_id = "None"
for value in Mammals_orthos["Human Gene ID"]:
    if value != ref_id:
        ref_id = value
        retrieve_length_of_sequence(value, DATASETS_PATH+"Primates/Human/human.sqlite")
        #human_lengths.append(seq_len)
    else:
        continue
        #human_lengths.append(seq_len)
"""
#Put rest of species intron lengths
species_lengths = []
for i,row in Mammals_orthos.iterrows():
    if i == 0:
        continue
    else:
        name = row["Species"].replace(" ", "_")
        db_path = DATASETS_PATH+"Primates/"+name+"/"+name.lower()+".sqlite"
        retrieve_length_of_sequence(row["Species Gene ID"], db_path)
    #species_lengths.append(seq_len)


#Mammals_orthos["Human Gene Length"] = human_lengths
#Mammals_orthos["Species Gene Length"] = species_lengths
#Mammals_orthos.to_csv(OUTPATH, sep="\t")
