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


#VARIABLES

PATH="/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/"
OUTPATH= PATH + "Primates/FirsIntronLengths.tsv"
"""
First we open the ortholog list of the species groups we work with
"""
cols = pd.read_csv(PATH + "Primates/ortho1to1.tsv", sep='\t', nrows=1).columns
Primate_orthos = pd.read_csv(PATH + "Primates/ortho1to1.tsv", sep='\t', low_memory=False, usecols=cols[1:])


"""
API FUNCTION TO RETRIEVE SEQUENCE WITH INTRONS
"""
def retrieve_length_of_sequence(ID_seq):
    ensembl_id = ID_seq
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/" + ensembl_id + "?" + ";mask_feature=1"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    while not r.ok:
        sleep(5)
        r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    if r.ok:
        print("ok")
    sequence = r.text
    list = re.findall('[a-z]+', sequence)
    if list:
        intron_length = len(list[0])
    else:
        intron_length = "NULL"
    return(intron_length)


"""
Then, for each gene family we loop and print the lengths
"""

#Put human gene lengths
human_lengths = []
ref_id = "None"
for value in Primate_orthos["Human Gene ID"]:
    if value != ref_id:
        ref_id = value
        seq_len = retrieve_length_of_sequence(value)
        human_lengths.append(seq_len)
    else:
        human_lengths.append(seq_len)

#Put rest of species intron lengths
species_lengths = []
for value in Primate_orthos["Species Gene ID"]:
    seq_len = retrieve_length_of_sequence(value)
    species_lengths.append(seq_len)


Primate_orthos["Human First Intron Length"] = human_lengths
Primate_orthos["Species First Intron Length"] = species_lengths
Primate_orthos.to_csv(OUTPATH, sep="\t")
"""
