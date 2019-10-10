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
OUTPATH= PATH + "Mammals/GeneLengths.tsv"
"""
First we open the ortholog list of the species groups we work with
"""
cols = pd.read_csv(PATH + "Mammals/ortho1to1.tsv", sep='\t', nrows=1).columns
Mammals_orthos = pd.read_csv(PATH + "Mammals/ortho1to1.tsv", sep='\t', low_memory=False, usecols=cols[1:])


"""
API FUNCTION TO RETRIEVE SEQUENCE LENGTH
"""
def retrieve_length_of_sequence(ID_seq):
    ensembl_id = ID_seq
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/" + ensembl_id + "?"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    while not r.ok:
        sleep(8)
        r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    if r.ok:
        print("ok")
    sequence = r.text
    if sequence:
        gene_length = len(sequence)
    else:
        gene_length = 0
    return(gene_length)


"""
Then, for each gene family we loop and print the lengths
"""

#Put human gene lengths
human_lengths = []
ref_id = "None"
for value in Mammals_orthos["Human Gene ID"]:
    if value != ref_id:
        ref_id = value
        seq_len = retrieve_length_of_sequence(value)
        human_lengths.append(seq_len)
    else:
        human_lengths.append(seq_len)

#Put rest of species intron lengths
species_lengths = []
for value in Mammals_orthos["Species Gene ID"]:
    seq_len = retrieve_length_of_sequence(value)
    species_lengths.append(seq_len)


Mammals_orthos["Human Gene Length"] = human_lengths
Mammals_orthos["Species Gene Length"] = species_lengths
Mammals_orthos.to_csv(OUTPATH, sep="\t")

