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

PATH=sys.argv[1]
OUTPATH=sys.argv[2]
"""
First we open the ortholog list of the species groups we work with
"""
cols = pd.read_csv(PATH, sep='\t', nrows=1).columns
Primate_orthos = pd.read_csv(PATH, sep='\t', low_memory=False, usecols=cols[1:])
time=0

"""
API FUNCTION TO RETRIEVE SEQUENCE WITH INTRONS
"""
def retrieve_length_of_sequence(ID_seq, x):
    ensembl_id = ID_seq
    print(ensembl_id)
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/" + ensembl_id + "?" + ";mask_feature=1"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    while not r.ok:
        x+=6
        sleep(x)
        r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    if r.ok:
        print("ok")
        x=6
        sleep(x)
    sequence = r.text
    list = re.findall('[a-z]+', sequence)
    if list:
        intron_length = sum(len(el) for el in list)
    else:
        intron_length = 0
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
        seq_len = retrieve_length_of_sequence(value, time)
        human_lengths.append(seq_len)
    else:
        human_lengths.append(seq_len)

#Put rest of species intron lengths
species_lengths = []
for value in Primate_orthos["Species Gene ID"]:
    seq_len = retrieve_length_of_sequence(value, time)
    species_lengths.append(seq_len)


Primate_orthos["Human Total Intron Length"] = human_lengths
Primate_orthos["Species Total Intron Length"] = species_lengths
Primate_orthos.to_csv(OUTPATH, sep="\t")
