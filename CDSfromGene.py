#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from requests.exceptions import ConnectionError
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

#We are going to parallelize in subsets to spend less time doing the
#analysis with REST API from Ensembl
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
def retrieve_length_of_cds(ID_seq, x):
    all_transcripts_lengths = {}
    ensembl_id = ID_seq
    print(ensembl_id)
    server = "https://rest.ensembl.org"
    ext1 = "/overlap/id/" + ensembl_id + "?" + ";feature=transcript;"
    r = requests.get(server+ext1, headers={ "Content-Type" : "text/plain"})
    while not r.ok:
        x+=4
        sleep(x)
        r = requests.get(server+ext1, headers={ "Content-Type" : "text/plain"})
    if r.ok:
        print("ok")
        x=4
        sleep(x)
    data = r.text.split("id: ")
    transcripts = [field.split("\n")[0] for field in data if field.startswith("ENST")]
    for trans in transcripts:
        ext2 = "/sequence/id/" + trans + "?"
        t = requests.get(server+ext2, headers={ "Content-Type" : "text/plain"})
        while not r.ok:
            x+=4
            sleep(x)
            t = requests.get(server+ext2, headers={ "Content-Type" : "text/plain"})
        if r.ok:
            print("ok")
            x=4
            sleep(x)
        all_transcripts_lengths[trans] = t.text
    right_trans = max(all_transcripts_lengths, key=lambda key: len(all_transcripts_lengths[key]))
    ext3 = "/sequence/id/" + right_trans + "?" + ";type=cds"
    cds = requests.get(server+ext3, headers={ "Content-Type" : "text/plain"})
    while not r.ok:
        x+=4
        sleep(x)
        cds = requests.get(server+ext3, headers={ "Content-Type" : "text/plain"})
    if r.ok:
        print("ok")
        x=4
        sleep(x)
    top_CDS_len = len(cds.text)
    return(top_CDS_len)


"""
Then, for each gene family we loop and print the lengths
"""

#Put human gene lengths
human_lengths = []
ref_id = "None"
for value in Primate_orthos["Human Gene ID"]:
    if value != ref_id:
        ref_id = value
        seq_len = retrieve_length_of_cds(value, time)
        human_lengths.append(seq_len)
    else:
        human_lengths.append(seq_len)

#Put rest of species intron lengths
species_lengths = []
for value in Primate_orthos["Species Gene ID"]:
    seq_len = retrieve_length_of_cds(value, time)
    species_lengths.append(seq_len)


Primate_orthos["Human Ortho CDS Length"] = human_lengths
Primate_orthos["Species Ortho CDS Length"] = species_lengths
Primate_orthos.to_csv(OUTPATH, sep="\t")
