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

#We are going to parallelize in subsets to spend less time doing the
#analysis with REST API from Ensembl
#VARIABLES

PATH="/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/"
OUTPATH= PATH + "Primates/"
"""
First we open the ortholog list of the species groups we work with
"""
cols = pd.read_csv(PATH + "Primates/ortho1to1.tsv", sep='\t', nrows=1).columns
Primate_orthos = pd.read_csv(PATH + "Primates/ortho1to1.tsv", sep='\t', low_memory=False, usecols=cols[1:])
size = 13000 #number of ortholog lines

def split_dataframe_orthologs(db, size):
    list_of_dfs = (db.loc[i:i+size-1,:] for i in range(0, len(db), size))
    return(list_of_dfs)

subsets = split_dataframe_orthologs(Primate_orthos, size)
count = 0
for each in subsets:
    count +=1
    Primate_orthos.to_csv(OUTPATH + "ortho1to1_subset" + str(count) + ".tsv", sep="\t")
