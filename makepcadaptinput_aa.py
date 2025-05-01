from asyncio import protocols
from collections import defaultdict
import os
import sys
from unicodedata import name
import numpy as np
import pandas as pd
import pickle
from statistics import mean
"""
snape-output has the following fields (without the spectrum):
    0:chromosome 
  
    1:position along the chromosome
  
    2:reference
  
    3:number of reference nucleotides
  
    4:number of minor (alternative) nucleotides
  
    5:average quality of the reference nucleotides
  
    6:average quality of the alternative nucleotides
  
    7:first and second most frequent nucleotides in the pileup
  
    8:$1 - p (0)$ where $p (f)$ is the probability distribution function for
     the minor allele freqeuncy
  
    9:$p (1)$
  
    10:$E (f)$ mean value of $f$

"""
## load in the files created in snape-pooled.sh
folder = "./maf_aa/final_set/"
files =[ x for x in os.listdir(folder) if x.endswith(".sfs")]

snps = defaultdict(lambda : defaultdict(lambda: [None,None,None,None,None]))
with open("./aalo_snape_20kbwindows_finalset2.pcadapt","w") as out:
    for i,f in enumerate(files):
        with open(f"{folder}{f}") as snape:
            for line in snape:
                c, p, r, nr, na, rq,aq, bases,_,_,maf = line.strip().split("\t")
                
                alt = list(set(bases).difference(set(r)))
                alt = alt[0] if alt else None
                snps[c][p][i]=maf
               

class Window():
    def __init__(self,name,freqs):
        self.name = name
        self.freqs = freqs

windows = []
    
sorted_chrms = sorted(snps.keys(),key = lambda x :(int(x.split(".")[0][-3:])))
for chrm in sorted_chrms:
    snps_sorted = sorted(snps[chrm].keys(),key = lambda x : int(x))
    w = 0
    while True:
        
        snps_to_use = snps_sorted[w:w+20]
        if len(snps_to_use)<20:
            break
        w+=10
        try:
            windows.append(Window(f"{chrm}:{snps_to_use[0]}-{snps_to_use[-1]}",[mean([float(snps[chrm][x][i]) for x in snps_to_use]) for i in range(len(files))]))
        except:
            print(chrm,snps_to_use,[[snps[chrm][x][i] for x in snps_to_use] for i in range(len(files))])
            input()
with open("./aalo_snape_20snpswindows_finalset.pcadapt","w") as out:
    out.write("pops\t{}\n".format("\t".join([x.name for x in windows])))
    for i in range(len(files)):
        out.write("{}\t{}\n".format(files[i].split(".")[0],"\t".join([str(x.freqs[i]) for x in windows])))