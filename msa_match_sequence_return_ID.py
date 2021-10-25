#!/usr/bin/env python
# coding: utf-8

# uso: python3 script.py seq.fasta output.csv
# the precise reference sequence is on the first line

from Bio import AlignIO
from Bio import SeqIO
import pandas as pd
import numpy as np
import sys


file = AlignIO.read(sys.argv[1],'fasta')
matrix = np.array(file)
rows, cols = matrix.shape



l_seq = []
l_id = []
for i in range(rows):
    l_seq.append(str(file[i].seq))
    l_id.append(file[i].id)
l_seq = pd.Series(l_seq)
l_id = pd.Series(l_id)



list_row = []
for j in range(rows):
    if l_seq[0] == l_seq[j]:
        list_row.append(l_id[j])
    df_match = pd.DataFrame()
    df_match['ID'] = list_row


df_match.to_csv(sys.argv[2], index=False)

