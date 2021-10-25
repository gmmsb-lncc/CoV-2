# uso: python3 script.py ID.csv spike_ref_seq.fasta name_output.fasta
import pandas as pd
import numpy as np
import sys



# reading reference file
df = pd.read_csv(sys.argv[1])

# turning column reference into list
compare = list(df['col1'])

# reading line by line from fasta file
file = open(sys.argv[2])
lines = file.readlines()
file.close()


inp = sys.argv[1].split('.')
name_trim = inp[0]


# declaring vectors for the iteration
identifiers = []
lengths = []
sequence = []

# even line = header, odd line = sequence
for i in range(len(lines)):
    if i%2==0:                     # verificando se é par
        identifiers.append(lines[i])
    else:
        sequence.append(lines[i][:-1])
        lengths.append(len(sequence[-1]))

# breaking the line at position "ID"
identifiers_tag = [x.split('|')[3] for x in identifiers]


#print(ident)


# searching for identifiers referring to item "i" in the reference list (compare)
l_id = []
for i in range(len(compare)):
	try:
		l_id.append(identifiers_tag.index(compare[i]))
	except:
		print(compare[i],' Não encontrado')
# adding the string "i" to the reference list
l_seq = []
for i in l_id:
    l_seq.append(sequence[i])

# generating data frame
df = pd.DataFrame()
df['ID'] = compare
df['SEQ'] = l_id
df['SEQUENCE'] = l_seq

# export in fast format
ofile = open(sys.argv[3], 'w')

'''
# PRINT ID + INFO
for i in l_id:

    ofile.write(str(identifiers[i]) + str(sequence[i])+'\n')

ofile.close()

'''
# PRINT ID ONLY
for i in l_id:

    ofile.write(str(identifiers[i]) + str(sequence[i])+'\n')

ofile.close()


# exporting in csv
#df.to_csv('./df_compare3.csv', index=False)
