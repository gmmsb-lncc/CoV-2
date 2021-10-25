# uso: python3 script.py ref.csv comp.csv output_intersec.csv

import pandas as pd
import sys

csv_ref = pd.read_csv(sys.argv[1], low_memory=False)
csv_comp = pd.read_csv(sys.argv[2], low_memory=False) 


df_ref = pd.DataFrame(csv_ref)
df_comp = pd.DataFrame(csv_comp)


df_ref.rename(columns={df_ref.columns[0]:'col1'}, inplace = True)
df_comp.rename(columns={df_comp.columns[0]:'col2'}, inplace = True)


column_ref = df_ref.iloc[:,0]
column_comp = df_comp.iloc[:,0]

df_all = pd.concat([column_ref, column_comp], axis=1)

# df_all = arquivo.csv contendo 2 colunas
# col1 = referencia
# col2 = amostra à comparar 

df_all.columns = ["col1","col2"]

col1 = df_all.col1
col2 = df_all.col2

#intersec
intersec = df_all.loc[
	col1.isin(col2), ['col1']
	]

intersec.to_csv(sys.argv[3],index = False)
