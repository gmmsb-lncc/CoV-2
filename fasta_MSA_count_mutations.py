# uso: python3 script.py msa_input.fasta  output_count.csv

import sys
import numpy as np
import pandas as pd
from Bio import AlignIO

def variations_count(alignment, result_variations):
    
    matrix = np.array(alignment)
    rows, cols = matrix.shape

    list_mat = []
    for j in range(cols):
        list_row = []
        list_row.append(matrix[0][j])
        for i in range(1,rows):
            if matrix[i][j] != matrix[0][j]:
                list_row.append(matrix[i][j])
            else:
                list_row.append('0')


        list_mat.append(list_row)

    matrix_temp = np.array(list_mat)
    matrix_temp = matrix_temp.transpose()
    df = pd.DataFrame(matrix_temp)

    result = []
    for col in df.columns:
        result.append(df[col].value_counts())

    for row in result:
        print(row)

    # Create the new table
    rows, cols = matrix_temp.shape
    l_changes = []
    l_qtd = []
    l_row_changes = []


    for j in range(cols):
        for i in range(1,rows):
            if (matrix_temp[i][j] != 'X') and (matrix_temp[i][j] != '0') and (matrix_temp[i][j] != '-'): # para desconsiderar o gap '-'
                str_change = matrix_temp[0][j]+str((j+1))+matrix_temp[i][j]
                if str_change not in l_changes:
                    l_changes.append(str_change)
                    l_qtd.append(1)
                    l_row_changes.append(str(i+1))
                else:
                    idx = l_changes.index(str_change)
                    l_qtd[idx] = l_qtd[idx]+1
                    l_row_changes[idx] = l_row_changes[idx]+','+str(i+1)

    df_result = pd.DataFrame()
    df_result['Changes'] = l_changes
    df_result['Times'] = l_qtd
    df_result['Lines'] = l_row_changes
    print(df_result)
    
    return df_result.to_csv(result_variations, index=False)
    

if __name__=="__main__":
    alignment = AlignIO.read(sys.argv[1], 'fasta')
    result_variations = sys.argv[2]
    
    variations_count(alignment, result_variations)
