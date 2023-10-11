import sys
import numpy as np
import pandas as pd
from Bio import AlignIO

class MutationCount:
    
    def __init__(self, alignment):
        self.matrix = np.array(alignment)
        
    def _get_variation_matrix(self):
        rows, cols = self.matrix.shape

        list_mat = []
        for j in range(cols):
            list_row = [self.matrix[0][j]]
            for i in range(1, rows):
                list_row.append(self.matrix[i][j] if self.matrix[i][j] != self.matrix[0][j] else '0')
            list_mat.append(list_row)

        return np.array(list_mat).transpose()

    def variations_count(self):
        matrix_temp = self._get_variation_matrix()
        df = pd.DataFrame(matrix_temp)

        rows, cols = matrix_temp.shape
        l_changes, l_qtd, l_row_changes = [], [], []

        for j in range(cols):
            for i in range(1, rows):
                if matrix_temp[i][j] not in ('X', '0', '-'):  # Ignoring 'X', '0', and gaps '-'
                    str_change = matrix_temp[0][j] + str(j + 1) + matrix_temp[i][j]
                    if str_change not in l_changes:
                        l_changes.append(str_change)
                        l_qtd.append(1)
                        l_row_changes.append(str(i + 1))
                    else:
                        idx = l_changes.index(str_change)
                        l_qtd[idx] += 1
                        l_row_changes[idx] += ',' + str(i + 1)

        df_result = pd.DataFrame({
            'Changes': l_changes,
            'Times': l_qtd,
            'Lines': l_row_changes
        })

        print(df_result)
        return df_result

if __name__ == "__main__":
    alignment = AlignIO.read(sys.argv[1], 'fasta')
    mutation_counter = MutationCount(alignment)
    result = mutation_counter.variations_count()
    result.to_csv(sys.argv[2], index=False)
