#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Align import AlignInfo


class Consensus:
    
    def __init__(self, alignment_file):
        self.alignment = AlignIO.read(alignment_file, "fasta")
        self.matrix, self.rows = self.create_matrix()

    def create_matrix(self):
        matrix = np.array(self.alignment)
        rows, cols = matrix.shape
        list_mat = []
        for j in range(cols):
            list_row = [matrix[0][j]]
            for i in range(1, rows):
                list_row.append(matrix[i][j] if matrix[i][j] != matrix[0][j] else matrix[0][j])
            list_mat.append(list_row)
        matrix_temp = np.array(list_mat).transpose()
        return matrix_temp, rows

    def canonical_consensus(self, output_file_name):
        summary_align = AlignInfo.SummaryInfo(self.alignment)
        consensus = summary_align.dumb_consensus(float(0.5), require_multiple=1)
        with open(output_file_name, "w") as f:
            f.write(f"> {output_file_name}\n")
            f.write(str(consensus) + "\n")

    def verify_domain_scores(self):
        df = pd.DataFrame(self.matrix)
        result = []
        for col in df.columns:
            result.append(df[col].value_counts())

        l_idx_con = []
        l_con = []
        norm = [291,5,14,58,146,5,19,4,51,40,21,619]

        for row in range(1,len(df)):
            if row%100 == 0:
                print(row)

            score = [0]*12
            for col in df.columns:                
                flag = -1        
                if (col>=12 and col<=302): #NTD
                    flag = 0
                elif (col>=400 and col<=404): #RGD:
                    flag = 1
                elif (col>=472 and col<=485): #HIGH-CFR1
                    flag = 2
                elif (col>=436 and col<=507): #RBM
                    flag = 3
                elif (col>=318 and col<=540): #RBD
                    flag = 4
                elif (col>=680 and col<=684): #S1/S2
                    flag = 5
                elif (col>=787 and col<=805): #FP
                    flag = 6
                elif (col>=890 and col<=893): #HIGH-CFR2
                    flag = 7
                elif (col>=919 and col<=969): #HR1
                    flag = 8
                elif (col>=1162 and col<=1201): #HR2
                    flag = 9
                elif (col>=1213 and col<=1233): #TM
                    flag = 10
                else:
                    flag = 11           
            
                residuo = df.loc[row, col]
                freq = list(result[col].index).index(residuo)
                score[flag] += result[col][freq]/sum(result[col])
            
            l_idx_con.append(row)
            l_con.append(np.array(score)/np.array(norm))

        m_score = np.transpose(np.array(l_con))
        return m_score, l_idx_con


    def find_seq_best_score(self):
        domain_scores, idx_line = self.verify_domain_scores()
        df_score = pd.DataFrame()
        df_score['ROW'] = idx_line_
        df_score['NTD'] = domain_score_[0]
        df_score['RGD'] = domain_score_[1]
        df_score['HIGH-CFR1'] = domain_score_[2]
        df_score['RBM'] = domain_score_[3]
        df_score['RBD'] = domain_score_[4]
        df_score['S1/S2'] = domain_score_[5]
        df_score['FP'] = domain_score_[6]
        df_score['HIGH-CFR2'] = domain_score_[7]
        df_score['HR1'] = domain_score_[8]
        df_score['HR2'] = domain_score_[9]
        df_score['TM'] = domain_score_[10]
        df_score['OTHERS'] = domain_score_[11]
        
        df_temp = df_score.sort_values(['GENERAL'], ascending=False).copy()
        df_temp.index = range(len(df_temp))

        l_scores = list(df_temp.drop_duplicates(['GENERAL'])['GENERAL'])
        df_scorefinal = pd.DataFrame(columns=['SEQ','SCORE','FREQ'])
        
        # loop for populating df_scorefinal
        for i in range(len(l_scores)):
            rows = list(df_temp[df_temp['GENERAL'] == l_scores[i]]['ROW'])
            l_seq = []
            l_country = []
            for row in rows:
                l_seq.append(str(alignment_[row].seq))
                l_country.append(alignment_[row].id.split('/')[1])
            l_seq = pd.Series(l_seq)
            l_country = pd.Series(l_country)
            
            df_seq_country = pd.DataFrame()
            df_seq_country['SEQ'] = l_seq
            df_seq_country['OCCUR'] = l_country
            
            df_temp_group = df_seq_country.groupby(['SEQ'], as_index=False).count()
            
            if len(df_temp_group) > 1:
                l_country = []
                l_qtd_country = []
                for k in range(len(df_temp_group)):
                    l_country.append(list(df_seq_country[df_seq_country['SEQ'] == df_temp_group.loc[k,'SEQ']]['OCCUR'].value_counts().index))
                    l_qtd_country.append(list(df_seq_country[df_seq_country['SEQ'] == df_temp_group.loc[k,'SEQ']]['OCCUR'].value_counts()))
                df_temp_group['COUNTRIES'] = l_country
                df_temp_group['FREQ'] = l_qtd_country
                df_temp_group['SCORE'] = l_scores[i]
                df_scorefinal = df_scorefinal.append(df_temp_group).copy()
            else:
                df_scorefinal_temp = pd.DataFrame(columns=['SEQ','SCORE','OCCUR','COUNTRIES','FREQ'])
                df_scorefinal_temp['SEQ'] = list(l_seq.value_counts().index)
                df_scorefinal_temp['SCORE'] = l_scores[i]
                df_scorefinal_temp['OCCUR'] = list(l_seq.value_counts())
                df_scorefinal_temp['COUNTRIES'] = [list(l_country.value_counts().index)]
                df_scorefinal_temp['FREQ'] = [list(l_country.value_counts())]
                df_scorefinal = df_scorefinal.append(df_scorefinal_temp).copy()

        
        df_scorefinal.index = range(len(df_scorefinal))
        return df_scorefinal

    def save_seq_scores(self, csv_filename):
        seq_score = self.find_seq_best_score()
        seq_score[seq_score['OCCUR'] >= (self.rows / 100)].to_csv(csv_filename, index=False)


if __name__ == "__main__":
    consensus_obj = Consensus(sys.argv[1])
    consensus_obj.save_seq_scores(sys.argv[2])
    consensus_obj.canonical_consensus(sys.argv[3])
