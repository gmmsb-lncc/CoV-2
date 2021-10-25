# python3 script.py sequence.fasta

import matplotlib.pyplot as plt
import pandas as pd
import sys


def main():

	input_01 = sys.argv[1]

	# lendo linhas - OK 
	rows = read_rows_from_file(input_01)

	# removendo a extenção do nome do input que será utilizado pelas funções - OK 
	inp = sys.argv[1].split('.')
	input_name_trim = inp[0]

	# formatar arquivo: linha par = ID, linha impar = sequencia - OK
	file_id_seq_format = format_row_id_seq(rows, input_name_trim)

	# verificando frequencia dos paises - OK
	frequency = frequency_filter(input_name_trim)



def read_rows_from_file(fasta_file):

	file = open(fasta_file)
	lines = file.readlines()
	file.close()

	return lines


def format_row_id_seq(rows, input_name_trim):

	# criando novo arquivo
	ofile = open(input_name_trim +'_formatado'+'.fasta', 'w')

	for i in range(len(rows)):
		if rows[i][0] == '>':
			ofile.write('\n'+rows[i])
		else:
			ofile.write(rows[i][:-1])

	ofile.close()

	line_del = 0
	remove_line =  remove_first_line(input_name_trim, line_del)

def remove_first_line(input_name_trim, line_del):

	# deletando primeira linha
	filename = input_name_trim +'_formatado'+'.fasta'
	line_to_delete = line_del
	initial_line = 0
	file_lines = {}
	
	content = read_rows_from_file(filename)

	for line in content:
	    file_lines[initial_line] = line.strip()
	    initial_line += 1

	file = open(filename, "w")
	for line_number, line_content in file_lines.items():
	    if line_number != line_to_delete:
	    	file.write('{}\n'.format(line_content))

	file.close()


def frequency_filter(input_name_trim):

	txt_spike_prot = read_rows_from_file(input_name_trim +'_formatado'+'.fasta')  

	# separar cabeçalho em colunas
	spike_prot = [line.split('/') for line in txt_spike_prot if '>' in line]

	# converte para data frame
	df_spike_prot = pd.DataFrame(spike_prot)

	# pega coluna de interesse
	df_columns_spike_prot = df_spike_prot.iloc[:,1] # spike_ref.fast

	df = pd.DataFrame(df_columns_spike_prot)

	df.rename(columns={ df.columns[0]: 'col1'}, inplace = True)

	freq = df['col1'].value_counts()

	print(freq)
	print('Frequencias extraidas com sucesso!')

	freq.to_csv(input_name_trim +'_frequency_country.csv', index = True)





main()
