import pandas as pd
import sys

def find_intersection(df_ref, df_comp, output_filename):
    """
    Find the intersection between two CSV files.
    
    Parameters:
    - df_ref: DataFrame containing the reference data.
    - df_comp: DataFrame to compare with the reference.
    - output_filename: Name of the output CSV file.
    
    Returns:
    None
    """
    
    # Ensure that the first column from both dataframes are named consistently
    df_ref.rename(columns={df_ref.columns[0]:'col1'}, inplace=True)
    df_comp.rename(columns={df_comp.columns[0]:'col2'}, inplace=True)
    
    # Find the intersection
    intersection = df_ref[df_ref['col1'].isin(df_comp['col2'])]
    
    # Write the intersection to the output file
    intersection.to_csv(output_filename, columns=['col1'], index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 script.py ref.csv comp.csv output_intersec.csv")
        sys.exit(1)
    
    csv_ref_filename = sys.argv[1]
    csv_comp_filename = sys.argv[2]
    output_filename = sys.argv[3]
    
    # Load the data
    df_ref = pd.read_csv(csv_ref_filename, low_memory=False)
    df_comp = pd.read_csv(csv_comp_filename, low_memory=False)
    
    find_intersection(df_ref, df_comp, output_filename)
