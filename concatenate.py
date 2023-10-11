import shutil

# List of file names to concatenate
file_names = [
    'B.1.1.7-UK-pt1_MSA.afa',
    'B.1.1.7-UK-pt2_MSA.afa',
    'B.1.1.7-UK-pt3_MSA.afa',
    'B.1.1.7-UK-pt4_MSA.afa',
    'B.1.1.7-UK-pt4_MSA.afa',
    'B.1.1.7-UK-pt6_MSA.afa',
    'B.1.1.7-UK-pt7_MSA.afa',
] 

# Open the destination file in write binary mode
with open('B.1.1.7.afa', 'wb') as output_file:
    # Iterate through each file in the list
    for file_name in file_names:
        # Open the source file in read binary mode
        with open(file_name, 'rb') as input_file:
            # Copy the contents of the source file to the destination file
            shutil.copyfileobj(input_file, output_file)
