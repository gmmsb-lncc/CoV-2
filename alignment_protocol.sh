#!/bin/bash
# execute as: bash modified_alignment_protocol.sh [CSV_REF] [CSV_COMP] [CSV_OUTPUT] [MSA_INPUT] [FASTA_ALIGNMENT] [CSV_MUTATION_OUTPUT] [CONSENSUS_OUTPUT]
# This script follows the SARS-CoV-2 multiple sequence alignment protocol as described in the README.

# Define input variables
CSV_REF="$1"
CSV_COMP="$2"
CSV_OUTPUT="$3"
MSA_INPUT="$4"
FASTA_ALIGNMENT="$5"
CSV_MUTATION_OUTPUT="$6"
CONSENSUS_OUTPUT="$7"

# Step 1: Extract Sequences from the intersection of studied sets and the Spike reference set used by GISAID
echo "Executing Step 1: Extracting Sequences..."
python3 csv_extract_columns_find_intersec.py "$CSV_REF" "$CSV_COMP" "$CSV_OUTPUT"

# Note: The next step would be the MSA using MAFFT, which is an external tool and not a script provided.
# It's assumed that the user would execute this step separately if needed.
# Execute MAFFT for multiple sequence alignment
echo "Executing MAFFT for multiple sequence alignment..."
mafft "$MSA_INPUT" > "$FASTA_ALIGNMENT"
if [ $? -ne 0 ]; then
    echo "Error executing MAFFT. Exiting."
    exit 1
fi

# Step 2: Process the MSA to quantify the variations between the amino acid residues
echo "Executing Step 2: Processing MSA..."
python3 fasta_MSA_count_mutations.py "$FASTA_ALIGNMENT" "$CSV_MUTATION_OUTPUT"

# Step 3: Identify the consensus sequences
echo "Executing Step 3: Identifying Consensus Sequences..."
python3 consensus_msa.py "$FASTA_ALIGNMENT" "$CONSENSUS_OUTPUT"

# Step 4: Visualize the heatmap corresponding to the RBD and RBM regions
echo "Executing Step 4: Generating Heatmap Visualization..."
Rscript heatmap_msa.R

echo "All steps completed!"
