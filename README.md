# SARS-CoV-2 multiple sequece alignment protocol

The procedure from obtaining the sequences, building the datasets, describing the impacts of non-synonymous variations and applying multiple alignment, will be described. 

The following image displays the workflow undertaken, acquisition of sequences, and construction of datasets for multiple alignment. After identifying the relevant features contained in the multiple alignment, structural modeling will be responsible for creating three-dimensional models for the subsequent application of virtual screening for the repurposing of anti-SARS-CoV-2 drugs.

![Alt text da image](https://github.com/gmmsb-lncc/CoV-2/blob/main/workflow.png)

# Sequence Acquisition

Currently, the largest database containing SARS-CoV-2 genomic sequences is **[GISAID](https://gisaid.org) (Global Initiative on Sharing All Influenza Data)**, an international initiative where researchers from all over the world submit biologically obtained sequences following pre-established criteria.

From this database, 785,158 sequences were obtained containing the variants (Alpha(α), Beta(α), Delta(α), Gamma(α)) and 821,793 sequences containing the mutant trio (K417N, E484K, N501Y), totaling 1,606,951 sequences. These sets are crucial for identifying relevant mutations present in the receptor-binding domain (RBD). The data was collected between the dates 05/01/2020 and 06/01/2021.

It is of utmost importance to note that the announcement of the new Omicron (o) variant occurred on November 26, 2021, that is, after the data collection period. However, due to the severity of the current situation, exceptionally for this case, the information on this variant contained in the literature was considered, along with the information provided by the UK Health Security Agency (UKHSA) for structural modeling and conducting other analyses for this strain. The procedures carried out will be detailed below.

Following is an illustration of the ectodomain of the spike protein. It consists of the S1 and S2 domains. The S1 domain contains the Receptor Binding Domain (RBD) responsible for recognizing and binding to the host cell receptor. The S2 domain is responsible for fusion and contains the putative fusion peptide (FP, in turquoise) and the heptad repeat HR1 (orange) and HR2 (brown), TM is the transmembrane domain represented in purple.

![Alt text image](https://github.com/gmmsb-lncc/CoV-2/blob/main/spike_sub_units.png)

# Dataset Construction

A set composed of sequences identified in GISAID as associated with variants Alpha(α), Beta(α), Delta(α), Gamma(α) containing mutations (K417N, N439K, L452R, F456L, G476S, T478K, V483A, E484K, N501Y) was initially chosen, as these are the Variants of Concern (VoC).

The mutant trio set (K417N, E484K, N501Y) was selected because the combination of these three mutations is present in both the VoC and the Variants of Interest, and their presence results in a change in affinity with the extracellular receptor ACE2 (Angiotensin Converting Enzyme 2). This can influence the reduction of the neutralizing capacity of vaccines and/or convalescent plasma. Therefore, some sequences in the mutant trio may also be present in the variant set.

# Addressing the issue of downloading large numbers of sequences from GISAID

It's crucial to note that GISAID has a limit of 10,000 sequences for download at once; however, the file containing only the access codes (EPI_ISL) allows downloading about 270,000 access codes at once (each of the sequences has its EPI_ISL code).

It's clear that for a dataset containing 1,606,951 sequences, downloading every 10,000 sequences is an unfeasible task. Due to this characteristic of the database, the access codes to the sequences of the intersection between the previously mentioned sets and the reference file of Spike amino acid sequences used by GISAID were sought. The Python algorithm, whose function is to extract sequences that intersect between the studied sets and the Spike reference set used by GISAID (csv_extract_columns_find_intersec), is in the current directory.

# Identification of VNS through Multiple Sequence Alignment

A multiple sequence alignment (MSA) was performed using the MAFFT software (Multiple Alignment using Fast Fourier Transform) based on the previously mentioned sets, and the Spike sequence from Wuhan (NCBI YP009724390.1) was used as a reference. An illustrative visualization of the msa is shown below.

![Alt text image](https://github.com/gmmsb-lncc/CoV-2/blob/main/msa.png)

After processing the MSA, the variations between the amino acid residues for each of the subsets described here were quantified (fasta_MSA_count_mutations). For each reference residue contained in the study sequences, the physicochemical properties of the amino acids concerning the reference sequence Wuhan-Hu-1 (NCBI YP009724390.1) were annotated. Subsequently, the most frequent residues in each set were identified, and from this, the consensus sequence was sought. The Python algorithm, whose function is the identification of the consensus sequences (consensus_msa), is present in the current directory.

## The analysis of 1,606,951 sequences revealed that there are a total of 8,324 VNSs in the Spike compared to the Wuhan-Hu reference sequence (NCBI YP009724390.1) within the analyzed set. All selected variations for Spike are contained within the RBD (residues 319-541) or RBM (residues 437-508).

For the visualization of the heatmap corresponding to the RBD and RBM regions, a function was built using the gmmsa package through the R language; the script (heatmap) is in the current directory. The heatmap figure is shown below. The image displays that: A) Correlation between the consensus sequences, where the variation in shade relative to the vertical reference indicates a change in the amino acid. B) Graphical representation of the consensus sequences, the gap indicates that the residue at this position underwent some variation, and the intensity of the blue color represents high prevalence of the amino acid. The less intense the blue, the lower the prevalence of the amino acid.

![Alt text image](https://github.com/sulfierry/CoV-2/blob/main/heatmap.png)

It is important to note that consensus sequences generated from multiple alignments that have a high number of sequences (ranging from hundreds of thousands to millions) tend to be reliable when checking the prevalence of residues within the studied set. The previous figure containing the heatmap displays different visualization methods stemming from the multiple alignment for consensus sequences of variants whose binding domain comprises the RBD and RBM, particularly between residues 415 and 502. This range was selected for visualization because it is where the majority of the VNSs are concentrated in the respective binding domains. The regions between residues 319 to 414 and 503 to 541 did not show a high prevalence of VNSs in the consensus sequence set related to the α, β, δ, γ variants.

Except for the gamma variant, all the others showed mutations that do not preserve physicochemical properties. For example, the alpha variant displays variations in the polarity of residues within the RBM (S438E, N440R), while the beta variant shows more drastic polarity changes among residues in the RBD (S359K, S359E, V395K, V395T, V395E). Similarly, the mutations in the delta variant alter the polarity (K462Q, K462F, F464Y, F464R) of the RBM region. Through the approach described in this section, it was observed that the multiple alignment of Spike sequences allowed for a comprehensive overview of the variations that occur in domains critical to its function.

**Description of the `alignment_protocol.sh` Script**

**Input Variables Definition:**
   - The script takes in several arguments, including paths for:
     - Reference CSV (`CSV_REF`).
     - Comparison CSV (`CSV_COMP`).
     - Output CSV (`CSV_OUTPUT`).
     - Multiple Sequence Alignment (MSA) input (`MSA_INPUT`).
     - Alignment output (`MSA_OUTPUT`).
     - Mutation output CSV (`CSV_MUTATION_OUTPUT`).
     - Consensus output (`CONSENSUS_OUTPUT`).

**Step 1: Sequence Extraction**
   - The script utilizes the Python script `csv_extract_columns_find_intersec.py` to extract sequences from the intersection of the studied sets and the Spike reference set used by GISAID.

**Multiple Sequence Alignment (MSA) using MAFFT:**
   - The script notes that the next step would be the MSA using MAFFT, an external tool.
   - It assumes the user would execute this step separately if needed.
   - The script then attempts to run MAFFT.

**Step 2: MSA Processing**
   - The script uses the Python script `fasta_MSA_count_mutations.py` to process the MSA and quantify the variations between amino acid residues.

**Step 3: Consensus Sequence Identification**
   - The script uses the Python script `consensus_msa.py` to identify consensus sequences.

**Step 4: Heatmap Visualization**
   - The script runs the R script `heatmap_msa.R` to generate a heatmap visualization corresponding to the RBD and RBM regions.
