import pandas as pd 
from pandas import DataFrame
import Bio
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections
from Bio import AlignIO
import random
import numpy as np
import argparse
import re
import os 
import multiprocessing
import subprocess
import itertools


# Dictionaries

## IUPAC dictionary
IUPAC_CODES = {
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],  
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"]
}

## Ligation junction dictionary
ligation_junctions_dict = {'TA': 'preferred',
                        'TA': 'preferred',
                        'GA': 'preferred',
                        'AG': 'preferred',
                        'TT': 'neutral',
                        'CT': 'neutral',
                        'CA': 'neutral',
                        'TC': 'neutral',
                        'AC': 'neutral',
                        'CC': 'neutral',
                        'TG': 'neutral',
                        'AA': 'neutral', 
                        'CG': 'non-preferred', 
                        'GT': 'non-preferred',
                        'GG': 'non-preferred',
                        'GC': 'non-preferred'}


# Functions

## Evaluate IUPAC mismatches format:
def parse_iupac_mismatches(mismatch_str):
    """
    Parses a string of mismatches formatted as "pos:base,pos:base" into a list of tuples.
    
    Args:
        mismatch_str (str): Mismatch input string (e.g., "5:R,10:G")
    
    Returns:
        list: A list of (position, base) tuples, e.g., [(5, 'R'), (10, 'G')].
    """
    mismatches = []
    try:
        for pair in mismatch_str.split(","):
            pos, base = pair.split(":")
            pos = int(pos.strip())  # Convert position to integer
            base = base.strip().upper()  # Ensure base is uppercase
            mismatches.append((pos, base))

    except (ValueError, IndexError):
        raise InputValueError("Invalid format for --iupac_mismatches. Use 'pos:base,pos:base', e.g., '5:R,10:G'.", 
                              field="iupac_mismatches", code="invalid_mismatch_format")
    
    return mismatches
## Evaluate ligation junction functions:
def evaluate_ligation_junction(targets, iupac_mismatches=None, plp_length=30, num_probes=10):
    """
    Evaluates the ligation junction of a probe and introduces mismatches if needed.

    Args:
        probe_seq (str): The probe sequence.
        iupac_mismatches (list of tuples): List of positions and IUPAC codes or bases to introduce mismatches. Note that the number of mismatches should be less than or equal to 2.
                                            Example: 5:R,10:G
        plp_length (int): Probe length (default: 30).

    Returns:
        tuple: (updated probe sequence, ligation junction category)
    """
    # Add Ligation junction column to the DataFrame
    if 'Ligation junction' not in targets.columns:
        targets['Ligation junction'] = 'non-preferred'


    # Extract the ligation junction (2 bases around the center of the probe)
    junction_position = int((plp_length / 2) - 1)
    new_rows = []
    for idx in targets.index:

        probe_seq = targets.loc[idx]['Sequence']
        ligation_junction = probe_seq[junction_position] + probe_seq[junction_position + 2]

        # Determine the ligation status
        ligation_status = ligation_junctions_dict.get(ligation_junction, "non-preferred")
        targets.loc[idx, 'Ligation junction'] = ligation_status

        # Apply IUPAC mismatches if provided
        if iupac_mismatches is not None:
            if isinstance(iupac_mismatches, str):  # Only parse if it's a string
                iupac_mismatches = parse_iupac_mismatches(iupac_mismatches)

            num_mismatches = len(iupac_mismatches)
#s            print("Num mismatch:", num_mismatches)
            if num_mismatches <= 2: # Limit to 2 mismatches
                for r in range(1, num_mismatches + 1): 
                    for subset in itertools.combinations(range(num_mismatches), r): # Generate all possible combinations of mismatches
                        selected_mismatches = [iupac_mismatches[i] for i in subset]
                        replacement_options = [IUPAC_CODES[symbol] for pos, symbol in selected_mismatches]
#                        print('Replacement options:', replacement_options)

                        # Generate all possible combinations of replacements
                        for replacement in itertools.product(*replacement_options):
                            new_probe_seq = list(probe_seq)
                            new_id_suffix = []
                            for (pos, iupac_symbol), new_base in zip(selected_mismatches, replacement):
                                new_probe_seq[pos] = new_base
                                new_id_suffix.append(f"{pos}_{iupac_symbol}_{new_base}")
                            new_probe_seq = "".join(new_probe_seq)
                            new_probe_id = f"{targets.loc[idx, 'Probe_id']}|{'_'.join(new_id_suffix)}"
                            # Revaluate the ligation junction
                            new_ligation_junction = new_probe_seq[junction_position] + new_probe_seq[junction_position + 2]
                            new_ligation_status = ligation_junctions_dict.get(new_ligation_junction, "non-preferred")
                            new_row = targets.loc[idx].copy()
                            new_row['Sequence'] = new_probe_seq
                            new_row['Ligation junction'] = new_ligation_status
                            new_row['Probe_id'] = new_probe_id
                            new_rows.append((new_probe_id, new_row))
            else:
                raise InputValueError("The number of mismatches should be less than or equal to 2", field="iupac_mismatches", code="mismatches_exceed_limit")
   
    # Append the new rows to the DataFrame
    new_rows_df = pd.DataFrame([row[1] for row in new_rows], index=[row[0] for row in new_rows])
    targets = pd.concat([targets, new_rows_df])
    return targets



## Custom exception classes:
class InputValueError(ValueError):
    """
    Custom exception class for invalid
    input values.
    """
    def __init__(self, message: str, field: str, code: str):
        super().__init__(message)
        self.field = field
        self.code = code


## Extract features functions:

def parse_gtf_to_dataframe(gtf_path: str) -> pd.DataFrame:
    """
    Parses a GTF (Gene Transfer Format) file into a structured Pandas DataFrame.

    Args:
        gtf_path (str): Path to the GTF file.

    Returns:
        pd.DataFrame: A DataFrame containing parsed GTF data.
    
    Reference:
        Adapted from Ricardo Filipe dos Santos script:
        https://gist.github.com/rf-santos/22f521c62ca2f85ac9582bf0d91e4054
    """
    print('Loading GTF file....')
    # Define column names
    col_names = [
        'seqname', 'source', 'feature', 'start', 'end', 'score', 
        'strand', 'frame', 'attribute'
    ]
    
    # Read the GTF file, skipping comments
    df = pd.read_csv(gtf_path, sep='\t', header=None, names=col_names, comment='#', dtype={'score': str})

    # Define attributes to extract
    attributes = [
        'gene_id', 'transcript_id', 'exon_number', 'gene_name', 'gene_source', 
        'gene_biotype', 'transcript_name', 'transcript_source', 'transcript_biotype', 
        'protein_id', 'exon_id', 'tag'
    ]

    # Create regex patterns for attribute extraction
    attr_patterns = {attr: rf'{attr} "([^"]*)"' for attr in attributes}

    # Extract attributes using vectorized operations
    for attr, pattern in attr_patterns.items():
        df[attr] = df['attribute'].str.extract(pattern)

    # Drop the original attribute column
    df.drop(columns=['attribute'], inplace=True)
    
    return df

def parse_gtf(gtf_file, genes_str=None, identifier_type='gene_id', gene_feature='CDS'):
    """
    Parses a GTF file and yields data only for the specified genes_of_interest.

    Args:
        gtf_file (str): Path to the GTF file.
        genes_of_interest (set or None): A set of gene IDs or names to parse.
                                         If None, parse all genes in the GTF.
        identifier_type (str): Type of identifier provided ('gene_id' or 'gene_name').

    Returns:
        pd.DataFrame: A DataFrame containing parsed GTF data.
    """
    print('Parsing GTF file....')
    # Read the GTF file into a DataFrame
    gtf_df = parse_gtf_to_dataframe(gtf_file)

    # Convert gene names to lowercase for case-insensitive matching
    if genes_str:
            genes_of_interest = set([g.strip().lower() for g in genes_str.split(",")])
            print(f"Processing genes: {', '.join(genes_of_interest)}")
    else:
        genes_of_interest = None
        raise InputValueError("No gene list provided. Processing all genes. $genes_of_interest", field="genes_of_interest", code="no_gene_list_provided")

    # Check if a gene list is provided and filter accordingly
    
    if genes_of_interest:
        if identifier_type == 'gene_id':
            gtf_df = gtf_df[gtf_df['gene_id'].str.lower().isin(genes_of_interest) & (gtf_df['feature'] == gene_feature)]
        elif identifier_type == 'gene_name':
            gtf_df = gtf_df[gtf_df['gene_name'].str.lower().isin(genes_of_interest)& (gtf_df['feature'] == gene_feature)]
        else:
            raise InputValueError("Gene identifier type must be 'gene_id' or 'gene_name'", field="identifier_type", code="no_identifier_type_provided")
            

    if len(gtf_df) == 0:
        raise InputValueError("No matching genes found in the GTF file. This can be due to either incomplete gtf file or errors in gene identifications.", field="genes_of_interest", code="no_matching_genes_found")

    return gtf_df, genes_of_interest


def merge_regions_and_coverage(genes_of_interest, gtf_df):
    """
    Merge CDS regions and computes average coverage for the specified genes.

    Args:
        genes_of_interest (set): A set of gene names to process.
        gtf_df (pd.DataFrame): A DataFrame containing parsed GTF data.

    Returns:
        pd.DataFrame: A DataFrame containing merged CDS regions and average coverage.
    """
    print('Merge regions and calculating coverage....')
    # Initialize an empty list to store results
    merged_regions = []

    for gn in genes_of_interest:
        # Subset the DataFrame for the given gene_name
        isoforms = gtf_df[gtf_df['gene_name'].str.lower() == gn].copy()
        
        # Sort by chromosome, strand, and start position
        isoforms = isoforms.sort_values(by=['seqname', 'strand', 'start'])

        # Determine the full genomic range
        min_start = isoforms['start'].min()
        max_end = isoforms['end'].max()


        # Create a NumPy array to track coverage over this genomic range
        coverage_array = np.zeros(max_end - min_start + 1, dtype=int)

        # Dictionary to track merged CDS regions
        merged = []
        
        for _, row in isoforms.iterrows():
            if not merged:
                merged.append(row.to_dict())  
            else:
                prev = merged[-1]

                if row['start'] <= prev['end']:  
                    merged[-1]['end'] = max(prev['end'], row['end'])  
                    merged[-1]['transcript_id'] += ";" + row['transcript_id']  
                else:
                    merged.append(row.to_dict())  

            # Update coverage using NumPy slicing (avoids looping over each position)
            coverage_array[row['start'] - min_start : row['end'] - min_start + 1] += 1

        # Convert merged results to a DataFrame
        merged_df = pd.DataFrame(merged)

        # Compute the average coverage for each merged region efficiently
        coverage_values = []
        for _, region in merged_df.iterrows():
            region_slice = coverage_array[region['start'] - min_start : region['end'] - min_start + 1]
            avg_coverage = np.mean(region_slice)  # Vectorized mean calculation
            coverage_values.append(avg_coverage)

        # Add coverage column to DataFrame
        merged_df['coverage'] = coverage_values

        merged_regions.append(merged_df)

    # Concatenate all gene-specific DataFrames into a final result
    final_df = pd.concat(merged_regions, ignore_index=True)

    # Add region_id column for easier downstream analysis
    final_df['region'] = final_df['seqname'].astype(str) + ":" + final_df['start'].astype(str) + "-" + final_df['end'].astype(str)
    
    # Keep only relevant columns
    final_df = final_df[['seqname', 'start', 'end', 'strand', 'gene_name', 'transcript_id', 'coverage', 'region']]

    # Display the DataFrame with coverage
    return(final_df)

## Extract sequences functions:
def save_regions_for_faidx(df, output_file, plp_length=30, identifier_type = 'gene_name'):
    """
    Saves genomic regions in a format compatible with `samtools faidx`.

    Args:
        df (pd.DataFrame): DataFrame with 'seqname', 'start', 'end' columns.
        output_file (str): Path to save the region file.
    """
    print(f"Saving regions to {output_file} for `samtools faidx`...")

    # Setting the datatype
    df['seqname'] = df['seqname'].astype(str)
    df['start'] = df['start'].astype(int)  
    df['end'] = df['end'].astype(int) 
    plp_length = int(plp_length)


    # Filter regions based on length
    df = df[(df['end'] - df['start']) >= (plp_length + plp_length / 2)]
    
    # Export regions in the format expected by `samtools faidx`
    df[['region']].to_csv(output_file + ".txt", sep = '\t', header=False, index=False)

def check_fasta_index(fasta_file):
    """
    Checks if a FASTA index file (.fai) exists.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        str: Path to the valid FASTA index file.
    """
    fai_file_1 = fasta_file + ".fai"  # Example: reference.fa.fai
    fai_file_2 = os.path.splitext(fasta_file)[0] + ".fai"  # Example: reference.fai

    if os.path.exists(fai_file_1):
        return fai_file_1
    elif os.path.exists(fai_file_2):
        return fai_file_2
    else:
        print(f"⚠️ FASTA index file not found for {fasta_file}!")
        print("Creating the index now... This may take a while for large genomes.")
        subprocess.run(["samtools", "faidx", fasta_file])
        return fasta_file + ".fai"  # Assume samtools follows the convention


def extract_sequences(fasta_file, regions_file, output_fasta, gtf_df):
    """
    Extracts CDS sequences using `samtools faidx`, adjusts for strand orientation, 
    and includes gene identifiers in the FASTA headers.

    Args:
        fasta_file (str): Path to the indexed FASTA file.
        regions_file (str): File with genomic regions (one per line).
        output_fasta (str): Path to save the extracted sequences.
        gtf_df (pd.DataFrame): DataFrame containing 'seqname', 'start', 'end', 'strand', 'gene_name' for strand correction.
    """
    print(f"Extracting sequences from {fasta_file} using `samtools faidx`...")

    # Get number of CPUs for multi-threading
    num_cpus = multiprocessing.cpu_count()

    # Temporary output file before strand correction
    temp_fasta = "temp_extracted.fa"

    # Run samtools faidx to extract sequences
    command = [
        "samtools", "faidx", "-@", str(num_cpus), 
        "-r", regions_file, 
        "-o", temp_fasta, 
        fasta_file
    ]
    subprocess.run(command, check=True)

    print("✅ Sequences extracted. Now adjusting for strand orientation and updating headers...")

    # Load region-to-gene mapping from gtf_df
    feature_dict = dict(zip(
        gtf_df['region'],
        zip(gtf_df['gene_name'], gtf_df['strand'])  
    ))

    # Read extracted sequences and apply strand correction
    updated_sequences = []
    seq_dict = SeqIO.to_dict(SeqIO.parse(temp_fasta, "fasta"))  

    for region, record in seq_dict.items():
        gene_name, strand = feature_dict.get(region, ("UNKNOWN", "+"))  
        sequence = record.seq

        # Reverse-complement if the gene is on the negative strand
        if strand == "-":
            sequence = sequence.reverse_complement()

        # Update header: >gene_name|region
        record.id = f"{gene_name}|{region}"
        record.description = ""  # Remove extra description
        record.seq = sequence
        updated_sequences.append(record)

    # Write the updated FASTA
    SeqIO.write(updated_sequences, output_fasta, "fasta")

    print(f"✅ Final sequences saved to {output_fasta}, with correct strand orientation and gene names.")


## Find target functions:
def find_targets(selected_features, fasta_file, plp_length, min_coverage, output_file='Candidate_probes.txt', gc_min=50, gc_max=65, num_probes=10, iupac_mismatches=None):
    """
    Function to extract target sequences fulfilling the following criteria:
    Adapted from sequence developed by Sergio 
    Args:
        selected_features (str): Path to the selected features file (TSV format).
        fasta_file (str): Path to the indexed FASTA file.
        output_file (str): Path to the output file.
        plp_length (int): Minimum probe length.
        min_coverage (int): Minimum coverage of the region.
        gc_min (int): Minimum GC content (default: 50).
        gc_max (int): Maximum GC content (default: 65).
        num_probes (int): Number of probes to select per gene (default: 10).
        iupac_mismatches (list of tuples): List of positions and IUPAC codes to introduce mismatches. Note that the number of mismatches should be less than or equal to 2.
                                            Example: 5:R,10:G
    Returns:
        DataFrame: DataFrame with extracted probe sequences.
    """
    # Create a dataframe to store the targets
    targets = pd.DataFrame(columns=['Probe_id', 'Gene', 'Region', 'Sequence', 'GC', 'Coverage'])


    # check if size of the probe is an even number
    if plp_length % 2 != 0:
        raise InputValueError("The size of the probe should be an even number", field="plp_length", code="odd_value_for_plp_length_provided")
    
        
    # load the fasta file
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))  

    # load the selected features table with gene name and region
    selected_features= pd.read_csv(selected_features, sep='\t')

    # Convert selected_features into a dictionary for fast lookup
    coverage_dict = dict(zip(selected_features['region'], selected_features['coverage']))


    # Initialize DataFrame for results
    targets = []

    # loop through the fasta file
    for keys in seq_dict:
        # extract gene and region of the sequence
        gene, region = keys.split("|")

        # Extract chromosome and start position from region
        chr_name, start_end = region.split(":")
        initial_start = int(start_end.split("-")[0])

        # Extract the sequence
        seq = seq_dict[keys].seq
        seq_len = len(seq)-(plp_length-1) # Limit the sequence length to extract probes

        # Keep track of the coverage value for the region and check the coverage to be greater than the minimum coverage
        coverage_value = coverage_dict.get(region, 0)
        if coverage_value < min_coverage:
            continue 

        # Loop through the sequence
        for i in range(0, seq_len):
            # Extract the probe sequence
            tmp_seq = seq[i:i+plp_length]
            # Calculate the gc content and check if it is within the range
            gc_content = gc_fraction(tmp_seq)*100

            if gc_content < gc_min or gc_content > gc_max:
                continue

            if not any(nucleotide * 3 in tmp_seq for nucleotide in "ACGT"):
                continue

            start = initial_start + i
            end = start + plp_length -1

            # save the target to the dataframe
            targets.append({
                "Probe_id": f"{gene}|{start}-{end}",
                "Gene": gene,
                "Region": f"{chr_name}:{start}-{end}",
                "Sequence": str(tmp_seq),
                "GC": gc_content,
                "Coverage": coverage_value
            })
    targets_df = pd.DataFrame(targets)
    targets_df = evaluate_ligation_junction(targets_df, iupac_mismatches=iupac_mismatches, plp_length=plp_length)
    create_fasta(targets_df, output_file)
    print(f"✅ Final sequences saved to {output_file} and fasta file {output_file}.fa .")
    return targets_df

def create_fasta(targets_df, output_file):
    """
    Create a FASTA file from the DataFrame with probe sequences.

    Args:
        targets_df (DataFrame): DataFrame with probe sequences.
        output_file (str): Path to the output file.
    """
    output_file = output_file + ".fa"
    with open(output_file, "w") as f:
        for _,row in targets_df.iterrows():
            f.write(f">{row['Probe_id']}\n{row['Sequence']}\n")

def select_top_probes(df, num_probes):
    """
    Select the top probes based on coverage and ligation junction preferences.

    Args:
        df (DataFrame): DataFrame with probe sequences.
        num_probes (int): Number of probes to select.

    Returns:
        DataFrame: DataFrame with the top probes.
    """

    # Sort the dataframe: Highest coverage first, then ligation junction order
    df_sorted = df.sort_values(by=["Coverage", "Ligation junction"], 
                               ascending=[False, True], 
                               key=lambda col: col.map({'preferred': 0, 'neutral': 1, 'non-preferred': 2}).fillna(3))
    
    # Select top `num_probes`
    final_df = df_sorted.groupby("Gene").head(num_probes)
    
    return final_df

## Check specificity functions:
