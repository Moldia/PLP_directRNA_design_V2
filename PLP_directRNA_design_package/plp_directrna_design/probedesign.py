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
from tqdm import tqdm  
import matplotlib.pyplot as plt
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt

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
def evaluate_ligation_junction(targets, iupac_mismatches=None, plp_length=30):
    """
    Evaluates the ligation junction of a probe and introduces mismatches if needed.

    Args:
        targets (pd.DataFrame): DataFrame containing probe sequences in the 'Sequence' column.
        iupac_mismatches (str or list of tuples): Mismatch instructions in the form "5:R,10:G" or [(5, 'R'), (10, 'G')].
        plp_length (int): Length of the probe (default: 30).

    Returns:
        pd.DataFrame: Updated DataFrame with modified probes and ligation junction statuses.
    """
    # Ensure the DataFrame has a 'Ligation junction' column
    if 'Ligation junction' not in targets.columns:
        targets['Ligation junction'] = 'non-preferred'

    # Determine junction position (using the original probe indexing)
    junction_position = int((plp_length / 2) - 1)
    new_rows = []

    for idx in targets.index:
        probe_seq = targets.loc[idx]['Sequence']
        ligation_junction = probe_seq[junction_position] + probe_seq[junction_position + 2]
        ligation_status = ligation_junctions_dict.get(ligation_junction, "non-preferred")
        targets.loc[idx, 'Ligation junction'] = ligation_status

        if iupac_mismatches is not None:
            # If mismatches are provided as a string, parse them into a list of (pos, symbol) tuples.
            if isinstance(iupac_mismatches, str):
                iupac_mismatches = parse_iupac_mismatches(iupac_mismatches)
                
            # Limit to 2 mismatches
            if len(iupac_mismatches) > 2:
                raise InputValueError("The number of mismatches should be less than or equal to 2",
                                      field="iupac_mismatches", code="mismatches_exceed_limit")

            # Convert user-provided (1-indexed) positions to 0-indexed for internal use,
            # while preserving the original 1-indexed value for the probe ID suffix.
            # Each tuple becomes (adjusted_pos, original_pos, iupac_symbol)
            mismatches_converted = [(pos - 1, pos, symbol) for pos, symbol in iupac_mismatches]

            # Try all combinations of 1 or 2 mismatches
            for r in range(1, len(mismatches_converted) + 1):
                for subset in itertools.combinations(range(len(mismatches_converted)), r):
                    selected_mismatches = [mismatches_converted[i] for i in subset]
                    # For each mismatch, retrieve the possible replacement bases from IUPAC_CODES.
                    replacement_options = [IUPAC_CODES[symbol] for _, _, symbol in selected_mismatches]

                    # Generate all possible replacement combinations.
                    for replacement in itertools.product(*replacement_options):
                        original_seq_list = list(probe_seq)
                        modified_seq = original_seq_list.copy()
                        new_id_suffix = []
                        changes_made = False

                        for (adj_pos, orig_pos, iupac_symbol), new_base in zip(selected_mismatches, replacement):
                            original_base = original_seq_list[adj_pos]
                            # Apply the replacement only if it results in an actual change.
                            if original_base != new_base:
                                modified_seq[adj_pos] = new_base
                                new_id_suffix.append(f"{orig_pos}_{original_base}_{new_base}")
                                changes_made = True

                        # Only add a new probe row if at least one change occurred.
                        if changes_made:
                            new_probe_seq = "".join(modified_seq)
                            new_probe_id = f"{targets.loc[idx, 'Probe_id']}|{'_'.join(new_id_suffix)}"
                            
                            # Re-evaluate the ligation junction with the new sequence.
                            new_ligation_junction = new_probe_seq[junction_position] + new_probe_seq[junction_position + 2]
                            new_ligation_status = ligation_junctions_dict.get(new_ligation_junction, "non-preferred")

                            new_row = targets.loc[idx].copy()
                            new_row['Sequence'] = new_probe_seq
                            new_row['Ligation junction'] = new_ligation_status
                            new_row['Probe_id'] = new_probe_id

                            new_rows.append((new_probe_id, new_row))

    # Append the new rows to the DataFrame, if any.
    if new_rows:
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

## Extract mrna function:
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

def parse_attributes(attr_str):
    """
    Parse the attributes column from a GTF file.
    Example: 'transcript_id "TX1"; gene_id "G1"; ...'
    Returns a dictionary mapping keys to values.
    """
    attrs = {}
    for part in attr_str.strip().split(';'):
        part = part.strip()
        if not part:
            continue
        m = re.match(r'(\S+)\s+"(.+)"', part)
        if m:
            key, value = m.groups()
            attrs[key] = value
        else:
            pieces = part.split()
            if len(pieces) >= 2:
                key = pieces[0]
                value = pieces[1].strip('"')
                attrs[key] = value
    return attrs

def extract_mrna_sequences(fasta_file, gtf_file, output_file=None,
                           plus_strand_only=False, revcomp=False,
                           translate=False, codon_table=1,
                           alternative_start_codon=False,
                           clean_final_stop=False, clean_internal_stop=False,
                           verbose=False):
    """
    Extracts mRNA sequences from a FASTA file using exon records from a GTF file.
    
    The function assumes that the GTF file uses 1-indexed, inclusive coordinates.
    Each exon is extracted as: [start-1:end] (Python slicing).
    
    For both plus and negative strands, exons are first sorted in ascending order.
    For negative strand transcripts the merged sequence is then reverse complemented
    so that the output is in the 5'->3' orientation.
    
    Upstream/downstream extractions have been removed.
    
    Parameters:
      fasta_file (str): Path to the reference FASTA file.
      gtf_file (str): Path to the GTF file with exon annotations.
      output_file (str, optional): If provided, the output FASTA will be written here.
      plus_strand_only (bool): If True, output sequence in plus strand orientation.
      revcomp (bool): If True, force reverse complement of the final sequence.
      translate (bool): If True, translate the nucleotide sequence.
      codon_table (int): NCBI codon table ID (default 1).
      alternative_start_codon (bool): If True, force a methionine (M) at the start if valid.
      clean_final_stop (bool): If True, remove a trailing stop codon after translation.
      clean_internal_stop (bool): If True, replace internal stop codons with 'X'.
      verbose (bool): If True, print progress messages.
      
    Returns:
      List of SeqRecord objects (one per transcript).
    """
    # Index the genome FASTA for quick lookup.
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    # Add error handling for missing FASTA and gtf entries.
    if not genome:
        raise InputValueError(f"FASTA file not found: {fasta_file}", field="fasta_file", code="fasta_file_not_found")
    if not os.path.isfile(gtf_file):
        raise InputValueError(f"GTF file not found: {gtf_file}", field="gtf_file", code="gtf_file_not_found")
    

    # Group exon features by transcript_id.
    transcripts = {}  # transcript_id -> {'chrom': ..., 'strand': ..., 'gene_id': ..., 'exons': [(start, end), ...]}
    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            if feature.lower() != "exon":
                continue
            start = int(start)
            end = int(end)
            attr_dict = parse_attributes(attributes)
            transcript_id = attr_dict.get("transcript_id")
            gene_id = attr_dict.get("gene_id", "")
            gene_name = attr_dict.get("gene_name", "")
            if transcript_id is None:
                continue
            if transcript_id not in transcripts:
                transcripts[transcript_id] = {"chrom": chrom, "strand": strand, "gene_id": gene_id, "exons": [], "gene_name": gene_name}
            transcripts[transcript_id]["exons"].append((start, end))
    
    records = []
    for transcript_id, info in transcripts.items():
        chrom = info["chrom"]
        strand = info["strand"]
        gene_id = info["gene_id"]
        exons = info["exons"]
        gene_name = info["gene_name"]
        if chrom not in genome:
            if verbose:
                print(f"Warning: Chromosome {chrom} not found in FASTA for transcript {transcript_id}.")
            continue
        chrom_seq = genome[chrom].seq

        # Always sort exons in ascending order.
        exons_sorted = sorted(exons, key=lambda x: x[0])
        # Extract each exon using 1-indexed, inclusive conversion:
        # Python slice [s-1:e] returns bases s to e (inclusive).
        exon_seqs = [chrom_seq[s-1:e] for s, e in exons_sorted]
        merged_seq = Seq("").join(exon_seqs)
        
        # For negative strand, reverse complement the merged sequence (unless forced to plus strand).
        if strand == "-" and not plus_strand_only:
            merged_seq = merged_seq.reverse_complement()
        # Additionally, if revcomp is set, always reverse complement.
        if revcomp:
            merged_seq = merged_seq.reverse_complement()
        
        final_seq = merged_seq
        
        record_description = f"gene={gene_id} gene_name={gene_name} seq_id={chrom} type=mrna"
        
        # Optional translation.
        if translate:
            prot_seq = final_seq.translate(table=codon_table, to_stop=False)
            if alternative_start_codon:
                table_obj = CodonTable.unambiguous_dna_by_id[codon_table]
                start_codon = str(final_seq[0:3])
                if start_codon in table_obj.start_codons and prot_seq[0] != "M":
                    prot_seq = "M" + str(prot_seq)[1:]
            if clean_final_stop and str(prot_seq).endswith("*"):
                prot_seq = prot_seq[:-1]
            if clean_internal_stop:
                if str(prot_seq).endswith("*"):
                    prot_seq = prot_seq[:-1].replace("*", "X") + "*"
                else:
                    prot_seq = str(prot_seq).replace("*", "X")
                prot_seq = Seq(prot_seq)
            final_seq = prot_seq
            record_description += " translated"
        
        record = SeqRecord(final_seq, id=transcript_id, description=record_description)
        records.append(record)
    
    if output_file:
        with open(output_file, "w") as out_handle:
            SeqIO.write(records, out_handle, "fasta")
        if verbose:
            print(f"Wrote {len(records)} transcripts to {output_file}")
    
    return records

## Find target functions:
def find_targets_deprecated(selected_features, fasta_file, plp_length, min_coverage, output_file='Candidate_probes.txt', gc_min=50, gc_max=65, num_probes=10, iupac_mismatches=None):
    """
    Function to extract target sequences fulfilling the following criteria:
    Adapted from sequence developed by Sergio 
    Args:
        selected_features (str): Path to the selected features file (TSV format).
        fasta_file (str): Path to the indexed FASTA file.
        output_file (str): Path to the output file.
        plp_length (int): Probe length (default: 30).
        min_coverage (int): Minimum coverage of the region.
        gc_min (int): Minimum GC content (default: 50).
        gc_max (int): Maximum GC content (default: 65).
        num_probes (int): Number of probes to select per gene (default: 10).
        iupac_mismatches (position:base): 
            List of positions (1-based) and IUPAC codes to introduce mismatches.        
            Recommended positions for mismatches are **15 and 16**.
            Example: 15:R,16:G
            **Note:** The total number of mismatches must be ≤2. Exceeding this limit may lead to unexpected behavior.

    Returns:
        DataFrame: DataFrame with extracted probe sequences.
    """
    # Create a dataframe to store the targets
    targets = pd.DataFrame(columns=['Probe_id', 'Gene', 'Region', 'Sequence', 'GC', 'Coverage', 'Transcript_id'])


    # check if size of the probe is an even number
    if plp_length % 2 != 0:
        raise InputValueError("The size of the probe should be an even number", field="plp_length", code="odd_value_for_plp_length_provided")
    
        
    # load the fasta file
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))  

    # load the selected features table with gene name and region
    selected_features= pd.read_csv(selected_features, sep='\t')
    selected_features.index = selected_features['region']
    
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
                "Coverage": coverage_value,
                "Transcript_id": selected_features.loc[region, 'transcript_id']
            })
    targets_df = pd.DataFrame(targets)
    targets_df = evaluate_ligation_junction(targets_df, iupac_mismatches=iupac_mismatches, plp_length=plp_length)

    # Remove non-preferred ligation junctions
    targets_df = select_top_probes(targets_df, num_probes)
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
    # Remove probes non-preferred ligation junctions
    df_sorted = df_sorted[df_sorted['Ligation junction'] != 'non-preferred']
    
    # Select top `num_probes`
    final_df = df_sorted.groupby("Gene").head(num_probes)
    
    return final_df

## Check specificity functions:
import sys
import os
import pandas as pd
import dnaio
from Bio import SeqIO
from Bio.Seq import Seq
from cutadapt.adapters import BackAdapter
from dataclasses import dataclass
from io import StringIO
import csv

@dataclass
class Alignment:
    query_name: str
    target_name: str
    target_start: int
    target_end: int
    mismatches: int
    query_sequence: str
    target_sequence: str

def find_probes_in_targets(targets_df, reference_fasta, max_errors=1, output_file=None):
    """
    Function to check specificity of extracted probes against a reference genome.
    Uses Cutadapt's aligner for finding target regions.

    Args:
        targets_df (DataFrame): DataFrame containing extracted probe sequences.
        reference_fasta (str): Path to the reference genome FASTA.
        max_errors (int): Maximum number of allowed mismatches.
        output_file (str, optional): Path to save the results as a CSV file.

    Returns:
        DataFrame: A DataFrame containing matched probe alignments.
    """
    results = []

    with dnaio.open(reference_fasta) as references:
        references = list(references)  # Load references into a list
        total_probes = len(targets_df)  # Get total probe count

        with tqdm(total=total_probes, desc="Testing probes for specificity", unit=" alignments") as pbar:
            for reference_record in references:
                ref_id = reference_record.id
                ref_seq = reference_record.sequence

                for _, row in targets_df.iterrows():
                    probe_id = row["Probe_id"]
                    probe_seq = row["Sequence"]

                    adapter = BackAdapter(probe_seq, max_errors=max_errors, min_overlap=len(probe_seq), indels=False)
                    aligner = adapter.aligner

                    # Forward strand search
                    for t_start, t_end, errors, target_seq in find_all(ref_seq, aligner):
                        results.append(Alignment(
                            query_name=probe_id,
                            target_name=ref_id,
                            target_start=t_start + 1,  # Convert to 1-based index
                            target_end=t_end,
                            mismatches=errors,
                            query_sequence=probe_seq,
                            target_sequence=target_seq
                        ))
                    pbar.update(1)  # Update progress bar

                    # Reverse complement search
                    rev_ref_seq = str(Seq(ref_seq).reverse_complement())
                    adapter = BackAdapter(probe_seq, max_errors=max_errors, min_overlap=len(probe_seq), indels=False)
                    aligner = adapter.aligner

                    for t_start, t_end, errors, target_seq in find_all(rev_ref_seq, aligner):
                        results.append(Alignment(
                            query_name=probe_id,
                            target_name=f"{ref_id}(reverse)",
                            target_start=t_start + 1,
                            target_end=t_end,
                            mismatches=errors,
                            query_sequence=probe_seq,
                            target_sequence=target_seq
                        ))
                    pbar.update(1)  # Update progress bar again for reverse search

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Save results to file if requested
    if output_file:
        results_df.to_csv(output_file, index=False)
        sys.stderr.write(f"Results saved to {output_file}\n")

    return results_df


def find_all(ref, aligner):
    """Find all occurrences of a probe in a reference sequence."""
    offset = 0
    while True:
        result = aligner.locate(ref)
        if result is None:
            break
        ref_start, ref_end, query_start, query_end, score, errors = result
        t_start = query_start + offset
        t_end = query_end + offset
        target_seq = ref[query_start:query_end]
        yield (t_start, t_end, errors, target_seq)
        offset += query_start + 1
        ref = ref[query_start + 1:]

def find_targets(selected_features, fasta_file, reference_fasta, plp_length=30, min_coverage=1,
                 output_file='Candidate_probes', gc_min=50, gc_max=65, num_probes=10,
                 iupac_mismatches=None, max_errors=1, check_specificity=False):
    """
    Extract target sequences based on defined probe criteria and optionally check specificity.

    Args:
        selected_features (str): Path to the selected features file (TSV format).
        fasta_file (str): Path to the indexed FASTA file.
        output_file (str): Path to the output file.
        plp_length (int): Probe length (default: 30).
        min_coverage (int): Minimum coverage of the region.
        gc_min (int): Minimum GC content (default: 50).
        gc_max (int): Maximum GC content (default: 65).
        num_probes (int): Number of probes to select per gene (default: 10).
        iupac_mismatches (position:base): 
            List of positions (1-based) and IUPAC codes to introduce mismatches.        
            Recommended positions for mismatches are **15 and 16**.
            Example: 15:R,16:G
            **Note:** The total number of mismatches must be ≤2. Exceeding this limit may lead to unexpected behavior.

        max_errors (int): Maximum mismatches allowed during specificity checking.
        check_specificity (bool): Whether to check probe specificity against reference.

    Returns:
        DataFrame: DataFrame with extracted probe sequences (and specificity results if enabled).
    """

    targets = []

    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    selected_features = pd.read_csv(selected_features, sep='\t')
    selected_features.index = selected_features['region']
    coverage_dict = dict(zip(selected_features['region'], selected_features['coverage']))

    if plp_length % 2 != 0:
        raise InputValueError("The size of the probe should be an even number", field="plp_length", code="odd_value_for_plp_length_provided")


    for keys in seq_dict:
        gene, region = keys.split("|")

        chr_name, start_end = region.split(":")
        initial_start = int(start_end.split("-")[0])
        seq = seq_dict[keys].seq
        seq_len = len(seq) - (plp_length - 1)

        coverage_value = coverage_dict.get(region, 0)
        if coverage_value < min_coverage:
            continue 

        for i in range(0, seq_len):
            tmp_seq = seq[i:i + plp_length]
            gc_content = (tmp_seq.count("G") + tmp_seq.count("C")) / len(tmp_seq) * 100

            if not (gc_min <= gc_content <= gc_max):
                continue
            if not any(nucleotide * 3 in tmp_seq for nucleotide in "ACGT"):
                continue

            start = initial_start + i
            end = start + plp_length - 1

            targets.append({
                "Probe_id": f"{gene}|{chr_name}:{start}-{end}",
                "Gene": gene,
                "Region": f"{chr_name}:{start}-{end}",
                "Sequence": str(tmp_seq),
                "GC": gc_content,
                "Coverage": coverage_value,
                "Transcript_id": selected_features.loc[region, 'transcript_id']
            })

    targets_df = pd.DataFrame(targets)

    # Introduction of IUPAC mismatches
    if iupac_mismatches:
        targets_df = evaluate_ligation_junction(targets_df, iupac_mismatches=iupac_mismatches, plp_length=plp_length)


    # Check probe specificity against reference genome if requested
    if check_specificity:
        specificity_results = find_probes_in_targets(targets_df, reference_fasta, max_errors, output_file=f"{output_file}_specificity.csv")

        # Identify probes that have **only one unique match**
        unique_queries = specificity_results.groupby('query_name').size().reset_index(name='counts')
        unique_queries = unique_queries[unique_queries['counts'] == 1]['query_name']

        # Select only probes with a unique match
        targets_df = targets_df[targets_df['Probe_id'].isin(unique_queries)]

        # Merge mismatches into targets_df
        targets_df = targets_df.merge(
            specificity_results[['query_name', 'mismatches']],
            left_on='Probe_id',
            right_on='query_name',
            how='left'
        ).drop(columns=['query_name'])  # Drop duplicate column

        print(f"✅ Specificity results saved to {output_file}_specificity.csv")


    # Debugging
#    print("Debugging:")
#    print(targets_df.shape[0], "probes extracted.")
#    targets_df = select_top_probes(targets_df, num_probes)
#    print(targets_df)
    return targets_df


def parse_specificity_results(specificity_results):
    """
    Parse the specificity results from the CSV file.

    Args:
        specificity_results (str): Path to the specificity results CSV file.

    Returns:
        DataFrame: DataFrame containing parsed specificity results.
    """
    specificity_results.group_by('Probe_id').size().reset_index(name='counts')
    return pd.read_csv(specificity_results)

def filter_probes_by_distance(group, min_dist_probes):
    """
    Filters probes based on a minimum distance between them.
    This function takes a DataFrame containing probe information and filters out probes
    that are too close to each other based on a specified minimum distance. The 'Region'
    column in the DataFrame should contain coordinates in the format 'chr:start-end'.
    Args:
        group (pd.DataFrame): A DataFrame containing probe information with a 'Region' column.
        min_dist_probes (int): The minimum distance required between probes.
    Returns:
        pd.DataFrame: A DataFrame containing the filtered probes that meet the minimum distance requirement.
    """

    # Extract start and end coordinates from the 'Region' column
    group['Start'] = group['Region'].apply(lambda x: int(x.split(':')[1].split('-')[0]))
    group['End'] = group['Region'].apply(lambda x: int(x.split(':')[1].split('-')[1]))
    
    # Sort by start coordinate
    group = group.sort_values(by='Start')
    
    # Filter probes based on the minimum distance
    filtered_probes = []
    last_end = -min_dist_probes  # Initialize with a value that ensures the first probe is included
    
    for idx, row in group.iterrows():
        if row['Start'] - last_end >= min_dist_probes:
            filtered_probes.append(row)
            last_end = row['End']
    filtered_probes = pd.DataFrame(filtered_probes)        
    filtered_probes = filtered_probes.drop(columns=['Start', 'End'])
    return pd.DataFrame(filtered_probes)


## Calculating weighted melting temperature functions:
def weighted_tm_gc_scoring(sequence, Tm_oligo, Tm_min=55.0, Tm_opt=60.0, Tm_max=65.0, 
                           GC_min=40.0, GC_opt=50.0, GC_max=60.0, w_Tm=1.0, w_GC=1.0):
    """Compute the weighted score for a given oligo based on melting temperature and GC content."""
    GC_oligo = gc_fraction(sequence) * 100  # Convert fraction to percentage
    Tm_opt = (Tm_max + Tm_min) / 2
    # Compute Tm deviation score
    if Tm_oligo >= Tm_opt:
        score_Tm = abs(Tm_oligo - Tm_opt) / (Tm_max - Tm_opt)
    else:
        score_Tm = abs(Tm_oligo - Tm_opt) / (Tm_opt - Tm_min)
    
    # Compute GC deviation score
    if GC_oligo >= GC_opt:
        score_GC = abs(GC_oligo - GC_opt) / (GC_max - GC_opt)
    else:
        score_GC = abs(GC_oligo - GC_opt) / (GC_opt - GC_min)
    
    return w_Tm * score_Tm + w_GC * score_GC  # Lower score is better


def analyze_scores(scores, percentile=5):
    """Analyze the score distribution and suggest a cutoff."""
    suggested_cutoff = np.percentile(scores, percentile)  # Get the threshold for top X%
    return suggested_cutoff


def score_padlock_probe(sequence, Tm_min=55.0, Tm_opt=60.0, Tm_max=65.0, 
                        GC_min=40.0, GC_opt=50.0, GC_max=60.0, w_Tm=1.0, w_GC=1.0, percentile=5):
    """Score a padlock probe by splitting it into two arms and analyzing distribution."""
    mid = len(sequence) // 2
    left_arm, right_arm = sequence[:mid], sequence[mid:]

    # Compute Tm for each arm
    Tm_left = mt.Tm_NN(left_arm)
    Tm_right = mt.Tm_NN(right_arm)

    # Score each arm separately
    score_left = weighted_tm_gc_scoring(left_arm, Tm_left, Tm_min, Tm_opt, Tm_max, GC_min, GC_opt, GC_max, w_Tm, w_GC)
    score_right = weighted_tm_gc_scoring(right_arm, Tm_right, Tm_min, Tm_opt, Tm_max, GC_min, GC_opt, GC_max, w_Tm, w_GC)

    # Combine scores (you can take the average, max, or another approach)
    final_score = (score_left + score_right) / 2  # Averaging the scores

    return final_score
def visualize_score_distribution(scores, cutoff):
    """Plot the score distribution and cutoff."""
    plt.figure(figsize=(8,5))
    plt.hist(scores, bins=30, alpha=0.7, color='blue', edgecolor='black')
    plt.axvline(cutoff, color='red', linestyle='dashed', linewidth=2, label=f'Cutoff ({cutoff:.2f})')
    plt.xlabel("Score")
    plt.ylabel("Frequency")
    plt.title("Distribution of Probe Scores")
    plt.legend()
    plt.show()
