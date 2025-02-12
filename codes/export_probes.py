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

def export_probes(selected_features, fasta_file, plp_length, min_coverage, output_file='Candidate_probes.txt', gc_min=50, gc_max=65, num_probes=10):
    """
    Function to extract probe sequences fulfilling the following criteria:
    Adapted from sequence developed by Sergio 
    Args:
        selected_features (str): Path to the selected features file (TSV format).
        fasta_file (str): Path to the indexed FASTA file.
        output_file (str): Path to the output file.
        plp_length (int): Minimum probe length.
        min_coverage (int): Minimum coverage of the region.
        gc_min (int): Minimum GC content (default: 50).
        gc_max (int): Maximum GC content (default: 65).
    Returns:
        DataFrame: DataFrame with extracted probe sequences.
    """
    # Create a dataframe to store the targets
    targets = pd.DataFrame(columns=['Probe_id', 'Gene', 'Region', 'Sequence', 'GC', 'Ligation junction', 'Coverage'])

    # check if size of the probe is an even number
    if plp_length % 2 != 0:
        raise ValueError("The size of the probe should be an even number")
    
    # Ligation junction dictionary
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
                # extract the ligations junction
            ligation_junction = tmp_seq[int((plp_length/2)-1)] + tmp_seq[int((plp_length/2)+1)]
            ligation_status = ligation_junctions_dict.get(ligation_junction, "non-preferred")
            
            if ligation_status == "non-preferred":
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
                "Ligation junction": ligation_status,
                "Coverage": coverage_value
            })
    targets_df = pd.DataFrame(targets)
    targets_df = select_top_probes(targets_df, num_probes)
    targets_df.to_csv(output_file, sep='\t', index=False)
    create_fasta(targets_df, output_file)
    print(f"✅ Final sequences saved to {output_file} and fasta file {output_file}.fa .")

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

def main(selected_features, fasta_file, output_file, plp_length, min_coverage, gc_min=50, gc_max=65, num_probes=10):
    """
    Main function for probe extraction.
    """
    print(f"🔹 Loading selected features from {selected_features}...")
#    selected_features = pd.read_csv(selected_features, sep='\t')

    export_probes(selected_features, fasta_file, plp_length, min_coverage, output_file, gc_min, gc_max, num_probes)


if __name__ == "__main__":
    import argparse

parser = argparse.ArgumentParser(
    description=(
        "Extracts probe sequences fulfilling the following criteria:\n"
        "\n"
        "- GC content between 50-65% (default).\n"
        "- Ligation junctions must be 'preferred' or 'neutral' (not 'non-preferred').\n"
        "  See Xenium Custom Panel Design guide:\n"
        "  https://cdn.10xgenomics.com/image/upload/v1716400584/support-documents/CG000683_TechNote_Xenium_Custom_Panel_Design_RevD.pdf\n"
        "\n"
        "  Ligation junction preferences:\n"
        "    Preferred:      AT, TA, GA, AG\n"
        "    Neutral:        TT, CT, CA, TC, AC, CC, TG, AA\n"
        "    Non-Preferred:  CG, GT, GG, GC (filtered out)\n"
        "\n"
        "- No homopolymers of length 3 or more (e.g., AAA, TTT, GGG, CCC).\n"
        "- Minimum coverage of the region is met (default: 1, based on CDS/exon overlap).\n"
        "- The probe size must be an even number (default: 30)."
        ),
        formatter_class=argparse.RawTextHelpFormatter  # Ensures multiline formatting
    )

#parser.parse_args(["--help"])  # Simulate --help call for testing

parser.add_argument("--selected_features", required=True, help="Path to the selected features file (TSV format)")
parser.add_argument("--fasta_file", required=True, help="Path to the extracted sequences file (CDS/exons)")
parser.add_argument("--output_file", required=True, help="Path to the output file")
parser.add_argument("--plp_length", default=30, type=int, help="Minimum probe length")
parser.add_argument("--min_coverage", default=1, type=int, help="Minimum coverage of the region")
parser.add_argument("--gc_min", default=50, type=int, help="Minimum GC content")
parser.add_argument("--gc_max", default=65, type=int, help="Maximum GC content")
parser.add_argument("--num_probes", default=10, type=int, help="Number of probes to select")
args = parser.parse_args()
main(args.selected_features, args.fasta_file, args.output_file, args.plp_length, 
     args.min_coverage, args.gc_min, args.gc_max, args.num_probes)



