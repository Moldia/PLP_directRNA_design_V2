import os
import pandas as pd
import argparse
import subprocess
import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq

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

def main(gtf_output, fasta_file, output_fasta, plp_length, identifier_type):
    """
    Main function to extract CDS sequences from a pre-processed GTF output.

    Args:
        gtf_output (str): Path to the GTF-derived DataFrame (TSV format).
        fasta_file (str): Path to the indexed FASTA file.
        output_fasta (str): Path to the output FASTA file.
    """
    print(f"🔹 Loading GTF output from {gtf_output}...")
    df = pd.read_csv(gtf_output, sep='\t')

    # Ensure FASTA index exists
    check_fasta_index(fasta_file)

    # Save regions for fast retrieval
    regions_file = "regions"
    save_regions_for_faidx(df, regions_file, plp_length, identifier_type=identifier_type)

    # Extract sequences
    extract_sequences(fasta_file, regions_file + ".txt", output_fasta, df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract CDS sequences from an indexed FASTA file")

    parser.add_argument("--gtf_output", required=True, help="Path to the GTF-derived TSV file (CDS regions)")
    parser.add_argument("--fasta", required=True, help="Path to the indexed FASTA file")
    parser.add_argument("--output_fasta", required=True, help="Path to the output FASTA file")
    parser.add_argument("--identifier_type", default='gene_id', choices=['gene_id', 'gene_name'], help="Type of identifier provided ('gene_id' or 'gene_name')")
    parser.add_argument("--plp_length", default=30, help="probe length")
    args = parser.parse_args()

    main(args.gtf_output, args.fasta, args.output_fasta, args.plp_length, args.identifier_type)