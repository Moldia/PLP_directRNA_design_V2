from plp_directrna_design import probedesign as plp
import pandas as pd
import argparse

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
    plp.check_fasta_index(fasta_file)

    # Save regions for fast retrieval
    regions_file = "regions"
    plp.save_regions_for_faidx(df, regions_file, plp_length, identifier_type=identifier_type)

    # Extract sequences
    plp.extract_sequences(fasta_file, regions_file + ".txt", output_fasta, df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract CDS sequences from an indexed FASTA file")

    parser.add_argument("--gtf_output", required=True, help="Path to the GTF-derived TSV file (CDS regions)")
    parser.add_argument("--fasta", required=True, help="Path to the indexed FASTA file")
    parser.add_argument("--output_fasta", required=True, help="Path to the output FASTA file")
    parser.add_argument("--identifier_type", default='gene_id', choices=['gene_id', 'gene_name'], help="Type of identifier provided ('gene_id' or 'gene_name')")
    parser.add_argument("--plp_length", default=30, help="probe length")
    args = parser.parse_args()

    main(args.gtf_output, args.fasta, args.output_fasta, args.plp_length, args.identifier_type)