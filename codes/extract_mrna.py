import sys
import argparse
from plp_directrna_design import probedesign as plp
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description="Extract mRNA sequences using a GTF and FASTA file."
    )
    parser.add_argument("--gtf", required=True, help="Path to the GTF file")
    parser.add_argument("--fasta", required=True, help="Path to the FASTA file")
    parser.add_argument("--output_file", required=True, help="Path to the output FASTA file")
    args = parser.parse_args()
    
    # Call the extraction function from your package.
    # You can adjust the options below as needed.
    records = plp.extract_mrna_sequences(
        fasta_file=args.fasta,
        gtf_file=args.gtf,
        output_file=args.output_file,
        plus_strand_only=False,
        revcomp=False,
        translate=False,
        codon_table=1,
        alternative_start_codon=True,
        clean_final_stop=True,
        clean_internal_stop=False,
        verbose=False
    )

    print(f"Extracted mRNA sequences have been saved to {args.output_file}")
    
if __name__ == "__main__":
    main()
