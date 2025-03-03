
from plp_directrna_design import probedesign as plp
import pandas as pd
import argparse

if __name__ == "__main__":
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

def main(selected_features, fasta_file, output_file, plp_length, min_coverage, gc_min=50, gc_max=65, num_probes=10, iupac_mismatches=None):
    """
    Main function for probe extraction.
    """
    print(f"🔹 Loading selected features from {selected_features}...")


    targets_df = plp.find_targets(selected_features = selected_features, fasta_file = fasta_file, 
                                 plp_length = plp_length, min_coverage = min_coverage, output_file=output_file, 
                                 gc_min=gc_min, gc_max=gc_max, num_probes=num_probes, iupac_mismatches=iupac_mismatches)
    targets_df.to_csv(output_file, sep='\t', index=False)

    

#parser.parse_args(["--help"])  # Simulate --help call for testing

parser.add_argument("--selected_features", required=True, help="Path to the selected features file (TSV format)")
parser.add_argument("--fasta_file", required=True, help="Path to the extracted sequences file (CDS/exons)")
parser.add_argument("--output_file", required=True, help="Path to the output file")
parser.add_argument("--plp_length", default=30, type=int, help="Minimum probe length")
parser.add_argument("--min_coverage", default=1, type=int, help="Minimum coverage of the region")
parser.add_argument("--gc_min", default=50, type=int, help="Minimum GC content")
parser.add_argument("--gc_max", default=65, type=int, help="Maximum GC content")
parser.add_argument("--num_probes", default=10, type=int, help="Number of probes to select")
parser.add_argument("--iupac_mismatches", default=None, help="IUPAC mismatches to consider. Note that the number of mismatches should be less than or equal to 2. Example: 5:R,6:A")
args = parser.parse_args()
main(args.selected_features, args.fasta_file, args.output_file, args.plp_length, 
     args.min_coverage, args.gc_min, args.gc_max, args.num_probes, args.iupac_mismatches)



