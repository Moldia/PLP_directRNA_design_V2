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

def main(selected_features, fasta_file, output_file, reference_fasta, min_coverage, 
         gc_min=50, gc_max=65, num_probes=10, iupac_mismatches=None, max_errors = 1, 
         check_specificity = False, plp_length=30, Tm_min=55, Tm_max=65, 
         lowest_percentile_Tm_score_cutoff=5, min_dist_probes=10,
         filter_ligation_junction=True):
    """
    Main function for probe extraction.
    """
    print(f"🔹 Loading selected features from {selected_features}...")

    print("DEBUG:", plp_length)
    targets_df = plp.find_targets(selected_features = selected_features, fasta_file = fasta_file, reference_fasta = reference_fasta,
                                 plp_length = plp_length, min_coverage = min_coverage, output_file=output_file, 
                                 gc_min=gc_min, gc_max=gc_max, num_probes=num_probes, iupac_mismatches=iupac_mismatches,
                                 max_errors=max_errors, check_specificity=check_specificity)
    # Calculate the melting temperature scores
    sequences = targets_df['Sequence']
    scores = [plp.score_padlock_probe(seq, Tm_min = Tm_min, Tm_max= Tm_max) for seq in sequences]
    targets_df['Melt_Tm_scores'] = scores
    # Calculate the suggested cutoff based on the 5th percentile
    suggested_cutoff = plp.analyze_scores(scores, percentile=lowest_percentile_Tm_score_cutoff)
    # Filter the targets based on the suggested cutoff
    targets_df = targets_df[targets_df['Melt_Tm_scores'] <= suggested_cutoff]
    # Filteer the probes based on the minimum distance between probes
    targets_df = plp.filter_probes_by_distance(targets_df, min_dist_probes=min_dist_probes)
    # filter the probes based on the ligation junction preferences
    if filter_ligation_junction:
        targets_df = targets_df[targets_df['Ligation junction'] != 'non-preferred']
    # Save the output    
    targets_df.to_csv(output_file, sep='\t', index=False)

parser.add_argument("--selected_features", required=True, help="Path to the selected features file (TSV format)")
parser.add_argument("--fasta_file", required=True, help="Path to the extracted sequences file (CDS/exons)")
parser.add_argument("--output_file", required=True, help="Path to the output file")
parser.add_argument("--reference_fasta", required=True, help="Path to the reference (genome/transcriptome) FASTA file; Please note that with genome reference, the probe design will be performed on the whole genome which may take a long time.")
parser.add_argument("--min_coverage", default=1, type=int, help="Minimum coverage of the region")
parser.add_argument("--gc_min", default=50, type=int, help="Minimum GC content")
parser.add_argument("--gc_max", default=65, type=int, help="Maximum GC content")
parser.add_argument("--num_probes", default=10, type=int, help="Number of probes to select")
parser.add_argument("--iupac_mismatches", default=None, help="IUPAC mismatches to consider. Note that the number of mismatches should be less than or equal to 2. Example: 5:R,6:A")
parser.add_argument("--max_errors", default=1, type=float, help="Maximum error rate (or number of errors if an integer, recommended range is 0-6)")
parser.add_argument("--check_specificity", action="store_true", help="Check probe specificity")
parser.add_argument("--plp_length", default=30, type=int, help="Minimum probe length")
parser.add_argument("--Tm_min", default = 55.0, type=float, help="Minimum melting temperature for each arm")
parser.add_argument("--Tm_max", default = 65.0, type=float, help="Maximum melting temperature for each arm")
parser.add_argument("--lowest_percentile_Tm_score_cutoff", default = 5, type=int, help="The lowest percentile of the melting temperature score cutoff to filter")
parser.add_argument("--min_dist_probes", default = 10, type=int, help="Minimum distance between probes")
parser.add_argument("--filter_ligation_junction", action="store_true", help="Filter probes based on ligation junction preferences; exclude non-preferred")

args = parser.parse_args()

main(args.selected_features, args.fasta_file, args.output_file, args.reference_fasta, 
     args.min_coverage, args.gc_min, args.gc_max, args.num_probes, 
     args.iupac_mismatches, args.max_errors, args.check_specificity, args.plp_length)
