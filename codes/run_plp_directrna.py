from plp_directrna_design import cli_utils as cli
import argparse

def main():
    # Build a master parser that includes all sub-steps
    parser = master_parser()
    args = parser.parse_args()
    
    # Step 1: Extract features
    print("EXTRACTING FEATURES")
    cli.extract_features(
        args.gtf, args.features_output, args.genes, args.identifier_type, args.gene_feature, 
    )

    # Step 2: Extract mRNA
    print("EXTRACTING mRNA")
    cli.extract_mrna(args.fasta, args.gtf, args.transcriptome_output)

    # Step 3: Extract sequences
    print("EXTRACTING SEQUENCES")
    cli.extract_sequences(
        args.features_output,
        args.fasta,
        args.sequences_output,
        args.plp_length,
        args.identifier_type,
    )

    # # Step 4: Find targets
    print("FINDING TARGETS")
    if args.iupac_mismatches == "None":
        args.iupac_mismatches = None

    cli.find_target(
        selected_features=args.features_output,
        sequences_output=args.sequences_output,
        output_file=args.targets_output,              # ← probes table
        reference_fasta=args.sequences_output,        # ← transcriptome sequences
        min_coverage=args.min_coverage,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        num_probes=args.num_probes,
        iupac_mismatches=args.iupac_mismatches,
        max_errors=args.max_errors,
        check_specificity=args.check_specificity,
        plp_length=args.plp_length,
        Tm_min=args.Tm_min,
        Tm_max=args.Tm_max,
        lowest_percentile_Tm_score_cutoff=args.lowest_percentile_Tm_score_cutoff,
        min_dist_probes=args.min_dist_probes,
        filter_ligation_junction=args.filter_ligation_junction,
        off_target_output=args.off_target_output
    )


def master_parser():
    parser = argparse.ArgumentParser(description="PLP DirectRNA Design Workflow")
    
    # Shared and individual inputs
    parser.add_argument("--gtf", required=True, help="Path to the GTF file")
    parser.add_argument("--genes", type=cli.parse_genes, required=True, help="Comma-separated list of gene IDs or names to filter")
    parser.add_argument("--identifier_type", default="gene_name", choices=["gene_id", "gene_name"], help="Type of identifier provided")
    parser.add_argument("--gene_feature", default="CDS", help="Feature type to extract")
    parser.add_argument("--fasta", required=True, help="Path to the FASTA file")

    # Outputs
    parser.add_argument("--features_output", default="extract_features_output.txt", help="Path to features output file")
    parser.add_argument("--transcriptome_output", default="data/transcriptome_out.fa", help="Path to transcriptome output file")
    parser.add_argument("--sequences_output", default="extract_seqs_output.fa", help="Path to extracted sequences output file (transcriptome)")
    parser.add_argument("--targets_output", default="targets.txt", help="Path to final probes sequences table output file")

    # Find targets specific args
    parser.add_argument("--min_coverage", type=int, default=1, help="Minimum coverage threshold")
    parser.add_argument("--gc_min", type=int, default=50, help="Minimum GC content percentage")
    parser.add_argument("--gc_max", type=int, default=65, help="Maximum GC content percentage")
    parser.add_argument("--num_probes", type=int, default=10, help="Number of probes to generate")
    parser.add_argument("--iupac_mismatches", default=None, help="IUPAC mismatches parameter")
    parser.add_argument("--max_errors", type=int, default=1, help="Maximum number of errors allowed")
    parser.add_argument("--check_specificity", action="store_true", help="Enable specificity checking")
    parser.add_argument("--plp_length", type=int, default=30, help="PLP length for probe design")
    parser.add_argument("--off_target_output", action="store_true", help="Enable saving of off-target output")
    parser.add_argument("--Tm_min", type=int, default=50, help="Minimum Tm for probes")
    parser.add_argument("--Tm_max", type=int, default=70, help="Maximum Tm for probes")
    parser.add_argument("--lowest_percentile_Tm_score_cutoff", type=int, default=10, help="Lowest percentile Tm score cutoff")
    parser.add_argument("--min_dist_probes", type=int, default=100, help="Minimum distance between probes")
    parser.add_argument("--filter_ligation_junction", action="store_true", help="Enable filtering of ligation junctions")

    return parser

if __name__ == "__main__":
    main()
