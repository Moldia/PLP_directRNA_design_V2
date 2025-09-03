from plp_directrna_design import cli_utils as cli
import argparse


def main():
    # Build a master parser that includes all sub-steps
    parser = cli.run_plp_directrna_parser()
    args = parser.parse_args()
    cli.run_plp_directrna(
        args.gtf,
        args.features_output,
        args.genes,
        args.identifier_type,
        args.fasta, 
        args.transcriptome_output,
        args.sequences_output,
        args.plp_length,
        args.targets_output,
        args.min_coverage,
        args.gc_min,
        args.gc_max,
        args.num_probes,
        args.iupac_mismatches,
        args.max_errors,
        args.check_specificity,
        args.Tm_min,
        args.Tm_max,
        args.lowest_percentile_Tm_score_cutoff,
        args.min_dist_probes,
        args.filter_ligation_junction,
        args.off_target_output,
    )


if __name__ == "__main__":
    main()
