from plp_directrna_design import cli_utils as cli


if __name__ == "__main__":
    parser = cli.find_target_parser()
    args = parser.parse_args()
    cli.find_target(
        args.selected_features,
        args.fasta_file,
        args.output_file,
        args.reference_fasta,
        args.min_coverage,
        args.gc_min,
        args.gc_max,
        args.num_probes,
        args.iupac_mismatches,
        args.max_errors,
        args.check_specificity,
        args.plp_length,
    )
