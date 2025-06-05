from plp_directrna_design import cli_utils as cli


if __name__ == "__main__":
    parser = cli.extract_sequences_parser()
    args = parser.parse_args()
    cli.extract_sequences(
        args.gtf_output,
        args.fasta,
        args.output_fasta,
        args.plp_length,
        args.identifier_type,
    )
