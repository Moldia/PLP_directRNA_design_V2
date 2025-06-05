from plp_directrna_design import cli_utils as cli


if __name__ == "__main__":
    parser = cli.extract_mrna_parser()
    args = parser.parse_args()
    cli.extract_mrna(args.fasta, args.gtf, args.output_file)
