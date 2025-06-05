from plp_directrna_design import cli_utils as cli


if __name__ == "__main__":
    parser = cli.extract_features_parser()
    args = parser.parse_args()
    cli.extract_features(
        args.gtf, args.output, args.genes, args.identifier_type, args.gene_feature
    )
