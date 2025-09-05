from plp_directrna_design import cli_utils as cli


if __name__ == "__main__":
    parser = cli.extract_features_parser()
    args = parser.parse_args()
    cli.extract_features(
        **vars(args)
    )
