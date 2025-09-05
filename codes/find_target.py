from plp_directrna_design import cli_utils as cli


if __name__ == "__main__":
    parser = cli.find_targets_parser()
    args = parser.parse_args()
    cli.find_targets(
        **vars(args)
    )
