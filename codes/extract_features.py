
from plp_directrna_design import probedesign as plp

import argparse

def main(gtf_file, output_file, genes_str=None, identifier_type='gene_id', gene_feature='CDS'):
 
    # Parse the GTF file and filter by gene list
    gtf_df, genes_of_interest = plp.parse_gtf(gtf_file, genes_str, identifier_type, gene_feature)

    # Merge regions and calculate coverage
    merged_cov_df = plp.merge_regions_and_coverage(genes_of_interest, gtf_df)

    # Write the merged results to an output file
    merged_cov_df.to_csv(output_file, sep='\t', index=False)

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse GTF file and calculate coverage")

    parser.add_argument("--gtf", required=True, help="Path to the GTF file")
    parser.add_argument("--output", required=True, help="Path to the output file")
    parser.add_argument("--genes", default=None, help="Comma-separated list of gene IDs or names to filter")
    parser.add_argument("--identifier_type", default='gene_id', choices=['gene_id', 'gene_name'], help="Type of identifier provided ('gene_id' or 'gene_name')")
    parser.add_argument("--gene_feature", default='CDS', help="Feature type to extract (default = 'CDS')")
    args = parser.parse_args()
    main(args.gtf, args.output, args.genes, args.identifier_type, args.gene_feature)