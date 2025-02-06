import pandas as pd
import numpy as np
import argparse
import re

def parse_gtf_to_dataframe(gtf_path: str) -> pd.DataFrame:
    """
    Parses a GTF (Gene Transfer Format) file into a structured Pandas DataFrame.

    Args:
        gtf_path (str): Path to the GTF file.

    Returns:
        pd.DataFrame: A DataFrame containing parsed GTF data.
    
    Reference:
        Adapted from Ricardo Filipe dos Santos script:
        https://gist.github.com/rf-santos/22f521c62ca2f85ac9582bf0d91e4054
    """
    print('Loading GTF file....')
    # Define column names
    col_names = [
        'seqname', 'source', 'feature', 'start', 'end', 'score', 
        'strand', 'frame', 'attribute'
    ]
    
    # Read the GTF file, skipping comments
    df = pd.read_csv(gtf_path, sep='\t', header=None, names=col_names, comment='#', dtype={'score': str})

    # Define attributes to extract
    attributes = [
        'gene_id', 'transcript_id', 'exon_number', 'gene_name', 'gene_source', 
        'gene_biotype', 'transcript_name', 'transcript_source', 'transcript_biotype', 
        'protein_id', 'exon_id', 'tag'
    ]

    # Create regex patterns for attribute extraction
    attr_patterns = {attr: rf'{attr} "([^"]*)"' for attr in attributes}

    # Extract attributes using vectorized operations
    for attr, pattern in attr_patterns.items():
        df[attr] = df['attribute'].str.extract(pattern)

    # Drop the original attribute column
    df.drop(columns=['attribute'], inplace=True)

    return df


def parse_gtf(gtf_file, genes_of_interest=None, identifier_type='gene_id', gene_feature='CDS'):
    """
    Parses a GTF file and yields data only for the specified genes_of_interest.

    Args:
        gtf_file (str): Path to the GTF file.
        genes_of_interest (set or None): A set of gene IDs or names to parse.
                                         If None, parse all genes in the GTF.
        identifier_type (str): Type of identifier provided ('gene_id' or 'gene_name').

    Returns:
        pd.DataFrame: A DataFrame containing parsed GTF data.
    """
    print('Parsing GTF file....')
    # Read the GTF file into a DataFrame
    gtf_df = parse_gtf_to_dataframe(gtf_file)

    # Check if a gene list is provided and filter accordingly
    if genes_of_interest:
        if identifier_type == 'gene_id':
            gtf_df = gtf_df[gtf_df['gene_id'].str.lower().isin(genes_of_interest) & (gtf_df['feature'] == gene_feature)]
        elif identifier_type == 'gene_name':
            gtf_df = gtf_df[gtf_df['gene_name'].str.lower().isin(genes_of_interest)& (gtf_df['feature'] == gene_feature)]
        else:
            raise ValueError("identifier_type must be 'gene_id' or 'gene_name'")

    return gtf_df


def merge_regions_and_coverage(genes_of_interest, gtf_df):
    """
    Merge CDS regions and computes average coverage for the specified genes.

    Args:
        genes_of_interest (set): A set of gene names to process.
        gtf_df (pd.DataFrame): A DataFrame containing parsed GTF data.

    Returns:
        pd.DataFrame: A DataFrame containing merged CDS regions and average coverage.
    """
    print('Merge regions and calculating coverage....')
    # Initialize an empty list to store results
    merged_regions = []

    for gn in genes_of_interest:
        # Subset the DataFrame for the given gene_name
        isoforms = gtf_df[gtf_df['gene_name'].str.lower() == gn].copy()
        
        # Sort by chromosome, strand, and start position
        isoforms = isoforms.sort_values(by=['seqname', 'strand', 'start'])

        # Determine the full genomic range
        min_start = isoforms['start'].min()
        max_end = isoforms['end'].max()

        # Create a NumPy array to track coverage over this genomic range
        coverage_array = np.zeros(max_end - min_start + 1, dtype=int)

        # Dictionary to track merged CDS regions
        merged = []
        
        for _, row in isoforms.iterrows():
            if not merged:
                merged.append(row.to_dict())  # First entry, add as-is
            else:
                prev = merged[-1]

                if row['start'] <= prev['end']:  # Overlapping region
                    merged[-1]['end'] = max(prev['end'], row['end'])  # Expand region
                    merged[-1]['transcript_id'] += ";" + row['transcript_id']  # Append transcript_id
                else:
                    merged.append(row.to_dict())  # Start new region

            # Update coverage using NumPy slicing (avoids looping over each position)
            coverage_array[row['start'] - min_start : row['end'] - min_start + 1] += 1

        # Convert merged results to a DataFrame
        merged_df = pd.DataFrame(merged)

        # Compute the average coverage for each merged region efficiently
        coverage_values = []
        for _, region in merged_df.iterrows():
            region_slice = coverage_array[region['start'] - min_start : region['end'] - min_start + 1]
            avg_coverage = np.mean(region_slice)  # Vectorized mean calculation
            coverage_values.append(avg_coverage)

        # Add coverage column to DataFrame
        merged_df['coverage'] = coverage_values

        merged_regions.append(merged_df)

    # Concatenate all gene-specific DataFrames into a final result
    final_df = pd.concat(merged_regions, ignore_index=True)

    # Keep only relevant columns
    final_df = final_df[['seqname', 'start', 'end', 'strand', 'gene_name', 'transcript_id', 'coverage']]

    # Display the DataFrame with coverage
    return(final_df)


def main(gtf_file, output_file, genes_str=None, identifier_type='gene_id', gene_feature='CDS'):
    """
    Main function to parse a GTF file (optionally filtered by a gene list)
    and print matched gene information, also merge regions and calculate coverage.

    Args:
        gtf_file (str): Path to the GTF file.
        output_file (str): Path to the output text file.
        genes_str (str): Comma-separated string of gene IDs or names to filter.
        identifier_type (str): Type of identifier provided ('gene_id' or 'gene_name').
        gene_feature (str): Feature type to extract (default = 'CDS').
    """
    # Convert the comma-separated string into a set of gene IDs or names (if provided)
    if genes_str:
        genes_of_interest = set([g.strip().lower() for g in genes_str.split(",")])
        print(f"Processing genes: {', '.join(genes_of_interest)}")
    else:
        genes_of_interest = None
        print("No gene list provided. Processing all genes.")

    # Parse the GTF file and filter by gene list
    gtf_df = parse_gtf(gtf_file, genes_of_interest, identifier_type, gene_feature)

    # Merge regions and calculate coverage
    merged_cov_df = merge_regions_and_coverage(genes_of_interest, gtf_df)

    # Write the merged results to an output file
    merged_cov_df.to_csv(output_file, sep='\t', index=False)

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse GTF file and calculate coverage")

    parser.add_argument("--gtf", required=True, help="Path to the GTF file")
    parser.add_argument("--output", required=True, help="Path to the output file")
    parser.add_argument("--genes", default=None, help="Comma-separated list of gene IDs or names to filter")
    parser.add_argument("--identifier_type", default='gene_id', help="Type of identifier provided ('gene_id' or 'gene_name')")
    parser.add_argument("--gene_feature", default='CDS', help="Feature type to extract (default = 'CDS')")
    args = parser.parse_args()
    main(args.gtf, args.output, args.genes, args.identifier_type, args.gene_feature)