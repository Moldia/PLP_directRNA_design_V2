from typing import Tuple
import argparse
from . import probedesign as plp, InputValueError
import pandas as pd
from typing import Literal, Optional
import logging
from typing import Optional, Tuple, List

logger = logging.getLogger(__name__)


# Commands
def find_targets(
    selected_features: str,
    sequences_output: str,
    output_file: str,
    reference_fasta: str,
    min_coverage: int,
    gc_min: int,
    gc_max: int,
    num_probes: int,
    iupac_mismatches: list[Tuple[int, str]],
    max_errors: int,
    check_specificity: bool,
    plp_length: int,
    min_dist_probes: int,
    filter_ligation_junction: bool,
    off_target_output: bool,
    Tm_min: Optional[int] = None,
    Tm_max: Optional[int] = None,
    lowest_percentile_Tm_score_cutoff: Optional[int] = None,
):
    """
    Main function for probe extraction.
    """
    logger.info(f"🔹 Loading selected features from {selected_features}...")

    if off_target_output:
        logger.info(
            "🔹 Off-target output is enabled. Off-target information will be saved to a separate file."
        )
        targets_df, off_target_info = plp.find_targets(
            selected_features=selected_features,
            sequences_output=sequences_output,
            reference_fasta=reference_fasta,
            plp_length=plp_length,
            min_coverage=min_coverage,
            output_file=output_file,
            gc_min=gc_min,
            gc_max=gc_max,
            num_probes=num_probes,
            iupac_mismatches=iupac_mismatches,
            max_errors=max_errors,
            check_specificity=check_specificity,
            off_target_output=True,
        )
        # Save the off-target information
        if off_target_info is not None:
            off_target_info.to_csv(
                f'{output_file}_off_target.csv', sep=",", index=False
            )
        else:
            logger.warning("No off-target information to save (off_target_info is None).")
    else:
        targets_df, off_target_info = plp.find_targets(
            selected_features=selected_features,
            sequences_output=sequences_output,
            reference_fasta=reference_fasta,
            plp_length=plp_length,
            min_coverage=min_coverage,
            output_file=output_file,
            gc_min=gc_min,
            gc_max=gc_max,
            num_probes=num_probes,
            iupac_mismatches=iupac_mismatches,
            max_errors=max_errors,
            check_specificity=check_specificity,
            off_target_output=False,
        )

    # Calculate the melting temperature scores
    if Tm_min is None or Tm_max is None or lowest_percentile_Tm_score_cutoff is None:
        suggested_cutoff = 0
        scores = 'NA' 
        targets_df["Melt_Tm_scores"] = scores
    else:
        sequences = targets_df["Sequence"]
        scores = [
            plp.score_padlock_probe(seq, Tm_min=Tm_min, Tm_max=Tm_max) for seq in sequences
        ]
        # Calculate the suggested cutoff based on the 5th percentile
        suggested_cutoff = plp.analyze_scores(
            scores, percentile=lowest_percentile_Tm_score_cutoff
        )
        # Filter the targets based on the suggested cutoff
        targets_df["Melt_Tm_scores"] = scores
        targets_df = targets_df[targets_df["Melt_Tm_scores"] <= suggested_cutoff]
    # Filter the probes based on the minimum distance between probes
    targets_df = plp.filter_probes_by_distance(
        targets_df, min_dist_probes=min_dist_probes
    )
    # filter the probes based on the ligation junction preferences
    if filter_ligation_junction:
        targets_df = targets_df[targets_df["Ligation junction"] != "non-preferred"]

    targets_df = plp.select_top_probes(targets_df, num_probes)
    # Save the output
    targets_df.to_csv(output_file, sep="\t", index=False)
    return targets_df


def extract_features(
    gtf: str,
    output: Optional[str] = None,
    genes: Optional[set[str]] = None,
    identifier_type: Literal["gene_id", "gene_name"] = "gene_id",
):
    # Parse the GTF file and filter by gene list
    gtf_df = plp.parse_gtf(gtf, genes, identifier_type)

    logger.info(
        f"🔹 Extracted {len(gtf_df)} features from GTF file."
    )

    # Merge regions and calculate coverage
    merged_cov_df = plp.merge_regions_and_coverage(genes, gtf_df)

    if output is not None:
        # Write the merged results to an output file
        merged_cov_df.to_csv(output, sep="\t", index=False)

        logger.info(f"Results saved to {output}")

    return merged_cov_df


def extract_mrna(fasta_file: str, gtf_file: str, output_file: str):
    # Call the extraction function from your package.
    # You can adjust the options below as needed.
    records = plp.extract_mrna_sequences(
        fasta_file=fasta_file,
        gtf_file=gtf_file,
        output_file=output_file,
        plus_strand_only=False,
        revcomp=False,
        translate=False,
        codon_table=1,
        alternative_start_codon=True,
        clean_final_stop=True,
        clean_internal_stop=False,
        verbose=False,
    )

    logger.info(f"Extracted mRNA sequences have been saved to {output_file}")
    return records


def extract_sequences(
    gtf_output: str,
    fasta: str,
    output_fasta: str,
    plp_length: int,
    identifier_type: Literal["gene_name", "gene_id"],
    regions_file: str = "regions"
):
    """
    Main function to extract CDS sequences from a pre-processed GTF output.

    Args:
        gtf_output (str): Path to the GTF-derived DataFrame (TSV format).
        fasta (str): Path to the indexed FASTA file.
        output_fasta (str): Path to the output FASTA file.
    """
    logger.info(f"🔹 Loading GTF output from {gtf_output}...")
    df = pd.read_csv(gtf_output, sep="\t")

    # Ensure FASTA index exists
    plp.check_fasta_index(fasta)

    # Save regions for fast retrieval
    plp.save_regions_for_faidx(
        df, regions_file, plp_length, identifier_type=identifier_type
    )

    # Extract sequences
    plp.extract_sequences(fasta, regions_file + ".txt", output_fasta, df)


def run_plp_directrna(
    gtf: str,
    features_output: Optional[str],
    genes: Optional[set[str]],
    identifier_type: Literal["gene_id", "gene_name"],
    fasta: str,
    transcriptome_output: str,
    sequences_output: str,
    plp_length: int,
    targets_output: str,
    min_coverage: int,
    gc_min: int,
    gc_max: int,
    num_probes: int,
    iupac_mismatches: list[Tuple[int, str]],
    max_errors: int,
    check_specificity: bool,
    Tm_min: int,
    Tm_max: int,
    lowest_percentile_Tm_score_cutoff: int,
    min_dist_probes: int,
    filter_ligation_junction: bool,
    off_target_output: bool,
):
    # Step 1: Extract features
    logger.info("EXTRACTING FEATURES")
    extract_features(
        gtf, features_output, genes, identifier_type
    )

    # Step 2: Extract mRNA
    logger.info("EXTRACTING mRNA")
    extract_mrna(fasta, gtf, transcriptome_output)

    # Step 3: Extract sequences
    logger.info("EXTRACTING SEQUENCES")
    extract_sequences(
        features_output,
        fasta,
        sequences_output,
        plp_length,
        identifier_type,
    )

    # Step 4: Find targets
    logger.info("FINDING TARGETS")
    return find_targets(
        selected_features=features_output,
        sequences_output=sequences_output,
        output_file=targets_output,              # ← probes table
        reference_fasta=sequences_output,        # ← transcriptome sequences
        min_coverage=min_coverage,
        gc_min=gc_min,
        gc_max=gc_max,
        num_probes=num_probes,
        iupac_mismatches=iupac_mismatches,
        max_errors=max_errors,
        check_specificity=check_specificity,
        plp_length=plp_length,
        Tm_min=Tm_min,
        Tm_max=Tm_max,
        lowest_percentile_Tm_score_cutoff=lowest_percentile_Tm_score_cutoff,
        min_dist_probes=min_dist_probes,
        filter_ligation_junction=filter_ligation_junction,
        off_target_output=off_target_output
    )



def none_or_int(s):
    return None if s in ("None", "none", "") else int(s)

# Argument parsers

def run_plp_directrna_parser():
    parser = argparse.ArgumentParser(description="PLP DirectRNA Design Workflow")

    # Shared and individual inputs
    parser.add_argument("--gtf", required=True, help="Path to the GTF file")
    parser.add_argument("--genes", type=parse_genes, required=True, help="Comma-separated list of gene IDs or names to filter")
    parser.add_argument("--identifier_type", default="gene_name", choices=["gene_id", "gene_name"], help="Type of identifier provided")
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
    parser.add_argument("--iupac_mismatches", type=parse_iupac_mismatches, default=None, help="IUPAC mismatches parameter")
    parser.add_argument("--max_errors", type=int, default=1, help="Maximum number of errors allowed")
    parser.add_argument("--check_specificity", action="store_true", help="Enable specificity checking")
    parser.add_argument("--plp_length", type=int, default=30, help="PLP length for probe design")
    parser.add_argument("--off_target_output", action="store_true", help="Enable saving of off-target output")
    parser.add_argument("--Tm_min", type=none_or_int, help="Minimum Tm for probes")
    parser.add_argument("--Tm_max", type=none_or_int, help="Maximum Tm for probes")
    parser.add_argument("--lowest_percentile_Tm_score_cutoff", type=none_or_int, help="Lowest percentile Tm score cutoff")
    parser.add_argument("--min_dist_probes", type=int, default=100, help="Minimum distance between probes")
    parser.add_argument("--filter_ligation_junction", action="store_true", help="Enable filtering of ligation junctions")

    return parser


def extract_features_parser():
    parser = argparse.ArgumentParser(
        description="Parse GTF file and calculate coverage"
    )
    parser.add_argument("--gtf", required=True, help="Path to the GTF file")
    parser.add_argument("--output", required=True, help="Path to the output file")
    parser.add_argument(
        "--genes",
        type=parse_genes,
        default=None,
        help="Comma-separated list of gene IDs or names to filter",
    )
    parser.add_argument(
        "--identifier_type",
        default="gene_id",
        choices=["gene_id", "gene_name"],
        help="Type of identifier provided ('gene_id' or 'gene_name')",
    )
    return parser


def extract_mrna_parser():
    parser = argparse.ArgumentParser(
        description="Extract mRNA sequences using a GTF and FASTA file."
    )
    parser.add_argument("--gtf", required=True, help="Path to the GTF file")
    parser.add_argument("--fasta", required=True, help="Path to the FASTA file")
    parser.add_argument(
        "--output_file", required=True, help="Path to the output FASTA file"
    )
    return parser


def find_targets_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Extracts probe sequences fulfilling the following criteria:\n"
            "\n"
            "- GC content between 50-65% (default).\n"
            "- Ligation junctions must be 'preferred' or 'neutral' (not 'non-preferred').\n"
            "  See Xenium Custom Panel Design guide:\n"
            "  https://cdn.10xgenomics.com/image/upload/v1716400584/support-documents/CG000683_TechNote_Xenium_Custom_Panel_Design_RevD.pdf\n"
            "\n"
            "  Ligation junction preferences:\n"
            "    Preferred:      AT, TA, GA, AG\n"
            "    Neutral:        TT, CT, CA, TC, AC, CC, TG, AA\n"
            "    Non-Preferred:  CG, GT, GG, GC (filtered out)\n"
            "\n"
            "- No homopolymers of length 3 or more (e.g., AAA, TTT, GGG, CCC).\n"
            "- Minimum coverage of the region is met (default: 1, based on CDS/exon overlap).\n"
            "- The probe size must be an even number (default: 30)."
        ),
        formatter_class=argparse.RawTextHelpFormatter,  # Ensures multiline formatting
    )
    parser.add_argument(
        "--selected_features",
        required=True,
        help="Path to the selected features file (TSV format)",
    )
    parser.add_argument(
        "--sequences_output",
        required=True,
        help="Path to the extracted sequences file (CDS/exons)",
    )
    parser.add_argument("--output_file", required=True, help="Path to the output file")
    parser.add_argument(
        "--reference_fasta",
        required=True,
        help="Path to the reference (genome/transcriptome) FASTA file; Please note that with genome reference, the probe design will be performed on the whole genome which may take a long time.",
    )
    parser.add_argument(
        "--min_coverage", default=1, type=int, help="Minimum coverage of the region"
    )
    parser.add_argument("--gc_min", default=50, type=int, help="Minimum GC content")
    parser.add_argument("--gc_max", default=65, type=int, help="Maximum GC content")
    parser.add_argument(
        "--num_probes", default=10, type=int, help="Number of probes to select"
    )
    parser.add_argument(
        "--iupac_mismatches",
        type=parse_iupac_mismatches,
        default=None,
        help="IUPAC mismatches to consider. Note that the number of mismatches should be less than or equal to 2. Example: 5:R,6:A",
    )
    parser.add_argument(
        "--max_errors",
        default=1,
        type=float,
        help="Maximum error rate (or number of errors if an integer, recommended range is 0-6)",
    )
    parser.add_argument(
        "--check_specificity", action="store_true", help="Check probe specificity"
    )
    parser.add_argument(
        "--plp_length", default=30, type=int, help="Minimum probe length"
    )
    parser.add_argument(
        "--Tm_min",
        default=55.0,
        type=float,
        help="Minimum melting temperature for each arm",
    )
    parser.add_argument(
        "--Tm_max",
        default=65.0,
        type=float,
        help="Maximum melting temperature for each arm",
    )
    parser.add_argument(
        "--lowest_percentile_Tm_score_cutoff",
        default=5,
        type=int,
        help="The lowest percentile of the melting temperature score cutoff to filter",
    )
    parser.add_argument(
        "--min_dist_probes",
        default=10,
        type=int,
        help="Minimum distance between probes",
    )
    parser.add_argument(
        "--filter_ligation_junction",
        action="store_true",
        help="Filter probes based on ligation junction preferences; exclude non-preferred",
    )
    parser.add_argument(
        "--off_target_output", action="store_true", help="Output off-target information"
    )
    return parser


def extract_sequences_parser():
    parser = argparse.ArgumentParser(
        description="Extract CDS sequences from an indexed FASTA file"
    )
    parser.add_argument(
        "--gtf_output",
        required=True,
        help="Path to the GTF-derived TSV file (CDS regions)",
    )
    parser.add_argument("--fasta", required=True, help="Path to the indexed FASTA file")
    parser.add_argument(
        "--output_fasta", required=True, help="Path to the output FASTA file"
    )
    parser.add_argument(
        "--identifier_type",
        default="gene_id",
        choices=["gene_id", "gene_name"],
        help="Type of identifier provided ('gene_id' or 'gene_name')",
    )
    parser.add_argument("--plp_length", default=30, help="probe length")
    return parser


# Input parsers
def parse_iupac_mismatches(mismatch_str: str) -> list[Tuple[int, str]]:
    """
    Parses a string of mismatches formatted as "pos:base,pos:base" into a list of tuples.

    Args:
        mismatch_str (str): Mismatch input string (e.g., "5:R,10:G")

    Returns:
        list: A list of (position, base) tuples, e.g., [(5, 'R'), (10, 'G')].
    """
    if mismatch_str == "None":
        return None

    mismatches = []
    try:
        for pair in mismatch_str.split(","):
            pos, base = pair.split(":")
            pos = int(pos.strip())  # Convert position to integer
            base = base.strip().upper()  # Ensure base is uppercase
            mismatches.append((pos, base))

    except (ValueError, IndexError):
        raise InputValueError(
            "Invalid format for iupac mismatches. Use 'pos:base,pos:base', e.g., '5:R,10:G'.",
            field="iupac_mismatches",
            code="invalid_mismatch_format"
        )

    return mismatches


def parse_genes(genes_str: str) -> set[str]:
    logger.info(f"Parsing genes: {genes_str}")
    selected_genes = set([g.strip().lower() for g in genes_str.split(",")])
    return selected_genes
