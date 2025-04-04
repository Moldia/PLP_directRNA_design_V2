import argparse
import subprocess
import threading

def run_extract_features(args):
    cmd = [
        "python3", "codes/extract_features.py",
        "--gtf", args.extract_features_gtf,
        "--genes", args.extract_features_genes,
        "--identifier_type", args.extract_features_identifier_type,
        "--gene_feature", args.extract_features_gene_feature,
        "--output", args.extract_features_output
    ]
    subprocess.run(cmd, check=True)

def run_extract_transcriptome(args):
    cmd = [
        "python3", "codes/extract_mrna.py",
        "--gtf", args.extract_transcriptome_gtf,
        "--fasta", args.extract_transcriptome_fasta,
        "--output_file", args.extract_transcriptome_output_file
    ]
    subprocess.run(cmd, check=True)

def run_extract_sequences(args):
    cmd = [
        "python3", "codes/extract_sequences.py",
        "--fasta", args.extract_sequences_fasta,
        "--output_fasta", args.extract_sequences_output_fasta,
        "--identifier_type", args.extract_sequences_identifier_type,
        "--plp_length", str(args.extract_sequences_plp_length),
        "--gtf_output", args.extract_sequences_gtf_output
    ]
    subprocess.run(cmd, check=True)

def run_find_targets(args):
    cmd = [
        "python3", "codes/find_target.py",
        "--selected_features", args.find_target_selected_features,
        "--fasta_file", args.find_target_fasta_file,
        "--output_file", args.find_target_output_file,
        "--iupac_mismatches", args.find_target_iupac_mismatches,
        "--reference_fasta", args.find_target_reference_fasta,
        "--max_errors", str(args.find_target_max_errors),
        "--Tm_min", str(args.find_target_Tm_min),
        "--Tm_max", str(args.find_target_Tm_max),
        "--lowest_percentile_Tm_score_cutoff", str(args.find_target_lowest_percentile_Tm_score_cutoff),
        "--min_dist_probes", str(args.find_target_min_dist_probes),
        "--num_probes", args.find_target_num_probes
    ]
    if args.find_target_off_target_output:
        cmd.append("--off_target_output")

        
    if args.find_target_filter_ligation_junction:
        cmd.append("--filter_ligation_junction")
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Master script to run the workflow: Extract transcriptome, Extract features, Extract sequences, and Find targets."
    )

    # Group: Extract Features
    parser.add_argument("--extract_features_gtf", default="data/tmp.gtf", help="Path to the GTF file for extract_features.py")
    parser.add_argument("--extract_features_genes", default="Grik2", help="Gene names for extract_features.py")
    parser.add_argument("--extract_features_identifier_type", default="gene_name", help="Identifier type for extract_features.py")
    parser.add_argument("--extract_features_gene_feature", default="CDS", help="Gene feature for extract_features.py")
    parser.add_argument("--extract_features_output", default="extract_features_output.txt", help="Output file for extract_features.py")

    # Group: Extract Transcriptome
    parser.add_argument("--extract_transcriptome_gtf", default="data/tmp.gtf", help="Path to the GTF file for extract_mrna.py")
    parser.add_argument("--extract_transcriptome_fasta", default="data/Mus.fa", help="Path to the FASTA file for extract_mrna.py")
    parser.add_argument("--extract_transcriptome_output_file", default="data/transcriptome_out.fa", help="Output file for extract_mrna.py")

    # Group: Extract Sequences
    parser.add_argument("--extract_sequences_fasta", default="data/Mus.fa", help="Path to the FASTA file for extract_sequences.py")
    parser.add_argument("--extract_sequences_output_fasta", default="extract_seqs_output.fa", help="Output FASTA file for extract_sequences.py")
    parser.add_argument("--extract_sequences_identifier_type", default="gene_name", help="Identifier type for extract_sequences.py")
    parser.add_argument("--extract_sequences_plp_length", type=int, default=30, help="PLP length for extract_sequences.py")
    parser.add_argument("--extract_sequences_gtf_output", default="extract_features_output.txt", help="Output GTF from extract_features.py for extract_sequences.py")

    # Group: Find Targets
    parser.add_argument("--find_target_selected_features", default="extract_features_output.txt", help="Selected features file for find_target.py")
    parser.add_argument("--find_target_fasta_file", default="extract_seqs_output.fa", help="FASTA file for find_target.py")
    parser.add_argument("--find_target_output_file", default="targets.txt", help="Output file for find_target.py")
    parser.add_argument("--find_target_iupac_mismatches", default="5:R,10:G", help="IUPAC mismatches for find_target.py")
    parser.add_argument("--find_target_reference_fasta", default="data/transcriptome_out.fa", help="Reference FASTA for find_target.py")
    parser.add_argument("--find_target_max_errors", type=int, default=4, help="Max errors for find_target.py")
    parser.add_argument("--find_target_Tm_min", type=int, default=58, help="Tm_min for find_target.py")
    parser.add_argument("--find_target_Tm_max", type=int, default=62, help="Tm_max for find_target.py")
    parser.add_argument("--find_target_lowest_percentile_Tm_score_cutoff", type=int, default=5, help="Lowest percentile Tm score cutoff for find_target.py")
    parser.add_argument("--find_target_min_dist_probes", type=int, default=8, help="Minimum distance between probes for find_target.py")
    parser.add_argument("--find_target_filter_ligation_junction", action="store_true", help="Include this flag to filter ligation junction for find_target.py")
    parser.add_argument("--find_target_num_probes", default="10", help="Number of probes to select for find_target.py")
    parser.add_argument("--find_target_off_target_output", action="store_true", help="Include this flag to output off-target information for find_target.py")

    args = parser.parse_args()

    # Run extract_features and extract_transcriptome in parallel
    thread_features = threading.Thread(target=run_extract_features, args=(args,))
    thread_transcriptome = threading.Thread(target=run_extract_transcriptome, args=(args,))
    thread_features.start()
    thread_transcriptome.start()

    # Wait for both parallel tasks to finish
    thread_features.join()
    thread_transcriptome.join()

    # Run extract_sequences (depends on output from extract_features)
    run_extract_sequences(args)

    # Run find_targets (requires outputs from previous steps)
    run_find_targets(args)

    print("Workflow completed successfully.")
