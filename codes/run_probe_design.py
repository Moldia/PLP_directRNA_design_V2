import subprocess
import argparse

def run_command(command):
    """
    Runs a shell command and prints output/errors in real-time.

    Args:
        command (list): Command to execute as a list of strings.
    """
    print(f"\n🚀 Running: {' '.join(command)}\n")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    for line in process.stdout:
        print(line, end="")  

    stderr_output = process.stderr.read()
    if stderr_output:
        print(f"\n❌ Error:\n{stderr_output}")

    process.wait()
    if process.returncode != 0:
        print(f"\n❌ Command failed: {' '.join(command)}")
        exit(1)

def main():
    parser = argparse.ArgumentParser(description="Run probe design pipeline")

    # Arguments for extract_features.py
    parser.add_argument("--gtf", required=True, help="Path to GTF file")
    parser.add_argument("--output", default="output_sequences.txt", help="Output file for extracted features")
    parser.add_argument("--genes", required=True, help="Comma-separated list of gene names or IDs")
    parser.add_argument("--gene_feature", default="CDS", help="Feature type to extract (default: CDS)")
    parser.add_argument("--identifier_type", default="gene_name", choices=["gene_id", "gene_name"], help="Identifier type")

    # Arguments for extract_sequences.py
    parser.add_argument("--fasta", required=True, help="Path to the reference FASTA file")
    parser.add_argument("--output_fasta", default="regions.fa", help="Output FASTA file with extracted sequences")
    parser.add_argument("--plp_length", default=30, type=int, help="Minimum probe length")

    args = parser.parse_args()

    # Step 1: Run extract_features.py
    extract_features_cmd = [
        "python", "extract_features.py",
        "--gtf", args.gtf,
        "--output", args.output,
        "--genes", args.genes,
        "--gene_feature", args.gene_feature,
        "--identifier_type", args.identifier_type
    ]
    run_command(extract_features_cmd)

    # Step 2: Run extract_sequences.py
    extract_sequences_cmd = [
        "python", "extract_sequences.py",
        "--gtf_output", args.output,  # Uses output from extract_features.py
        "--fasta", args.fasta,
        "--output_fasta", args.output_fasta,
        "--plp_length", str(args.plp_length),
        "--identifier_type", args.identifier_type
    ]
    run_command(extract_sequences_cmd)

    print("\n✅ All steps completed successfully!")

if __name__ == "__main__":
    main()
