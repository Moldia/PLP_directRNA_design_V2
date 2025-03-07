#!/usr/bin/env python3
#
# Run this script with 'pipx run findprobes.py'
#
# /// script
# dependencies = ["dnaio", "cutadapt"]
# ///
# Developed by:
#   Marcel Martin
# Modifed by:
#   Nima Rafati

# DONE
# - handle reverse complements

import sys
import csv
from argparse import ArgumentParser
from dataclasses import dataclass
import dnaio
from cutadapt.align import Aligner
from cutadapt.adapters import BackAdapter
from Bio.Seq import Seq

@dataclass
class Alignment:
    query_name: str         # Probe ID
    target_name: str        # Reference ID
    target_start: int       # Start position (1-based index)
    target_end: int         # End position
    mismatches: int         # Number of mismatches
    query_sequence: str     # The probe sequence
    target_sequence: str    # The aligned segment of the reference

def parse_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description="Find probe sequences in a reference genome.")
    parser.add_argument(
        "--errors", "-e",
        dest="max_errors",
        type=float,
        default=1,
        help="Maximum error rate (or number of errors if an integer >=1)"
    )
    parser.add_argument(
        "--output_file", "-o",
        type=str,
        default=None,
        help="Output CSV file to save results (default: None, prints to console)"
    )
    parser.add_argument(
        "--revcomp", "-rc",
        action="store_true",
        help="Also search for reverse complement matches."
    )
    parser.add_argument(
        "--reference_path",
        metavar="reference",
        help="Reference FASTA (can be compressed)"
    )
    parser.add_argument(
        "--probe_path",
        metavar="probes",
        help="Probe sequences (FASTA, FASTQ, can be compressed)"
    )
    return parser.parse_args()

def load_probes(probe_path):
    """Load probes from file; return list of (probe_id, probe_seq)."""
    probes = []
    with dnaio.open(probe_path) as f:
        for record in f:
            probes.append((record.id, record.sequence))
    return probes

def find_all(ref, aligner):
    """Yield alignment positions and mismatches for a given reference sequence."""
    offset = 0
    while True:
        result = aligner.locate(ref)
        if result is None:
            break
        ref_start, ref_end, query_start, query_end, score, errors = result
        t_start = query_start + offset
        t_end = query_end + offset
        target_seq = ref[query_start:query_end]
        yield (t_start, t_end, errors, target_seq)
        offset += query_start + 1
        ref = ref[query_start + 1:]

def find_alignments(probes, reference_record, max_errors, revcomp):
    """
    Finds alignments of all probes in a single reference sequence.
    """
    results = []
    ref_id = reference_record.id
    ref_seq = reference_record.sequence

    for probe_id, probe_seq in probes:
        # Forward search
        adapter = BackAdapter(probe_seq, max_errors=max_errors, min_overlap=len(probe_seq), indels=False)
        aligner = adapter.aligner
        for t_start, t_end, errors, target_seq in find_all(ref_seq, aligner):
            results.append(Alignment(
                query_name=probe_id,
                target_name=ref_id,
                target_start=t_start + 1,  # Convert to 1-based index
                target_end=t_end,
                mismatches=errors,
                query_sequence=probe_seq,
                target_sequence=target_seq
            ))
        # Reverse complement search if enabled
        if revcomp:
            rev_ref_seq = str(Seq(ref_seq).reverse_complement())
            adapter = BackAdapter(probe_seq, max_errors=max_errors, min_overlap=len(probe_seq), indels=False)
            aligner = adapter.aligner
            for t_start, t_end, errors, target_seq in find_all(rev_ref_seq, aligner):
                results.append(Alignment(
                    query_name=probe_id,
                    target_name=f"{ref_id}(reverse)",
                    target_start=t_start + 1,
                    target_end=t_end,
                    mismatches=errors,
                    query_sequence=probe_seq,
                    target_sequence=target_seq
                ))
    return results

def process_all(probes, reference_path, max_errors, revcomp, output_file):
    """
    Processes all reference sequences and finds alignments for all probes.
    Saves the output to a CSV file if specified or prints to the console.
    Also returns the results as a list of Alignment objects.
    """
    results = []
    with dnaio.open(reference_path) as references:
        for reference_record in references:
            alignments = find_alignments(probes, reference_record, max_errors, revcomp)
            results.extend(alignments)

    # Print or save results
    if output_file:
        with open(output_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["query_name", "target_name", "target_start", "target_end", "mismatches", "query_sequence", "target_sequence"])
            for aln in results:
                writer.writerow([aln.query_name, aln.target_name, aln.target_start, aln.target_end, aln.mismatches, aln.query_sequence, aln.target_sequence])
        sys.stderr.write(f"Results saved to {output_file}\n")
    else:
        print("query_name,target_name,target_start,target_end,mismatches,query_sequence,target_sequence")
        for aln in results:
            print(f"{aln.query_name},{aln.target_name},{aln.target_start},{aln.target_end},{aln.mismatches},{aln.query_sequence},{aln.target_sequence}")

    return results

def main():
    """Main function to execute the probe search."""
    args = parse_args()
    probes = load_probes(args.probe_path)
    
    if not probes:
        sys.stderr.write("No probes loaded; exiting.\n")
        return

    sys.stderr.write(f"Loaded {len(probes)} probes.\n")
    sys.stderr.write(f"Processing reference sequences from {args.reference_path}...\n")

    # Run the probe search and store the results
    result_table = process_all(probes, args.reference_path, args.max_errors, args.revcomp, args.output_file)

    sys.stderr.write("Processing completed.\n")
    return result_table  # Returning the object for further use

if __name__ == "__main__":
    main()
