#!/usr/bin/env python3
from argparse import ArgumentParser
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
import sys
import os

import dnaio
from cutadapt.adapters import BackAdapter
from Bio.Seq import Seq

@dataclass
class Alignment:
    query_name: str         # probe ID
    target_name: str        # reference ID (with note if reverse)
    target_start: int       # starting position (1-indexed in output)
    target_end: int         # ending position
    mismatches: int         # number of mismatches (errors)
    query_sequence: str     # the probe sequence
    target_sequence: str    # the aligned segment of the reference

def parse_args():
    parser = ArgumentParser(description="Find probe sequences in a reference genome.")
    parser.add_argument(
        "--errors", "-e",
        dest="max_errors",
        type=float,
        default=1,
        help="Maximum error rate (or number of errors if an integer >=1)"
    )
    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=os.cpu_count(),
        nargs="?",
        help="Number of parallel processes to use (default: all available CPUs)"
    )
    parser.add_argument(
        "--reverse-complement", "-rc",
        action="store_true",
        help="Also search for reverse complement matches."
    )
    parser.add_argument(
        "reference_path",
        metavar="reference",
        help="Reference FASTA (can be compressed)"
    )
    parser.add_argument(
        "probe_path",
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
    """
    Yield (target_start, target_end, errors, target_seq) for each alignment in ref.
    Note: The aligner (from BackAdapter) returns coordinates relative to the current
    reference substring, so we adjust them based on how much of ref was trimmed.
    """
    offset = 0
    while True:
        result = aligner.locate(ref)
        if result is None:
            break
        # Unpack result; note: roles of reference and query are reversed compared to Cutadapt terminology.
        ref_start, ref_end, query_start, query_end, score, errors = result
        # Adjust coordinates to the original reference
        t_start = query_start + offset
        t_end = query_end + offset
        target_seq = ref[query_start:query_end]
        yield (t_start, t_end, errors, target_seq)
        # Move past the first matched position in the current substring
        offset += query_start + 1
        ref = ref[query_start + 1:]

def find_alignments_for_reference(probes, reference_record, max_errors, reverse_complement):
    """
    For a given reference record, go over each probe and try to find all alignments.
    Recreate the adapter/aligner for each probe inside the worker.
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
                target_start=t_start + 1,  # convert to 1-indexed
                target_end=t_end,
                mismatches=errors,
                query_sequence=probe_seq,
                target_sequence=target_seq
            ))
        # Reverse complement search if enabled
        if reverse_complement:
            # Use the same adapter (its aligner) on the reverse complement of the reference.
            rev_ref = str(Seq(ref_seq).reverse_complement())
            adapter = BackAdapter(probe_seq, max_errors=max_errors, min_overlap=len(probe_seq), indels=False)
            aligner = adapter.aligner
            for t_start, t_end, errors, target_seq in find_all(rev_ref, aligner):
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

def process_references(probes, reference_path, max_errors, threads, reverse_complement):
    """Process each reference record in parallel and output CSV rows."""
    # Print CSV header
    print("query_name,target_name,target_start,target_end,mismatches,query_sequence,target_sequence")
    with dnaio.open(reference_path) as references:
        records = list(references)  # load all records to avoid IO bottlenecks in workers
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for record in records:
            futures.append(
                executor.submit(find_alignments_for_reference, probes, record, max_errors, reverse_complement)
            )
        for future in as_completed(futures):
            try:
                alignments = future.result()
                for aln in alignments:
                    print(f"{aln.query_name},{aln.target_name},{aln.target_start},{aln.target_end},{aln.mismatches},{aln.query_sequence},{aln.target_sequence}")
            except Exception as e:
                sys.stderr.write(f"Error processing a reference: {e}\n")

def main():
    args = parse_args()
    probes = load_probes(args.probe_path)
    if not probes:
        sys.stderr.write("No probes loaded; exiting.\n")
        return
    sys.stderr.write(f"Loaded {len(probes)} probes.\n")
    sys.stderr.write(f"Processing references from {args.reference_path} using {args.threads} threads...\n")
    process_references(probes, args.reference_path, args.max_errors, args.threads, args.reverse_complement)
    sys.stderr.write("Processing completed.\n")

if __name__ == "__main__":
    main()
