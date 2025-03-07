#!/usr/bin/env python3
from argparse import ArgumentParser
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import sys

import dnaio
from cutadapt.adapters import BackAdapter
from Bio.Seq import Seq

@dataclass
class Alignment:
    query_name: str         # probe ID
    target_name: str        # reference ID (or note if reverse)
    target_start: int       # starting position (1-indexed in output)
    target_end: int         # ending position
    mismatches: int         # number of mismatches (errors)
    query_sequence: str     # the probe sequence
    target_sequence: str    # the aligned segment of the reference

def parse_args():
    parser = ArgumentParser(description="Find probe sequences in a reference genome by splitting the probe list.")
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
        "--chunk-size", "-c",
        type=int,
        default=10,
        help="Number of probes per chunk (default: 10)"
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

def load_references(reference_path):
    """Load reference records as a list of (ref_id, ref_seq) tuples."""
    references = []
    with dnaio.open(reference_path) as f:
        for record in f:
            references.append((record.id, record.sequence))
    return references

def split_list(lst, chunk_size):
    """Yield successive chunk-sized lists from lst."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def find_all_alignments(ref, aligner):
    """
    Generator that yields (t_start, t_end, errors, target_seq)
    for each alignment in the reference string using the given aligner.
    Coordinates are adjusted to the original reference.
    """
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

def process_probe_chunk(probes_chunk, references, max_errors, reverse_complement):
    """
    For a given chunk of probes, iterate over all references.
    For each probe, recreate the adapter/aligner and find alignments.
    Returns a list of Alignment objects.
    """
    results = []
    for probe_id, probe_seq in probes_chunk:
        # Process each reference record
        for ref_id, ref_seq in references:
            # Forward search
            adapter = BackAdapter(probe_seq, max_errors=max_errors, min_overlap=len(probe_seq), indels=False)
            aligner = adapter.aligner
            for t_start, t_end, errors, target_seq in find_all_alignments(ref_seq, aligner):
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
                rev_ref_seq = str(Seq(ref_seq).reverse_complement())
                adapter = BackAdapter(probe_seq, max_errors=max_errors, min_overlap=len(probe_seq), indels=False)
                aligner = adapter.aligner
                for t_start, t_end, errors, target_seq in find_all_alignments(rev_ref_seq, aligner):
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

def process_all(probes, references, max_errors, threads, chunk_size, reverse_complement):
    """
    Splits the probe list into chunks, processes each chunk in parallel,
    prints the results in real time, and returns the complete table as a list.
    """
    table = []
    probe_chunks = list(split_list(probes, chunk_size))
    
    # Print CSV header
    header = "query_name,target_name,target_start,target_end,mismatches,query_sequence,target_sequence"
    print(header)

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for chunk in probe_chunks:
            futures.append(executor.submit(process_probe_chunk, chunk, references, max_errors, reverse_complement))
        
        for future in as_completed(futures):
            try:
                alignments = future.result()
                for aln in alignments:
                    line = f"{aln.query_name},{aln.target_name},{aln.target_start},{aln.target_end},{aln.mismatches},{aln.query_sequence},{aln.target_sequence}"
                    print(line)
                    table.append(aln)
            except Exception as e:
                sys.stderr.write(f"Error processing a probe chunk: {e}\n")
    return table

def main():
    args = parse_args()
    probes = load_probes(args.probe_path)
    if not probes:
        sys.stderr.write("No probes loaded; exiting.\n")
        return
    sys.stderr.write(f"Loaded {len(probes)} probes.\n")
    
    references = load_references(args.reference_path)
    if not references:
        sys.stderr.write("No references loaded; exiting.\n")
        return
    sys.stderr.write(f"Loaded {len(references)} reference records.\n")
    
    sys.stderr.write(f"Processing using {args.threads} threads and chunk size {args.chunk_size}...\n")
    # Process the probes in chunks over all references; the result table is returned
    result_table = process_all(probes, references, args.max_errors, args.threads, args.chunk_size, args.reverse_complement)
    sys.stderr.write("Processing completed.\n")
    return result_table

if __name__ == "__main__":
    # The table is both printed (real-time) and returned as an object.
    main()
