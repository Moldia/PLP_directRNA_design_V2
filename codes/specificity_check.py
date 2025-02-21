from argparse import ArgumentParser
from dataclasses import dataclass

import dnaio
from cutadapt.align import Aligner
from cutadapt.adapters import BackAdapter


@dataclass
class Alignment:
    start: int
    end: int
    score: int
    errors: int


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "--errors",
        "-e",
        dest="max_errors",
        type=float,
        default=1,
        help="Maximum error rate (or number of errors if an integer >=1)",
    )
    parser.add_argument(
        "reference_path",
        metavar="reference",
        help="Reference FASTA (can be compressed)",
    )
    parser.add_argument(
        "primer_path",
        metavar="primers",
        help="Primer sequences (FASTA, FASTQ, can be compressed)",
    )
    args = parser.parse_args()
    run(**vars(args))


def run(reference_path, primer_path, max_errors):
    with dnaio.open(primer_path) as f:
        records = list(f)
    adapters = [
        BackAdapter(record.sequence, max_errors=max_errors, min_overlap=len(record.sequence), indels=False)
        for record in records
    ]
    aligners = [
        (record.id, adapter.aligner) for record, adapter in zip(records, adapters)
    ]
    del adapters

    with dnaio.open(reference_path) as references:
        for record in references:
            for name, aligner in aligners:
                for alignment in find_all(record.sequence, aligner):
                    print(f"Found {name} in {record.id} starting at {alignment.start+1} with {alignment.errors} errors")


def find_all(ref, aligner):
    start = 0
    while (alignment := aligner.locate(ref)) is not None:
        reference_start, reference_end, query_start, query_end, score, errors = alignment
        # Note roles of reference and query are reversed compared to the
        # terminology used in Cutadapt
        alignment = Alignment(query_start + start, query_end + start, score, errors)
        yield alignment
        start = query_start + 1
        ref = ref[start:]


if __name__ == "__main__":
    main()