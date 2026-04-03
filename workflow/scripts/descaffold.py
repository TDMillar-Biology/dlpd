#!/usr/bin/env python3

import argparse
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def split_sequence(seq, threshold):
    """
    Split sequence on runs of N/n >= threshold.
    Returns list of sequence strings.
    """
    pattern = re.compile(f"[Nn]{{{threshold},}}")
    fragments = pattern.split(str(seq))

    # remove empty fragments
    return [frag for frag in fragments if frag]

def descaffold_record(record, threshold):
    """
    Split a SeqRecord into multiple SeqRecords.
    """
    fragments = split_sequence(record.seq, threshold)

    new_records = []
    base_name = record.id

    for i, frag in enumerate(fragments, start=1):
        new_id = f"{base_name}_{i}"

        new_record = SeqRecord(
            seq=frag,
            id=new_id,
            description=""
        )

        new_records.append(new_record)

    return new_records

def process_fasta(in_fasta, out_fasta, threshold):
    """
    Read input FASTA, descaffold all records, write output FASTA.
    """
    total_in = 0
    total_out = 0

    with open(out_fasta, "w") as out_handle:
        for record in SeqIO.parse(in_fasta, "fasta"):
            total_in += 1

            new_records = descaffold_record(record, threshold)
            total_out += len(new_records)

            SeqIO.write(new_records, out_handle, "fasta")

    return total_in, total_out

def main():
    parser = argparse.ArgumentParser(
        description="Descaffold FASTA by splitting on poly-N runs"
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA file"
    )

    parser.add_argument(
        "--threshold",
        type=int,
        default=10,
        help="Minimum N-run length to split (default: 10)"
    )

    args = parser.parse_args()

    total_in, total_out = process_fasta(
        args.fasta,
        args.output,
        args.threshold
    )

    print(f"Processed {total_in} sequences")
    print(f"Produced {total_out} sequences")


if __name__ == "__main__":
    main()
