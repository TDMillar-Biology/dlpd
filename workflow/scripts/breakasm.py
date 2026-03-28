'''
Trevor Millar
A simple tool for breaking genome assemblies as a means of curation
'''

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os

def parse_breakfile(filepath):
    try:
        df = pd.read_csv(filepath, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError(f"Breakfile '{filepath}' not found.")

    # ensure idempotency by sorting
    df = df.sort_values(["contig", "break"])

    # Validate columns
    required_cols = {"contig", "break"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"{filepath} must contain columns: contig\\tbreak"
        )

    # Drop completely empty rows
    df = df.dropna(subset=["contig", "break"])

    # Validate break column
    try:
        df["break"] = pd.to_numeric(df["break"], errors="raise").astype(int)
    except Exception:
        bad_rows = df[pd.to_numeric(df["break"], errors="coerce").isna()]
        raise ValueError(f"Invalid break coordinates in:\n{bad_rows}")

    # Build dictionary
    contig_breaks = (
        df.groupby("contig")["break"]
        .apply(list)
        .to_dict()
    )

    return contig_breaks

def break_sequences(seqs, contig_breaks):
    new_records = []
    log_entries = []

    for contig, seq_record in seqs.items():
        if contig not in contig_breaks:
            new_records.append(seq_record)
            continue

        breaks = contig_breaks[contig]

        # Validate breakpoints
        for b in breaks:
            if b < 1 or b >= len(seq_record):
                raise ValueError(
                    f"Break coordinate {b} out of range for contig '{contig}' "
                    f"(length {len(seq_record)})"
                )

        starts = [0] + breaks
        ends = breaks + [len(seq_record)]

        fragment_ids = []

        for i, (start, end) in enumerate(zip(starts, ends), 1):
            fragment = seq_record.seq[start:end]
            new_id = f"{contig}_{i}"

            new_record = SeqRecord(fragment, id=new_id, description="")
            new_records.append(new_record)
            fragment_ids.append(new_id)

        log_entries.append(
            f"{contig} broken into: {', '.join(fragment_ids)} "
            f"at positions: {', '.join(map(str, breaks))}"
        )

    return new_records, log_entries

def read_fasta_as_dict(fasta_path):
    try:
        return SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    except FileNotFoundError:
        raise FileNotFoundError(f"Fasta file '{fasta_path}' not found.")
    
def write_outputs(records, log_entries, fasta_path, output_fasta, output_log):
    # Write FASTA
    SeqIO.write(records, output_fasta, "fasta-2line")

    # Write log
    with open(output_log, "w") as log:
        log.write(f"Input FASTA:   {fasta_path}\n")
        log.write(f"Output FASTA:  {output_fasta}\n\n")
        log.write("\n".join(log_entries) + "\n")

    print(f"Wrote broken assembly to '{output_fasta}'.")
    print(f"Wrote break log to '{output_log}'.")

def parse_cli():
    parser = argparse.ArgumentParser(
        description="Break contigs at specified coordinates."
    )

    parser.add_argument('fasta', help="Input assembly fasta file")

    parser.add_argument(
        '--breakfile',
        required=True,
        help="TSV file with columns: contig, break"
    )

    # Explicit outputs (preferred)
    parser.add_argument(
        '--out-fasta',
        help="Output FASTA path"
    )

    parser.add_argument(
        '--out-log',
        help="Output log path"
    )

    # Convenience mode (fallback)
    parser.add_argument(
        '--prefix',
        default="assembly",
        help="Output file prefix (used if explicit outputs not provided)"
    )

    parser.add_argument(
        '--out_dir',
        default=".",
        help="Output directory for prefix-based outputs"
    )

    args = parser.parse_args()

    return args

def main():
    args = parse_cli()

    contig_breaks = parse_breakfile(args.breakfile)
    seqs = read_fasta_as_dict(args.fasta)

    new_records, log_entries = break_sequences(
        seqs=seqs,
        contig_breaks=contig_breaks
    )

    if not contig_breaks:
        log_entries.append("No breakpoints provided; assembly unchanged.")

    # Resolve outputs
    if args.out_fasta and args.out_log:
        output_fasta = args.out_fasta
        output_log = args.out_log

    elif args.out_fasta or args.out_log:
        raise ValueError("Provide both --out-fasta and --out-log, or neither.")

    else:
        # fallback to prefix behavior
        output_prefix = os.path.join(args.out_dir, args.prefix)
        output_fasta = f"{output_prefix}.curated.fasta"
        output_log = f"{output_prefix}.curated.log"

    # ensure directories exist
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    os.makedirs(os.path.dirname(output_log), exist_ok=True)

    write_outputs(
        records=new_records,
        log_entries=log_entries,
        fasta_path=args.fasta,
        output_fasta=output_fasta,
        output_log=output_log
    )
    
if __name__ == "__main__":
    main()
