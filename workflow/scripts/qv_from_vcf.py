'''
Trevor D. Millar and gpt

error_rate = (# homozygous variant bases) / (callable bases)
QV = -10 * log10(error_rate)
'''

#!/usr/bin/env python3

import argparse
import math
import pysam

def load_variants(vcf_path):
    """Open VCF using pysam."""
    return pysam.VariantFile(vcf_path)

def count_errors(vcf, min_gq=None):
    """
    Count homozygous ALT variants as errors.

    Returns:
        total_errors (int)
        n_snps (int)
        n_indel_bases (int)
    """
    total_errors = 0
    n_snps = 0
    n_indel_bases = 0

    sample = list(vcf.header.samples)[0] # single sample assumption fine

    for rec in vcf:
        # Filter: PASS only
        if rec.filter.keys() != {"PASS"}: # skip filtered records
            continue

        sample_data = rec.samples[sample] 

        # Genotype
        gt = sample_data.get("GT")
        if gt != (1, 1): # skip anything other than homozygous variants
            continue

        # Optional genotype quality filter
        if min_gq is not None:
            gq = sample_data.get("GQ")
            if gq is None or gq < min_gq: # gq is None raises potential problem, monitor
                continue

        ref = rec.ref
        alt = rec.alts[0]

        # SNP
        if len(ref) == 1 and len(alt) == 1:
            total_errors += 1
            n_snps += 1

        # INDEL
        else:
            indel_len = abs(len(ref) - len(alt))
            total_errors += indel_len
            n_indel_bases += indel_len

    return total_errors, n_snps, n_indel_bases

def get_genome_size(fasta_path):
    """
    Compute total assembly size (excluding Ns).
    """
    fasta = pysam.FastaFile(fasta_path)
    total = 0

    for contig in fasta.references:
        seq = fasta.fetch(contig)
        total += len(seq.replace("N", "").replace("n", "")) # dont count gaps

    return total

def compute_qv(total_errors, genome_size):
    """
    Compute QV from error count.
    """
    if total_errors == 0:
        return float("inf")

    error_rate = total_errors / genome_size
    qv = -10 * math.log10(error_rate)

    return qv

def parse_args():
    parser = argparse.ArgumentParser(description="Compute QV from VCF + assembly")

    parser.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
    parser.add_argument("--fasta", required=True, help="Assembly FASTA")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("--min_gq", type=int, default=None, help="Minimum GQ filter")

    args = parser.parse_args()

    return args

def main():
    # cli parse
    args = parse_args()

    # Load
    vcf = load_variants(args.vcf)

    # Count errors
    total_errors, n_snps, n_indel_bases = count_errors(
        vcf, min_gq=args.min_gq
    )

    # Genome size
    genome_size = get_genome_size(args.fasta)

    # Compute QV
    qv = compute_qv(total_errors, genome_size)

    # Write output
    with open(args.output, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"genome_size\t{genome_size}\n")
        out.write(f"total_errors\t{total_errors}\n")
        out.write(f"snp_errors\t{n_snps}\n")
        out.write(f"indel_error_bases\t{n_indel_bases}\n")
        out.write(f"qv\t{qv:.4f}\n")


if __name__ == "__main__":
    main()
