#! /usr/bin/env python
"""
Date: 2023/03/27
Author: Ziwei Pan

Extract per-read CpG information from Guppy modified-base BAM

Input: Guppy modified-base BAM is 0-based file. 

Output: Per-read file:
"""
import argparse
import csv
import os

import pysam
from Bio.Seq import reverse_complement
from modbampy import ModBam


def process_ref(ref_fasta):
    chr_list = []
    # Read reference
    reference_fasta = pysam.FastaFile(ref_fasta)

    # modbam needs to get all contig lengths first
    tmp_bam = pysam.AlignmentFile(input_bam, "rb")
    for i in range(tmp_bam.nreferences):
        chr_name = tmp_bam.get_reference_name(i)
        chr_size = tmp_bam.get_reference_length(chr_name)
        chr_list.append([chr_name, chr_size])
    return reference_fasta, chr_list


def process_modbam2bed(input_bam, chr_list, reference_fasta, output_bed, canon_threshold, mod_threshold):
    with open(output_bed, 'w') as output_bed:
        output_bed = csv.writer(output_bed, delimiter='\t', quoting=csv.QUOTE_NONE)
        # output header
        output_bed.writerow(['ID', 'Chr', 'Pos', 'Strand', 'Prediction', 'Prob_methylation'])
        # https://github.com/epi2me-labs/modbam2bed
        for chr_name, chr_size in chr_list:
            with ModBam(input_bam) as bam:
                # iterates the BAM by read
                for read in bam.reads(chr_name, 0, int(chr_size)):
                    read_data = []
                    if read.is_unmapped or read.is_secondary:
                        continue
                    else:
                        for pos_mod in read.mod_sites:
                            # Check modbampy part in https://github.com/epi2me-labs/modbam2bed
                            read_id, ref_pos, read_pos, ref_strand, mod_strand, canon_base, mod_base, mod_score = pos_mod
                            prob_mod = round(mod_score / 255, 3)

                            ## Ignore the read that doesn't have lignment (ref_pos=-1)
                            # https://github.com/epi2me-labs/modbam2bed/issues/16
                            if ref_pos == -1:
                                continue
                            ### Filter out CpG sites
                            if ref_strand == "+":
                                try:
                                    reference_base = reference_fasta.fetch(chr_name, ref_pos, ref_pos + 2).upper()
                                except:
                                    continue
                            elif ref_strand == "-":
                                # For "-" strand, get the reverse complement of the reference sequence
                                try:
                                    reference_base = reverse_complement(
                                        reference_fasta.fetch(chr_name, ref_pos - 1, ref_pos + 1)).upper()
                                except:
                                    continue
                            else:
                                reference_base = "NANA"
                                print("ref_strand is not correct!")

                            if reference_base in ['CG' or 'GC']:
                                if prob_mod < canon_threshold:  # Extract canon_read and label as 0
                                    # read_data.append([chr_name, ref_pos, ref_pos + 1, read_id, prob_mod, ref_strand, 0])
                                    meth_state = 0
                                elif prob_mod >= mod_threshold:  # Extract mod_read and label as 1
                                    # read_data.append([chr_name, ref_pos, ref_pos + 1, read_id, prob_mod, ref_strand, 1])
                                    meth_state = 1
                                else:
                                    continue
                                read_data.append([read_id, chr_name, ref_pos + 1, ref_strand, meth_state, prob_mod])
                            else:
                                continue

                    for read in read_data:
                        output_bed.writerow(read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert Guppy modified-base BAM to per-read results')
    parser.add_argument("-i", "--input_file", action="store", type=str,
                        help='input BAM file after Guppy modified base calling', required=True)
    parser.add_argument('-o', '--output_bed_file', action="store", type=str,
                        help='output BED file', required=True)
    parser.add_argument('-r', '--reference_fasta', action="store", type=str,
                        help='reference fasta', required=True)
    parser.add_argument('-a', '--canon_threshold', action="store", type=float, default=0.33,
                        help='Bases with mod. probability < THRESHOLD are counted as canonical (default 0.33)',
                        required=False)
    parser.add_argument('-b', '--mod_threshold', action="store", type=float, default=0.66,
                        help='Bases with mod. probability >= THRESHOLD are counted as modified (default 0.66)',
                        required=False)
    args = parser.parse_args()

    input_bam = os.path.abspath(args.input_file)
    output_bed = os.path.abspath(args.output_bed_file)
    ref_fasta = os.path.abspath(args.reference_fasta)
    canon_threshold = args.canon_threshold
    mod_threshold = args.mod_threshold

    print("Loading reference genome...")
    chr_list = process_ref(ref_fasta)
    reference_fasta, chr_list = process_ref(ref_fasta)

    print("Process modbam...")
    process_modbam2bed(input_bam, chr_list, reference_fasta, output_bed, canon_threshold, mod_threshold)
