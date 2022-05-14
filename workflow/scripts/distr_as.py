import argparse
import HTSeq
import pysam
import pandas as pd
import numpy as np
import random
from collections import defaultdict
from copy import deepcopy
import matplotlib.pyplot as plt


def reverse_complement(seq):
    seq_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    seq.rstrip()
    return "".join([seq_dict[base] for base in reversed(seq)])


def get_distribution(genome, sam):

    sequence_iterator = HTSeq.FastaReader(genome)

    bam_reader = HTSeq.BAM_Reader(sam)

    sequences = dict((s.name, s) for s in HTSeq.FastaReader(genome))

    distr_as = dict((i, 0) for i in range(0, 12))
    num_seqs = 0
    for line in bam_reader:
        read = None
        read = line.iv.copy()
        chrom_len = len(sequences[read.chrom].seq)
        if abs(read.end - read.start) < 15:
            continue
        if read.strand == "+":
            read.strand = "-"
            read.end = read.start + 1
            read.start = read.start - 10
            if read.start < 0:
                read.start = 0
            if read.end > chrom_len:
                read.end = chrom_len
        elif read.strand == "-":
            read.strand = "+"
            read.start = read.end + 1
            read.end = read.end + 12
            if read.start > chrom_len:
                read.start = chrom_len
            if read.end > chrom_len:
                read.end = chrom_len

        fasta = pysam.faidx(
            sequence_iterator.fos, "%s:%d-%d" % (read.chrom, read.start, read.end - 1)
        )
        seq = fasta.splitlines()[1]
        if read.strand == "-":
            seq = reverse_complement(seq)

        distr_as[seq.count("A")] += 1
        with open(f"{sam}.seqs", "a") as fout:
            fout.write(seq)
        num_seqs += 1
    for key in distr_as.keys():
        distr_as[key] = (distr_as[key] / num_seqs) * 100
    df_tmp = pd.DataFrame(
        np.array([list(distr_as.values())]), index=[sam], columns=distr_as.keys()
    )
    return df_tmp


def get_random_transcript_position(transcripts):

    sequence_iterator = HTSeq.FastaReader(transcripts)

    sequences = dict((s.name, s) for s in HTSeq.FastaReader(transcripts))
    distr_as = dict((i, 0) for i in range(0, 12))
    num_seqs = 0
    for i in range(0, 1000000):
        r = random.choice(list(sequences.values()))

        pos = random.randint(1, len(r) - 11)

        fasta = pysam.faidx(sequence_iterator.fos, "%s:%d-%d" % (r.name, pos, pos + 10))
        seq = fasta.splitlines()[1]

        distr_as[seq.count("A")] += 1

        num_seqs += 1

    for key in distr_as.keys():
        distr_as[key] = (distr_as[key] / num_seqs) * 100
    df_tmp = pd.DataFrame(
        np.array([list(distr_as.values())]), index=["random"], columns=distr_as.keys()
    )
    return df_tmp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome")
    ap.add_argument("--sam")
    ap.add_argument("--transcripts")
    ap.add_argument("--output")

    args = ap.parse_args()

    sam_list = args.sam.split(" ")
    df = pd.DataFrame()
    for sam in sam_list:
        df = df.append(get_distribution(args.genome, sam))
        print(f"Esse é o df {df}")

    df = df.append(get_random_transcript_position(args.transcripts))
    df.T.to_csv(args.output, sep=",")
    # df.T.plot()
    # plt.savefig(f"{args.output}.png")


if __name__ == "__main__":
    main()
