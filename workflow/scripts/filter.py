import argparse
import HTSeq
import pysam
from collections import defaultdict
from copy import deepcopy


def reverse_complement(seq):
    seq_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    seq.rstrip()
    return "".join([seq_dict[base] for base in reversed(seq)])


def filter(genome, sam, out):

    sequence_iterator = HTSeq.FastaReader(genome)

    bam_reader = HTSeq.BAM_Reader(sam)

    bam_writer = HTSeq.BAM_Writer.from_BAM_Reader(out, bam_reader)

    sequences = dict((s.name, s) for s in HTSeq.FastaReader(genome))

    for line in bam_reader:
        read = None
        read = line.iv.copy()
        chrom_len = len(sequences[read.chrom].seq)
        if abs(read.end - read.start) < 15:
            continue
        if read.strand == "+":
            read.strand = "-"
            read.end = read.start + 1
            read.start = read.start - 40
            if read.start < 0:
                read.start = 0
            if read.end > chrom_len:
                read.end = chrom_len
        elif read.strand == "-":
            read.strand = "+"
            read.start = read.end + 1
            read.end = read.end + 42
            if read.start > chrom_len:
                read.start = chrom_len
            if read.end > chrom_len:
                read.end = chrom_len

        t_count = int(line.read.name.split("Ts=")[1])

        fasta = pysam.faidx(
            sequence_iterator.fos, "%s:%d-%d" % (read.chrom, read.start, read.end - 1)
        )
        seq = fasta.splitlines()[1]
        if read.strand == "-":
            seq = reverse_complement(seq)
        # print(f"{line}\t{seq}\t{seq[0:(t_count)]}\t{seq[(11-t_count):11]}\t{t_count}")
        if seq[0:t_count].count("A") < t_count - 1:
            bam_writer.write(line)
    bam_writer.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome")
    ap.add_argument("--sam")
    ap.add_argument("--output")

    args = ap.parse_args()

    filter(args.genome, args.sam, args.output)


if __name__ == "__main__":
    main()
