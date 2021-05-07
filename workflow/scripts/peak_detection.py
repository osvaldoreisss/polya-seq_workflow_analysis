import argparse
import HTSeq
import pandas as pd
import numpy as np
from scipy.signal import find_peaks

def invert_read_strand(read):
    if read.iv.strand == '+':
        read.iv.strand = '-'
        read.iv.end = read.iv.start+1
    elif read.iv.strand == '-':
        read.iv.strand = '+'
        read.iv.start = read.iv.end
        read.iv.end = read.iv.end+1

    return read

def peak_detection(sam,chrom_len):

    bam_reader = HTSeq.BAM_Reader(sam)
    coverage = HTSeq.GenomicArray("auto", stranded=True, typecode="i")
    chrom_array = {chrom: list() for chrom in chrom_len.index}

    for read in bam_reader:
        read = invert_read_strand(read)
        coverage[read.iv] += 1

    for chrom in chrom_len.index:
        for strand in ['+', '-']:
            chrom_array[chrom].append((find_peaks(np.array([coverage.chrom_vectors[chrom][strand][i] for i in range(0, chrom_len.loc[chrom].length)]), height=5, distance=25)))
   
    return chrom_array


def write_to_bed(chrom_array, bedfile, chrom_len):
    with open(bedfile, 'w') as fout:
        for chrom in chrom_len.index:
            peaks, properties = chrom_array[chrom][0]
            peaks = peaks.tolist()
            peak_heights = properties['peak_heights'].tolist()

            for i in range(0,len(peaks)):
                fout.write(f"{chrom}\t{(peaks[i]-24)}\t{(peaks[i]+25)}\tpeak_{i}\t{peak_heights[i]}\t+\n")

            peaks, properties = chrom_array[chrom][1]
            peaks = peaks.tolist()
            peak_heights = properties['peak_heights'].tolist()

            for i in range(0,len(peaks)):
                fout.write(f"{chrom}\t{(peaks[i]-24)}\t{(peaks[i]+25)}\tpeak_{i}\t{peak_heights[i]}\t-\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sam')
    ap.add_argument('--chrom_len_file')
    ap.add_argument('--output')

    args = ap.parse_args()

    chrom_len = pd.read_csv(args.chrom_len_file, sep=",").set_index("chr", drop=False)

    chrom_array = peak_detection(args.sam, chrom_len)

    write_to_bed(chrom_array, args.output, chrom_len)

if __name__ == "__main__":
    main()

