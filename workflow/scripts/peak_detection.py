import argparse
import HTSeq
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import poly_A_classify
from collections import defaultdict

def invert_read_strand(read):
    if read.iv.strand == '+':
        read.iv.strand = '-'
        read.iv.end = read.iv.start+1
    elif read.iv.strand == '-':
        read.iv.strand = '+'
        read.iv.start = read.iv.end
        read.iv.end = read.iv.end+1

    return read

def peak_detection(sam,chrom_len,output_bam):

    bam_reader = HTSeq.BAM_Reader(sam)
    bam_writer = HTSeq.BAM_Writer.from_BAM_Reader( output_bam, bam_reader )
    coverage = HTSeq.GenomicArray("auto", stranded=True, typecode="i")
    chrom_array = {chrom: list() for chrom in chrom_len.index}

    for line in bam_reader:
        read = invert_read_strand(line)
        bam_writer.write(read)
        coverage[read.iv] += 1

    for chrom in chrom_len.index:
        for strand in ['+', '-']:
            chrom_array[chrom].append((find_peaks(np.array([coverage.chrom_vectors[chrom][strand][i] for i in range(0, chrom_len.loc[chrom].length)]), height=10, distance=25)))
   
    return chrom_array


def annotate_peak(chrom_array, chrom_len, gaos):
    peaks_dict = defaultdict()
    for chrom in chrom_len.index:
        for s in range(0,2):
            if s == 0:
                strand = '+'
            else:
                strand = '-'
            peaks, properties = chrom_array[chrom][s]
            peaks = peaks.tolist()
            peak_heights = properties['peak_heights'].tolist()

            for i in range(0,len(peaks)):
                iv = HTSeq.GenomicInterval( chrom, peaks[i], peaks[i]+1, strand)
                feature, gene = poly_A_classify.feature_classify(iv, gaos)
                if feature:
                    #peaks_dict[chrom][strand][gene][feature][peaks[i]]=peak_heights[i]
                    peaks_dict[(chrom, strand, gene, feature, peaks[i])]=peak_heights[i]
    return peaks_dict

def write_to_bed(peaks_dict, bedfile, region_len):
    i = 1
    with open(bedfile, 'w') as fout:
        for key in peaks_dict.keys():
            (chrom, strand, gene, feature, peak) = key
            fout.write(f"{chrom}\t{(peak-region_len)}\t{(peak+region_len)}\tpeak|{gene}|{feature}|{i}\t{peaks_dict[key]}\t{strand}\n")
            i+=1

def get_expression_by_feature(peaks_dict, feature_name):
    gene_expression = defaultdict()
    for key in peaks_dict.keys():
        (chrom, strand, gene, feature, peak) = key
        if feature == feature_name:
            if gene in gene_expression.keys():
                gene_expression[gene]+=peaks_dict[key]
            else:
                gene_expression[gene]=peaks_dict[key]
    return gene_expression

    """ for chrom in peaks_dict.keys():
        for strand in peaks_dict[chrom].keys():
            for gene in peaks_dict[chrom][strand].keys():
                for feature in peaks_dict[chrom][strand][gene].keys():
                    for peak in peaks_dict[chrom][strand][gene][feature].keys():
                        if feature == feature_name:
                            if gene in gene_expression.keys():
                                gene_expression[gene]+=peaks_dict[chrom][strand][gene][feature][peak]
                            else:
                                gene_expression[gene]=peaks_dict[chrom][strand][gene][feature][peak]
    return gene_expression """

def write_bed_by_feature(peaks_dict, bedfile, feature_name, region_len):
    i = 1
    gene_expression=get_expression_by_feature(peaks_dict, feature_name)
    with open(bedfile, 'w') as fout:
        for key in peaks_dict.keys():
            (chrom, strand, gene, feature, peak) = key
            if feature == feature_name:
                relative_peak_expression = float(int(peaks_dict[key])/gene_expression[gene])
                fout.write(f"{chrom}\t{(peak-region_len)}\t{(peak+region_len)}\tpeak|{gene}|{feature}|{i}\t{relative_peak_expression}\t{strand}\n")
            i+=1

        """ for chrom in peaks_dict.keys():
            for strand in peaks_dict[chrom].keys():
                for gene in peaks_dict[chrom][strand].keys():
                    for feature in peaks_dict[chrom][strand][gene].keys():
                        for peak in peaks_dict[chrom][strand][gene][feature].keys(): """""" 
                            relative_peak_expression = float(int(peaks_dict[chrom][strand][gene][feature][peak])/gene_expression[gene])
                            fout.write(f"{chrom}\t{(peak-24)}\t{(peak+25)}\tpeak|{gene}|{feature}|{i}\t{relative_peak_expression}\t{strand}\n")
                            i+=1 """

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sam')
    ap.add_argument('--gtf')
    ap.add_argument('--chrom_len_file')
    ap.add_argument('--output')
    ap.add_argument('--output_three')
    ap.add_argument('--output_bam')
    ap.add_argument('--region_len')

    args = ap.parse_args()

    chrom_len = pd.read_csv(args.chrom_len_file, sep=",").set_index("chr", drop=False)

    chrom_array = peak_detection(args.sam, chrom_len, args.output_bam)

    gaos = poly_A_classify.create_utr(args.gtf, 300, 600)

    peaks_dict= annotate_peak(chrom_array, chrom_len, gaos)

    write_to_bed(peaks_dict, args.output, int(args.region_len))

    write_bed_by_feature(peaks_dict, args.output_three, 'three_prime_utr', int(args.region_len))

if __name__ == "__main__":
    main()

