import argparse
import HTSeq
import pandas as pd
import numpy as np
import random

def add_genes(gtf, chrom_len):
    gtf = HTSeq.GFF_Reader(gtf)
    print(chrom_len)
    chrom_dict = chrom_len.to_dict('dict')['length']
    coverage = HTSeq.GenomicArray(chrom_dict, stranded=True, typecode="i")
    for feature in gtf:
        if feature.type == "gene":
            ch_len = chrom_len.loc[feature.iv.chrom].length
            #print(feature.iv, ch_len)
            if feature.iv.strand == "+":
                feature.iv.start = feature.iv.start - 300 
                feature.iv.end = feature.iv.end + 600 
            elif feature.iv.strand == "-":
                feature.iv.start = feature.iv.start - 600
                feature.iv.end = feature.iv.end + 300
            if feature.iv.start < 0:
                    feature.iv.start = 0
            if feature.iv.end > ch_len:
                    feature.iv.end = ch_len
            #print(feature.iv)
            coverage[feature.iv] = 1
    return coverage

def subtract_bed(coverage, bed_list, chrom_len, region_len):
    for file in bed_list:
        bed = HTSeq.BED_Reader(file)
        for feature in bed:
            coverage[feature.iv] = 0

    chrom_dict = {chrom: {strand: list() for strand in ['+', '-']} for chrom in chrom_len.index}
    
    for chrom in chrom_len.index:
        print(chrom)
        for strand in ['+', '-']:
            i = 0
            while(i < chrom_len.loc[chrom].length-region_len):
            #for i in range(0, chrom_len.loc[chrom].length):
                print(chrom, strand, i, coverage.chrom_vectors[chrom][strand][i])
                region_sum = sum(coverage.chrom_vectors[chrom][strand][i:i+region_len])
                if region_sum == region_len:
                    #print(type(coverage.chrom_vectors[chrom][strand][i:i+region_len]))
                    #print(coverage.chrom_vectors[chrom][strand][i:i+region_len])
                    #coverage.chrom_vectors[chrom][strand][i:i+150] = [0]*150
                    #coverage.chrom_vectors[chrom][strand][i] = 0
                    for j in range(i, i+region_len):
                        coverage.chrom_vectors[chrom][strand][j] = 0 
                    i+=region_len
                    chrom_dict[chrom][strand].append(i)
                else:
                    i+=1

    return chrom_dict

def sample_regions(chrom_dict, number_regions, output, region_len):
    with open(output, 'w') as fout:
        for i in range(int(number_regions)):
            print(i)
            chrom = random.choice(list(chrom_dict.keys()))
            strand = random.choice(list(chrom_dict[chrom].keys()))
            if chrom_dict[chrom][strand]:
                pos = random.choice(chrom_dict[chrom][strand])
                chrom_dict[chrom][strand].remove(pos)
                print(chrom, strand, pos, len(chrom_dict[chrom][strand]))
                fout.write(f"{chrom}\t{pos}\t{(pos+region_len)}\tregion_{i}\t1000\t{strand}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf')
    ap.add_argument('--chrom_len_file')
    ap.add_argument('--bed')
    ap.add_argument('--number_regions')
    ap.add_argument('--region_len')
    ap.add_argument('--output')

    args = ap.parse_args()

    chrom_len = pd.read_csv(args.chrom_len_file, sep=",").set_index("chr", drop=False)
    bed_list = args.bed.split(',')

    coverage = add_genes(args.gtf, chrom_len)
    chrom_dict = subtract_bed(coverage, bed_list, chrom_len, int(args.region_len))
    sample_regions(chrom_dict, args.number_regions, args.output, int(args.region_len))

if __name__ == "__main__":
    main()
