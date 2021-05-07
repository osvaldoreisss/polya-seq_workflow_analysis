import argparse
import HTSeq
import re
from collections import defaultdict

def create_utr(gtf, five_length, three_length):
    gtf = HTSeq.GFF_Reader(gtf)
    gaos = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    for feature in gtf:
        if feature.type == "start_codon":
            five_utr_feature = feature
            five_utr_feature.type = "five_prime_utr"
            if five_utr_feature.iv.strand == '+':
                five_utr_feature.iv.end = five_utr_feature.iv.start - 1
                five_utr_feature.iv.start = five_utr_feature.iv.start - five_length -1
            elif five_utr_feature.iv.strand == '-':
                five_utr_feature.iv.start = five_utr_feature.iv.end + 1
                five_utr_feature.iv.end = five_utr_feature.iv.end + five_length + 1
            try:
                gaos[five_utr_feature.iv] += str(feature)
            except IndexError:
                pass
            
        if feature.type == "stop_codon":
            three_utr_feature = feature
            three_utr_feature.type = "three_prime_utr"
            if three_utr_feature.iv.strand == '+':
                three_utr_feature.iv.start = three_utr_feature.iv.end + 1
                three_utr_feature.iv.end = three_utr_feature.iv.end + three_length + 1
            elif three_utr_feature.iv.strand == '-':
                three_utr_feature.iv.end = three_utr_feature.iv.start - 1
                three_utr_feature.iv.start = three_utr_feature.iv.start - three_length - 1
            #print(feature, three_utr_feature)
            try:
                gaos[three_utr_feature.iv] += str(feature)
            except IndexError:
                pass
        #print(feature)
        try:
            gaos[feature.iv] += str(feature)
        except IndexError:
            pass

    return gaos

def invert_read_strand(read):
    if read.iv.strand == '+':
        read.iv.strand = '-'
        read.iv.end = read.iv.start+1
    elif read.iv.strand == '-':
        read.iv.strand = '+'
        read.iv.start = read.iv.end
        read.iv.end = read.iv.end+1

    return read

def feature_classify(read, gaos):
    step_set = [step_set for i, step_set in gaos[read.iv].steps()][0]
    
    feature_genes = set(re.findall(r"(\w+[prime|CDS]\w+\s'\w+[-]*\w+')", str(step_set)))
    feat_gene_dict = defaultdict(dict)
    gene = ""
    for f in feature_genes:
        feat, gene = f.split(' ')
        feat_gene_dict[feat] = gene.strip('\'')
    #if not feature_genes.intersection(genes):
    #    return ""
    classification = ""
    if step_set:
        features = str(step_set)
        #print(features)
        if 'three_prime_utr' in features:
            classification = 'three_prime_utr'
            gene = feat_gene_dict['three_prime_utr']
        elif 'five_prime_utr' in features:
            classification = 'five_prime_utr'
            gene = feat_gene_dict['five_prime_utr']
        elif 'CDS' in features:
            classification = 'CDS'
            gene = feat_gene_dict['CDS']
    #print(classification, gene, str(feature_genes), str(step_set))
    if classification:
        return classification, gene
    else:
        return "", ""

def classify(gtf, bam):
    
    gaos = create_utr(gtf, 300, 600)
    bam = HTSeq.BAM_Reader(bam)
    classification_result = defaultdict(dict)
    classification_result_gene = defaultdict(dict)
    reads_unclassified_tmp = []
    reads_unclassified = []
    for read in bam:    
        read = invert_read_strand(read)
        feature, gene = feature_classify(read, gaos)
        #print(f"-->>{feature} {gene} {classification_result}\n")
        if feature:
            try:
                classification_result[feature]+=1
            except (KeyError, TypeError) as _:
                classification_result[feature]=1
            try:
                classification_result_gene[feature][gene]+=1
            except(KeyError, TypeError) as _:
                classification_result_gene[feature][gene]=1
        else:
            reads_unclassified_tmp.append(read)
    
    gaos = create_utr(gtf, 300, 1200)
    for read in reads_unclassified_tmp:
        feature, gene = feature_classify(read, gaos)
        if feature:
            try:
                classification_result[feature]+=1
            except (KeyError, TypeError) as _:
                classification_result[feature]=1
            try:
                classification_result_gene[feature][gene]+=1
            except(KeyError, TypeError) as _:
                classification_result_gene[feature][gene]=1
        else:
            reads_unclassified.append(read)
    #print(f"-->>{classification_result}")
    return classification_result, classification_result_gene

def get_gene_set(gtf):
    gtf = HTSeq.GFF_Reader(gtf)
    gene_set = {feature.name for feature in gtf if feature.type == "gene"}
    return gene_set

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf')
    ap.add_argument('--bam')
    ap.add_argument('--prefix')
    #ap.add_argument('--genes')
    args = ap.parse_args()
    cds_out = open(f"{args.prefix}_cds_out.txt", 'w')
    five_out = open(f"{args.prefix}_five_out.txt", 'w')
    three_out = open(f"{args.prefix}_three_out.txt", 'w')
    
    classification, genes = classify(args.gtf, args.bam)
    print(dict(classification))

    gene_set = get_gene_set(args.gtf)

    for gene in gene_set:
        try:
            three_out.write(f"{gene}\t{genes['three_prime_utr'][gene]}\n")
        except (KeyError, TypeError) as _:
            genes['three_prime_utr'][gene] = 0
            three_out.write(f"{gene}\t{genes['three_prime_utr'][gene]}\n")
        try:
            five_out.write(f"{gene}\t{genes['five_prime_utr'][gene]}\n")
        except (KeyError, TypeError) as _:
            genes['five_prime_utr'][gene] = 0
            five_out.write(f"{gene}\t{genes['five_prime_utr'][gene]}\n")
        try:
            cds_out.write(f"{gene}\t{genes['CDS'][gene]}\n")
        except (KeyError, TypeError) as _:
            genes['CDS'][gene] = 0
            cds_out.write(f"{gene}\t{genes['CDS'][gene]}\n")

if __name__ == "__main__":
    main()
