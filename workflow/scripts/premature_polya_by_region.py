import argparse
import HTSeq
import math

def create_cds(gtf):
    gtf = HTSeq.GFF_Reader(gtf)
    gaos = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    gene_length = dict()
    for feature in gtf:
        if feature.type == "CDS":
            try:
                gene_length[feature.name] += abs(feature.iv.end - feature.iv.start +1)
            except:
                gene_length[feature.name] = abs(feature.iv.end - feature.iv.start +1)
            step_list = [step[0] for i, step in gaos[feature.iv].steps() if step]
            #print(step_list)
            step_list.append(feature)
            gaos[feature.iv] = step_list
            #except:
            #    print("Aqui")
            #    gaos[feature.iv] = feature
    return gaos, gene_length

def invert_read_strand(read):
    if read.iv.strand == '+':
        read.iv.strand = '-'
        read.iv.end = read.iv.start+1
    elif read.iv.strand == '-':
        read.iv.strand = '+'
        read.iv.start = read.iv.end
        read.iv.end = read.iv.end+1

    return read

def region_classify(read, gaos, gene_length):
    step_set = [step_set for i, step_set in gaos[read.iv].steps()][0]
    if not step_set:
        return -1
    #print(step_set)
    for step in step_set:
        #print(step_set)
        #print(step)
        if not step or isinstance(step, list):
            continue
        if step.type == 'CDS':
            
            if read.iv.strand == '+':
                #print(read.iv.strand, read.iv.start, step.iv.start, gene_length[step.name])
                r = math.ceil(((read.iv.start - step.iv.start+0.000001)/gene_length[step.name])*5)
            elif read.iv.strand == '-':
                #print(read.iv.strand, read.iv.start, step.iv.end, gene_length[step.name])
                r = math.ceil(((step.iv.end - read.iv.start+0.000001)/gene_length[step.name])*5)
    return str(r)


def classify(gtf, bam):
    
    gaos, gene_length = create_cds(gtf)
    bam = HTSeq.BAM_Reader(bam)
    
    region = {'1':0, '2':0, '3':0, '4':0, '5':0}
    #print(region)
    for read in bam:
        read = invert_read_strand(read)
        r = region_classify(read, gaos, gene_length)
        if r == -1:
            continue
        region[r]+=1
    
    return region
        


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf')
    ap.add_argument('--bam')
    ap.add_argument('--output')
    args = ap.parse_args()

    regions = classify(args.gtf, args.bam)

    with open(args.output, 'w') as out:
        out.write(str(regions))


if __name__ == "__main__":
    main()
