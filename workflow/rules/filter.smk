rule filter_uniq:
    input: 
        "results/mapped/{sample}.Aligned.sortedByCoord.out.bam"
    output: 
        "results/mapped_uniq/{sample}.bam"
    conda:
        "../envs/filter.yaml"
    shell:
        "samtools view -@ {threads} -h {input} | awk 'BEGIN{{OFS=\"\t\"}}{{if($1 ~ /^@/){{print}}else if($12==\"NH:i:1\"){{print}}}}' | samtools view -@ 1 -Sb - > {output};samtools index {output}"

rule filter_non_specific_annealing:
    input: 
        config["reference"]["genome"],
        "results/mapped_uniq/{sample}.bam"
    output:
        "results/filter_non_specific_annealing/{sample}.bam"
    conda:
        "../envs/filter.yaml"
    log:
        "results/logs/filter_non_specific_annealing/{sample}.log"
    shell:
        "python workflow/scripts/filter.py --genome {input[0]} --sam {input[1]} --output {output} 2> {log}; samtools index {output}"

rule classify_polya:
    input: 
        "results/filter_bam_by_peak/{sample}.bam",
        #"results/expressed/most_expressed.txt"
    output:
        "results/classify/{sample}.tsv",
        "results/quantification/{sample}"
    conda:
        "../envs/filter.yaml"
    params:
        gtf=config["reference"]["gtf"]
    shell:
        "python workflow/scripts/poly_A_classify.py --gtf {params.gtf} --bam {input[0]} --prefix {output[1]} > {output[0]}; touch {output[1]}"

rule classify_polya_old:
    input: 
        "results/filter_non_specific_annealing/{sample}.bam",
        #"results/expressed/most_expressed.txt"
    output:
        "results/classify_old/{sample}.tsv",
        "results/quantification_old/{sample}"
    conda:
        "../envs/filter.yaml"
    params:
        gtf=config["reference"]["gtf"]
    shell:
        "python workflow/scripts/poly_A_classify.py --gtf {params.gtf} --bam {input[0]} --prefix {output[1]} > {output[0]}; touch {output[1]}"

rule feauture_distribution:
    input:
        expand("results/classify/{sample}.tsv", sample=samples['sample'])
    output:
        "results/classify/feature_distribution.csv"
    run:
        import json
        import re
        import pandas as pd
        import matplotlib.pyplot as plt

        df = pd.DataFrame()
        
        for file in input:
            with open(file) as json_file:
                data = json_file.read().replace("\'", "\"")
                res = json.loads(data)
                #features = ['five_prime_utr', 'three_prime_utr', 'cds', 'ncrna', 'trna', 'snorna', 'snrna', 'pseudogene', 'transposable_element']
                features = ['five_prime_utr', 'three_prime_utr', 'CDS', 'ncRNA',  'tRNA', 'snoRNA', 'snRNA', 'pseudogene', 'transposable_element']
                out = list()
                for feature in features:
                    if feature in res:
                        out.append(res[feature])
                    else:
                        out.append(0)
                #res = [res['five_prime_utr'], res['three_prime_utr'],res['CDS'], res['ncRNA'], res['tRNA'], res['snoRNA'], res['snRNA'], res['pseudogene'], res['transposable_element']]
                df_tmp = pd.DataFrame(out).T
                df_tmp.index=[file.split("/")[2].split('.')[0]]
                df = df.append(df_tmp)

        df.columns = ['five_prime_utr', 'three_prime_utr', 'cds', 'ncrna', 'trna', 'snorna', 'snrna', 'pseudogene', 'transposable_element']
        df.loc[:,"five_prime_utr":"transposable_element"] = df.loc[:,"five_prime_utr":"transposable_element"].div(df.sum(axis=1), axis=0)
        data = df.sort_index()
        data = data * 100
        data.to_csv(output[0])

rule feauture_distribution_old:
    input:
        expand("results/classify_old/{sample}.tsv", sample=samples['sample'])
    output:
        "results/classify_old/feature_distribution.csv"
    run:
        import json
        import re
        import pandas as pd
        import matplotlib.pyplot as plt

        df = pd.DataFrame()
        
        for file in input:
            with open(file) as json_file:
                data = json_file.read().replace("\'", "\"")
                res = json.loads(data)
                #features = ['five_prime_utr', 'three_prime_utr', 'cds', 'ncrna', 'trna', 'snorna', 'snrna', 'pseudogene', 'transposable_element']
                features = ['five_prime_utr', 'three_prime_utr', 'CDS', 'ncRNA',  'tRNA', 'snoRNA', 'snRNA', 'pseudogene', 'transposable_element']
                out = list()
                for feature in features:
                    if feature in res:
                        out.append(res[feature])
                    else:
                        out.append(0)
                #res = [res['five_prime_utr'], res['three_prime_utr'],res['CDS'], res['ncRNA'], res['tRNA'], res['snoRNA'], res['snRNA'], res['pseudogene'], res['transposable_element']]
                df_tmp = pd.DataFrame(out).T
                df_tmp.index=[file.split("/")[2].split('.')[0]]
                df = df.append(df_tmp)

        df.columns = ['five_prime_utr', 'three_prime_utr', 'cds', 'ncrna', 'trna', 'snorna', 'snrna', 'pseudogene', 'transposable_element']
        df.loc[:,"five_prime_utr":"transposable_element"] = df.loc[:,"five_prime_utr":"transposable_element"].div(df.sum(axis=1), axis=0)
        data = df.sort_index()
        data = data * 100
        data.to_csv(output[0])
        
rule peak_detection:
    input: 
        "results/filter_non_specific_annealing/{sample}.bam",
    output:
        "results/peak_detection/{sample}.bed",
        "results/peak_detection/{sample}_three.bed",
        "results/peak_detection/{sample}.bam"
    conda:
        "../envs/peak_detection.yaml"
    params:
        chr_len=config["reference"]["chrom_len"],
        gtf=config["reference"]["gtf"],
        **config["params"]
    shell:
        "python workflow/scripts/peak_detection.py --sam {input} --chrom_len_file {params.chr_len} --gtf {params.gtf} --output {output[0]} --output_three {output[1]} --output_bam {output[2]} {params.peak_detection}"

rule filter_bam_by_peak:
    input:
        "results/peak_detection/{sample}.bed",
        "results/peak_detection/{sample}.bam"
    output:
        "results/filter_bam_by_peak/{sample}_read_names.txt"
    conda:
        "../envs/filter.yaml"
    shell:
        "bedtools intersect -wa -s -a {input[1]} -b {input[0]} | samtools view /dev/stdin | awk '{{print $1}}' > {output}"

rule select_peak_reads_by_name:
    input:
        "results/filter_non_specific_annealing/{sample}.bam",
        "results/filter_bam_by_peak/{sample}_read_names.txt"
    output:
        "results/filter_bam_by_peak/{sample}.bam",
    params:
        genome_fasta=config["reference"]["genome"]
    run:
        import HTSeq
        import pysam

        bam_reader = HTSeq.BAM_Reader(input[0])

        bam_writer = HTSeq.BAM_Writer.from_BAM_Reader( output[0], bam_reader )

        read_set = set()
        with open(input[1]) as reads:
            read_set = {read.rstrip() for read in reads}
        
        for line in bam_reader:
            if line.read.name in read_set:
                bam_writer.write( line )
        bam_writer.close()

rule plot_polyadenylation_by_feature:
    input:
        "results/classify/SRR11849623.tsv",
        "results/classify/SRR11849624.tsv"
    output:
        "results/plots/fig_1A.svg"
    run:
        import json
        import re
        import pandas as pd
        import matplotlib.pyplot as pl
        import numpy as np

        df = pd.DataFrame()
        for file in input:
            with open(file) as json_file:
                data = json_file.read().replace("\'", "\"")
                res = json.loads(data)
                res = [res['five_prime_utr'], res['three_prime_utr'],res['CDS']]
                df_tmp = pd.DataFrame(res).T
                df_tmp.index=[file.split("/")[2].split('.')[0]]
                df = df.append(df_tmp)
        
        df.columns = ['five_prime_utr', 'three_prime_utr', 'cds']
        df.loc[:,"five_prime_utr":"cds"] = df.loc[:,"five_prime_utr":"cds"].div(df.sum(axis=1), axis=0)
        data = df.sort_index()
        data = data * 100

        ax = data.T.plot(kind='bar', y=['SRR11849623', 'SRR11849624'])

        ax.bar_label(ax.containers[0], fmt='%.1f%%')
        ax.bar_label(ax.containers[1], fmt='%.1f%%')
        ax.legend(labels=['WT_Rep1','WT_Rep2'])
        ax.figure.savefig(output[0])


rule plot_premature_polyadenylation_x_expression:
    input:
        get_quantification
    output:
        "results/plots/fig_1B_{condition}.svg"
    run:
        from sklearn.metrics import r2_score
        import pandas as pd
        import matplotlib.pyplot as pl

        df_cds = pd.DataFrame()
        df_three = pd.DataFrame()

        cds_file_count = 0
        three_file_count = 0

        for file in input:
            df_tmp = pd.DataFrame()
            if 'cds_out' in file:
                df_tmp = pd.read_csv(file, sep='\t', header=None, index_col=0).sort_index()
                if not df_cds.empty:
                    df_cds = df_cds + df_tmp
                else:
                    df_cds = df_tmp
                cds_file_count+=1
            elif 'three_out' in file:
                df_tmp = pd.read_csv(file, sep='\t', header=None, index_col=0).sort_index()
                if not df_three.empty:
                    df_three = df_three + df_tmp
                else:
                    df_three = df_tmp.copy()
                three_file_count+=1

        df_cds = df_cds / cds_file_count
        df_three = df_three / three_file_count


        fig = pl.figure()
        ax = pl.gca()
        ax.scatter(df_three, df_cds , c='blue')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0, 1000000)
        ax.set_ylim(0, 1000000)
        ax.set_xlabel("Counts 3' UTR")
        ax.set_ylabel("Counts CDS")
        ax.set_title(wildcards.condition + ' R2: ' + "{:.4f}".format(r2_score(df_three, df_cds)))
        ax.figure.savefig(output[0])
