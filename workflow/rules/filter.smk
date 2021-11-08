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
        "results/filter_non_specific_annealing/{sample}.bam",
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
                res = [res['five_prime_utr'], res['three_prime_utr'],res['CDS'], res['ncRNA'], res['tRNA'], res['snoRNA'], res['snRNA'], res['pseudogene'], res['transposable_element']]
                df_tmp = pd.DataFrame(res).T
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
        "results/peak_detection/{sample}.bed"
    conda:
        "../envs/peak_detection.yaml"
    params:
        chr_len=config["reference"]["chrom_len"]
    shell:
        "python workflow/scripts/peak_detection.py --sam {input} --chrom_len_file {params.chr_len} --output {output}"

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