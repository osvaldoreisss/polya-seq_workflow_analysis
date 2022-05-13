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

rule classify_polya_peak:
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

rule classify_polya:
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

rule feauture_distribution_peak:
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

rule feauture_distribution:
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

rule premature_polya_by_region:
    input: 
        config["reference"]["gtf"],
        "results/filter_non_specific_annealing/{sample}.bam"
    output:
        "results/premature_polya_by_region/{sample}.txt",
        "results/premature_polya_by_region/{sample}_stop_in_frame.txt"
    conda:
        "../envs/filter.yaml"
    log:
        "results/logs/premature_polya_by_region/{sample}.log"
    shell:
        "python workflow/scripts/premature_polya_by_region.py --gtf {input[0]} --bam {input[1]} --output {output[0]} --out_not_stop_in_frame {output[1]} 2> {log};"

rule plot_polyadenylation_by_feature:
    input:
        "results/classify_old/SRR11849623.tsv",
        "results/classify_old/SRR11849624.tsv"
    output:
        multiext("results/plots/fig_1A", ".pdf", ".svg", ".png")
    run:
        import json
        import re
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np

        df = pd.DataFrame()
        for file in input:
            with open(file) as json_file:
                data = json_file.read().replace("\'", "\"")
                res = json.loads(data)
                res = [res['five_prime_utr'], res['CDS'],res['three_prime_utr']]
                df_tmp = pd.DataFrame(res).T
                df_tmp.index=[file.split("/")[2].split('.')[0]]
                df = df.append(df_tmp)

        df.columns = ['five_prime_utr', 'cds', 'three_prime_utr']
        df.loc[:,"five_prime_utr":"three_prime_utr"] = df.loc[:,"five_prime_utr":"three_prime_utr"].div(df.sum(axis=1), axis=0)
        data = df.sort_index()
        data = data * 100

        labels = ['five_prime_utr', 'cds', 'three_prime_utr']
        x_pos = np.arange(len(labels))
        means = list(data.mean())
        stds = list(data.std())
        
        fig = plt.figure()
        ax = plt.gca()

        bars = ax.bar(x_pos, means,
            yerr=stds,
            align='center',
            alpha=0.5,
            ecolor='black',
            capsize=10,
            color='#0072b2')
        ax.bar_label(bars, fmt='%.1f%%')
        ax.set_ylabel('% of transcripts')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels)
        ax.set_title('Distribution of clivages sites across features')
        #ax.plt.tight_layout()
        ax.figure.savefig(output[0], transparent=True)
        ax.figure.savefig(output[1], transparent=True)
        ax.figure.savefig(output[2], transparent=True)


rule plot_premature_polyadenylation_x_expression:
    input:
        get_quantification
    output:
        multiext("results/plots/fig_1B_{condition}", ".pdf", ".svg", ".png")
    run:
            from sklearn.linear_model import LinearRegression
            import pandas as pd
            import matplotlib.pyplot as pl
            

            df_cds = pd.DataFrame()
            df_three = pd.DataFrame()

            cds_file_count = 0
            three_file_count = 0

            files = [f"{prefix}_{sufix}" for prefix in input for sufix in ['three_out.txt', 'cds_out.txt']]

            for file in files:
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

            model = LinearRegression()
            model.fit(df_three, df_cds)
            r_squared = model.score(df_three, df_cds)

            fig = pl.figure()
            ax = pl.gca()
            ax.scatter(df_three, df_cds , alpha=0.5, color='#0072b2')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(0, 1000000)
            ax.set_ylim(0, 1000000)
            ax.set_xlabel("Counts 3' UTR")
            ax.set_ylabel("Counts CDS")
            ax.set_title(wildcards.condition + ' R²: ' + "{:.2f}".format(r_squared))
            ax.figure.savefig(output[0], transparent=True)
            ax.figure.savefig(output[1], transparent=True)
            ax.figure.savefig(output[2], transparent=True)

rule plot_polyadenylation_in_cds_by_region:
    input:
        get_polya_by_region
    output:
        multiext("results/plots/fig_1C_{condition}", ".pdf", ".svg", ".png")
    run:
        import json
        import re
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np
        from collections import Counter

        data = pd.DataFrame()
        for file in input:
            with open(file) as json_file:
                d = json_file.read().replace("\'", "\"")
                res = json.loads(d)
                data = data.append(pd.DataFrame(res, index=[0]))


        data.loc[:,"1":"5"] = data.loc[:,"1":"5"].div(data.sum(axis=1), axis=0)
        #data = data.sort_index()
        data = data * 100

        labels = [str(int((int(x)/len(data.columns))*100))+'%' for x in data.columns]
        x_pos = np.arange(len(labels))
        means = list(data.mean())
        stds = list(data.std())

        plt.rcParams['figure.figsize'] = [6.8, 6.8]
        plt.rcParams['figure.autolayout'] = True
        print(plt.rcParams['figure.figsize'])
        

        fig = plt.figure()
        ax = plt.gca()

        bars = ax.bar(x_pos, means,
            yerr=stds,
            align='center',
            alpha=0.5,
            ecolor='black',
            capsize=10,
            color='#0072b2')
        ax.bar_label(bars, fmt='%.1f%%')
        ax.set_ylabel('% of premature transcripts')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels)
        ax.set_title(f'Premature Polyadenylation by CDS region ({wildcards.condition})')
        plt.tight_layout()
        ax.figure.savefig(output[0], transparent=True)
        ax.figure.savefig(output[1], transparent=True)
        ax.figure.savefig(output[2], transparent=True)

rule plot_polyadenylation_in_cds_with_stop_in_frame:
    input:
        expand("results/premature_polya_by_region/{sample}_stop_in_frame.txt", sample=samples['sample'])
    output:
        multiext("results/plots/fig_1D", ".pdf", ".svg", ".png")
    run:
        import json
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np

        cond_sample_dict = dict()

        for line in input:
            sample = line.split('/')[2].split('_')[0]
            x = samples.loc[(sample, slice(None)),'condition'].dropna()
            cond_sample_dict[sample] = x[sample][1]

        #print(cond_sample_dict)

        cond_dict_out_frame = dict()
        cond_dict_in_frame = dict()
        for file in input:
            with open(file) as json_file:
                data = json_file.read().replace("\'", "\"")
                res = json.loads(data)
                sample = file.split("/")[2].split('_')[0]
                cond = cond_sample_dict[sample]
                if cond in cond_dict_out_frame.keys():
                    cond_dict_out_frame[cond].append(res['out_frame'])
                else:
                    cond_dict_out_frame[cond] = []
                    cond_dict_out_frame[cond].append(res['out_frame'])
                if cond in cond_dict_in_frame.keys():
                    cond_dict_in_frame[cond].append(res['in_frame'])
                else:
                    cond_dict_in_frame[cond] = []
                    cond_dict_in_frame[cond].append(res['in_frame'])
        

        for cond in cond_dict_in_frame.keys():
            for i in range(len(cond_dict_in_frame[cond])):
                total_in_frame = cond_dict_in_frame[cond][i]
                total_out_frame = cond_dict_out_frame[cond][i]

                total = total_in_frame + total_out_frame
                cond_dict_in_frame[cond][i] = (cond_dict_in_frame[cond][i]/total)*100
                cond_dict_out_frame[cond][i] = (cond_dict_out_frame[cond][i]/total)*100

        data = pd.DataFrame.from_dict(cond_dict_out_frame)

        data = data[[
        'YPD_log_phase_yeast_cells', 
        'YPD_diauxic_yeast_cells', 
        'YPGal_yeast_cells', 
        'Minimal_yeast_cells', 
        'Sorbitol_yeast_cells', 
        'Rpb1_H1085Q_slower_yeast_cells',
        'Rpb1_F1086S_slow_yeast_cells',
        'Rpb1_L1101S_fast_yeast_cells',
        'Rpb1_E1103G_faster_yeast_cells',
        'spt4_yeast_cells',
        'hpr1_yeast_cells']]

        labels = list(data.columns)
        x_pos = np.arange(len(labels))
        means = list(data.mean())
        stds = list(data.std())
        print(len(labels))
        print(len(means))
        print(stds)


        plt.rcParams['figure.figsize'] = [7.8, 6.8]
        plt.rcParams['figure.autolayout'] = True
        print(plt.rcParams['figure.figsize'])
        

        fig = plt.figure()
        ax = plt.gca()

        bars = ax.bar(x_pos, means,
            yerr=stds,
            align='center',
            alpha=0.5,
                ecolor='black',
                capsize=10,
                color='#0072b2')
        ax.bar_label(bars, fmt='%.1f%%')
        ax.set_ylabel('% of transcripts')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels)
        ax.set_title('Premature Polyadenylation without STOP Codon in frame in different conditions')
        plt.xticks(rotation=90)
        plt.tight_layout()
        ax.figure.savefig(output[0], transparent=True)
        ax.figure.savefig(output[1], transparent=True)
        ax.figure.savefig(output[2], transparent=True)

rule plot_polyadenylation_in_cds_by_condition:
    input:
        expand("results/classify_old/{sample}.tsv", sample=samples['sample'])
    output:
        multiext("results/plots/fig_3A", ".pdf", ".svg", ".png")
    run:
        import json
        import re
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np

        cond_sample_dict = dict()

        for line in input:
            sample = line.split('/')[2].split('.')[0]
            x = samples.loc[(sample, slice(None)),'condition'].dropna()
            cond_sample_dict[sample] = x[sample][1]

        #print(cond_sample_dict)

        cond_dict_cds = dict()
        cond_dict_three = dict()
        cond_dict_five = dict()
        for file in input:
            with open(file) as json_file:
                data = json_file.read().replace("\'", "\"")
                res = json.loads(data)
                sample = file.split("/")[2].split('.')[0]
                cond = cond_sample_dict[sample]
                if cond in cond_dict_cds.keys():
                    cond_dict_cds[cond].append(res['CDS'])
                else:
                    cond_dict_cds[cond] = []
                    cond_dict_cds[cond].append(res['CDS'])
                if cond in cond_dict_five.keys():
                    cond_dict_five[cond].append(res['five_prime_utr'])
                else:
                    cond_dict_five[cond] = []
                    cond_dict_five[cond].append(res['five_prime_utr'])
                if cond in cond_dict_three.keys():
                    cond_dict_three[cond].append(res['three_prime_utr'])
                else:
                    cond_dict_three[cond] = []
                    cond_dict_three[cond].append(res['three_prime_utr'])
        

        for cond in cond_dict_cds.keys():
            for i in range(len(cond_dict_cds[cond])):
                total_cds = cond_dict_cds[cond][i]
                total_five = cond_dict_five[cond][i]
                total_three = cond_dict_three[cond][i]

                total = total_cds + total_five + total_three
                cond_dict_cds[cond][i] = (cond_dict_cds[cond][i]/total)*100
                cond_dict_five[cond][i] = (cond_dict_five[cond][i]/total)*100
                cond_dict_three[cond][i] = (cond_dict_three[cond][i]/total)*100

                

        data = pd.DataFrame.from_dict(cond_dict_cds)
        data = data[[
        'YPD_log_phase_yeast_cells', 
        'YPD_diauxic_yeast_cells', 
        'YPGal_yeast_cells', 
        'Minimal_yeast_cells', 
        'Sorbitol_yeast_cells', 
        'Rpb1_H1085Q_slower_yeast_cells',
        'Rpb1_F1086S_slow_yeast_cells',
        'Rpb1_L1101S_fast_yeast_cells',
        'Rpb1_E1103G_faster_yeast_cells',
        'spt4_yeast_cells',
        'hpr1_yeast_cells']]

        labels = list(data.columns)
        x_pos = np.arange(len(labels))
        means = list(data.mean())
        stds = list(data.std())
        print(len(labels))
        print(len(means))
        print(stds)


        plt.rcParams['figure.figsize'] = [6.8, 6.8]
        plt.rcParams['figure.autolayout'] = True
        print(plt.rcParams['figure.figsize'])
        

        fig = plt.figure()
        ax = plt.gca()

        bars = ax.bar(x_pos, means,
            yerr=stds,
            align='center',
            alpha=0.5,
            ecolor='black',
            capsize=10,
            color='#0072b2')
        ax.bar_label(bars, fmt='%.1f%%')
        ax.set_ylabel('% of transcripts')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels)
        ax.set_title('Premature Polyadenylation in different conditions')
        plt.xticks(rotation=90)
        plt.tight_layout()
        ax.figure.savefig(output[0], transparent=True)
        ax.figure.savefig(output[1], transparent=True)
        ax.figure.savefig(output[2], transparent=True)
        