rule quantify_cds_expression:
    input:
        "results/filter_non_specific_annealing/{sample}.bam",
    output:
        "results/htseq/{sample}.tsv",
    conda:
        "../envs/filter.yaml"
    threads: threads
    params:
        gtf=config["reference"]["gtf"],
    log:
        "logs/htseq/{sample}.log",
    shell:
        """
        htseq-count -f bam -s reverse -t CDS {input} {params.gtf} > {output}
        """


rule select_most_expressed:
    input:
        expand("results/htseq/{sample}.tsv", sample=samples["sample"]),
    output:
        "results/expressed/most_expressed.txt",
    params:
        cutoff=config["expression"]["cutoff"],
    run:
        import pandas as pd

        df = pd.DataFrame()

        for file in input:
            data = pd.read_csv(file, header=None, index_col=0, sep="\t")
            data = data.loc[~data.index.str.startswith("__")]
            df = df.append(data.T)

        mean_expression = df.sum().sort_values(ascending=False)
        selected = mean_expression
        with open(str(output), "w") as select_out:
            select_out.writelines("\n".join(list(selected.index)))
