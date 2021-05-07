rule trim_galore:
    input: 
        get_fastq
    output: 
        directory("results/quality_analysis/{sample}.{run}")
    conda:
        "../envs/trim_galore.yaml"
    threads: threads
    params:
        **config["params"]
    log: 
        "results/logs/trim_galore/{sample}-{run}.log"
    shell:
        "trim_galore {params.trim} --cores {threads} --output_dir {output} {input} 2> {log}"

rule concatenate_rv_comp:
    input: 
        lambda wildcards: \
            [f"results/quality_analysis/{wildcards.sample}.{run}" \
                for run in samples.loc[(wildcards.sample), ["run"]]['run'].dropna()
            ]
    output: 
        fq1="results/concatenated/{sample}.fq.gz"
    threads: 1
    log:
        "results/logs/concatenate/{sample}.log"
    shell:
        """
        set +e
        FQ1=`echo {input} | awk '{{for(X=1;X<=NF;X++){{OUT=OUT $X"/*_trimmed.fq.gz "}}}}END{{print OUT}}'`
        echo $FQ1
        cat $FQ1 > {output.fq1} 2> {log}
        """

rule trim_5_prime_Ts:
    input: 
        "results/concatenated/{sample}.fq.gz"
    output:
        "results/trim_5_prime_Ts/{sample}.fastq.gz"
    shell:
        """
        set +e
        zcat {input} | sed -r 's/(^T*)(.*)/\\2/'  | \
        awk -v mac="^@" '{{if($1 ~ mac){{flag=1;header=$1}}else if(flag==1){{len=length($1);seq=$1;flag=0}}else if($1 == "+"){{headerqual=$1}} else {{quallength=length($1);diff=quallength-len;text=substr($1,diff+1);if(length(seq)>=20 && diff>=2){{print header "Ts=" diff "\\n"seq"\\n"headerqual "\\n"text}}}}}}' | gzip -9 > {output}; 
        """
