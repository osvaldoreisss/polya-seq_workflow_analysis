rule generate_genome_index:
    input: 
        config["reference"]["genome"]
    output:
        directory("index/Saccharomyces_cerevisiae.R64")
    conda:
        "../envs/mapping.yaml"
    threads: threads
    params:
        gtf=config["reference"]["gtf"]
    log:
        "logs/STAR/generate_index.log"
    shell:
        "mkdir index/Saccharomyces_cerevisiae.R64; STAR --runMode genomeGenerate --runThreadN {threads} --genomeFastaFiles {input} --genomeDir {output} --sjdbGTFfile {params.gtf} --sjdbOverhang 99 --genomeSAindexNbases 11 2> {log}" 

rule map:
    input: 
        "results/trim_5_prime_Ts/{sample}.fastq.gz",
        "index/Saccharomyces_cerevisiae.R64"
    output:
        "results/mapped/{sample}.",
        "results/mapped/{sample}.Aligned.sortedByCoord.out.bam"
    conda:
        "../envs/mapping.yaml"
    threads: threads
    params:
        **config["params"]
    log:
        "logs/STAR/{sample}.log"
    shell:
        "STAR --runThreadN {threads} --genomeDir {input[1]} --readFilesIn {input[0]} --outFileNamePrefix {output[0]} {params.star} 2> {log}; touch {output[0]}"