from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep=",").set_index(["sample","run","condition"], drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

threads=config['threads']

def get_fastq(wildcards):
    fastqs = samples.loc[(wildcards.sample, int(wildcards.run)), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return f"libs/{fastqs.fq1}", f"libs/{fastqs.fq2}"
    return f"libs/{fastqs.fq1}"

def get_quantification(wildcards):
    s = samples.loc[(slice(None), slice(None), wildcards.condition),'sample'].dropna()
    files = list()
    for sample in s:
        files.append(f'results/quantification/{sample}_three_out.txt')
        files.append(f'results/quantification/{sample}_cds_out.txt')
    return files