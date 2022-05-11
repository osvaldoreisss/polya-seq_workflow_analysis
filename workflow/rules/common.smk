from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep=",").set_index(["sample", "run"], drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

threads=config['threads']

def get_fastq(wildcards):
    fastqs = samples.loc[(wildcards.sample, int(wildcards.run)), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return f"libs/{fastqs.fq1}", f"libs/{fastqs.fq2}"
    return f"libs/{fastqs.fq1}"

def get_quantification(wildcards):
    #print(wildcards.condition)
    samples_copy = samples.copy()
    s = samples_copy.set_index(["sample","run","condition"], drop=False).loc[(slice(None), slice(None), wildcards.condition),'sample'].dropna()
    result = expand( "results/quantification_old/{x}", x = s )
    return result

def get_polya_by_region(wildcards):
    samples_copy = samples.copy()
    s = samples_copy.set_index(["sample","run","condition"], drop=False).loc[(slice(None), slice(None), wildcards.condition),'sample'].dropna()
    result = expand( "results/premature_polya_by_region/{x}.txt", x = s )
    return result