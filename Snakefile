# Snakefile
# This is designed to use the cluster, so it may not be the completely portable
# But it provides an idea and a learning place for me to learn Snakemake :)
# Ugalde 2019

from collections import defaultdict

#FOLDERS]
FASTQ = "data/fastq/"
PROCESS = "data/process/"
REFERENCES = "data/references/"



# Get samples names from the metadata
# This should go maybe on the config file?

SAMPLES = ["Sample1", "Sample2"]

rule all:
    input:
        PROCESS + "megahit_assembly/BT1.contigs.fa",
        expand(PROCESS + "{sample}_spades/contigs.fasta", sample=SAMPLES)

rule run_bbduk_qc:
    input:
        R1 = FASTQ + "{sample}_1.fq.gz",
        R2 = FASTQ + "{sample}_2.fq.gz",

    output:
        R1 = PROCESS + "clean_reads/{sample}.trimmed_1.fastq.gz",
        R2 = PROCESS + "clean_reads/{sample}.trimmed_2.fastq.gz"

    params:
        adaptors = REFERENCES + "adapters.fa"

    threads: 10

    shell:
        "bbduk.sh in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} threads={threads} ref={params.adaptors} qtrim=rl trimq=20 ktrim=r k=23 mink=11 hdist=1 tpe tbo"

rule run_bbmap:
    input:
        R1 = PROCESS + "clean_reads/{sample}.trimmed_1.fastq.gz",
        R2 = PROCESS + "clean_reads/{sample}.trimmed_2.fastq.gz"

    output:
        R1 = PROCESS + "clean_reads/{sample}.clean_1.fastq.gz",
        R2 = PROCESS + "clean_reads/{sample}.clean_2.fastq.gz"

    params:
        filter_db = REFERENCES

    threads: 10

    shell:
        "bbmap.sh path={params.filter_db} minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 in={input.R1} in2={input.R2} outu={output.R1} outu2={output.R2} threads={threads}"

rule run_megahit_assembly:
    input:
        R1 = expand(PROCESS + "clean_reads/{sample}.clean_1.fastq.gz", sample=SAMPLES),
        R2 = expand(PROCESS + "clean_reads/{sample}.clean_2.fastq.gz", sample=SAMPLES)

    output:
        PROCESS + "megahit_assembly/BT1.contigs.fa"

    params:
        folder = PROCESS + "megahit_assembly",
        prefix = "BT1",
        R1 = expand("-1 " + PROCESS + "clean_reads/{sample}.clean_1.fastq.gz", sample=SAMPLES),
        R2 = expand("-2 " + PROCESS + "clean_reads/{sample}.clean_2.fastq.gz", sample=SAMPLES)

    threads: 20

    shell:
        """
        rm -rf {params.folder}
        megahit {params.R1} {params.R2} --kmin-1pass -m 0.5 --k-list 27,37,47,57,67,77,87 --min-contig-len 300 -t {threads} -o {params.folder} --out-prefix {params.prefix}
        """
rule run_spades_assembly:
    input:
        R1 = PROCESS + "clean_reads/{sample}.trimmed_1.fastq.gz",
        R2 = PROCESS + "clean_reads/{sample}.trimmed_2.fastq.gz"

    output:
        PROCESS + "{sample}_spades/contigs.fasta"

    params:
        folder = PROCESS + "{sample}_spades"

    threads: 10

    shell:
        """
        spades.py -1 {input.R1} -2 {input.R2} -o {params.folder} --meta --threads {threads} --memory 90
        """
