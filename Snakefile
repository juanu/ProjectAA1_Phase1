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
        PROCESS + "anvio_data/BT1_final.contigs.fa"

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
        R1 = expand(PROCESS + "clean_reads/{sample}.clean_1.fastq.gz", sample=SAMPLES),
        R2 = expand(PROCESS + "clean_reads/{sample}.clean_2.fastq.gz", sample=SAMPLES)

    output:
        conc_R1 = "clean_reads/conc_R1.fastq.gz",
        conc_R2 = "clean_reads/conc_R2.fastq.gz",
        output_contigs = PROCESS + "spades_assembly/contigs.fasta"

    params:
        folder = PROCESS + "spades_assembly"

    threads: 22

    shell:
        """
        cat {input.R1} > {output.conc_R1}
        cat {input.R2} > {output.conc_R2}
        spades.py -1 {output.conc_R1} -2 {output.conc_R2} -o {params.folder} --meta --threads {threads} --memory 90
        """

rules process_contigs:
    input:
        contigs_megahit = PROCESS + "megahit_assembly/BT1.contigs.fa",
        contigs_spades = PROCESS + "spades_assembly/contigs.fasta"

    output:
        filtered_contigs_megahit = PROCESS + "anvio_data/BT1_megahit.contigs.fa",
        filtered_contigs_spade = PROCESS + "anvio_data/BT1_spades.contig.fa"

    params:
        megahit_name_file = PROCESS + "anvio_data/megahit_name_conversions.txt",
        spades_name_file = PROCESS + "anvio_data/spades_name_conversions.txt"

    shell:
        """
        anvi-script-reformat-fasta {input.contigs_megahit} -o {output.filtered_contigs_megahit} --min-len 2500 --simplify-names --report {params.megahit_name_file}
        anvi-script-reformat-fasta {input.contigs_spades} -o {output.filtered_contigs_spade} --min-len 2500 --simplify-names --report {params.spades_name_file}
        """

#rules map_reads:
#    input:
#        contigs = PROCESS + "anvio_data/BT1_final.contigs.fa"
