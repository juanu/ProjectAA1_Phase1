# Snakefile
# This is designed to use the cluster, so it may not be the completely portable
# But it provides an idea and a learning place for me to learn Snakemake :)
# Ugalde 2019

from collections import defaultdict

#FOLDERS]
FASTQ = "data/fastq/"
PROCESS = "data/process/"
REFERENCES = "data/references/"
RESULTS = "results/"


# Get samples names from the metadata
# This should go maybe on the config file?

SAMPLES = ["Sample1", "Sample2"]

rule all:
    input:
        PROCESS + "anvio_data/BT1_megahit.contigs.fa",
        PROCESS + "anvio_data/BT1_spades.contig.fa",
        "metabat2_done.check",
        #expand("PhyloFlash-{sample}.phyloFlash.html", sample=SAMPLES),
        expand(RESULTS + "humann2/{sample}_genefamilies.tsv", sample=SAMPLES),
        RESULTS + "checkm_results/CheckM.txt"


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

rule process_contigs:
    input:
        contigs_megahit = PROCESS + "megahit_assembly/BT1.contigs.fa",
        contigs_spades = PROCESS + "spades_assembly/contigs.fasta"

    output:
        filtered_contigs_megahit = PROCESS + "anvio_data/BT1_megahit.contigs.fa",
        filtered_contigs_spade = PROCESS + "anvio_data/BT1_spades.contig.fa"

    params:
        megahit_name_file = PROCESS + "anvio_data/megahit_name_conversions.txt",
        spades_name_file = PROCESS + "anvio_data/spades_name_conversions.txt"

    conda:
        "anvio.yml"

    shell:
        """
        anvi-script-reformat-fasta {input.contigs_megahit} -o {output.filtered_contigs_megahit} --min-len 2500 --simplify-names --report {params.megahit_name_file}
        anvi-script-reformat-fasta {input.contigs_spades} -o {output.filtered_contigs_spade} --min-len 2500 --simplify-names --report {params.spades_name_file}
        """

rule make_bw_index:
    input:
        contigs = PROCESS + "anvio_data/BT1_spades.contig.fa"

    output:
        db = PROCESS + "anvio_data/BT1_spades.bowtie2.1.bt2"

    shell:
        """
        bowtie2-build {input.contigs} {output.db}
        """

rule map_reads:
    input:
        R1 = PROCESS + "clean_reads/{sample}.clean_1.fastq.gz",
        R2 = PROCESS + "clean_reads/{sample}.clean_2.fastq.gz",
        db = PROCESS + "anvio_data/BT1_spades.bowtie2.1.bt2"

    output:
        mapped = PROCESS + "anvio_data/{sample}.sorted.bam"

    params:
        out_sam = PROCESS + "anvio_data/{sample}.sam",
        temp_bam = PROCESS + "anvio_data/{sample}-RAW.bam",
        db = PROCESS + "anvio_data/BT1_spades.bowtie2"

    threads: 20

    shell:
        """
        bowtie2 --threads {threads} -x {params.db} -1 {input.R1} -2 {input.R2} -S {params.out_sam}
        samtools view -F 4 -bS {params.out_sam} > {params.temp_bam}
        samtools sort -o {output.mapped} {params.temp_bam}
        samtools index {output.mapped}
        rm {params.out_sam}
        rm {params.temp_bam}
        """

#change for phyloflash
rule run_phyloflash:
    input:
        R1 = PROCESS + "clean_reads/{sample}.clean_1.fastq.gz",
        R2 = PROCESS + "clean_reads/{sample}.clean_2.fastq.gz"

    output:
        "PhyloFlash-{sample}.phyloFlash.html"

    params:
        db = "/hpcudd/home/jugalde/storage/databases/phyloflash/132",
        lib = "PhyloFlash-{sample}"

    threads:20

    conda:
        "phyloflash.yml"

    shell:
        """
        unset MAFFT_BINARIES
        export MAFFT_BINARIES=~/miniconda2/envs/phyloflash/bin/mafft
        phyloFlash.pl -dbhome {params.db} -lib {params.lib} -read1 {input.R1} -read2 {input.R2} -CPUS {threads} -log -emirge -poscov
        """

# Running Metabat2
rule run_metabat2:
    input:
        bam_files = expand(PROCESS + "anvio_data/{sample}.sorted.bam", sample=SAMPLES),
        assembly =  PROCESS + "anvio_data/BT1_spades.contig.fa"

    output:
        "metabat2_done.check"

    shell:
        """
        runMetaBat.sh {input.assembly} {input.bam_files}
        touch metabat2_done.check
        """

# Running this manually
rule run_checkm:
    input:
        "metabat2_done.check"

    output:
        summary = RESULTS + "checkm_results/CheckM.txt"

    conda:
        "checkm.yml"

    params:
        bins = "BT1_spades.contig.fa.metabat-bins",
        output_folder = RESULTS + "checkm_results"

    threads: 20

    shell:
        """
        checkm data setRoot /hpcudd/home/jugalde/storage/databases/checkm
        checkm lineage_wf -f {output.summary} -t {threads} -x fa {params.bins} {params.output_folder}
        """

# Run Humman2

rule concatenate_reads:
    input:
        R1 = PROCESS + "clean_reads/{sample}.clean_1.fastq.gz",
        R2 = PROCESS + "clean_reads/{sample}.clean_2.fastq.gz"

    output:
        R1R2 = PROCESS + "clean_reads/{sample}.clean_R1R2.fastq.gz"

    shell:
        """
        cat {input.R1} {input.R2} > {output.R1R2}
        """

rule run_humann2:
    input:
        R1R2 = PROCESS + "clean_reads/{sample}.clean_R1R2.fastq.gz"

    output:
        RESULTS + "humann2/{sample}_genefamilies.tsv"

    params:
        output_folder = RESULTS + "humann2",
        sample_name = "{sample}",
        nuc_db = "/hpcudd/home/jugalde/storage/databases/humann2_dbs/chocophlan",
        prot_db = "/hpcudd/home/jugalde/storage/databases/humann2_dbs/uniref",
        metaphlan_pkl = "/hpcudd/home/jugalde/storage/databases/humann2_dbs/metaphlan2/mpa_v20_m200.pkl",
        metaphlan_bwt = "/hpcudd/home/jugalde/storage/databases/humann2_dbs/metaphlan2/"

    threads: 20

    conda:
        "humann2.yml"

    shell:
        """
        humann2 --input {input.R1R2} --output {params.output_folder} --output-basename {params.sample_name} --nucleotide-database {params.nuc_db} --protein-database {params.prot_db} --threads {threads} --metaphlan-options "--mpa_pkl {params.metaphlan_pkl} --bowtie2db {params.metaphlan_bwt}"
        """
