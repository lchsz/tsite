import os
import stat

import pandas as pd


configfile: "data/config/config.yml"


def prepare_dir():
    data_dir = "data"
    temp_dir = os.path.join(data_dir,"temp")
    trim_dir = os.path.join(temp_dir,"trim_reads")
    aln_dir = os.path.join(temp_dir,"aln")
    load_aln_dir = os.path.join(temp_dir,"load_aln")
    filter_aln_dir = os.path.join(temp_dir,"filter_aln")
    filter_fq_dir = os.path.join(temp_dir,"filter_fq")
    align_vector_dir = os.path.join(temp_dir,"align_vector")
    load_vector_aln_dir = os.path.join(temp_dir,"load_vector_aln")
    align_host_dir = os.path.join(temp_dir,"align_host")
    load_host_aln_dir = os.path.join(temp_dir,"load_host_aln")
    bam_dir = os.path.join(temp_dir,"bam")
    result_dir = os.path.join(data_dir,"result")

    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    if not os.path.exists(trim_dir):
        os.mkdir(trim_dir)
    if not os.path.exists(aln_dir):
        os.mkdir(aln_dir)
    if not os.path.exists(load_aln_dir):
        os.mkdir(load_aln_dir)
    if not os.path.exists(filter_aln_dir):
        os.mkdir(filter_aln_dir)
    if not os.path.exists(filter_fq_dir):
        os.mkdir(filter_fq_dir)
    if not os.path.exists(align_vector_dir):
        os.mkdir(align_vector_dir)
    if not os.path.exists(load_vector_aln_dir):
        os.mkdir(load_vector_aln_dir)
    if not os.path.exists(align_host_dir):
        os.mkdir(align_host_dir)
    if not os.path.exists(load_host_aln_dir):
        os.mkdir(load_host_aln_dir)
    if not os.path.exists(bam_dir):
        os.mkdir(bam_dir)
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    os.chmod(result_dir,stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
    os.chmod(temp_dir,stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)


def get_vector(wildcards):
    return os.path.join("data/resource",samples['vector'][wildcards.sample])


def get_host(wildcards):
    return os.path.join("data/resource",samples['host'][wildcards.sample])


def index_ref(ref):
    indexed_files = [f"data/resource/{ref}.{ext}" for ext in ["amb", "ann", "bwt", "pac", "sa"]]
    indexed = True
    for index_file in indexed_files:
        if not os.path.exists(index_file):
            indexed = False
            break

    if not indexed:
        os.system(f"bwa index data/resource/{ref}")


def start_on():
    prepare_dir()
    for sample in samples.index:
        index_ref(samples['host'][sample])
        index_ref(samples['vector'][sample])


samples = pd.read_table(os.path.join("data/config",config["samples"])).set_index("sample",drop=False)
start_on()


rule all:
    input:
        expand("data/result/{sample}_site.html",sample=samples.index)


rule trim_reads:
    input:
        "data/fastq/{sample}_R1.fastq.gz",
        "data/fastq/{sample}_R2.fastq.gz"
    output:
        "data/temp/trim_reads/{sample}_R1.fastq.gz",
        "data/temp/trim_reads/{sample}_R2.fastq.gz",
        "data/temp/trim_reads/{sample}_fastp.html",
        "data/temp/trim_reads/{sample}_fastp.json"
    shell:
        "fastp --dedup -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} "
        "-h {output[2]} -j {output[3]}"


rule align:
    input:
        get_vector,
        "data/temp/trim_reads/{sample}_R1.fastq.gz",
        "data/temp/trim_reads/{sample}_R2.fastq.gz",
    output:
        "data/temp/aln/{sample}_R1.bam",
        "data/temp/aln/{sample}_R2.bam"
    threads:
        config["bwa"]["threads"]
    shell:
        "bwa mem -t {threads} {input[0]} {input[1]} |samtools view -bo {output[0]} - && "
        "bwa mem -t {threads} {input[0]} {input[2]} |samtools view -bo {output[1]} -"


rule load_aln:
    input:
        "data/temp/aln/{sample}_R1.bam",
        "data/temp/aln/{sample}_R2.bam"
    output:
        "data/temp/load_aln/{sample}_R1.aln",
        "data/temp/load_aln/{sample}_R2.aln"
    shell:
        "java -jar script/hts.jar load_aln {input[0]} {output[0]} && "
        "java -jar script/hts.jar load_aln {input[1]} {output[1]}"


rule filter_aln:
    input:
        "data/temp/load_aln/{sample}_R1.aln",
        "data/temp/load_aln/{sample}_R2.aln"
    output:
        "data/temp/filter_aln/{sample}_read_name.ids"
    script:
        "script/filter_vector_aln.R"


rule filter_fq:
    input:
        "data/temp/filter_aln/{sample}_read_name.ids",
        "data/temp/trim_reads/{sample}_R1.fastq.gz",
        "data/temp/trim_reads/{sample}_R2.fastq.gz",
    output:
        "data/temp/filter_fq/{sample}_R1.fastq",
        "data/temp/filter_fq/{sample}_R2.fastq"
    shell:
        "java -jar script/hts.jar filter_fq {input[0]} {input[1]} {output[0]} &&"
        "java -jar script/hts.jar filter_fq {input[0]} {input[2]} {output[1]}"


rule align_vector:
    input:
        get_vector,
        "data/temp/filter_fq/{sample}_R1.fastq",
        "data/temp/filter_fq/{sample}_R2.fastq",
    output:
        "data/temp/align_vector/{sample}_R1.sam",
        "data/temp/align_vector/{sample}_R2.sam"
    threads:
        config["bwa"]["threads"]
    shell:
        "bwa mem -t {threads} -o {output[0]} {input[0]} {input[1]} && "
        "bwa mem -t {threads} -o {output[1]} {input[0]} {input[2]}"


rule load_vector_aln:
    input:
        "data/temp/align_vector/{sample}_R1.sam",
        "data/temp/align_vector/{sample}_R2.sam"
    output:
        "data/temp/load_vector_aln/{sample}_R1.aln",
        "data/temp/load_vector_aln/{sample}_R2.aln"
    shell:
        "java -jar script/hts.jar load_aln {input[0]} {output[0]} && "
        "java -jar script/hts.jar load_aln {input[1]} {output[1]}"


rule align_host:
    input:
        get_host,
        "data/temp/filter_fq/{sample}_R1.fastq",
        "data/temp/filter_fq/{sample}_R2.fastq",
    output:
        "data/temp/align_host/{sample}_R1.sam",
        "data/temp/align_host/{sample}_R2.sam"
    threads:
        config["bwa"]["threads"]
    shell:
        "bwa mem -t {threads} -o {output[0]} {input[0]} {input[1]} && "
        "bwa mem -t {threads} -o {output[1]} {input[0]} {input[2]}"


rule load_host_aln:
    input:
        "data/temp/align_host/{sample}_R1.sam",
        "data/temp/align_host/{sample}_R2.sam"
    output:
        "data/temp/load_host_aln/{sample}_R1.aln",
        "data/temp/load_host_aln/{sample}_R2.aln"
    shell:
        "java -jar script/hts.jar load_aln {input[0]} {output[0]} && "
        "java -jar script/hts.jar load_aln {input[1]} {output[1]}"


rule filter_site_aln:
    input:
        "data/temp/load_vector_aln/{sample}_R1.aln",
        "data/temp/load_vector_aln/{sample}_R2.aln",
        "data/temp/load_host_aln/{sample}_R1.aln",
        "data/temp/load_host_aln/{sample}_R2.aln"
    output:
        "data/result/{sample}_site_aln.tsv"
    script:
        "script/filter_site_aln.R"


rule detect_site:
    input:
        "data/result/{sample}_site_aln.tsv"
    output:
        "data/result/{sample}_site.tsv"
    script:
        "script/detect_site.py"


rule filter_sam:
    input:
        "data/result/{sample}_site.tsv",
        "data/result/{sample}_site_aln.tsv",
        "data/temp/align_vector/{sample}_R1.sam",
        "data/temp/align_vector/{sample}_R2.sam",
        "data/temp/align_host/{sample}_R1.sam",
        "data/temp/align_host/{sample}_R2.sam"
    output:
        directory("data/temp/bam/{sample}")
    script:
        "igv/filter_sam.py"


rule report_site:
    input:
        site_file="data/result/{sample}_site.tsv",
        site_aln_file="data/result/{sample}_site_aln.tsv",
        host_file=get_host,
        vector_file=get_vector,
        bam_dir="data/temp/bam/{sample}"
    output:
        "data/result/{sample}_site.html"
    script:
        "igv/report_site.py"
