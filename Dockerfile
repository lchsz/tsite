FROM lchsz/ubuntu

WORKDIR /tsite
COPY igv ./igv
COPY script ./script
COPY Snakefile ./

# python
RUN pip3 install pandas snakemake pysam Jinja2 \
        -i https://pypi.tuna.tsinghua.edu.cn/simple \
        --default-timeout=100

# fastp
RUN wget -O fastp http://opengene.org/fastp/fastp.0.23.4 \
    && chmod +x fastp \
    && mv fastp /usr/local/bin/

# bwa, samtools
RUN apt -y --no-install-recommends install bwa samtools

CMD ["snakemake", "--core", "4"]