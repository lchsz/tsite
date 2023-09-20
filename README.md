# TSite: A Dockerized, End-to-end Pipeline for Fast Identifying T-DNA Integration Sites in Plant Genomes

The following will introduce the TSite configuration and usage using the data directory in the code source repository.

## Datasets
A data directory should be presented as the following structure: raw reads, host and vector references, and configure files.

```sh
data
├── config
│   ├── config.yml
│   └── samples.tsv
├── fastq
│   ├── SRR16230416_R1.fastq.gz
│   ├── SRR16230416_R2.fastq.gz
│   ├── SRR16230426_R1.fastq.gz
│   └── SRR16230426_R2.fastq.gz
└── resource
    ├── 46bg.fasta
    ├── 49btm.fasta
    └── osa.fasta
```

### Read files
* `fastq/SRR16230416_R1.fastq.gz` and `SRR16230416_R2.fastq.gz` are raw reads from NGS sequencing for sample `SRR16230416`.
Input files must be in FASTQ format, with ".fastq.gz" extensions supported.
* Host and vector reference FASTA files are placed under the `resource` subfold.

### Configure file
* `config/config.yml` contains processing parameters:
    + `samples: samples.tsv` specify the sample information file. 
    + `max_filler_len: 10` maximum length of small DNA fragments insertion at T-DNA integration sites.
    + `max_microhomology_len: 5` maximum number of microhomology at the integration site.
    + `min_local_identity: 0.9` is the minimum percent identity of matched bases in aligned regions.
    + `min_reads_num: 3` minimum number of supporting informative split reads.
    + `cluster_dist: 5` maximum distance used to cluster sites into the same integration site.
    + `bwa:` 
        - `threads: 8` number of BWA threads.

* `config/samples.tsv` contains sample metadata:
    ```tsv
    sample  host    vector
    SRR16230416 osa.fasta   46bg.fasta
    SRR16230426 osa.fasta   49btm.fasta
    ```

## Running TSite
TSite is built as a Docker container image and shared through Docker Hub (https://hub.docker.com/repository/docker/lchsz/tsite/general). The container can run on various operating systems, such as Linux, Mac, and Windows, all supported by the Docker engine. 

### Requirements
Running the pipeline requires a recent version of the docker client, found on the docker website.
Follow the instructions to [download and install Docker](https://docs.docker.com/get-docker/)

### Run a TSite container
The data required by running TSite are stored on the host's file system. 
When running a container of TSite, we must use a bind mount to mount the data directory into the container.
We can configure and place the reads files into a data directory outside the container, which can watch the changes you make as soon as you save a file.

Before running the application, specify the parameter and place all the read files in their directories as the instructions mentioned.

```sh
docker run -v data-dir:/tsite/data lchsz/tsite
```

The –v option tells Docker to create a bind mount, where `data-dir` is the directory on your host machine, and `/tsite/data` would appear inside the container.


