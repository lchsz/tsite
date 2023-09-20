suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(tibble))

script_dir <- "script"
source(file.path(script_dir, "util.R"))

r1_vector_aln_file <- snakemake@input[[1]]
r2_vector_aln_file <- snakemake@input[[2]]
read_names_file <- snakemake@output[[1]]

min_local_identity <- snakemake@config[["min_local_identity"]]

min_anchor_len <- 10

r1_vector_aln <- read_tsv(r1_vector_aln_file) %>%
    filter_aln(min_anchor_len, min_local_identity)
r2_vector_aln <- read_tsv(r2_vector_aln_file) %>%
    filter_aln(min_anchor_len, min_local_identity)
read_names <- union(r1_vector_aln$q_name, r2_vector_aln$q_name)

write_tsv(as_tibble(read_names), file = read_names_file, col_names = FALSE)
