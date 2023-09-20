suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))

source(file.path("script", "util.R"))

r1_vector_aln_file <- snakemake@input[[1]]
r2_vector_aln_file <- snakemake@input[[2]]
r1_host_aln_file <- snakemake@input[[3]]
r2_host_aln_file <- snakemake@input[[4]]
site_read_file <- snakemake@output[[1]]

max_microhomology_len <- snakemake@config[["max_microhomology_len"]]
max_filler_len <- snakemake@config[["max_filler_len"]]
min_local_identity <- snakemake@config[["min_local_identity"]]

min_global_identity <- 0.9
min_anchor_len <- 10

r1_vector_aln <- read_tsv(r1_vector_aln_file)
r2_vector_aln <- read_tsv(r2_vector_aln_file)
r1_host_aln <- read_tsv(r1_host_aln_file)
r2_host_aln <- read_tsv(r2_host_aln_file)

r1_host_full_aln <-
  filter(r1_host_aln, global_identity > min_global_identity)
r2_host_full_aln <-
  filter(r2_host_aln, global_identity > min_global_identity)

r1_host_aln <- filter_aln(r1_host_aln, min_anchor_len, min_local_identity)
r2_host_aln <- filter_aln(r2_host_aln, min_anchor_len, min_local_identity)

r1_vector_aln <- filter_aln(r1_vector_aln, min_anchor_len, min_local_identity)
r2_vector_aln <- filter_aln(r2_vector_aln, min_anchor_len, min_local_identity)

sites1 <- find_site(
  r1_vector_aln, r1_host_aln, r1_host_full_aln,
  max_microhomology_len, max_filler_len
)
sites2 <- find_site(
  r2_vector_aln, r2_host_aln, r2_host_full_aln,
  max_microhomology_len, max_filler_len
)

sites1$read_type <- "R1"
sites2$read_type <- "R2"
sites <- bind_rows(sites1, sites2)
sites$pos_h <- ifelse(
  (sites$strand_h == "+" &
    sites$strand_v == "+" &
    sites$q_start_h > sites$q_start_v) |
    (sites$strand_h == "-" &
      sites$strand_v == "-" &
      sites$q_start_h < sites$q_start_v) |
    (sites$strand_h == "+" &
      sites$strand_v == "-" &
      sites$q_start_h > sites$q_start_v) |
    (sites$strand_h == "-" &
      sites$strand_v == "+" &
      sites$q_start_h < sites$q_start_v),
  sites$t_start_h, sites$t_end_h
)

sites$pos_v <- ifelse(
  (sites$strand_h == "+" &
    sites$strand_v == "+" &
    sites$q_start_h < sites$q_start_v) |
    (sites$strand_h == "-" &
      sites$strand_v == "-" &
      sites$q_start_h > sites$q_start_v) |
    (sites$strand_h == "+" &
      sites$strand_v == "-" &
      sites$q_start_h > sites$q_start_v) |
    (sites$strand_h == "-" &
      sites$strand_v == "+" &
      sites$q_start_h < sites$q_start_v),
  sites$t_start_v, sites$t_end_v
)
sites <- arrange(sites, t_name_h, pos_h, t_name_v, pos_v, type)
sites <- cbind(`id` = seq_len(nrow(sites)) - 1, sites)

write_tsv(sites, file = site_read_file, col_names = TRUE)
