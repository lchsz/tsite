suppressMessages(library(dplyr))
suppressMessages(library(purrr))

filter_aln <- function(alns,
                       min_anchor_len = 10,
                       min_local_identity = 0.9,
                       min_align_start = 10) {
  end_mis <- 5
  alns %>% filter(
    (q_end - q_start + 1) > min_anchor_len &
      t_start > min_align_start &
      (t_size - t_end + 1) > min_align_start &
      local_identity > min_local_identity
  ) %>% filter(
    (q_start > min_anchor_len & (q_size - q_end) < end_mis) |
      (q_size - q_end > min_anchor_len & q_start < end_mis)
  )
}


find_site <- function(r_vector_aln,
                      r_host_aln,
                      r_host_full_aln,
                      max_microhomology_len,
                      max_filler_len) {
  alns <- anti_join(r_vector_aln, r_host_full_aln, by = "q_name")
  alns <- left_join(alns, r_host_aln,
    by = "q_name",
    relationship = "many-to-many",
    suffix = c("_v", "_h")
  )

  alns1 <- alns %>%
    filter(
      (q_end_h - q_start_v >= 0 & q_end_h - q_start_v < max_microhomology_len) |
        (q_end_v - q_start_h >= 0 & q_end_v - q_start_h < max_microhomology_len)
    ) %>%
    mutate(type = "microhomology", insert_len = if_else(
      q_end_h - q_start_v >= 0 & q_end_h - q_start_v < max_microhomology_len,
      q_start_v - q_end_h, q_start_h - q_end_v
    ))

  alns2 <- alns %>%
    filter(
      (q_start_v - q_end_h >= 0 & q_start_v - q_end_h < max_filler_len) |
        (q_start_h - q_end_v >= 0 & q_start_h - q_end_v < max_filler_len)
    ) %>%
    mutate(type = "filler", insert_len = if_else(
      q_start_v - q_end_h >= 0 & q_start_v - q_end_h < max_filler_len,
      q_start_v - q_end_h, q_start_h - q_end_v
    ))

  bind_rows(alns1, alns2)
}
