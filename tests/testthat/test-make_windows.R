# Helpers ---
read_annot <- function(file_path) {
  df <- read.delim(file_path, header = F)
  names(df) <- c("seqnames", "start", "end", "gene_id", "score", "strand")
  return(df)
}

defactorize <- function(df) {
  df %>% dplyr::mutate(across(where(is.factor), as.character))
}

test_that("minimal arguments with gtf", {
  expect_equal(make_windows(input_file="data/example.gtf",
                            upstream = 1000, downstream = 2000,
                            feature = NULL, biotype = NULL)$result %>% defactorize(),
               read_annot("out_ref/minimal_gtf.bed"))
  expect_equal(make_windows(input_file="data/example_annotation_1000up_2000down.bed",
                            upstream = 1000, downstream = 2000,
                            feature = NULL, biotype = NULL)$result %>% defactorize(),
               read_annot("out_ref/minimal_bed.bed"))
})

test_that("---filter wrt feature=\"gene\" biotype=\"lncRNA\"", {
  expect_equal(make_windows(input_file="data/example.gtf",
                            upstream = 1000, downstream = 2000,
                            feature = "gene", biotype = "lncRNA")$result %>% defactorize(),
               read_annot("out_ref/example_gtf_gene_lncRNA.bed"))
})

test_that("---edge case with edge make_windows", {
  expect_equal(make_windows(input_file="data/example_edges.gtf",
                            upstream = 1000, downstream = 2000,
                            feature = "gene", biotype = "lncRNA")$result %>% defactorize(),
               read_annot("out_ref/example_gtf_gene_lncRNA_edges.bed"))
})

test_that("---with position=start", {
  expect_equal(make_windows(input_file="data/example.gtf",
                            upstream = 1000, downstream = 2000,
                            feature = "gene", biotype = "lncRNA",
                            position = "start")$result %>% defactorize(),
               read_annot("out_ref/example_gtf_gene_lncRNA_position=start.bed"))
  expect_equal(make_windows(input_file="data/example.gtf",
                            upstream = 1000, downstream = 2000,
                            feature = "gene", biotype = "lncRNA",
                            position = "end")$result %>% defactorize(),
               read_annot("out_ref/example_gtf_gene_lncRNA_position=end.bed"))
})


