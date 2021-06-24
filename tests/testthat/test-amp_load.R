test_that("loading example data with default options returns an ampvis2 class object", {
  d <- amp_load(
    example_otutable,
    example_metadata
  )
  expect_s3_class(d, class = "ampvis2", exact = TRUE)
})

test_that("loading example data with separate taxonomy returns identical object", {
  expect_identical(
    amp_load(
      otutable = example_otutable,
      metadata = example_metadata,
      taxonomy = example_taxonomy
    ),
    amp_load(
      otutable = example_otutable,
      metadat = example_metadata
    )
  )
})

test_that("loading otutable with both 'OTU' and 'ASV' columns fails", {
  otutable <- example_otutable
  otutable$ASV <- rownames(example_otutable)
  expect_error(
    amp_load(
      otutable = otutable,
      metadata = example_metadata,
      taxonomy = NULL,
      fasta = NULL,
      tree = NULL,
      pruneSingletons = FALSE
    ),
    regexp = "More than one column in otutable is named OTU/ASV, don't know which one to use."
  )
})

test_that("loading otutable with 'ASV' column return an ampvis2 class object", {
  otutable <- example_otutable
  colnames(otutable)[16] <- "ASV"
  d <- amp_load(
    otutable = otutable,
    metadata = example_metadata,
    taxonomy = NULL,
    fasta = NULL,
    tree = NULL,
    pruneSingletons = FALSE
  )
  expect_s3_class(d, class = "ampvis2", exact = TRUE)
})

test_that("loading otutable with '#OTU ID' column return an ampvis2 class object", {
  otutable <- example_otutable
  colnames(otutable)[16] <- "#OTU ID"
  d <- amp_load(
    otutable = otutable,
    metadata = example_metadata,
    taxonomy = NULL,
    fasta = NULL,
    tree = NULL,
    pruneSingletons = FALSE
  )
  expect_s3_class(d, class = "ampvis2", exact = TRUE)
})

test_that("warning if no OTU column in otutable", {
  testthat::expect_warning(
    amp_load(
      otutable = example_otutable[,-16],
      metadata = example_metadata,
      taxonomy = NULL,
      fasta = NULL,
      tree = NULL,
      pruneSingletons = FALSE
    ),
    regexp = "Could not find a column named OTU/ASV in otutable, using rownames as OTU ID's")
})

test_that("loading separate taxonomy table returns an ampvis2 class object", {
  otutable <- example_otutable[,1:8, drop = FALSE]
  otutable$OTU <- rownames(otutable)
  taxonomy <- example_otutable[,9:16, drop = FALSE]
  d <- amp_load(
    otutable = otutable,
    metadata = example_metadata,
    taxonomy = taxonomy,
    fasta = NULL,
    tree = NULL,
    pruneSingletons = FALSE
  )
  expect_s3_class(d, class = "ampvis2", exact = TRUE)
})

test_that("loading taxonomy with less than 7 levels returns an ampvis2 class object", {
  otutable <- example_otutable[,-15]
  d <- amp_load(
    otutable = otutable,
    metadata = example_metadata,
    taxonomy = NULL,
    fasta = NULL,
    tree = NULL,
    pruneSingletons = FALSE
  )
  expect_s3_class(d, class = "ampvis2", exact = TRUE)
})

test_that("loading data with no taxonomy creates dummy taxonomy table and throws warning", {
  otutable <- example_otutable[,1:8]
  otutable$OTU <- rownames(otutable)
  expect_warning(
    amp_load(
      otutable = otutable,
      metadata = example_metadata,
      taxonomy = NULL,
      fasta = NULL,
      tree = NULL,
      pruneSingletons = FALSE
    ),
    regexp = "Could not find or parse taxonomy, creating a dummy taxonomy table with only OTUs"
  )
})

test_that("loading otutable with more samples than metadata returns a warning", {
  expect_warning(
    amp_load(
      otutable = example_otutable,
      metadata = example_metadata[1:5,],
      taxonomy = NULL,
      fasta = NULL,
      tree = NULL,
      pruneSingletons = FALSE
    ),
    regexp = "Only 5 of 8 unique sample names match between metadata and otutable"
  )
})

test_that("loading data with only rownames containing OTUs", {
  otutable <- example_otutable[,1:8, drop = FALSE]
  taxonomy <- example_otutable[,9:15, drop = FALSE]
  rownames(otutable) <- NULL
  rownames(taxonomy) <- NULL
  suppressWarnings(
    d <- amp_load(
      otutable = otutable,
      metadata = example_metadata,
      taxonomy = taxonomy
    )
  )
  
  expect_setequal(d$tax$OTU, rownames(d$abund))
  expect_true(nrow(d$tax) == nrow(d$abund))
  expect_s3_class(
    d,
    class = "ampvis2",
    exact = TRUE
  )
})

test_that("loading data with more OTUs in taxonomy than in otutable", {
  otutable <- example_otutable[1:5,c(1:8, 16), drop = FALSE]
  taxonomy <- example_otutable[,9:16, drop = FALSE]
  expect_warning(
    amp_load(
      otutable = otutable, 
      taxonomy = taxonomy, 
      metadata = example_metadata
    ),
    regexp = "The OTU's between otutable and taxonomy do not match exactly. 5 OTU's in taxonomy not present in otutable have been removed from taxonomy."
  )
})

test_that("loading data with more OTUs in otutable than in taxonomy", {
  otutable <- example_otutable[,c(1:8,16)]
  taxonomy <- example_otutable[5:10,9:16, drop = FALSE]
  expect_warning(
    amp_load(
      otutable = otutable, 
      taxonomy = taxonomy, 
      metadata = example_metadata
    ),
    regexp = "The OTU's between otutable and taxonomy do not match exactly. 4 OTU's are missing "
  )
})

test_that("loading data with unique OTU's in both otutable and taxonomy", {
  otutable <- example_otutable[1:7, c(1:8, 16), drop = FALSE]
  taxonomy <- example_otutable[5:10, 9:16, drop = FALSE]
  expect_warning(
    amp_load(
      otutable = otutable, 
      taxonomy = taxonomy, 
      metadata = example_metadata
    ),
    regexp = "The OTU's between otutable and taxonomy do not match exactly. 3 OTU's in taxonomy not present in otutable have been removed from taxonomy. 4 OTU's are missing "
  )
})

test_that("loading a (non-hdf5 format) BIOM file with no taxonomy fails", {
  expect_error(
    amp_load("../testdata/min_sparse_otu_table.biom")
  )
})

test_that("loading a (hdf5 format) BIOM file with no taxonomy fails", {
  expect_error(
    amp_load("../testdata/min_sparse_otu_table_hdf5.biom")
  )
})

test_that("loading a (non-hdf5 format) BIOM file returns an ampvis2 class object", {
  expect_s3_class(
    suppressWarnings(amp_load("../testdata/rich_sparse_otu_table.biom")),
    class = "ampvis2",
    exact = TRUE
  )
})

test_that("loading a (hdf5 format) BIOM file returns an ampvis2 class object", {
  expect_s3_class(
    suppressWarnings(amp_load("../testdata/rich_sparse_otu_table_hdf5.biom")),
    class = "ampvis2",
    exact = TRUE
  )
})

test_that("loading hdf5 and non-hdf5 format BIOM file returns identical ampvis2 class objects", {
  biom_nonhdf5 <- suppressWarnings(amp_load("../testdata/rich_sparse_otu_table.biom"))
  biom_hdf5 <- suppressWarnings(amp_load("../testdata/rich_sparse_otu_table_hdf5.biom"))
  expect_identical(biom_hdf5, biom_nonhdf5)
})

test_that("loading sintax format taxonomy returns an ampvis2 class object", {
  expect_s3_class(
    suppressWarnings(
      amp_load("tests/testdata/ASVtable.tsv", taxonomy = "tests/testdata/ASVs.sintax")
    ),
    class = "ampvis2",
    exact = TRUE
  )
})
