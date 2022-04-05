test_that("Merging a single object returns an ampvis2 class object", {
  expect_s3_class(
    amp_merge_ampvis2(
      AalborgWWTPs
    ),
    class = "ampvis2",
    exact = TRUE
  )
})

test_that("Merging two ampvis2-class objects originating from the same data set returns an ampvis2 class object", {
  suppressMessages({
    obj1 <- amp_filter_samples(
      AalborgWWTPs,
      Year %in% "2014"
    )
    obj2 <- amp_filter_samples(
      AalborgWWTPs,
      Year %in% "2015"
    )
  })

  res <- expect_s3_class(
    amp_merge_ampvis2(
      obj1,
      obj2,
      by_refseq = FALSE
    ),
    class = "ampvis2",
    exact = TRUE
  )

  suppressMessages({
    compare <- amp_filter_samples(AalborgWWTPs, Year %in% c("2014", "2015"))
  })
  expect_equal(
    nrow(res$abund),
    nrow(compare$abund)
  )

  expect_equal(
    sort(colnames(res[["abund"]])),
    sort(colnames(compare[["abund"]]))
  )

  expect_equal(
    sort(rownames(res[["abund"]])),
    sort(rownames(compare[["abund"]]))
  )

  expect_equal(
    sort(rownames(res[["tax"]])),
    sort(rownames(compare[["tax"]]))
  )

  expect_equal(
    sort(res$metadata[[1]]),
    sort(compare$metadata[[1]])
  )
})

test_that("Merging two ampvis2-class objects with contradicting taxonomy errors", {
  suppressMessages({
    obj1 <- amp_filter_samples(
      AalborgWWTPs,
      Year %in% "2014"
    )
    obj2 <- amp_filter_samples(
      AalborgWWTPs,
      Year %in% "2015"
    )
  })

  obj2[["tax"]]["OTU_1", "Species"] <- "thisisnotthesame as in obj1"

  res <- expect_warning(
    amp_merge_ampvis2(
      obj1,
      obj2
    ),
    regexp = ".*conflicting taxonomy.*"
  )
})

test_that("Merging three ampvis2-class objects with some non-unique samples errors", {
  suppressMessages({
    MiDAS_2011 <- amp_filter_samples(
      MiDAS,
      Year %in% "2011"
    )
    MiDAS_2011_2012 <- amp_filter_samples(
      MiDAS,
      Year %in% c("2011", "2012")
    )
    MiDAS_2013 <- amp_filter_samples(
      MiDAS,
      Year %in% "2013"
    )
  })

  expect_error(
    amp_merge_ampvis2(
      MiDAS_2011,
      MiDAS_2011_2012,
      MiDAS_2013
    ),
    regexp = "^One or more samples occurs more than once between the objects.*$",
    fixed = FALSE
  )
})

test_that("Merging two objects where one is normalised and the other is not errors", {
  suppressMessages({
    obj1 <- amp_filter_samples(
      AalborgWWTPs,
      Year %in% "2014",
      normalise = FALSE
    )
    obj2 <- amp_filter_samples(
      AalborgWWTPs,
      Year %in% "2015",
      normalise = TRUE
    )
  })

  expect_error(
    amp_merge_ampvis2(
      obj1,
      obj2
    ),
    regexp = "^All objects must be either normalised or not, not mixed$",
    fixed = FALSE
  )
})

test_that("Merging non-ampvis2 class objects errors", {
  expect_error(
    amp_merge_ampvis2(
      AalborgWWTPs,
      iris
    ),
    regexp = "^One or more objects is not an ampvis2-class object$",
    fixed = FALSE
  )
})
