test_that("mixed_by does not mutate caller-owned data or external_df", {
  data_df <- data.frame(
    row_id = seq_len(8L),
    id = rep(c(2L, 1L, 4L, 3L), each = 2L),
    split = rep(c("b", "a"), 4L),
    y = seq_len(8L)
  )
  external_dt <- data.table::data.table(
    id = c(4L, 2L, 3L, 1L),
    external_value = c("d", "b", "c", "a")
  )
  data_before <- data_df
  external_before <- data.table::copy(external_dt)

  suppressMessages(mixed_by(
    data_df,
    outcomes = "y",
    rhs_model_formulae = list(m1 = ~ external_value),
    split_on = "split",
    external_df = external_dt,
    external_merge_by = "id",
    calculate = character(0),
    padjust_by = NULL
  ))

  expect_identical(data_df, data_before)
  expect_false(data.table::is.data.table(data_df))
  expect_identical(external_dt, external_before)
  expect_null(data.table::key(external_dt))
})

test_that("mixed_by does not sort or key a caller-owned data.table", {
  data_dt <- data.table::data.table(
    row_id = seq_len(8L),
    id = rep(c(2L, 1L, 4L, 3L), each = 2L),
    split = rep(c("b", "a"), 4L),
    y = seq_len(8L)
  )
  data_before <- data.table::copy(data_dt)

  suppressMessages(mixed_by(
    data_dt,
    outcomes = "y",
    rhs_model_formulae = list(m1 = ~ 1),
    split_on = "split",
    calculate = character(0),
    padjust_by = NULL
  ))

  expect_identical(data_dt, data_before)
  expect_null(data.table::key(data_dt))
})
