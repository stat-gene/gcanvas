test_that("progress notes are messages on stderr", {
  expect_message(.gcanvas_note("gcanvas::test", "hello"), "^\\[gcanvas::test\\] hello")
  expect_message(.gcanvas_warn_msg("careful"), "^\\[WARN\\] careful")
  expect_message(.gcanvas_msg("a", 1, "b"), "^a1b")
  expect_output(suppressMessages(.gcanvas_note("gcanvas::test", "hello")), NA)
  expect_silent(.gcanvas_note("gcanvas::test", "hello", silent = TRUE))
})
