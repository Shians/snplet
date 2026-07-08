# ==============================================================================
# Test Suite: Plot Themes
# Description: Tests for the shared ggplot2 theme helpers
# ==============================================================================

library(testthat)
library(ggplot2)

# ==============================================================================

test_that("theme_density() returns a ggplot theme object", {
    result <- theme_density()

    # Verify theme_density returns a ggplot theme
    expect_s3_class(result, "theme")
})

test_that("theme_density() can be added to a ggplot without error", {
    p <- ggplot(data.frame(x = 1:3), aes(x = x)) +
        geom_density() +
        theme_density()

    # Verify the theme applies without error when building the plot
    expect_no_error(ggplot_build(p))
})
