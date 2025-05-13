#' Create a consistent theme for density plots
#'
#' This function provides a minimalist ggplot2 theme optimized for density plots.
#' It removes background elements, grid lines, and y-axis text while keeping
#' essential elements like the axis lines and titles.
#'
#' @return A ggplot2 theme object that can be added to any ggplot
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(data, aes(x = value)) +
#'   geom_density() +
#'   theme_density()
#' }
#'
#' @internal
theme_density <- function() {
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16)
    )
}
