#' @importFrom Matrix Matrix rowSums colSums readMM writeMM rowMeans
#' @importFrom dplyr mutate inner_join filter left_join select rename distinct
#' @importFrom dplyr first summarise if_else count arrange desc slice any_of
#' @importFrom tidyr pivot_longer contains
#' @importFrom tibble as_tibble
#' @importFrom rlang enquos get_expr
#' @importFrom scales percent label_comma label_number cut_short_scale
#' @importFrom logger log_info log_success log_warn
#' @importFrom glue glue
#' @importFrom stringr str_remove str_to_title
#' @importFrom methods is new setClass setGeneric setMethod
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_density scale_y_continuous scale_x_log10
#' @importFrom ggplot2 geom_segment scale_x_continuous theme_void theme
#' @importFrom ggplot2 element_blank element_line element_text margin labs
#' @importFrom ggplot2 geom_hline geom_point facet_grid vars expansion
#' @importFrom ggplot2 coord_cartesian theme_classic theme_bw
#' @importFrom ggrepel geom_text_repel
#' @importFrom readr read_tsv read_csv write_tsv write_csv cols col_character col_integer
#' @importFrom plyranges as_granges join_overlap_left
#' @importFrom fs path
#' @importFrom furrr future_map2_dbl
#' @importFrom R.utils gzip
#' @importFrom grid unit
NULL
