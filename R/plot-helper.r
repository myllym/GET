############################################################################
# Plotting.
############################################################################

# Inner function for setting the position of the legend for 1d envelope plots
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
set_envelope_legend_position <- function() {
  theme(legend.position = 'bottom', legend.key.height = grid::unit(0, "inches"))
}
