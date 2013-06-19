############################################################################
# Plotting.
############################################################################

#' A plain ggplot2 theme.
#'
#' @param base_size base font size
#' @param base_family base font family
ThemePlain <- function(base_size=15, base_family='') {
    # Starts with theme_bw and then modify some parts
    theme_grey(base_size=base_size, base_family=base_family) %+replace%
            theme(
                    panel.background=element_blank(),
                    panel.grid.major=element_line(colour='grey90', size=rel(0.2)),
                    panel.grid.minor=element_line(colour='grey95', size=rel(0.5)),
                    axis.text.x=element_text(lineheight=0.9, vjust=0.5, angle=-90,
                            hjust=0),
                    axis.ticks=element_line(size=rel(0.2), colour='grey90'),
                    legend.key=element_blank(),
                    legend.position='bottom',
                    legend.key.height=grid::unit(0, "inches"),
                    legend.margin=grid::unit(0, "inches"),
                    plot.margin=grid::unit(c(0.01,0.01,0,0), "inches")
            )
}
