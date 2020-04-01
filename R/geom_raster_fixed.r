# A modified ggplot2::geom_raster to allow setting pixel size using the arguments width and height
#' @importFrom ggplot2 layer
geom_raster_fixed <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        hjust = 0.5,
                        vjust = 0.5,
                        width,
                        height,
                        interpolate = FALSE,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE)
{
  if (!(is.numeric(hjust) && length(hjust) == 1)) stop("`hjust` must be a numeric scalar")
  if (!(is.numeric(vjust) && length(vjust) == 1)) stop("`vjust` must be a numeric scalar")
  if (!(is.numeric(width) && length(width) == 1)) stop("`width` must be a numeric scalar")
  if (!(is.numeric(height) && length(height) == 1)) stop("`height` must be a numeric scalar")

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRasterFixed,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      hjust = hjust,
      vjust = vjust,
      width = width,
      height = height,
      interpolate = interpolate,
      na.rm = na.rm,
      ...
    )
  )
}

# This function is similar to GeomRaster, but uses the given w and h
#' @importFrom ggplot2 ggproto draw_key_rect Geom alpha
#' @importFrom grid rasterGrob
GeomRasterFixed <- ggproto("GeomRasterFixed", Geom,
  default_aes = aes(fill = "grey20", alpha = NA),
  non_missing_aes = c("fill", "xmin", "xmax", "ymin", "ymax"),
  required_aes = c("x", "y"),

  setup_data = function(data, params) {
    hjust <- params$hjust
    vjust <- params$vjust

    w <- params$width
    h <- params$height
    data$xmin <- data$x - w * (1 - hjust)
    data$xmax <- data$x + w * hjust
    data$ymin <- data$y - h * (1 - vjust)
    data$ymax <- data$y + h * vjust
    data
  },

  draw_panel = function(data, panel_params, coord, interpolate = FALSE,
                        hjust = 0.5, vjust = 0.5, height, width) {
    if (!inherits(coord, "CoordCartesian")) {
      stop("geom_raster only works with Cartesian coordinates")
    }
    # Use round instead of as.integer because round off errors
    x_pos <- round((data$x - min(data$x)) / width)
    y_pos <- round((data$y - min(data$y)) / height)
    data <- coord$transform(data, panel_params)

    # Convert vector of data to raster

    nrow <- max(y_pos) + 1
    ncol <- max(x_pos) + 1

    raster <- matrix(NA_character_, nrow = nrow, ncol = ncol)
    raster[cbind(nrow - y_pos, x_pos + 1)] <- alpha(data$fill, data$alpha)

    # Figure out dimensions of raster on plot
    x_rng <- c(min(data$xmin, na.rm = TRUE), max(data$xmax, na.rm = TRUE))
    y_rng <- c(min(data$ymin, na.rm = TRUE), max(data$ymax, na.rm = TRUE))

    rasterGrob(raster,
               x = mean(x_rng), y = mean(y_rng),
               width = diff(x_rng), height = diff(y_rng),
               default.units = "native", interpolate = interpolate
    )
  },
  draw_key = draw_key_rect
)
