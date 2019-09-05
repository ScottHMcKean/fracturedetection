#' Function to process dxf files for analysis
#'
#' @description Basic function to load a .dxf file
#' @param dxf_file filename, direct path preferred
#' @return list of dxf object including sf object, bounding box, and convex hull
#' @export
process_dxf <- function(dxf_file){
  sf_df <- readOGR(dsn = dxf_file) %>%
    st_as_sf(ogr)

  bound_box <- st_bbox(sf_df) %>% round(-1)

  convex_hull <- st_combine(sf_df) %>% st_convex_hull()

  list(sf_df = sf_df, bbox = bound_box, conv_hull = convex_hull)
}

#' Run digital scanlines
#'
#' @description Function to generate digital scanlines, get intersections,
#' and return coordinates for further analysis (statistics)
#'
#' @param dxf_sf output from process_dxf
#' @param scanline_direction should be 'horiz' or 'vert'
#' @param pixel_distance pixel distance to space scanlines
#' @return list of scanlines results including origins, spacing, intersections, and coordinates
#' @export
run_scanlines <- function(dxf_sf, scanline_direction = 'horiz',
                          pixel_distance = 100){

  if (!scanline_direction %in% c('horiz','vert')){
    warning('Scanline direction not specified as horizontal or vertical, setting to horizontal')
    scanline_direction = 'horiz'
  }

  # generate scanlines
  if (scanline_direction == 'horiz'){
    scanline_origins <- seq(dxf_sf$bbox[2], dxf_sf$bbox[4], pixel_distance)

    f_linestring <- function(y){
      st_linestring(matrix(c(dxf_sf$bbox[1],dxf_sf$bbox[3], y, y), ncol = 2), dim = 'XY')
    }

    scanlines <- lapply(scanline_origins,f_linestring) %>%
      st_sfc() %>%
      st_sf()

  } else {
    scanline_origins <- seq(dxf_sf$bbox[1], dxf_sf$bbox[3], pixel_distance)

    f_linestring <- function(x){
      st_linestring(matrix(c(x, x, dxf_sf$bbox[2],dxf_sf$bbox[4]), ncol = 2), dim = 'XY')
    }

    scanlines <- lapply(scanline_origins,f_linestring) %>%
      st_sfc() %>%
      st_sf()
  }

  # get rows of intersecting fractures with scanlines
  f_intersects <- function(x) st_intersects(x, dxf_sf$sf_df$geometry)[[1]]

  scanline_intersections <- map(scanlines$geometry, f_intersects)

  # get coordinates of intersections
  f_coords <- function(i) {
    st_intersection(scanlines[i,'geometry'], dxf_sf$sf_df$geometry) %>%
      st_coordinates()
  }

  scanline_coords <- map(1:nrow(scanlines), f_coords)

  list(distance = pixel_distance,
       direction = scanline_direction,
       origins = scanline_origins,
       intersections = scanline_intersections,
       coords = scanline_coords)
}

#' Get scanlines stats
#'
#' @description Calculates a dataframe with statistics on each scanline,
#' including count, line length, spacing, p10, and p21
#'
#' @param scanline_results output from calc_scanlines
#' @return dataframe of scanline stats
#' @export
scanline_stats <- function(scanline_results){

  n = length(scanline_results$coords)

  if (scanline_results$direction == 'horiz'){
    coords <- map(1:n, function(x) {scanline_results$coords[[x]][,1]})
  } else {
    coords <- map(1:n, function(x) {scanline_results$coords[[x]][,2]})
  }

  raw_frac_count <- sapply(coords, length)
  censored_frac_count <- ifelse(raw_frac_count < 2, 0, raw_frac_count - 2)

  # line length = max - min coords
  line_length <- function(coord){
    if (length(coord) == 0) return(0)
    min = min(coord, na.rm = TRUE)
    max = max(coord, na.rm = TRUE)
    max - min
  }

  lengths <- sapply(coords, line_length) %>%
    ifelse(is.finite(.), ., NA)

  # get P10 for each scanline
  p10 = censored_frac_count / lengths

  # get frac spacing using diff()
  diff_spacing <- function(coord){
    if (length(coord) < 3) return(0)

    sort(coord)[2:(length(coord)-1)] %>%
      diff() %>%
      abs() %>%
      mean()
  }

  diff_spacings <- sapply(coords, diff_spacing)

  calc_p21 <- function(spacing, length){
    if (!is.finite(spacing)) return(NA)
    if (!is.finite(length)) return(NA)
    spacing / length
  }

  p21 = map2_dbl(diff_spacings, lengths, calc_p21)

  data.frame('scanline' = 1:n,
             'coord' = scanline_results$origins,
             raw_frac_count,
             censored_frac_count,
             'scanline_length' = lengths,
             diff_spacings,
             p10,
             p21)
}
