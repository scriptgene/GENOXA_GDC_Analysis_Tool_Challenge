#' Theme for single plots
#' @noRd
extra_theme <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = element_rect(
        fill = NA,
        colour = "grey40",
        size = 0.5
      ),
      legend.title = element_blank(),
      axis.ticks = element_line(colour = "grey40"),
      axis.line = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    )
}

# Theme for faceted plots
extra_facet_theme <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = element_rect(
        fill = NA,
        colour = "grey40",
        size = 0.5
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_line(colour = "grey40"),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        colour = "grey20",
        hjust = 0.01,
        face = "bold",
        vjust = 0
      )
    )
}

# Pick viridis colors
set_col <- function(palette,
                    n,
                    direction = 1) {
  # Assign palette
  if (palette == "inferno") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.6,
      end = 1,
      direction = direction,
      option = "B"
    )(n)
  } else if (palette == "magma") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.4,
      end = 1,
      direction = direction,
      option = "A"
    )(n)
  } else if (palette == "plasma") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.55,
      end = 1,
      direction = direction,
      option = "C"
    )(n)
  } else if (palette == "cividis") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.2,
      end = 1,
      direction = direction,
      option = "E"
    )(n)
  } else if (palette == "rocket") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.3,
      end = 1,
      direction = direction,
      option = "F"
    )(n)
  } else if (palette == "mako") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.5,
      end = 1,
      direction = direction,
      option = "G"
    )(n)
  } else if (palette == "turbo") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.3,
      end = 1,
      direction = direction,
      option = "H"
    )(n)
  } else {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.5,
      end = 1,
      direction = direction,
      option = "D"
    )(n)
  }
  return(pal_col)
}



# Volcano Plot
# Load necessary library
library(grDevices)  # For colorRampPalette

set_col2 <- function(palette,
                    n,
                    direction = 1) {
  # Assign palette
  if (palette == "dark_brown") {
    pal_col <- colorRampPalette(c("#3B2C1C"))(n)
  } else if (palette == "green") {
    pal_col <- colorRampPalette(c("#2F4F2F"))(n)
  } else if (palette == "red") {
    pal_col <- colorRampPalette(c("#4B2C2F"))(n)
  } else if (palette == "gray") {
    pal_col <- colorRampPalette(c("#4C4C4C"))(n)
  } else if (palette == "teal") {
    pal_col <- colorRampPalette(c("#004d4d"))(n)
  } else if (palette == "purple") {
    pal_col <- colorRampPalette(c("#3d0a3b"))(n)
  } else if (palette == "orange") {
    pal_col <- colorRampPalette(c("#4b2c1c"))(n)
  } else if (palette == "blue") {
    pal_col <- colorRampPalette(c("#27408B"))(n)
  } else if (palette == "magenta") {
    pal_col <- colorRampPalette(c("#4b003d"))(n)
  } else {
    stop("Unknown palette name. Please choose from: dark_brown, dark_green, dark_red, dark_gray, dark_teal, dark_purple, dark_orange, dark_blue, dark_magenta.")
  }
  return(pal_col)
}


# Heatmap
# Function to generate color palette using hcl.colors
set_col3 <- function(palette, n) {
  # List of valid palettes
  valid_palettes <- c("BluYl", "RdYlBu", "Greens", "Oranges", "Reds", "YlOrBr", "YlOrRd", "BuPu", "GnBu", "OrRd")
  
  # Check if the palette is valid
  if (!(palette %in% valid_palettes)) {
    stop("Invalid palette name. Choose from: ", paste(valid_palettes, collapse = ", "))
  }
  
  # Generate and return the color palette
  hcl.colors(n, palette)
}
