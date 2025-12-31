# ---------------------------
# PORTFOLIO / ANONYMIZED TEMPLATE
# ---------------------------
# Portfolio-safe version of a real workflow.
# Client-specific identifiers (paths, sample IDs, locations, coordinates) were removed or generalized.
# Data files are NOT included. Configure paths via environment variables or edit the CONFIG section.
#
# Recommended usage:
#   PROJECT_DIR=/path/to/project Rscript <this_script>.R
# ---------------------------

project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
setwd(project_dir)

# Fig. 1. Spatial distribution of sampling sites (portfolio template).
# Reviewer fix: "All figures and fonts are dull. Please improve all figures with good quality sharp font (do not make it bold)."
# Strategy: increase contrast, slightly thicken key lines, increase font sizes, prefer vector export (PDF) + keep 600 dpi PNG.

library(sf)
library(dplyr)
library(ggplot2)
library(scales)
library(grid)

pdf_fonts <- names(grDevices::pdfFonts())
font_family <- if ("Times New Roman" %in% pdf_fonts) {
  "Times New Roman"
} else if ("Times" %in% pdf_fonts) {
  "Times"
} else {
  "serif"
}
base_text   <- element_text(family = font_family, face = "plain")

# ==== 1. Farm coordinates (with C as BOTH open & closed) ====
# ==== 1. Example site coordinates (REPLACE with your own) ====
# Portfolio version uses placeholder sites (not real client locations).
farms <- tibble::tribble(
  ~site_id, ~region,    ~system,    ~lat,  ~lon,
  "S1",     "Region 1", "System A", 14.50, 120.90,
  "S2",     "Region 2", "System B", 14.60, 121.00,
  "S3",     "Region 3", "System A", 14.70, 121.10,
  "S4",     "Region 4", "System B", 14.80, 121.20,
  "S5",     "Region 5", "System A", 14.90, 121.30
)

farms_sf <- st_as_sf(farms, coords = c("lon", "lat"), crs = 4326)

# Aggregate to 5 region sites, centering markers to cover paired locations
sites_sf <- farms_sf %>%
  group_by(region) %>%
  summarise(geometry = st_centroid(st_union(geometry)), .groups = "drop")

# ==== 2. Locate and read region shapefile ====
# ==== 2. Load admin boundaries shapefile (optional) ====
# Provide an external shapefile via env var ADMIN1_SHP, or place one at data/shapes/admin1.shp
shape_path <- Sys.getenv("ADMIN1_SHP", unset = "data/shapes/admin1.shp")
if (!file.exists(shape_path)) {
  warning("Admin boundary shapefile not found: ", shape_path, " (plot will show points only).")
}
th_prov <- if (file.exists(shape_path)) st_read(shape_path, quiet = TRUE) %>% st_make_valid() else NULL

central_names <- c("Region 1","Region 2","Region 3","Region 4","Region 5")

central_th <- if (!is.null(th_prov)) {
  th_prov %>%
    filter(ADM1_EN %in% central_names) %>%
    mutate(region = ADM1_EN)
} else {
  NULL
}

# ==== 3. Bounding box around farms ====
bbox <- st_bbox(farms_sf)
xpad <- (bbox$xmax - bbox$xmin) * 0.35
ypad <- (bbox$ymax - bbox$ymin) * 0.35
xlim <- c(bbox["xmin"] - xpad, bbox["xmax"] + xpad)
ylim <- c(bbox["ymin"] - ypad, bbox["ymax"] + ypad)

farms_sf$system <- factor(farms_sf$system, levels = c("Open cage", "Closed pond"))
farms_sf$region <- factor(farms_sf$region, levels = central_names)
sites_sf$region <- factor(sites_sf$region, levels = central_names)

site_palette <- c(
  "Region 1" = "#0f1c2d",
  "Region 2" = "#d1495b",
  "Region 3" = "#1f7a8c",
  "Region 4" = "#f0a202",
  "Region 5" = "#2d936c"
)

# Define key rivers to emphasize (case-insensitive match on OSM "name")
 to emphasize (case-insensitive match on OSM "name")
major_rivers <- tolower(c(
  # Add patterns here if you want named major rivers emphasized
))

# Read river network from downloaded HOTOSM shapefile
river_major <- NULL
river_minor <- NULL
water_path <- Sys.getenv("WATERWAYS_SHP", unset = "data/shapes/waterways.shp")

if (file.exists(water_path)) {
  river_sf <- st_read(water_path, quiet = TRUE) %>% st_make_valid()

  bbox_mat <- matrix(
    c(
      xlim[1], ylim[1],
      xlim[2], ylim[1],
      xlim[2], ylim[2],
      xlim[1], ylim[2],
      xlim[1], ylim[1]
    ),
    ncol = 2,
    byrow = TRUE
  )
  bbox_sfc <- st_sfc(st_polygon(list(bbox_mat)), crs = st_crs(river_sf))
  river_sf <- st_intersection(river_sf, bbox_sfc)

  # If no major river patterns provided, treat all as minor
  if (length(major_rivers) == 0) {
    river_minor <- river_sf
  } else {
    river_sf <- river_sf %>%
      mutate(
        name_combined = tolower(paste(
          dplyr::coalesce(name_en, ""),
          dplyr::coalesce(name, ""),
          dplyr::coalesce(name_th, ""),
          sep = " "
        )),
        is_major_name = purrr::reduce(
          lapply(major_rivers, function(pat) grepl(pat, name_combined, fixed = TRUE)),
          `|`,
          .init = FALSE
        )
      )

    river_major <- river_sf %>% filter(is_major_name)
    river_minor <- river_sf %>% filter(!is_major_name)
  }
}

# ==== 4. Optional helpers ====
add_spatial_decor <- function(p) {
  if (requireNamespace("ggspatial", quietly = TRUE)) {
    p +
      ggspatial::annotation_scale(location = "bl", width_hint = 0.35, text_cex = 0.75) +
      ggspatial::annotation_north_arrow(
        location = "bl",
        which_north = "true",
        height = unit(1.2, "cm"),
        width = unit(0.8, "cm"),
        pad_x = unit(0.8, "cm"),
        pad_y = unit(0.6, "cm"),
        style = ggspatial::north_arrow_fancy_orienteering
      )
  } else {
    p
  }
}

# ==== 5. Plot (Figure 1) ====
# Reviewer-facing adjustments:
# - Higher contrast map fill + white panel
# - Slightly thicker borders/lines (not bold fonts)
# - Larger, darker text for print readability
# - Crisp vector export option (PDF) + high-res PNG

base_plot <- ggplot()
if (!is.null(central_th)) {
  base_plot <- base_plot +
    geom_sf(
      data = central_th,
      fill = NA,
      color = "#0f1c2d",
      linewidth = 0.65
    )
}
base_plot <- base_plot +
{ if (!is.null(river_minor) && nrow(river_minor) > 0)
    geom_sf(
      data = river_minor,
      color = "#0B4F9C",
      alpha = 0.70,
      linewidth = 0.35,
      lineend = "round",
      show.legend = "none"
    ) else NULL } +
  { if (!is.null(river_major) && nrow(river_major) > 0)
    geom_sf(
      data = river_major,
      color = "#08306B",
      alpha = 1.00,
      linewidth = 0.50,
      lineend = "round",
      show.legend = "none"
    ) else NULL } +
  geom_sf(
    data = sites_sf,
    aes(fill = region),
    shape = 21,
    size = 6.2,
    stroke = 1.6,
    color = "#0F1C2D"
  ) +
  scale_fill_manual(
    values = site_palette,
    breaks = central_names,
    labels = central_names,
    name = "Region"
  ) +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 4.4, color = "#0F1C2D", stroke = 1.3)
  )) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Sampling sites across a study region",
    subtitle = "Admin boundaries (optional) and waterways (optional)"
  ) +
  theme_void() +
  theme(
    text            = base_text,
    plot.background  = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    panel.border     = element_rect(color = "#1C2C3A", fill = NA, linewidth = 0.8),

    legend.position  = "right",
    legend.box       = "vertical",
    legend.title     = element_text(color = "#0F1C2D", size = 12, family = font_family, face = "plain"),
    legend.text      = element_text(color = "#1F2D3D", size = 11, family = font_family, face = "plain"),

    plot.title       = element_text(color = "#0F1C2D", size = 19, family = font_family, face = "plain"),
    plot.subtitle    = element_text(color = "#2B3B4C", size = 13, margin = margin(b = 6), family = font_family, face = "plain"),

    axis.text        = element_text(color = "#1F2D3D", size = 10, family = font_family, face = "plain"),
    axis.title       = element_text(color = "#0F1C2D", size = 11, family = font_family, face = "plain"),

    plot.margin      = margin(10, 16, 12, 16)
  )

final_plot <- add_spatial_decor(base_plot)

dir.create("figures", showWarnings = FALSE)

# ---- Preferred (vector) export for "sharp font" ----
# If you have cairo installed, this yields crisp text in PDF.
if (capabilities("cairo")) {
  ggsave(
    filename = "figures/Figure1_map.pdf",
    plot = final_plot,
    width = 7.5,
    height = 7.5,
    units = "in",
    device = cairo_pdf
  )
  ggsave(
    filename = "figures/Figure1_map.tiff",
    plot = final_plot,
    width = 7.5,
    height = 7.5,
    dpi = 600,
    units = "in",
    device = "tiff",
    compression = "lzw",
    type = "cairo"
  )
} else {
  # Fallback PDF device
  ggsave(
    filename = "figures/Figure1_map.pdf",
    plot = final_plot,
    width = 7.5,
    height = 7.5,
    units = "in"
  )
}

# ---- High-res raster export (if journal wants PNG/TIFF) ----
ggsave(
  filename = "figures/Figure1_map.png",
  plot = final_plot,
  width = 7.5,
  height = 7.5,
  dpi = 600,
  units = "in",
  type = "cairo",
  bg = "white"
)

final_plot
