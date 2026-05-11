# ==========================================
# Global asthma burden maps (GBD 2019)
# Outputs required by Project.tex:
#   - map_asmr.pdf (facets: 1990, 2019)
#   - map_asdr.pdf (facets: 1990, 2019)
# ==========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(readr)
  library(stringr)
  library(rnaturalearth)
  library(cowplot)
  library(ggpubr)
  library(scales)
})

# Use the script directory as the working root to keep outputs in the current folder.
script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Read only heatmap.csv from the current folder.
input_csv <- file.path(script_dir, "heatmap.csv")
 if (!file.exists(input_csv)) {
  stop("Input CSV not found. Place heatmap.csv in the current folder.")
}

# Read and standardize core fields used for filtering/joining.
raw <- read_csv(input_csv, show_col_types = FALSE) %>%
  mutate(
    location = str_squish(as.character(location)),
    metric = str_squish(as.character(metric)),
    measure = str_squish(as.character(measure)),
    age = str_squish(as.character(age)),
    year = as.integer(year),
    val = as.numeric(val)
  )

# Keep age-standardized rates for 1990 and 2019 only.
core <- raw %>%
  filter(year %in% c(1990, 2019), metric == "Rate", age == "Age-standardized") %>%
  select(location, measure, year, val)

# World geometry from Natural Earth plus deterministic name harmonization for GBD labels.
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  mutate(
    location_match = case_when(
      name == "Vietnam" ~ "Viet Nam",
      name == "Venezuela" ~ "Venezuela (Bolivarian Republic of)",
      name == "Micronesia" ~ "Micronesia (Federated States of)",
      name == "Marshall Is." ~ "Marshall Islands",
      name == "N. Mariana Is." ~ "Northern Mariana Islands",
      name == "U.S. Virgin Is." ~ "United States Virgin Islands",
      name == "Tanzania" ~ "United Republic of Tanzania",
      name == "Taiwan" ~ "Taiwan (Province of China)",
      name == "Syria" ~ "Syrian Arab Republic",
      name == "eSwatini" ~ "Eswatini",
      name == "S. Sudan" ~ "South Sudan",
      name == "South Korea" ~ "Republic of Korea",
      name == "Solomon Is." ~ "Solomon Islands",
      str_detect(name, "Tom") & str_detect(name, "Principe") ~ "Sao Tome and Principe",
      name == "St. Vin. and Gren." ~ "Saint Vincent and the Grenadines",
      name == "St. Kitts and Nevis" ~ "Saint Kitts and Nevis",
      name == "Russia" ~ "Russian Federation",
      name == "North Korea" ~ "Democratic People's Republic of Korea",
      name == "Cook Is." ~ "Cook Islands",
      name == "Moldova" ~ "Republic of Moldova",
      name == "Laos" ~ "Lao People's Democratic Republic",
      name == "Iran" ~ "Iran (Islamic Republic of)",
      name == "Eq. Guinea" ~ "Equatorial Guinea",
      name == "Dominican Rep." ~ "Dominican Republic",
      name == "Dem. Rep. Congo" ~ "Democratic Republic of the Congo",
      name == "Central African Rep." ~ "Central African Republic",
      name == "Brunei" ~ "Brunei Darussalam",
      name == "Bosnia and Herz." ~ "Bosnia and Herzegovina",
      name == "Bolivia" ~ "Bolivia (Plurinational State of)",
      name == "Antigua and Barb." ~ "Antigua and Barbuda",
      name == "Ivory Coast" ~ "Cote d'Ivoire",
      TRUE ~ name
    )
  )

# Keep an explicit data quality check to surface unmatched country names during cleaning.
report_unmatched <- function(df_metric, metric_tag) {
  unmatched <- df_metric %>%
    distinct(location) %>%
    anti_join(world %>% st_drop_geometry() %>% distinct(location_match), by = c("location" = "location_match")) %>%
    arrange(location)

  if (nrow(unmatched) > 0) {
    message("[", metric_tag, "] Unmatched locations (not drawn): ", paste(unmatched$location, collapse = "; "))
  }
}

# Build one two-panel map (1990 and 2019) for a target metric.
make_metric_map <- function(df_metric, metric_label, out_pdf) {
  report_unmatched(df_metric, metric_label)

  # Build one year-specific panel and insert the Top 5 table into the map area.
  make_year_panel <- function(target_year) {
    year_df <- df_metric %>% filter(year == target_year)

    # Left join preserves world polygons; unmatched countries stay NA and are rendered in grey.
    map_sf <- world %>%
      left_join(year_df, by = c("location_match" = "location")) %>%
      mutate(year = target_year) %>%
      st_transform(crs = "+proj=robin")

    top5 <- year_df %>%
      filter(!is.na(val)) %>%
      arrange(desc(val)) %>%
      slice_head(n = 5) %>%
      transmute(Country = location, Value = round(val, 2))

    table_plot <- ggtexttable(
      top5,
      rows = NULL,
      theme = ttheme(
        "minimal",
        base_size = 8,
        colnames.style = colnames_style(face = "bold"),
        tbody.style = tbody_style(fill = c("#FFFFFF", "#F5F5F5"))
      )
    )

    base_map <- ggplot(map_sf) +
      geom_sf(aes(fill = val), color = "grey20", size = 0.08) +
      scale_fill_gradientn(
        colours = c("#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#fee090", "#fdae61", "#f46d43", "#d73027"),
        trans = "log10",
        na.value = "#eeeeee",
        labels = label_number(big.mark = ","),
        name = metric_label
      ) +
      theme_void(base_size = 11) +
      theme(
        legend.position = "bottom",
        legend.key.width = unit(1.6, "cm"),
        plot.margin = margin(0, 0, 0, 0)
      )

    # Reserve a left strip for the table to avoid occluding the world map.
    ggdraw() +
      draw_plot(base_map, x = 0.18, y = 0.00, width = 0.82, height = 1.00) +
      draw_plot(table_plot, x = 0.02, y = 0.57, width = 0.15, height = 0.34)
  }

  panel_1990 <- make_year_panel(1990)
  panel_2019 <- make_year_panel(2019)

  p <- plot_grid(panel_1990, panel_2019, ncol = 1, align = "v", axis = "lr")

  ggsave(
    filename = file.path(script_dir, out_pdf),
    plot = p,
    width = 14,
    height = 12,
    device = "pdf",
    useDingbats = FALSE
  )
}

# ASMR: deaths rate.
asmr_df <- core %>%
  filter(str_detect(str_to_lower(measure), "death")) %>%
  select(location, year, val)

# ASDR: DALY rate.
asdr_df <- core %>%
  filter(str_detect(str_to_lower(measure), "daly")) %>%
  select(location, year, val)

if (nrow(asmr_df) == 0) {
  stop("No ASMR data found. Check the 'measure' values for deaths rate rows.")
}

if (nrow(asdr_df) == 0) {
  stop("No ASDR data found. Check the 'measure' values for DALY rate rows.")
}

make_metric_map(asmr_df, "ASMR", "map_asmr.pdf")
make_metric_map(asdr_df, "ASDR", "map_asdr.pdf")

message("Done. Saved: map_asmr.pdf and map_asdr.pdf")