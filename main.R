 install.packages("pacman")
pacman::p_load(
    terra, elevatr, tidyverse, ggnewscale, ggspatial,
    geodata, sf, scales
)

# --- 1) Load your custom state boundary shapefile ---
base_path <- "C:/Users/Rahul Kumar/OneDrive/Documents/shapefile/"
states <- sf::st_read(paste0(base_path, "STATE_BOUNDARY_CLEAN.shp"), quiet = TRUE)
states_vect <- terra::vect(states)

# --- 2) DEM for India (SRTM 30 arc-sec, in WGS84) ---
dem_srtm <- geodata::elevation_30s(country = "IND", path = tempdir())

# --- 3) WorldPop 100m population count 2020 (local file, also WGS84) ---
pop_100m <- terra::rast("C:/Users/Rahul Kumar/Downloads/korea/NAROBI/2015/ind_ppp_2020_constrained.tif")

# --- 4) Reproject your state vector to match DEM CRS ---
# (WGS84, EPSG:4326)
states_vect_wgs84 <- terra::project(states_vect, dem_srtm)

# --- 5) Crop and mask DEM and population raster with states boundary ---
dem_masked <- terra::crop(dem_srtm, states_vect_wgs84) |> terra::mask(states_vect_wgs84)
pop_masked <- terra::crop(pop_100m, states_vect_wgs84) |> terra::mask(states_vect_wgs84)

# --- 6) Exaggerate DEM for hillshade ---
dem_exaggerated <- dem_masked * 1.3

# --- 7) Compute hillshade ---
slope <- terra::terrain(dem_exaggerated, v = "slope", unit = "radians")
aspect <- terra::terrain(dem_exaggerated, v = "aspect", unit = "radians")
hillshade_raw <- terra::shade(
    slope, aspect,
    angle = 40, direction = 225
)

# --- 8) Resample population raster to hillshade grid ---
pop_on_hillshade <- terra::resample(pop_masked, hillshade_raw, method = "bilinear")

# --- 9) Mask hillshade where population exists ---
hillshade_no_pop <- terra::ifel(is.na(pop_on_hillshade), hillshade_raw, NA)

# --- 10) Convert rasters to data frames for ggplot ---
hillshade_df <- terra::as.data.frame(hillshade_no_pop, xy = TRUE, na.rm = TRUE)
colnames(hillshade_df)[3] <- "hillshade"

pop_df <- terra::as.data.frame(pop_on_hillshade, xy = TRUE, na.rm = TRUE)
names(pop_df)[3] <- "ind_ppp_2020_constrained"
pop_df$ind_ppp_2020_constrained[pop_df$ind_ppp_2020_constrained <= 0.1] <- NA

# --- 11) Plot with ggplot2 ---
brks <- c(1, 10, 100, 1e3)

p <- ggplot() +
    # Hillshade background
    geom_raster(data = hillshade_df, aes(
        x, y, fill = hillshade
    )) +
    scale_fill_gradient(
        low = "grey70", high = "grey10",
        guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    # Population layer
    geom_raster(data = pop_df, aes(
        x, y, fill = ind_ppp_2020_constrained
    )) +
    scale_fill_viridis_c(
        name = "Population",
        option = "plasma",
        alpha = 1, begin = .2, end = 1,
        trans = "log10", breaks = brks,
        labels = scales::comma,
        guide = guide_colourbar(
            title.position = "top",
            barheight = unit(30, "mm"),
            barwidth = unit(2, "mm"),
            ticks.color = "grey10",
            frame.colour = "grey10"
        )
    ) +
    # State boundaries overlay
    geom_sf(
        data = sf::st_transform(states, 4326), fill = NA,
        color = "black", linewidth = .25
    ) +
    # Cartographic extras
    ggspatial::annotation_north_arrow(
        location = "tl", which_north = "true",
        height = unit(10, "mm"),
        width = unit(10, "mm"),
        style = north_arrow_orienteering
    ) +
    annotation_scale(
        location = "br", pad_y = unit(2, "mm"),
        height = unit(2, "mm")
    ) +
    coord_sf(expand = FALSE) +
    labs(
        title = "India · Population (2020)",
        subtitle = "WorldPop 100m constrained grid, cropped to custom state boundaries",
        caption = "Data: WorldPop · SRTM via geodata | Design: Rahul4planets Maps"
    ) +
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = .02),
        plot.subtitle = element_text(size = 14, hjust = .02),
        plot.caption = element_text(hjust = .5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.margin = margin(t = 0, r = 5, b = 0, l = 3),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    theme_void()

# --- 12) Save the figure ---
ggsave(
    "india_population_relief_custom_states.png",
    plot = p,
    width = 8, height = 5, dpi = 600,
    bg = "white"
)
