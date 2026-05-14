library(sf)
library(ggplot2)
library(dplyr)
library(terra)
library(cowplot)
library(geosphere)
library(scales)
library(maptiles)
library(tidyterra)

set.seed(1234)

out_dir <- "C:/Users/HP/Documents/GEO/eda_output"

county_sf <- st_read("C:/Users/HP/Documents/GEO/data/County.shp", quiet = TRUE)
county_sf$county <- county_sf$COUNTY
county_sf$county[county_sf$COUNTY == "Murang'a"] <- "Muranga"
county_sf$county[county_sf$COUNTY == "Keiyo-Marakwet"] <- "Elgeyo Marakwet"
county_sf$county[county_sf$COUNTY == "Tharaka"] <- "Tharaka Nithi"

pop_rast <- rast("C:/Users/HP/Documents/GEO/covariates/KEN_pd_2020_1km.tif")

vil_sf <- st_read("C:/Users/HP/Documents/GNSAMPLING/kenya_villages.shp", quiet = TRUE)
vil_df <- st_drop_geometry(vil_sf)
colnames(vil_df) <- tolower(colnames(vil_df))
vil_df <- vil_df[, c("name", "longdd", "latdd")]
colnames(vil_df)[2:3] <- c("lon", "lat")
vil_df <- vil_df[!duplicated(cbind(vil_df$lon, vil_df$lat)), ]
vil_df <- vil_df[complete.cases(vil_df), ]
vil_sf <- st_as_sf(vil_df, coords = c("lon", "lat"), crs = 4326)

wards <- st_read("C:/Users/HP/Documents/GNSAMPLING/wards.shp", quiet = TRUE)
wards <- st_make_valid(wards)
wards$county_w <- trimws(gsub("[-_]", " ", tolower(wards$county)))

vil_sf <- st_join(vil_sf, wards[, "county_w"], join = st_within)

# assign nearest ward to villages that didn't fall inside any ward polygon
no_ward <- is.na(vil_sf$county_w)
if (any(no_ward)) {
  nn <- st_nearest_feature(vil_sf[no_ward, ], wards)
  vil_sf$county_w[no_ward] <- wards$county_w[nn]
}

vil_km <- vil_sf[grepl("kirinyaga|murang", vil_sf$county_w, ignore.case = TRUE), ]
vil_km$county_label <- ifelse(grepl("murang", vil_km$county_w, ignore.case = TRUE), "Muranga", "Kirinyaga")

# worldpop density as sampling weight, minimum 0.5
pop_vals <- extract(pop_rast, project(vect(vil_km), crs(pop_rast)))[, 2]
pop_vals[is.na(pop_vals)] <- 0
vil_km$weight <- ifelse(pop_vals < 0.5, 0.5, pop_vals)


# draw primary sites one at a time, only keep candidates at least delta away from existing sites
select_primary <- function(coords, n, delta, weights) {
  chosen <- c()
  chosen_xy <- matrix(NA, nrow = n, ncol = 2)

  for (i in 1:n) {
    placed <- FALSE

    for (attempt in 1:1000) {
      pool <- setdiff(1:nrow(coords), chosen)
      w <- weights[pool] / sum(weights[pool])
      cand <- sample(pool, 1, prob = w)
      xy <- coords[cand, ]

      if (length(chosen) == 0) {
        chosen <- c(chosen, cand)
        chosen_xy[i, ] <- xy
        placed <- TRUE
        break
      }

      dists <- distHaversine(xy, chosen_xy[1:(i - 1), , drop = FALSE])
      if (min(dists) >= delta) {
        chosen <- c(chosen, cand)
        chosen_xy[i, ] <- xy
        placed <- TRUE
        break
      }
    }

    # if nothing worked after 1000 tries, take whichever village is furthest away
    if (!placed) {
      pool <- setdiff(1:nrow(coords), chosen)
      min_d <- sapply(pool, function(j)
        min(distHaversine(coords[j, ], chosen_xy[1:(i - 1), , drop = FALSE])))
      best <- pool[which.max(min_d)]
      chosen <- c(chosen, best)
      chosen_xy[i, ] <- coords[best, ]
    }
  }

  list(indices = chosen, coords = chosen_xy)
}

# for k of the primary sites, find a nearby companion village within zeta
add_close_pairs <- function(prim_idx, prim_xy, all_coords, all_idx, zeta, k) {
  hosts <- sort(sample(1:length(prim_idx), k))
  pairs <- c()

  for (j in 1:k) {
    used <- c(prim_idx, pairs)
    pool <- setdiff(all_idx, used)
    if (length(pool) == 0) next

    dists <- distHaversine(prim_xy[hosts[j], ], all_coords[pool, , drop = FALSE])
    in_range <- pool[dists <= zeta]

    if (length(in_range) > 0) {
      pairs <- c(pairs, sample(in_range, 1))
    } else {
      pairs <- c(pairs, pool[which.min(dists)])
    }
  }
  pairs
}

# count how many of the n primary slots fall back to forced placement at a given delta
count_fails <- function(coords, n, delta, weights) {
  fails <- 0
  chosen <- c()
  chosen_xy <- matrix(NA, nrow = n, ncol = 2)

  for (i in 1:n) {
    ok <- FALSE
    for (attempt in 1:500) {
      pool <- setdiff(1:nrow(coords), chosen)
      w <- weights[pool] / sum(weights[pool])
      cand <- sample(pool, 1, prob = w)
      xy <- coords[cand, ]

      if (length(chosen) == 0 ||
          min(distHaversine(xy, chosen_xy[1:length(chosen), , drop = FALSE])) >= delta) {
        chosen <- c(chosen, cand)
        chosen_xy[length(chosen), ] <- xy
        ok <- TRUE
        break
      }
    }
    if (!ok) fails <- fails + 1
  }
  fails
}


run_icp <- function(county_name) {
  df <- vil_km[vil_km$county_label == county_name, ]
  coords <- st_coordinates(df)
  weights <- df$weight
  idx <- 1:nrow(df)

  # binary search over 8 iterations to find largest delta with no forced placements
  lo <- 500
  hi <- 20000
  for (i in 1:8) {
    mid <- (lo + hi) / 2
    if (count_fails(coords, n_primary, mid, weights) == 0) lo <- mid else hi <- mid
  }
  delta <- floor(lo / 1000) * 1000
  zeta <- delta / 3

  prim <- select_primary(coords, n_primary, delta, weights)
  cp <- add_close_pairs(prim$indices, prim$coords, coords, idx, zeta, k)

  out <- df[c(prim$indices, cp), ]
  out$type <- c(rep("Primary", length(prim$indices)), rep("Close pair", length(cp)))
  out
}

n_primary <- 20
k <- 5

kir <- run_icp("Kirinyaga")
mur <- run_icp("Muranga")


pt_col <- c("Primary" = "forestgreen", "Close pair" = "red3")
cty_col <- c("Kirinyaga" = "darkorange", "Muranga" = "steelblue")

make_inset <- function(county_name) {
  cty <- county_sf[county_sf$county == county_name, ]
  ggplot() +
    geom_sf(data = county_sf, fill = "grey88", colour = "white", linewidth = 0.12) +
    geom_sf(data = cty, fill = cty_col[county_name], colour = "white", linewidth = 0.3) +
    coord_sf(expand = FALSE) +
    theme_void() +
    theme(
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.6),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

make_pop_map <- function(county_name) {
  cty <- county_sf[county_sf$county == county_name, ]
  pts <- if (county_name == "Kirinyaga") kir else mur
  bbox <- st_bbox(cty) + c(-0.04, -0.04, 0.04, 0.04)

  pop <- crop(pop_rast, project(vect(cty), crs(pop_rast)), mask = TRUE)
  pdat <- as.data.frame(pop, xy = TRUE)
  colnames(pdat)[3] <- "pop"
  pdat <- pdat[!is.na(pdat$pop), ]

  p <- ggplot() +
    geom_raster(data = pdat, aes(x = x, y = y, fill = pop)) +
    scale_fill_distiller(
      palette = "Purples", direction = 1,
      name = "Persons/km²", labels = scales::comma,
      limits = c(0, quantile(pdat$pop, 0.99, na.rm = TRUE)),
      oob = scales::squish
    ) +
    geom_sf(data = cty, fill = NA, colour = "grey40", linewidth = 0.7) +
    geom_sf(data = pts, colour = "white", size = 3.2, alpha = 0.5) +
    geom_sf(data = pts, aes(colour = type), size = 2.0) +
    scale_colour_manual(values = pt_col, name = "ICP cluster") +
    guides(
      fill = guide_colorbar(order = 1, barheight = unit(2.5, "cm")),
      colour = guide_legend(order = 2, override.aes = list(size = 3))
    ) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9, face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )

  ggdraw(p) + draw_plot(make_inset(county_name), x = 0.01, y = 0.70, width = 0.18, height = 0.26)
}

make_osm_map <- function(county_name) {
  cty <- county_sf[county_sf$county == county_name, ]
  pts <- if (county_name == "Kirinyaga") kir else mur
  cty_buf <- st_transform(st_buffer(st_transform(cty, 3857), 2000), 4326)
  tiles <- get_tiles(cty_buf, provider = "OpenStreetMap", zoom = 12, crop = TRUE)
  cty_3857 <- st_transform(cty, 3857)
  pts_3857 <- st_transform(pts, 3857)
  bbox <- st_bbox(cty_3857)

  p <- ggplot() +
    geom_spatraster_rgb(data = tiles) +
    geom_sf(data = cty_3857, fill = NA, colour = "black", linewidth = 0.9) +
    geom_sf(data = pts_3857, colour = "white", size = 3.2, alpha = 0.5) +
    geom_sf(data = pts_3857, aes(colour = type), size = 2.0) +
    scale_colour_manual(values = pt_col, name = "ICP cluster") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    coord_sf(
      crs = 3857,
      xlim = c(bbox["xmin"] - 1500, bbox["xmax"] + 1500),
      ylim = c(bbox["ymin"] - 1500, bbox["ymax"] + 1500),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9, face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )

  ggdraw(p) + draw_plot(make_inset(county_name), x = 0.01, y = 0.70, width = 0.18, height = 0.26)
}
