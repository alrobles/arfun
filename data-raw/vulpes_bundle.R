#!/usr/bin/env Rscript
# Rebuild inst/extdata/vulpes/ from the real data services.
#
#   tree        : frozen product of the paper pipeline (VertLife clade cut)
#   occurrences : GBIF parquet index, cleaned with the paper pipeline (KU HPC)
#   climate     : mean of yearly GeoTIFFs from https://bioclim.ecoseek.org
#                 (ERA5-Land BIO01-BIO19, 0.1 deg), cropped + aggregated
#   M polygons  : oneearthr::oe_layer("ecoregions") -- RESOLVE Ecoregions 2017
#                 containing each species' occurrences, dissolved per species
#
# Sections 1-2 need KU HPC paths and are skipped when the products already
# exist in OUT_DIR; sections 3-4 run anywhere with terra + oneearthr.
#
#   OUT_DIR=inst/extdata/vulpes Rscript data-raw/vulpes_bundle.R

suppressPackageStartupMessages({
  library(terra); library(ape); library(oneearthr)
})

OUT_DIR <- Sys.getenv("OUT_DIR", "inst/extdata/vulpes")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CLIMATE_YEARS <- 1991:2020           # standard normals period
VARS          <- c("bio01", "bio12")
BIOCLIM_URL   <- "https://bioclim.ecoseek.org/api/download"

## ---- 1. Tree -------------------------------------------------------------
nwk <- file.path(OUT_DIR, "vulpes_tree.nwk")
if (file.exists(nwk)) {
  cat("tree: keeping existing", nwk, "\n")
} else {
  prep <- readRDS("/beegfs/a474r867/arfun-paper/use_case/Vulpes/prepared_data.rds")
  write.tree(prep$tree, nwk)
}
sp <- gsub("_", " ", read.tree(nwk)$tip.label)

## ---- 2. Occurrences (GBIF index, paper cleaning) ---------------------------
csv <- file.path(OUT_DIR, "vulpes_occurrences.csv")
if (file.exists(csv)) {
  cat("occurrences: keeping existing", csv, "\n")
} else {
  library(arrow); library(dplyr)
  occ <- open_dataset("/home/a474r867/work/gbifdata/indexed/V") |>
    filter(species %in% sp) |>
    select(species, decimallongitude, decimallatitude) |>
    collect() |>
    filter(!is.na(decimallongitude), !is.na(decimallatitude),
           between(decimallongitude, -180, 180),
           between(decimallatitude, -90, 90))
  names(occ) <- c("species", "lon", "lat")
  occ <- occ |>
    group_by(species) |>
    group_modify(~ {
      pts <- .x
      if (nrow(pts) >= 4) {
        keep <- rep(TRUE, nrow(pts))
        for (cl in c("lon", "lat")) {
          q <- quantile(pts[[cl]], c(0.25, 0.75))
          f <- 1.5 * diff(q)
          keep <- keep & pts[[cl]] >= q[1] - f & pts[[cl]] <= q[2] + f
        }
        pts <- pts[keep, ]
      }
      cell <- paste(round(pts$lon / 0.1), round(pts$lat / 0.1))
      pts[!duplicated(cell), ]
    }) |> ungroup() |> as.data.frame()
  write.csv(occ, csv, row.names = FALSE)
}
occ <- read.csv(csv)
occ$species <- as.character(occ$species)

## ---- 3. M polygons: ecoregions containing each species' points ------------
pts <- vect(occ, geom = c("lon", "lat"), crs = "EPSG:4326")
eco <- oe_layer("ecoregions")                    # packaged geometries, offline
cat("ecoregions:", nrow(eco), "layer\n")

hit <- relate(eco, pts, "intersects")            # [n_eco x n_pts]
m_list <- lapply(sp, function(s) {
  rows <- which(rowSums(hit[, occ$species == s, drop = FALSE]) > 0)
  m <- aggregate(eco[rows])                      # dissolve -> one polygon
  m$species <- s
  m
})
msv <- do.call(rbind, m_list)
msv <- simplifyGeom(msv, tolerance = 0.05, preserveTopology = TRUE)
writeVector(msv, file.path(OUT_DIR, "vulpes_m_polygons.gpkg"),
            layer = "m", overwrite = TRUE)
cat("M polygons:", nrow(msv), "\n")

## ---- 4. Climate: bioclim.ecoseek.org yearly means -------------------------
cache <- file.path(tempdir(), "bioclim_ecoseek")
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
fetch <- function(var, year) {
  f <- file.path(cache, sprintf("%s_%d.tif", var, year))
  if (!file.exists(f))
    download.file(sprintf("%s/%d/%s_%d.tif", BIOCLIM_URL, year, var, year),
                  f, mode = "wb", quiet = TRUE)
  f
}
clim <- rast(lapply(VARS, function(v) {
  stk <- rast(vapply(CLIMATE_YEARS, function(y) fetch(v, y), character(1)))
  mean(stk, na.rm = TRUE)
}))
names(clim) <- VARS

bb <- ext(msv)
e <- ext(xmin(bb) - 5, xmax(bb) + 5,
         max(ymin(bb) - 5, -90), min(ymax(bb) + 5, 90))
clim_c <- crop(clim, e)
crs(clim_c) <- "EPSG:4326"
fact <- max(1, round(0.5 / res(clim_c)[1]))       # coarsen to ~0.5 deg
if (fact > 1)
  clim_c <- aggregate(clim_c, fact = fact, fun = mean, na.rm = TRUE)
cat("clim dims:", paste(dim(clim_c), collapse = "x"), "\n")
writeRaster(clim_c, file.path(OUT_DIR, "vulpes_climate.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))

cat("bundle written to", OUT_DIR, "\n")
print(file.info(list.files(OUT_DIR, full.names = TRUE))[, "size", drop = FALSE])
