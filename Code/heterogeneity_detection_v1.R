# ------------------------------------------------------------
# heterogeneity_detection_v1.R  (v1.)
# ------------------------------------------------------------
# Author: Janeth Alpala
# Date: 30‑Apr‑2025
# Description:
#   Script for heterogeneity detection in SAR images (ENVI format)
#   using Shannon or Rényi non‑parametric estimators (with or without
#   bootstrap resampling).
# 
# Quick-start:
#   1. Edit the `opt` list below with your own paths, entropy type, window size, etc.
#      Example:
#        opt <- list(
#          image     = "path/to/image.img",
#          header    = "path/to/image.hdr",
#          entropy   = "renyi",
#          looks     = 16,
#          bootstrap = TRUE,
#          B         = 10,
#          window    = 7,
#          renyi     = 0.9
#        )
#   2. Run the script in R:  source("heterogeneity_detection_v1.R")   
#   3. Watch the progress bar; when finished a 2×2 composite plot pops up:
#        • Equalised SAR image
#        • p-value map (grey)
#        • p-value map (colour + scale bar)
#        • Significance mask  (p < opt$p_thr)
# Outputs:
#   • On‑screen plots (e.g., from plot_composite and quick previews).
#   • High‑resolution PNGs saved to ./PNG/ (SAR image, statistics, p-values, threshold maps).
#   • An .RData file with img_mat, entropy_test, p_values, and opt saved in ./Data/.


# ----------------------------------------------------------------------
# Required packages
# ----------------------------------------------------------------------
req_pkgs <- c("viridisLite", "png", "progress", "fields")
for (p in req_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(viridisLite)
library(png)
library(progress)
library(fields)   # for image.plot legend

# ----------------------------------------------------------------------
# External functions 
# ----------------------------------------------------------------------
# 

source("./Code/al_omari_1_estimator.R")
source("./Code/bootstrap_al_omari_1_estimator.R")
source("./Code/renyi_entropy_estimator_v1.R")
source("./Code/bootstrap_renyi_entropy_estimator_v1.R")
source("./Code/read_ENVI_images.R")

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
shannon_theoretical <- function(L, mu) {
  log(mu) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L))
}
renyi_theoretical <- function(L, lambda, mu) {
  num <- lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
    (lambda * (L - 1) + 1) * log(lambda)
  (num / (lambda - 1)) + log(mu) - log(L)
}
calculate_p_values_matrix <- function(stat_mat) {
  mu  <- mean(stat_mat, na.rm = TRUE)
  sig <- sd(stat_mat,   na.rm = TRUE)
  eps <- (stat_mat - 0) / sig
  2 * pnorm(-abs(eps))
}

# ----------------------------------------------------------------------
# Manual parameters 
# ----------------------------------------------------------------------



#image         = "../../../Data/SAR/L24_envi_Foggia_size_600/Intensity_VV.img",
#header        = "../../../Data/SAR/L24_envi_Foggia_size_600/Intensity_VV.hdr",

opt <- list(
  image         = "./Data/SAR/L16_envi_dublin_size_600/Intensity_HH.img",
  header        = "./Data/SAR/L16_envi_dublin_size_600/Intensity_HH.hdr",
  entropy       = "renyi",   # "shannon" | "renyi"
  looks         = 16,        # number of looks L
  bootstrap     = TRUE,      # TRUE to enable bootstrap
  B             = 10,        # bootstrap replicates
  window        = 7,         # Window side length (results in 7x7 = 49 pixels per window)
  prefix        = "Foggia_size_600",# Prefix for saved files, e.g., "image_size_sar.png"
  renyi         = 0.9,        # Rényi parameter λ (only for "renyi"); recommend 0.9 if L > 1, 3 if L = 1
  progress_step = 10,         # progress bar refresh
  p_thr         = 0.05        # significance threshold 
)

# ----------------------------------------------------------------------
# Basic checks
# ----------------------------------------------------------------------
stopifnot(!is.null(opt$image), !is.null(opt$header))
if (opt$window < 1) {
  stop("Window size must be a positive integer")
}
if (opt$entropy == "renyi") {
  if (opt$renyi <= 0 || opt$renyi == 1) {
    stop("Rényi parameter λ must satisfy: λ > 0 and λ ≠ 1")
  }
}

# ----------------------------------------------------------------------
# Read image
# ----------------------------------------------------------------------
cat("Reading image", opt$image, "...\n")
img_mat <- myread.ENVI(file = opt$image, headerfile = opt$header)
rows <- nrow(img_mat); cols <- ncol(img_mat)
cat("Dimensions:", rows, "x", cols, "pixels\n")

# ----------------------------------------------------------------------
# Select estimator dynamically
# ----------------------------------------------------------------------

get_estimator <- function(type, L, bootstrap) {
  if (type == "shannon") {
    if (bootstrap && L > 1) return(function(z) bootstrap_al_omari_1_estimator(z, opt$B))
    else                     return(al_omari_1_estimator)
  }
  if (type == "renyi") {
    if (bootstrap && L > 1) return(function(z) bootstrap_renyi_entropy_estimator_v1(z, opt$B, opt$renyi))
    else                     return(function(z) renyi_entropy_estimator_v1(z, opt$renyi))
  }
  stop("Unknown entropy type")
}
entropy_estimator <- get_estimator(opt$entropy, opt$looks, opt$bootstrap)
get_theoretical   <- function(type, L, lambda, mu) {
  if (type == "shannon") return(shannon_theoretical(L, mu))
  if (type == "renyi")   return(renyi_theoretical(L, lambda, mu))
}

# ----------------------------------------------------------------------
# Sliding‑window processing
# ----------------------------------------------------------------------
cat("Processing", opt$window, "x", opt$window, "windows...\n")
start_time <- Sys.time()

out_nrow <- rows - opt$window + 1
out_ncol <- cols - opt$window + 1
entropy_test <- matrix(NA_real_, nrow = out_nrow, ncol = out_ncol)

pb <- progress_bar$new(total = out_nrow,
                       format = "[:bar] :current/:total (:percent) in :elapsed")
for (i in 1:out_nrow) {
  for (j in 1:out_ncol) {
    slice <- img_mat[i:(i + opt$window - 1), j:(j + opt$window - 1)]
    est   <- entropy_estimator(slice)
    mu    <- mean(slice)
    theo  <- get_theoretical(opt$entropy, opt$looks, opt$renyi, mu)
    entropy_test[i, j] <- est - theo
  }
  pb$tick()
}

# ----------------------------------------------------------------------
# p‑value matrix
# ----------------------------------------------------------------------
cat("Computing p‑values...\n")
p_values <- calculate_p_values_matrix(entropy_test)



# ----------------------------------------------------------------------
# Quick exploratory plots
# ----------------------------------------------------------------------

source("Code/imagematrix_visualizer.R")     # visual helper library
# dev.new()
# plot(imagematrix(equalize(img_mat)), main = "SAR image") # equalize SAR image gray scale
# plot(imagematrix(p_values),main = "p-values (greyscale)" ) # plot p-values gray scale
# plot(imagematrix(p_values< 0.05), ,main = "p-values < 0.05")
# previewImagematrix(imagematrix_color(p_values), palette_option = "viridis-H",main = "p-values (colour)")# color image with colorbar


# ----------------------------------------------------------------------
# Composite plot function
# ----------------------------------------------------------------------

plot_composite <- function(img_mat, p_values,
                           thr = 0.05,
                           palette_option = "viridis-H") {
  
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
  
  # SAR image
  plot(imagematrix(equalize(img_mat)), main = "SAR image")
  
  # grayscale p-values
  plot(imagematrix(p_values), main = "p-values (grayscale)")
 
  # colour p-values 
  
  previewImagematrixPanel(
    imagematrix_color(p_values),
    palette_option     = palette_option,
    significance_level = thr
    #main               = "p-values (colour)"
  )
  par(mar = c(1, 1, 3, 1))
  #  threshold map
  plot(imagematrix(p_values < thr), main = paste0("p-values < ", thr))
}

# ----------------------------------------------------------------------
# Show composite in RStudio panel
# ----------------------------------------------------------------------
cat("Plotting 2×2 composite figure...\n")
plot_composite(img_mat, p_values, thr = opt$p_thr)

#plot_composite(img_mat, p_values, thr = opt$p_thr, palette_option = "viridis-H")


# ----------------------------------------------------------------------
# Save Individual PNG exports (optional)
# ----------------------------------------------------------------------
cat("Saving PNG outputs...\n")
if (!dir.exists("./PNG")) dir.create("./PNG")

imagematrixPNG(
  imagematrix(equalize(img_mat)),
  file.path("./PNG", paste0(opt$prefix, "_sar.png"))
)

imagematrixPNG(
  imagematrix(equalize(entropy_test)),
  file.path("./PNG", paste0(opt$prefix, "_statistic.png"))
)

imagematrixPNG(
  imagematrix(p_values),
  file.path("./PNG", paste0(opt$prefix, "_pvals_gray.png"))
)

imagematrixPNG(
  imagematrix(p_values < opt$p_thr),
  file.path("./PNG", paste0(opt$prefix, "_pvals_thr.png"))
)

imagematrix_colorPNG(
  imagematrix_color(p_values),
  name            = file.path("./PNG", paste0(opt$prefix, "_pvals_color.png")),
  palette_option  = "viridis-H",
  legend_width_px = 540,
  scale_factor    = 2
)

# ----------------------------------------------------------------------
# Save workspace objects
# ----------------------------------------------------------------------
if (!dir.exists("./Data")) dir.create("./Data")
base_name <- sprintf("%s_%s_L%d%s_w%d%s%s",
                     opt$prefix,
                     ifelse(opt$entropy == "shannon", "shannon", "renyi"),
                     opt$looks,
                     ifelse(opt$bootstrap, paste0("_b", opt$B), ""),
                     opt$window,
                     ifelse(opt$entropy == "renyi", paste0("_lambda", gsub("\\.", "", sprintf("%.2f", opt$renyi))), ""),
                     format(Sys.Date(), "_%Y%m%d"))
save(file = file.path("./Data", paste0(base_name, ".RData")),
     list = c("img_mat", "entropy_test", "p_values", "opt"))

# -------------------- Execution time -----------------------------------
end_time <- Sys.time()
cat("Total run time:", round(difftime(end_time, start_time, units = "secs"),1), "seconds\n")
