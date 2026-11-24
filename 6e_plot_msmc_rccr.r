library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(patchwork)

# RCCR + Ne plots for BCCH / CACH / MOCH, styled after the Sulidae
# RCCR_BFPE_MANA_Temperature scripts (one stacked panel per species pair, each
# curve drawn for three generation times, glacial-maximum markers, optional
# temperature panel).
#
# Copy analyses/msmc_rccr/output (and optionally bootstrap/outputs) locally:
#   see 6_MSMC_RCCR.md, step 5.

setwd("~/Desktop/chickadee_msmc_rccr/output")

# SCRIPT PARAMETERS ===========================================
mu   <- 4.6e-9          # per site per generation; CHANGE (must match MU in 6c if MSMC-IM is used)
gens <- c(1, 2, 3)      # generation times (years) to draw; CHANGE. Middle value = solid line
xmax <- 1e6
boot_dir <- "../bootstrap/outputs"   # set to NA to skip bootstrap lines

pairs <- c(BCCH_CACH = "Black-capped vs. Carolina",
           BCCH_MOCH = "Black-capped vs. Mountain",
           CACH_MOCH = "Carolina vs. Mountain")
pair_color <- "#7570b3"
species_colors <- c(BCCH = "#1b9e77", CACH = "#d95f02", MOCH = "#7570b3")

PGM_time  <- 140000
PGM_color <- "#e78ac3"
LGM_time  <- 21000
LGM_color <- "#9d196aff"

temps_file <- "temperatures.csv"     # same file as the Sulidae figure; skipped if absent
# =============================================================

gen_levels <- as.character(gens)
lty_vals <- setNames(c("dashed", "solid", "dotted"), gen_levels)
lwd_vals <- setNames(c(0.5, 0.9, 0.5), gen_levels)

glacial_lines <- list(
  geom_vline(xintercept = LGM_time, linewidth = 2.5, color = LGM_color, alpha = 0.35),
  geom_vline(xintercept = PGM_time, linewidth = 2.5, color = PGM_color, alpha = 0.35)
)

# ----------------------------
# Helper to build RCCR df for a pair (one copy per generation time)
# ----------------------------
make_rccr_df <- function(file, mu) {
  data <- read.table(file, header = TRUE)

  df_base <- data %>%
    mutate(
      CCR_raw = 2 * lambda_01 / (lambda_00 + lambda_11),
      CCR = pmin(pmax(CCR_raw, 0), 1)
    )

  bind_rows(lapply(gens, function(g)
    df_base %>% mutate(gen = g, time_years = (left_time_boundary / mu) * g))) %>%
    mutate(gen = factor(gen, levels = gens))
}

# bootstrap replicates, drawn at the middle generation time only
read_boot_rccr <- function(pair_lc) {
  if (is.na(boot_dir) || !dir.exists(boot_dir)) return(NULL)
  files <- list.files(boot_dir, pattern = paste0("^combined\\.", pair_lc, "\\.final\\.txt$"),
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  bind_rows(lapply(files, function(f) read.table(f, header = TRUE) %>%
                     mutate(rep = basename(dirname(f))))) %>%
    mutate(CCR = pmin(pmax(2 * lambda_01 / (lambda_00 + lambda_11), 0), 1),
           time_years = (left_time_boundary / mu) * gens[2])
}

# ----------------------------
# One RCCR panel per pair
# ----------------------------
rccr_panels <- list()
split_summary <- list()
for (i in seq_along(pairs)) {
  pair <- names(pairs)[i]
  pair_lc <- tolower(pair)
  f <- paste0("combined.", pair_lc, ".final.txt")
  if (!file.exists(f)) { message("Missing ", f, "; skipping"); next }

  df <- make_rccr_df(f, mu)
  if (any(df$CCR_raw > 1.05, na.rm = TRUE))
    message(pair, ": raw RCCR exceeds 1 in some intervals (clamped for plotting)")

  boots <- read_boot_rccr(pair_lc)
  is_last <- i == length(pairs)

  p <- ggplot(df, aes(x = time_years, y = CCR)) +
    { if (!is.null(boots)) geom_step(data = boots, aes(group = rep),
                                     color = pair_color, alpha = 0.2, linewidth = 0.3) } +
    geom_step(aes(linetype = gen, linewidth = gen), color = pair_color) +
    scale_linetype_manual(values = lty_vals) +
    scale_linewidth_manual(values = lwd_vals) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50") +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
    labs(x = if (is_last) "Years before present" else NULL, y = "RCCR",
         linetype = "Generation time (years)", linewidth = "Generation time (years)",
         title = pairs[[pair]]) +
    glacial_lines +
    theme_classic() +
    theme(legend.justification = c(0, 0.5))
  if (i > 1) p <- p + guides(linetype = "none", linewidth = "none")
  rccr_panels[[pair]] <- p

  # split-time summary: first interval where RCCR < 0.5, per generation time
  split_summary[[pair]] <- df %>% group_by(gen) %>% arrange(time_years) %>%
    summarise(years_rccr_below_0.5 = time_years[which(CCR < 0.5)[1]], .groups = "drop") %>%
    mutate(pair = pair)
}
print(bind_rows(split_summary))

# ----------------------------
# Optional temperature panel (as in the Sulidae figure)
# ----------------------------
p_temp <- NULL
if (file.exists(temps_file)) {
  temps <- read_csv(temps_file, show_col_types = FALSE)
  time_col <- names(temps)[grepl("^Time", names(temps), ignore.case = TRUE)][1]
  temps <- temps %>% rename(Time_kyr_BP = all_of(time_col)) %>%
    mutate(time_years = Time_kyr_BP * 1000)
  q2_5 <- names(temps)[grepl("2.5", names(temps))]
  q50  <- names(temps)[grepl("^50$", names(temps))]
  q97_5 <- names(temps)[grepl("97.5", names(temps))]
  temps <- temps %>% transmute(time_years, t_lo = .data[[q2_5]], t_md = .data[[q50]], t_hi = .data[[q97_5]])
  temps_win <- temps %>% filter(time_years >= 0, time_years <= xmax)
  t_min <- min(temps_win$t_lo, na.rm = TRUE)
  t_max <- max(temps_win$t_hi, na.rm = TRUE)
  temp_breaks <- pretty(c(t_min, t_max), n = 5)

  p_temp <- ggplot(temps_win, aes(x = time_years)) +
    geom_ribbon(aes(ymin = (t_lo - t_min) / (t_max - t_min), ymax = (t_hi - t_min) / (t_max - t_min)),
                alpha = 0.25) +
    geom_line(aes(y = (t_md - t_min) / (t_max - t_min)), linewidth = 0.25) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
    scale_y_continuous(name = "Temp (°C)", breaks = (temp_breaks - t_min) / (t_max - t_min),
                       labels = round(temp_breaks, 1)) +
    labs(x = NULL, title = "Δ Global avg. surface temperature (relative to 0–5 ka)") +
    glacial_lines +
    theme_classic()
}

panels <- c(if (!is.null(p_temp)) list(p_temp), rccr_panels)
heights <- c(if (!is.null(p_temp)) 0.8, rep(1.1, length(rccr_panels)))
p_stack <- wrap_plots(panels, ncol = 1) + plot_layout(heights = heights)

ggsave("RCCR_BCCH_CACH_MOCH_stacked.png", p_stack, width = 7, height = 1.7 * sum(heights) + 1)
ggsave("RCCR_BCCH_CACH_MOCH_stacked.pdf", p_stack, width = 7, height = 1.7 * sum(heights) + 1)

# ----------------------------
# Ne per species (middle generation time), as in plotMSMC.allsamples
# ----------------------------
ne_files <- list.files(".", pattern = "^within\\..*\\.final\\.txt$")
ne_data <- bind_rows(lapply(ne_files, function(f) {
  dat <- read.table(f, header = TRUE)
  data.frame(time = dat$left_time_boundary / mu * gens[2],
             pop.size = (1 / dat$lambda) / (2 * mu),
             species = toupper(str_match(f, "^within\\.([a-z]+)\\.")[, 2]))
}))

p_ne <- ggplot(ne_data %>% filter(time > 0), aes(x = time, y = pop.size, color = species)) +
  geom_step(linewidth = 0.9) +
  scale_color_manual(values = species_colors) +
  scale_x_log10() +
  labs(x = "Years before present (log scale)", y = "Effective population size",
       title = paste0("MSMC2 per species (gen = ", gens[2], " yr)")) +
  theme_minimal(base_size = 14)

ggsave("MSMC_Ne_BCCH_CACH_MOCH.png", p_ne, width = 8, height = 6)

# ----------------------------
# Optional: MSMC-IM cumulative migration M(t), if 6c was run with RUN_IM=1
# ----------------------------
im_files <- list.files("msmcim", pattern = "MSMC_IM\\.estimates\\.txt$", full.names = TRUE)
if (length(im_files) > 0) {
  im_df <- bind_rows(lapply(im_files, function(f) read.delim(f) %>%
    mutate(pair = toupper(str_match(basename(f), "^([a-z]+_[a-z]+)\\.")[, 2])))) %>%
    mutate(time_years = left_time_boundary * gens[2])   # MSMC-IM times are in generations
  p_im <- ggplot(im_df %>% filter(time_years > 0), aes(x = time_years, y = M, color = pair)) +
    geom_step(linewidth = 0.9) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dotted", color = "grey50") +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
    labs(x = "Years before present", y = "Cumulative migration M(t)", color = "") +
    glacial_lines +
    theme_classic()
  ggsave("MSMCIM_M_BCCH_CACH_MOCH.png", p_im, width = 7, height = 3)
}
