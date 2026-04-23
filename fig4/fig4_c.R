# =============================================================================
# Script:  BUSCO Metrics Comparison — Global, Individual, and Stacked Analysis
# Author:  Yedomon Ange Bovys Zoclanclounon
#          Bioinformatics Scientist, Rothamsted Research
#          Twitter/X: @angeomics
# Date:    2026-04-02
# Purpose: Statistical comparison of BUSCO metrics across sequencing types.
#          Generates: 1) Individual static & interactive plots per metric.
#                     2) A publication-ready vertical stack using patchwork.
# =============================================================================

# ── 1. Working environment ────────────────────────────────────────────────────

# Clear workspace and free memory
remove(list = ls())
gc()

# Set random seed for reproducibility of jittered points
set.seed(1000)

# ── 2. Load required packages ─────────────────────────────────────────────────
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,        # Data wrangling
  ggplot2,      # Base plotting
  ggstatsplot,  # Statistical plotting
  patchwork,    # Combining plots (essential for wrap_plots)
  plotly,       # Interactive plotting
  htmlwidgets,  # Saving HTML widgets
  svglite       # SVG device support
)

# ── 3. Load and Prepare Data ──────────────────────────────────────────────────
master_data <- read.csv("master_data_filtered.csv", header = TRUE)

plot_data <- master_data %>%
  mutate(
    # Handle missing sequencing types and set factor levels
    Seq_type = ifelse(is.na(Seq_type) | Seq_type == "", "not_provided", Seq_type),
    Seq_type = factor(Seq_type, levels = c("not_provided", "long_read", "short_read"))
  ) %>%
  # Focus analysis on Long-read vs Short-read
  filter(Seq_type %in% c("long_read", "short_read"))

# ── 4. Global Configuration ───────────────────────────────────────────────────
custom_colors <- c(
  long_read    = "#755cc5ff",
  short_read   = "#da57a1ff"
)

# Map of column names to descriptive labels
busco_map <- list(
  "C" = "Complete BUSCO (%)",
  "S" = "Single-copy BUSCO (%)",
  "D" = "Duplicated BUSCO (%)",
  "F" = "Fragmented BUSCO (%)",
  "M" = "Missing BUSCO (%)"
)

# =============================================================================
# 5. GENERATE INDIVIDUAL PLOTS (LOOP)
# =============================================================================

# Initialize a list to store plots for the combined patchwork later
busco_plots_list <- list()

for (metric in names(busco_map)) {
  
  metric_label <- busco_map[[metric]]
  message(paste0("Processing metric: ", metric_label))
  
  # A. Create the ggbetweenstats object
  p <- ggbetweenstats(
    data  = plot_data,
    x     = Seq_type,
    y     = !!sym(metric),
    title = paste0("Distribution of ", metric_label),
    xlab  = "Sequencing Type",
    ylab  = metric_label,
    type  = "parametric",           # Adjust to "nonparametric" if data is skewed
    pairwise.comparisons = TRUE,
    messages = FALSE,
    ggtheme = theme_classic(base_size = 12)
  ) + 
    scale_color_manual(values = custom_colors)
  
  # Store plot in list for later
  busco_plots_list[[metric]] <- p
  
  # B. Save Individual Static Plots
  ggsave(
    filename = paste0("BUSCO_", metric, "_by_seqtype.pdf"),
    plot = p, device = cairo_pdf, dpi = 600, width = 8, height = 6
  )
  
  ggsave(
    filename = paste0("BUSCO_", metric, "_by_seqtype.svg"),
    plot = p, device = "svg", width = 8, height = 6
  )
  
  # C. Save Individual Interactive Plot
  # Note: ggplotly might omit the complex ggstatsplot subtitle/stats text
  p_inter <- ggplotly(p, tooltip = c("y", "x"))
  
  saveWidget(
    widget = p_inter,
    file = paste0("BUSCO_", metric, "_interactive.html"),
    selfcontained = TRUE
  )
}

# =============================================================================
# 6. COMBINE INTO ONE COLUMN (PUBLICATION FORMAT)
# =============================================================================

# Note: We use wrap_plots here as requested for maximum compatibility
# with the complex ggstatsplot objects stored in our list.

combined_plot <- wrap_plots(busco_plots_list, ncol = 1) +
  plot_annotation(
    title = "BUSCO metrics across sequencing technologies",
    subtitle = "Parametric comparison of assembly completeness and quality",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# Save the combined figure
# Setting height to 18-20 inches to accommodate the vertical stack of 5 plots
ggsave(
  filename = "BUSCO_all_metrics_combined_stack.pdf",
  plot = combined_plot,
  device = cairo_pdf,
  dpi = 600,
  width = 3,
  height = 18,
  units = "in"
)

ggsave(
  filename = "BUSCO_all_metrics_combined_stack.svg",
  plot = combined_plot,
  device = "svg",
  width = 3,
  height = 18,
  units = "in"
)

# =============================================================================
# End of script
# =============================================================================