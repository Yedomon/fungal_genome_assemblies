# =============================================================================
# Script:  Contig N50 vs Assembly Release Year — static (PDF/SVG) + interactive (HTML)
# Author:  Yedomon Ange Bovys Zoclanclounon
#          Bioinformatics Scientist, Rothamsted Research
#          Twitter/X: @angeomics
# Date:    2026-04-02 (Updated)
# Purpose: Visualise the trend of Contig N50 (Mbp) over assembly release years
#          with italicized species names in the interactive tooltips.
# =============================================================================

# ── 1. Working environment ────────────────────────────────────────────────────

# Clear workspace and free memory
remove(list = ls())
gc()

# Set random seed for reproducibility
set.seed(1000)

# Set working directory
setwd("U:/brainstorming_ange/Sent")

# ── 2. Load required packages ─────────────────────────────────────────────────
require(pacman)
pacman::p_load(
  dplyr,        # data wrangling
  ggplot2,      # static plotting
  scales,       # axis label formatting
  plotly,       # interactive plotting
  htmlwidgets   # saving interactive plots as self-contained HTML
)

# ── 3. Load data ──────────────────────────────────────────────────────────────
master_data <- read.csv("master_data_filtered.csv", h = TRUE)

# ── 4. Prepare data ───────────────────────────────────────────────────────────
plot_data <- master_data %>%
  mutate(
    # Fill missing or blank Seq_type with "not_provided"
    Seq_type = ifelse(
      is.na(Seq_type) | Seq_type == "",
      "not_provided",
      Seq_type
    ),
    # Convert bp → Mbp
    Contig_N50_Mbp = Assembly_Stats_Contig_N50 / 1e6,
    # Extract release year as numeric
    Release_Year   = as.numeric(Assembly_Release_Year),
    # Create italicized species name for the interactive tooltip using HTML tags
    Species_Italic = paste0("<i>", Current_name, "</i>")
  ) %>%
  # Keep only positive N50 values and assemblies from 2002 onwards
  filter(
    Contig_N50_Mbp > 0,
    Release_Year   >= 2002
  )

# ── 5. Colour palette (shared by static and interactive plots) ────────────────
custom_colors <- c(
  long_read    = "#755cc5ff",
  short_read   = "#da57a1ff",
  not_provided = "#a8c66c"
)

# =============================================================================
# 6. Static plot (ggplot2)
# =============================================================================
# We map 'Species_Italic' to the 'label' aesthetic so plotly can find it later.
# Since we aren't using geom_text(), these labels won't appear on the PDF/SVG.
scatter_plot <- ggplot(plot_data,
                       aes(x     = Release_Year,
                           y     = Contig_N50_Mbp,
                           color = Seq_type,
                           label = Species_Italic)) +
  
  # Points
  geom_point(size = 1.5, alpha = 0.8) +
  
  # Apply shared colour palette
  scale_color_manual(values = custom_colors) +
  
  # Y-axis label
  scale_y_continuous(name = "Contig N50 (Mbp)") +
  
  # X-axis: one break per year across the full range
  scale_x_continuous(
    name   = "Assembly release year",
    breaks = seq(min(plot_data$Release_Year), max(plot_data$Release_Year), by = 1)
  ) +
  
  # Plot labels
  labs(
    title = "Contig N50 (Mbp) vs Assembly Release Year",
    color = "Sequencing type"
  ) +
  
  theme_classic(base_size = 13)

# ── 7. Save static plots ──────────────────────────────────────────────────────

# PDF — high-resolution (600 DPI), vector, Cairo renderer for font support
ggsave(
  filename = "ContigN50_vs_year_species_names_v3.pdf",
  plot     = scatter_plot,
  device   = cairo_pdf,
  dpi      = 600,
  width    = 14,
  height   = 8,
  units    = "in"
)

# SVG — vector format, ideal for editing in Inkscape / Illustrator
ggsave(
  filename = "ContigN50_vs_year_species_names_v3.svg",
  plot     = scatter_plot,
  device   = "svg",
  width    = 14,
  height   = 8,
  units    = "in"
)

# =============================================================================
# 8. Interactive plot (plotly)
# =============================================================================
# ggplotly() converts the ggplot object; we explicitly call "label" in the tooltip
interactive_plot <- ggplotly(
  scatter_plot,
  tooltip = c("label", "x", "y", "color")
)

# ── 9. Save interactive plot ──────────────────────────────────────────────────
# selfcontained = TRUE bundles all JS/CSS into a single portable HTML file
saveWidget(
  widget        = interactive_plot,
  file          = "ContigN50_vs_year_species_names_interactive_v3.html",
  selfcontained = TRUE
)

# =============================================================================
# End of script
# =============================================================================