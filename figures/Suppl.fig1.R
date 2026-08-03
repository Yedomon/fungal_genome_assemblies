# =============================================================================
# Script:  Assembly Size vs Contig N50 — static (PDF/SVG) + interactive (HTML)
# Author:  Yedomon Ange Bovys Zoclanclounon
#          Bioinformatics Scientist, Rothamsted Research
#          Twitter/X: @angeomics
# Date:    2026-04-02
# Purpose: Visualise the relationship between assembly size and contig N50
#          across sequencing types using log-log axes; export static and
#          interactive versions of the figure.
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
  ggrepel,      # non-overlapping text labels
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
    Assembly_Size_Mbp = Assembly_Stats_Total_Sequence_Length / 1e6,
    Contig_N50_Mbp    = Assembly_Stats_Contig_N50 / 1e6
  ) %>%
  # Remove non-positive values (required for log10 axes)
  filter(Assembly_Size_Mbp > 0, Contig_N50_Mbp > 0)

# ── 5. Colour palette (shared by static and interactive plots) ────────────────
custom_colors <- c(
  not_provided = "#a8c66c",
  long_read    = "#755cc5ff",
  short_read   = "#da57a1ff"
)

# =============================================================================
# 6. Static plot (ggplot2)
# =============================================================================
plot1 <- ggplot(plot_data,
                aes(x     = Assembly_Size_Mbp,
                    y     = Contig_N50_Mbp,
                    color = Seq_type)) +
  
  # Points
  geom_point(size = 3, alpha = 0.85) +
  
  # Italic species-name labels only for extreme assembly sizes
  # (very small: < 0.3 Mbp  |  very large: >= 1000 Mbp)
  geom_text_repel(
    data          = plot_data %>%
      filter(Assembly_Size_Mbp < 0.3 | Assembly_Size_Mbp >= 1000),
    aes(label     = Current_name),
    size          = 3,
    fontface      = "italic",
    max.overlaps  = Inf,
    box.padding   = 0.4,
    point.padding = 0.3,
    segment.color = "grey70",
    show.legend   = FALSE
  ) +
  
  # Log10 axes with readable tick labels
  scale_x_log10(labels = label_number(accuracy = 0.1)) +
  scale_y_log10(labels = label_number(accuracy = 0.01)) +
  
  # Apply shared colour palette
  scale_color_manual(values = custom_colors) +
  
  # Axis / legend / title labels
  labs(
    x     = "Assembly size (Mbp, log10)",
    y     = "Contig N50 (Mbp, log10)",
    color = "Sequencing type",
    title = "Assembly size vs Contig N50"
  ) +
  
  theme_classic(base_size = 13)

# ── 7. Save static plots ──────────────────────────────────────────────────────

# PDF — high-resolution (600 DPI), vector, Cairo renderer for font support
ggsave(
  filename = "assembly_size_vs_contig_N50_loglog.pdf",
  plot     = plot1,
  device   = cairo_pdf,
  dpi      = 600,
  width    = 8,
  height   = 8,
  units    = "in"
)

# SVG — vector format, ideal for editing in Inkscape / Illustrator
ggsave(
  filename = "assembly_size_vs_contig_N50_loglog.svg",
  plot     = plot1,
  device   = "svg",
  width    = 8,
  height   = 8,
  units    = "in"
)

# =============================================================================
# 8. Interactive plot (plotly) — built natively for richer hover tooltips
# =============================================================================
plot_interactive <- plot_ly(
  data   = plot_data,
  x      = ~Assembly_Size_Mbp,
  y      = ~Contig_N50_Mbp,
  color  = ~Seq_type,
  colors = custom_colors,
  type   = "scatter",
  mode   = "markers",
  marker = list(size = 8, opacity = 0.85),
  
  # Custom hover tooltip: species name, metrics, sequencing type, accession
  text = ~paste0(
    "<b>", Current_name, "</b><br>",
    "Assembly size : ", round(Assembly_Size_Mbp, 3), " Mbp<br>",
    "Contig N50    : ", round(Contig_N50_Mbp,    3), " Mbp<br>",
    "Seq type      : ", Seq_type,                    "<br>",
    "Accession     : ", Assembly_Accession
  ),
  hoverinfo = "text"
  
) %>%
  layout(
    title      = "Assembly size vs Contig N50",
    xaxis      = list(title = "Assembly size (Mbp, log10)", type = "log"),
    yaxis      = list(title = "Contig N50 (Mbp, log10)",    type = "log"),
    legend     = list(title = list(text = "Sequencing type")),
    hoverlabel = list(bgcolor = "white")
  )

# ── 9. Save interactive plot ──────────────────────────────────────────────────
# selfcontained = TRUE bundles all JS/CSS into a single portable HTML file
saveWidget(
  widget        = plot_interactive,
  file          = "assembly_size_vs_contig_N50_interactive.html",
  selfcontained = TRUE
)

# =============================================================================
# End of script
# =============================================================================