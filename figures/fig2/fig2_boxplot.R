# =============================================================================
# Script:  Contig N50 Distribution by Sequencing Type — static (PDF/SVG) + interactive (HTML)
# Author:  Yedomon Ange Bovys Zoclanclounon
#          Bioinformatics Scientist, Rothamsted Research
#          Twitter/X: @angeomics
# Date:    2026-04-02
# Purpose: Visualise the distribution of Contig N50 (Mbp) across sequencing
#          types using ggstatsplot (with pairwise comparisons) for static
#          export and a native plotly build for interactive exploration.
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
  ggplot2,      # static plotting base
  ggstatsplot,  # statistical visualisation with annotations
  scales,       # axis label formatting
  plotly,       # interactive plotting
  htmlwidgets   # saving interactive plots as self-contained HTML
)

# ── 3. Load data ──────────────────────────────────────────────────────────────
master_data <- read.csv("master_data_filtered.csv", h = TRUE)

# ── 4. Prepare data ───────────────────────────────────────────────────────────
THEDATA <- master_data %>%
  mutate(
    # Fill missing or blank Seq_type with "not_provided"
    Seq_type = ifelse(
      is.na(Seq_type) | Seq_type == "",
      "not_provided",
      Seq_type
    ),
    # Convert bp → Mbp
    Contig_N50_Mbp = Assembly_Stats_Contig_N50 / 1e6
  ) %>%
  # Remove non-positive values
  filter(Contig_N50_Mbp > 0)

# ── 5. Colour palette (shared by static and interactive plots) ────────────────
seq_colors <- c(
  long_read    = "#755cc5ff",
  short_read   = "#da57a1ff"
)

# =============================================================================
# 6. Static plot (ggstatsplot)
# =============================================================================
plot_stats <- ggbetweenstats(
  data                 = THEDATA,
  x                    = Seq_type,
  y                    = Contig_N50_Mbp,
  title                = "Distribution of Contig N50 (Mbp) across Sequencing type",
  xlab                 = "Sequencing type",
  ylab                 = "Contig N50 (Mbp)",
  type                 = "parametric",       # one-way ANOVA + pairwise t-tests
  pairwise.comparisons = TRUE,
  messages             = FALSE
) +
  # Apply shared colour palette
  ggplot2::scale_color_manual(values = unname(seq_colors)) +
  # Custom y-axis breaks every 2.5 Mbp
  ggplot2::scale_y_continuous(
    breaks = seq(0, max(THEDATA$Contig_N50_Mbp, na.rm = TRUE), by = 2.5),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  theme_classic(base_size = 13)

# ── 7. Save static plots ──────────────────────────────────────────────────────

# PDF — high-resolution (600 DPI), vector, Cairo renderer for font support
ggsave(
  filename = "N50_distribution_by_seqtype.pdf",
  plot     = plot_stats,
  device   = cairo_pdf,
  dpi      = 600,
  width    = 8,
  height   = 8,
  units    = "in"
)

# SVG — vector format, ideal for editing in Inkscape / Illustrator
ggsave(
  filename = "N50_distribution_by_seqtype.svg",
  plot     = plot_stats,
  device   = "svg",
  width    = 8,
  height   = 8,
  units    = "in"
)

# =============================================================================
# 8. Interactive plot (plotly) — built natively for richer hover tooltips
# =============================================================================

# Pre-compute per-group summary statistics for annotations
group_summary <- THEDATA %>%
  group_by(Seq_type) %>%
  summarise(
    n      = n(),
    mean   = mean(Contig_N50_Mbp,   na.rm = TRUE),
    median = median(Contig_N50_Mbp, na.rm = TRUE),
    sd     = sd(Contig_N50_Mbp,     na.rm = TRUE),
    .groups = "drop"
  )

# Unique Seq_type levels — one trace pair (box + jitter) built per group
seq_types <- unique(THEDATA$Seq_type)

plot_interactive <- plot_ly()

for (grp in seq_types) {
  
  grp_data <- THEDATA        %>% filter(Seq_type == grp)
  grp_col  <- seq_colors[grp]
  grp_sum  <- group_summary  %>% filter(Seq_type == grp)
  
  # Boxplot layer — outlier dots hidden (replaced by jittered points below)
  plot_interactive <- plot_interactive %>%
    add_trace(
      data      = grp_data,
      x         = ~Seq_type,
      y         = ~Contig_N50_Mbp,
      type      = "box",
      name      = grp,
      marker    = list(color = grp_col, opacity = 0),  # hide default outlier dots
      line      = list(color = grp_col),
      fillcolor = adjustcolor(grp_col, alpha.f = 0.2), # light semi-transparent fill
      boxmean   = "sd",                                 # overlay mean ± SD marker
      showlegend = FALSE
    )
  
  # Jittered points layer — individual observations with rich tooltips
  plot_interactive <- plot_interactive %>%
    add_trace(
      data  = grp_data,
      x     = ~Seq_type,
      y     = ~Contig_N50_Mbp,
      type  = "scatter",
      mode  = "markers",
      name  = grp,
      marker = list(
        color   = grp_col,
        size    = 7,
        opacity = 0.75,
        line    = list(color = "white", width = 0.5)
      ),
      # Hover tooltip: species name, N50, sequencing type, accession
      text = ~paste0(
        "<b>", Current_name, "</b><br>",
        "Contig N50 : ", round(Contig_N50_Mbp, 3), " Mbp<br>",
        "Seq type   : ", Seq_type,                 "<br>",
        "Accession  : ", Assembly_Accession
      ),
      hoverinfo  = "text",
      # JS-side jitter to spread overlapping points horizontally
      transforms = list(
        list(type = "jitter", target = "x", strength = 0.3)
      )
    )
}

# Per-group n / mean / sd annotations placed just above the tallest point
annotations <- lapply(seq_along(seq_types), function(i) {
  grp     <- seq_types[i]
  grp_sum <- group_summary %>% filter(Seq_type == grp)
  list(
    x         = grp,
    y         = max(THEDATA$Contig_N50_Mbp, na.rm = TRUE) * 1.05,
    text      = paste0(
      "n = ",    grp_sum$n,                  "<br>",
      "mean = ", round(grp_sum$mean, 2),     "<br>",
      "sd = ",   round(grp_sum$sd,   2)
    ),
    showarrow = FALSE,
    font      = list(size = 11),
    align     = "center"
  )
})

plot_interactive <- plot_interactive %>%
  layout(
    title = list(
      text = "Distribution of Contig N50 (Mbp) across Sequencing type",
      font = list(size = 15)
    ),
    xaxis      = list(title = "Sequencing type"),
    yaxis      = list(
      title      = "Contig N50 (Mbp)",
      tickformat = ".1f",
      dtick      = 2.5,
      rangemode  = "tozero"
    ),
    boxmode     = "group",
    annotations = annotations,
    hoverlabel  = list(bgcolor = "white"),
    legend      = list(title = list(text = "Sequencing type"))
  )

# ── 9. Save interactive plot ──────────────────────────────────────────────────
# selfcontained = TRUE bundles all JS/CSS into a single portable HTML file
saveWidget(
  widget        = plot_interactive,
  file          = "N50_distribution_by_seqtype_interactive.html",
  selfcontained = TRUE
)

# =============================================================================
# End of script
# =============================================================================