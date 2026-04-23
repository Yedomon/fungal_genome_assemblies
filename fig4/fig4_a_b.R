# =============================================================================
# Script:  Contig N50 vs BUSCO Completeness — Global & High-Res Zoom (0.1 Mbp)
# Author:  Yedomon Ange Bovys Zoclanclounon
#          Bioinformatics Scientist, Rothamsted Research
#          Twitter/X: @angeomics
# Date:    2026-04-02
# Purpose: Compare assembly contiguity (N50) against gene completeness (BUSCO).
#          Outputs static (PDF/SVG) and interactive (HTML) versions for the
#          full dataset and a zoomed micro-scale (≤ 0.1 Mbp).
# =============================================================================

# ── 1. Working environment ────────────────────────────────────────────────────

# Clear workspace and free memory
remove(list = ls())
gc()

# Set random seed for consistent label placement in ggrepel
set.seed(1000)

# ── 2. Load required packages ─────────────────────────────────────────────────
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,        # data wrangling
  ggplot2,      # static plotting
  ggrepel,      # intelligent label placement
  scales,       # axis label formatting
  plotly,       # interactive plotting
  htmlwidgets,  # saving interactive plots
  svglite       # SVG device support
)

# ── 3. Load and Prepare Master Data ───────────────────────────────────────────
master_data <- read.csv("master_data_filtered.csv", h = TRUE)

# Shared data cleaning
plot_data_base <- master_data %>%
  mutate(
    # Clean Seq_type and set factor order
    Seq_type = ifelse(is.na(Seq_type) | Seq_type == "", "not_provided", Seq_type),
    Seq_type = factor(Seq_type, levels = c("not_provided", "long_read", "short_read")),
    
    # Unit conversion: bp -> Mbp
    Contig_N50_Mbp = Assembly_Stats_Contig_N50 / 1e6,
    
    # Create italicized label string for species names
    label_italic = paste0("italic('", Current_name, "')")
  ) %>%
  filter(
    Contig_N50_Mbp > 0, 
    C >= 0, 
    Seq_type %in% c("long_read", "short_read")
  )

# ── 4. Global Aesthetic Settings ──────────────────────────────────────────────
custom_colors <- c(
  long_read    = "#755cc5ff",
  short_read   = "#da57a1ff"
)

# =============================================================================
# 5. GLOBAL VIEW: FULL RANGE
# =============================================================================

plot_global <- ggplot(plot_data_base, 
                      aes(x = Contig_N50_Mbp, y = C, color = Seq_type, label = Current_name)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(
    size = 3, max.overlaps = 15, box.padding = 0.3, segment.color = "grey70"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(name = "Contig N50 (Mbp)") +
  scale_y_continuous(name = "BUSCO Complete (C, %)") +
  labs(
    title = "Contig N50 vs BUSCO Complete (Global View)",
    color = "Sequencing type"
  ) +
  theme_classic(base_size = 13)

# Save Global Static
ggsave("ContigN50_vs_BUSCO_C_Global.pdf", plot_global, device = cairo_pdf, dpi = 600, width = 8, height = 6)
ggsave("ContigN50_vs_BUSCO_C_Global.svg", plot_global, device = "svg", width = 8, height = 6)

# Save Global Interactive
saveWidget(
  widget = ggplotly(plot_global, tooltip = c("x", "y", "label", "color")), 
  file = "ContigN50_vs_BUSCO_C_Global_interactive.html",
  selfcontained = TRUE
)

# =============================================================================
# 6. HIGH-RESOLUTION ZOOM: 0 to 0.1 Mbp
# =============================================================================

# Filter specifically for the zoom-in range
zoom_data <- plot_data_base %>% filter(Contig_N50_Mbp <= 0.1)

plot_zoom <- ggplot(zoom_data, 
                    aes(x = Contig_N50_Mbp, y = C, color = Seq_type, label = Current_name)) +
  geom_point(size = 0.8, alpha = 0.8) +
  geom_text_repel(
    aes(label = label_italic), 
    parse = TRUE,            # Renders 'italic()' syntax
    size = 3, 
    max.overlaps = 20, 
    box.padding = 0.4, 
    point.padding = 0.3,
    segment.color = "grey70"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(
    name = "Contig N50 (Mbp)", 
    limits = c(0, 0.1), 
    breaks = seq(0, 0.1, by = 0.01),
    labels = label_number(accuracy = 0.001)
  ) +
  scale_y_continuous(name = "BUSCO Complete (C, %)") +
  labs(
    title = "Contig N50 vs BUSCO Complete (Zoom: ≤ 0.1 Mbp)",
    color = "Sequencing type"
  ) +
  theme_classic(base_size = 13)

# Save Zoom Static
ggsave("ContigN50_vs_BUSCO_C_Zoom_0.1Mbp.pdf", plot_zoom, device = cairo_pdf, dpi = 600, width = 8, height = 8)
ggsave("ContigN50_vs_BUSCO_C_Zoom_0.1Mbp.svg", plot_zoom, device = "svg", width = 8, height = 8)

# Save Zoom Interactive
saveWidget(
  widget = ggplotly(plot_zoom, tooltip = c("x", "y", "label", "color")), 
  file = "ContigN50_vs_BUSCO_C_Zoom_0.1Mbp_interactive.html",
  selfcontained = TRUE
)

# =============================================================================
# End of script
# =============================================================================