# =============================================================================
# Script:  Filtered Assembly Quality — Specific N50 & BUSCO Range
# Author:  Yedomon Ange Bovys Zoclanclounon
# Date:    2026-04-02
# Purpose: Isolate assemblies within precise quantitative boundaries:
#          - Contig N50: [1.05, 1.55] Mbp
#          - BUSCO C:    [48.5, 51.7] %
# =============================================================================

# ── 1. Environment & Packages ─────────────────────────────────────────────────
remove(list = ls())
gc()
set.seed(1000)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2, ggrepel, scales, plotly, htmlwidgets, svglite)

# ── 2. Load and Apply Quantitative Filters ────────────────────────────────────
master_data <- read.csv("master_data_filtered.csv", h = TRUE)

filtered_subset <- master_data %>%
  mutate(
    Seq_type = ifelse(is.na(Seq_type) | Seq_type == "", "not_provided", Seq_type),
    Seq_type = factor(Seq_type, levels = c("not_provided", "long_read", "short_read")),
    Contig_N50_Mbp = Assembly_Stats_Contig_N50 / 1e6,
    label_italic = paste0("italic('", Current_name, "')")
  ) %>%
  # Applying your specific quantitative filters
  filter(
    Contig_N50_Mbp >= 1.05 & Contig_N50_Mbp <= 1.55,
    C >= 48.5 & C <= 51.7,
    Seq_type %in% c("long_read", "short_read")
  )

# ── 3. Aesthetic Settings ─────────────────────────────────────────────────────
custom_colors <- c(
  long_read    = "#755cc5ff",
  short_read   = "#da57a1ff"
)

# =============================================================================
# 4. STATIC PLOT: Filtered Range
# =============================================================================

plot_subset <- ggplot(filtered_subset, 
                      aes(x = Contig_N50_Mbp, y = C, color = Seq_type)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(
    aes(label = label_italic), 
    parse = TRUE, 
    size = 3.5, 
    max.overlaps = Inf, 
    box.padding = 0.6, 
    point.padding = 0.5,
    segment.color = "grey70"
  ) +
  scale_color_manual(values = custom_colors) +
  # Fixed limits to reflect your filter range
  scale_x_continuous(
    name = "Contig N50 (Mbp)", 
    limits = c(1.04, 1.56),
    breaks = seq(1.05, 1.55, by = 0.05)
  ) +
  scale_y_continuous(
    name = "BUSCO Complete (C, %)", 
    limits = c(48, 52),
    breaks = seq(48.5, 51.5, by = 0.5)
  ) +
  labs(
    title = "Filtered Assemblies: Targeted Quality Window",
    subtitle = "N50: 1.05–1.55 Mbp | BUSCO C: 48.5–51.7%",
    color = "Sequencing type"
  ) +
  theme_classic(base_size = 13)

# Save Static
ggsave("Filtered_Range_N50_BUSCO.pdf", plot_subset, device = cairo_pdf, width = 8, height = 7)
ggsave("Filtered_Range_N50_BUSCO.svg", plot_subset, device = "svg", width = 8, height = 7)

# =============================================================================
# 5. INTERACTIVE PLOT: Filtered Range
# =============================================================================

interactive_subset <- ggplotly(
  plot_subset, 
  tooltip = c("x", "y", "color")
) %>%
  style(text = paste0("<b>Species: ", filtered_subset$Current_name, "</b>",
                      "<br>N50: ", round(filtered_subset$Contig_N50_Mbp, 3), " Mbp",
                      "<br>BUSCO: ", filtered_subset$C, "%",
                      "<br>Accession: ", filtered_subset$Assembly_Accession))

saveWidget(interactive_subset, "Filtered_Range_Interactive.html", selfcontained = TRUE)

# =============================================================================
# End of script
# =============================================================================
































# =============================================================================
# Script:  Contig N50 vs BUSCO — Specific Focus: Hanseniaspora
# Author:  Yedomon Ange Bovys Zoclanclounon
#          Bioinformatics Scientist, Rothamsted Research
# Date:    2026-04-02
# Purpose: Focused visualization of assembly quality for the genus Hanseniaspora.
#          Generates static (PDF/SVG) and interactive (HTML) versions.
# =============================================================================

# ── 1. Working environment ────────────────────────────────────────────────────
remove(list = ls())
gc()
set.seed(1000)

# ── 2. Load required packages ─────────────────────────────────────────────────
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2, ggrepel, scales, plotly, htmlwidgets, svglite)

# ── 3. Load and Prepare Master Data ───────────────────────────────────────────
master_data <- read.csv("master_data_filtered.csv", h = TRUE)

# Filter specifically for Hanseniaspora
# We use grepl to catch any "Hanseniaspora" in the Current_name or Organism_Name
hansen_data <- master_data %>%
  mutate(
    Seq_type = ifelse(is.na(Seq_type) | Seq_type == "", "not_provided", Seq_type),
    Seq_type = factor(Seq_type, levels = c("not_provided", "long_read", "short_read")),
    Contig_N50_Mbp = Assembly_Stats_Contig_N50 / 1e6,
    label_italic = paste0("italic('", Current_name, "')")
  ) %>%
  filter(
    grepl("Hanseniaspora", Current_name, ignore.case = TRUE),
    Contig_N50_Mbp > 0,
    C >= 0
  )

# ── 4. Aesthetic Settings ─────────────────────────────────────────────────────
custom_colors <- c(
  long_read    = "#755cc5ff",
  short_read   = "#da57a1ff",
  not_provided = "#A9A9A9"
)

# =============================================================================
# 5. STATIC PLOT: Hanseniaspora Focus
# =============================================================================

plot_hansen <- ggplot(hansen_data, 
                      aes(x = Contig_N50_Mbp, y = C, color = Seq_type)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    aes(label = label_italic), 
    parse = TRUE, 
    size = 3.5, 
    max.overlaps = Inf, 
    box.padding = 0.5, 
    segment.color = "grey70"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(name = "Contig N50 (Mbp)", labels = label_number(accuracy = 0.001)) +
  scale_y_continuous(name = "BUSCO Complete (C, %)") +
  labs(
    title = "Assembly Quality: Hanseniaspora genus",
    subtitle = "Contig N50 vs. BUSCO Completeness",
    color = "Sequencing type"
  ) +
  theme_classic(base_size = 13)

# Save Static Versions
ggsave("Hanseniaspora_N50_vs_BUSCO.pdf", plot_hansen, device = cairo_pdf, width = 8, height = 7)
ggsave("Hanseniaspora_N50_vs_BUSCO.svg", plot_hansen, device = "svg", width = 8, height = 7)

# =============================================================================
# 6. INTERACTIVE PLOT: Hanseniaspora Focus
# =============================================================================

interactive_hansen <- ggplotly(
  plot_hansen, 
  tooltip = c("x", "y", "color")
) %>%
  # Adding the species name to the hover info manually since parse=TRUE 
  # can sometimes mess up simple tooltips in plotly
  style(text = paste0("Species: ", hansen_data$Current_name, 
                      "<br>N50: ", round(hansen_data$Contig_N50_Mbp, 4), " Mbp",
                      "<br>BUSCO: ", hansen_data$C, "%"))

# Save Interactive Version
saveWidget(
  widget = interactive_hansen, 
  file = "Hanseniaspora_N50_vs_BUSCO_interactive.html",
  selfcontained = TRUE
)

# =============================================================================
# End of script
# =============================================================================