# =============================================================================
# Script:  Proportion of NCBI fungal assemblies with annotations
# Author:  Rowena Hill
#          Earlham Institute
# Date:    2026-04-27
# =============================================================================

library(tidyverse)    # v2.0.0
library(tgutil)       # v0.1.20

# ── Annotation barplot ────────────────────────────────────────────────────────

# Read in NCBI data
ncbi.df <- read.csv("NCBI_genomes_with_mycobank_taxonomy_seqtype-2025-10-27.csv")

ncbi.df <- ncbi.df %>%
  # Format date
  mutate(year=year(as.Date(Assembly_Release_Date, format="%d/%m/%Y")),
         # Detect annotation
         annotated=!is.na(Annotation_Count_Gene_Total))

# Plot barplot
gg.annotated <- ggplot(ncbi.df %>%
                         filter(!is.na(year)),
                       aes(x=factor(year), fill=annotated)) +
  # Separate long- and short-read assemblies
  facet_grid(.~Seq_type, space="free", scales="free",
             labeller=labeller(
               Seq_type=c(long_read="Long-read", short_read="Short-read"))
  ) +
  geom_bar(position="fill") +
  scale_fill_manual(values=c("dimgrey", "lightgrey"),
                    breaks=c("TRUE", "FALSE")) +
  scale_y_continuous(labels=scales::comma,
                     expand=c(0, 0)) +
  scale_x_discrete(expand=c(0, 0)) +
  labs(x=NULL, y="Proportion of\ngenomes in NCBI", fill="Annotated:") +
  theme_minimal() +
  theme(legend.position="top",
        legend.title=element_text(size=8, face="bold"),
        legend.text=element_text(size=7, margin=margin(0, 0, 0, 1)),
        legend.key.size=unit(8, "pt"),
        legend.key.spacing.x=unit(1.5, "pt"),
        legend.margin=margin(0, 0, -5, 0),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks.x=element_line(linewidth=0.3),
        axis.text.x=element_text(size=5, angle=90, vjust=0.5, hjust=4)) +
  ggpreview(width=3, height=2)

# Write to file
pdf(paste0("supp_fig_2-", Sys.Date(), ".pdf"), width=3, height=2)
gg.annotated
dev.off()

# ── Stats for text ────────────────────────────────────────────────────────────

# Proportion of assemblies with no annotation
ncbi.df %>%
  group_by(annotated) %>%
  summarise(prop=n()/nrow(ncbi.df))

# Proportion of species that have no annotation
ncbi.df %>%
  group_by(Current_name, annotated) %>%
  summarise(n()) %>%
  group_by(Current_name) %>%
  filter(all(annotated == FALSE)) %>%
  ungroup() %>%
  summarise(prop=n()/length(unique(ncbi.df$Current_name)))

# ──────────────────────────────────────────────────────────────────────────────