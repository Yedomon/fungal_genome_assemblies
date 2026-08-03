# =============================================================================
# Script:  Order-level taxonomic tree of fungi with associated assembly counts
#          and completeness
# Author:  Rowena Hill
#          Earlham Institute
# Date:    2026-04-27
# =============================================================================

library(tidyverse)    # v2.0.0
library(ape)          # v5.8-1
library(ggnewscale)   # v0.5.2
library(ggtree)       # v4.0.5
library(ggtreeExtra)  # v1.20.1
library(patchwork)    # v1.3.2
library(scales)       # v1.4.0
library(tgutil)       # v0.1.20

# ── Read in data ──────────────────────────────────────────────────────────────

tax <- read.csv("MycoBank_processed_taxonomy-2026-04-27.csv")
df <- read.csv("NCBI_and_JGI_genomes_updated_taxonomy-2026-04-27.csv")

# ── Plot taxonomy tree ────────────────────────────────────────────────────────

# Format to order-level
tax.2 <- tax %>%
  mutate(across(everything(), as_factor)) %>%
  select(Current_phylum, Current_class, Current_order) %>%
  unique()

# Generate tree object
tree <- as.phylo(~Current_phylum/Current_class/Current_order, data=tax.2)

# Generate dataframe for tip labels
labels.df <- data.frame(order=tree$tip.label, Label=tree$tip.label)
labels.df$Label[grep("Incertae", labels.df$Label)] <- "Incertae sedis"

# Plot base tree
gg.tree <- ggtree(tree, layout="circular", linewidth=NA) %<+% labels.df

# Generate dataframe for nodes (phylum-level)
phyla.nodes <- data.frame(
  phylum=unique(
    tax.2$Current_phylum[match(
      # Extract nodes from tree object
      gg.tree$data %>%
        filter(isTip) %>%
        arrange(y) %>%
        pull(label),
      tax.2$Current_order
    )]
  ),
  node=NA,
  # Generate field for alternating highlights
  box=rep(c(0, 1),
          # Total number of nodes
          length.out=length(
            unique(
              tax.2$Current_phylum[match(
                gg.tree$data %>%
                  filter(isTip) %>%
                  arrange(y) %>%
                  pull(label),
                tax.2$Current_order)]
            )
          )
  )
)

# Find most recent common ancestor for each phylum
for (i in 1:length(phyla.nodes$phylum)) {
  
  phyla.nodes$node[i] <-
    MRCA(tree, tax.2$Current_order[tax.2$Current_phylum == phyla.nodes$phylum[i]])
  
}

# Generate dataframe for nodes (class-level)
class.nodes <- data.frame(
  class=
    unique(
      tax.2$Current_class[match(
        gg.tree$data %>%
          arrange(y) %>%
          filter(isTip == "TRUE") %>%
          pull(label),
        tax.2$Current_order
      )]
    ),
  node=NA,
  box=0
)

# Add binary for alternating highlights
class.nodes$box[
  class.nodes$class %in% unique(
    tax.2$Current_class[
      tax.2$Current_phylum %in%
        names(which(
          table(unique(tax.2[
            c("Current_phylum", "Current_class")
          ])$Current_phylum) > 1
        ))
    ]
  )
] <-
  rep(c(0,1),
      length.out=
        length(
          class.nodes$box[
            class.nodes$class %in% unique(
              tax.2$Current_class[
                tax.2$Current_phylum %in%
                  names(which(
                    table(unique(tax.2[
                      c("Current_phylum", "Current_class")
                    ])$Current_phylum) > 1
                  ))
              ]
            )
          ]
        )
  )

# Find most recent common ancestor for each class
for (i in 1:length(class.nodes$class)) {
  
  class.nodes$node[i] <-
    MRCA(tree, tax.2$Current_order[tax.2$Current_class == class.nodes$class[i]])
  
}

# Abbreviate incertae sedis
class.nodes$class <- sub("Incertae sedis", "Inc. sed.", class.nodes$class)
class.nodes$class <- sub("Inc. sed. Incertae sedis", "Inc. sed.", class.nodes$class)

# Add to tree plot
gg.tree.2 <- gg.tree +
  xlim(-2, 10) +
  # Alternating highlights for classes
  geom_highlight(data=class.nodes,
                 aes(node=node, fill=as.factor(box)),
                 alpha=0.2,
                 extend=4.2,
                 show.legend=FALSE) +
  scale_fill_manual(values=c("white", "#BABABA")) +
  new_scale_fill() +
  # Alternating highlights for phyla
  geom_highlight(data=phyla.nodes,
                 aes(node=node, fill=as.factor(box)),
                 extend=6.9,
                 alpha=0.2,
                 show.legend=FALSE) +
  scale_fill_manual(values=c("#A3A3A3", "white")) +
  # Tree branches
  geom_tree(size=0.3) +
  # Tip labels
  geom_tiplab(aes(label=Label),
              size=0.8,
              offset=2.2,
              fontface="italic",
              align=TRUE,
              linesize=0.2) +
  # Labels for classes
  geom_cladelab(data=class.nodes,
                mapping=aes(node=node, label=class),
                angle="auto",
                horizontal=TRUE,
                fontsize=1,
                fontface="italic",
                barsize=0.2,
                #offset=4,
                offset=4.2,
                offset.text=0.05) +
  # Lines for phyla
  geom_cladelab(data=phyla.nodes,
                mapping=aes(node=node, label=phylum),
                angle="auto",
                fontsize=0,
                fontface="italic",
                barsize=0.3,
                offset=6.9,
                offset.text=0.05) +
  ggpreview(width=7, height=7)

# Summarise number of long- and short-read assemblies per class
df.counts <- df %>%
  group_by(Current_order, Seq_type) %>%
  summarise(num=n())

# Dummy plot to extract legend
gg.legend <- gg.tree.2 +
  new_scale_fill() +
  # Add assembly counts
  geom_fruit(data=df.counts,
             geom=geom_point,
             mapping=aes(y=Current_order, x=Seq_type, fill=Seq_type, size=num),
             shape=21,
             colour="black",
             offset=0.03,
             pwidth=0.11,
             alpha=0.5) +
  scale_fill_manual(values=c("#755cc5ff", "#da57a1ff"),
                    labels=c("long-read", "short-read")) +
  scale_size_continuous(trans="sqrt",
                        range=c(0.05, 5),
                        breaks=c(1, 10, 100, 1000),
                        labels=comma) +
  guides(size=guide_legend(title=NULL,
                           label.position = "bottom",
                           direction="horizontal"),
         fill=guide_legend(title="# genome\nassemblies",
                           title.position="top",
                           title.hjust=0.5,
                           direction="vertical",
                           order=1)) +
  theme(legend.direction="horizontal",
        legend.margin=margin(0, 0, 0, 0),
        legend.text=element_text(size=6, margin=margin(0, 0, 0, 0)),
        legend.title=element_text(size=6, face="bold"),
        legend.key.size=unit(8, "pt"),
        legend.key.spacing.x=unit(1.5, "pt"),
        legend.box.just="center",
        legend.spacing.y=unit(2, "pt"),
        legend.justification="center")

# Extract legend
size.legend.grob <- cowplot::get_legend(gg.legend)

# Add genome data to tree plot
gg.tree.3 <- gg.tree.2 +
  new_scale_fill() +
  # Add ring with number of assemblies
  geom_fruit(data=df.counts,
             geom=geom_point,
             mapping=aes(y=Current_order, x=Seq_type, fill=Seq_type, size=num),
             shape=21,
             colour="black",
             offset=0.03,
             pwidth=0.11,
             alpha=0.5,
             show.legend=FALSE) +
  # Add blank ring for boxplot background
  geom_fruit(data=tax.2,
             geom=geom_tile,
             mapping=aes(y=Current_order, x=1),
             fill="white",
             colour=NA,
             alpha=0.8,
             offset=0.05,
             pwidth=0.3,
             show.legend=FALSE) +
  # Add boxplot ring with BUSCO completeness
  geom_fruit(data=df,
             geom=geom_boxplot,
             mapping=aes(y=Current_order, x=C),
             outliers=FALSE,
             linewidth=0.1,
             offset=-0.15,
             pwidth=0.3,
             axis.params=list(
               axis="x",
               text.size=0.8,
               vjust=1.2
             ),
             grid.params=list(linewidth=0.1, color="grey")) +
  scale_fill_manual(values=c("#755cc5ff", "#da57a1ff"),
                    labels=c("long-read", "short-read")) +
  scale_size_continuous(trans="sqrt",
                        range=c(0.05, 5),
                        breaks=c(1, 10, 100, 1000),
                        labels=comma) +
  theme(panel.background=element_blank()) +
  # Insert legend in center
  inset_element(size.legend.grob, l=0.44, b=0.44, r=0.56, t=0.55, on_top=FALSE) +
  ggpreview(width=7, height=7)

# Write to file
pdf(paste0("supp_fig_1-", Sys.Date(), ".pdf"), width=7, height=7)
gg.tree.3
dev.off()

# ──────────────────────────────────────────────────────────────────────────────