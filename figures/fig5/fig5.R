# =============================================================================
# Script:  Class-level taxonomic tree of fungi with associated assembly counts
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

# ── MycoBank taxonomic framework ──────────────────────────────────────────────

# Read in MycoBank taxonomy (downloaded 20/10/2025) 
mb <- readxl::read_xlsx("MBList.xlsx")

# Check for contradictory taxonomies
mb.check <- mb %>%
  # Extract current genus and fill in missing rank
  mutate(Current_genus=str_extract(Synonymy,
                                   "(?<=Current name: )\\w+"),
         `Rank.Rank name`=ifelse(`Taxon name` == Current_genus &
                                   `Rank.Rank name` == "-",
                                 "gen.",
                                 `Rank.Rank name`)) %>%
  # Filter for legitimate genus-level fungal taxonomy
  filter(`Rank.Rank name` == "gen.",
         !`Name status` %in% c("Illegitimate", "Invalid"),
         !grepl("Protozoa|Chromista|Bacteria|Algae|Plantae|Fossil",
                Classification)) %>%
  # Count number of distinct taxonomies per current genus
  group_by(Current_genus) %>%
  mutate(num=n_distinct(Classification))

# Set aside genera with taxonomy consensus
mb.1 <- mb.check %>%
  filter(num == 1) %>%
  distinct(Classification, .keep_all=TRUE)

mb.check2 <- mb.check %>%
  # Filter for informative contradictory taxonomies
  filter(num > 1,
         !Classification %in% c("?", "-", "Fungi, Incertae sedis",
                                "Fungi, Fossil Fungi")) %>%
  distinct(Classification, .keep_all=TRUE) %>%
  arrange(Current_genus) %>%
  # Count number of taxonomic ranks assigned
  mutate(num_tax=str_count(Classification, "\\S+")) %>%
  select(Current_genus, Classification, num_tax, everything()) %>%
  # Filter for taxonomies where taxon name and current name are in agreement
  filter(if(any(`Current name.Taxon name` == first(Current_genus) &
                `Taxon name` == first(Current_genus))) {
    `Current name.Taxon name` == first(Current_genus) & 
      `Taxon name` == first(Current_genus)
  } else if (any(`Current name.Taxon name` == first(Current_genus))) {
    `Current name.Taxon name` == first(Current_genus)
  } else {
    TRUE
  })

# Set aside genera with taxonomy consensus
mb.2 <- mb.check2 %>%
  group_by(Current_genus) %>%
  filter(n_distinct(Classification) == 1)

# Filter genera with duplicate taxonomies, for manual checking
mb.check3 <- mb.check2 %>%
  group_by(Current_genus) %>%
  mutate(num=n()) %>%
  filter(n_distinct(Classification) > 1,
         # Where only two different taxonomies exist and the number of taxonomic
         # ranks is very different (diff >= 4), filter for the taxonomy with
         # more taxonomic ranks assigned
         if (first(num) == 2 & max(num_tax) - min(num_tax) >= 4) {
           num_tax == max(num_tax)
         } else {
           TRUE
         },
         `Name status` == "Legitimate",
         # Filter taxonomies ending on subfamily rank
         !grepl("ideae$", Classification)) %>%
  mutate(num=n())

# Read in researched fixes after manual appraisal
tax.fix <- read.csv("manual_curation.tsv", sep="\t")

# Set aside genera with taxonomy consensus
mb.3 <- tax.fix %>%
  bind_rows((mb.check3 %>%
               filter(num == 1))) %>%
  group_by(Current_genus) %>%
  distinct(Classification)

# Combine resolved taxonomies
tax <- bind_rows(mb.1, mb.2, mb.3) %>%
  select(Current_genus, Classification)

# Manually check genera with no classification
tax %>%
  filter(Classification %in% c("-", "?")) %>%
  print(n=50)

# Manual fixes:

# Fix missing Dikarya
tax$Classification[which(grepl("Fungi, Ascomycota", tax$Classification))] <-
  sub("Fungi, Ascomycota", "Fungi, Dikarya, Ascomycota",
      tax$Classification[which(grepl("Fungi, Ascomycota", tax$Classification))])

# Fix Microsporidia
tax$Classification[which(grepl("Microsporidiomycota", tax$Classification))] <-
  sub("Microsporidiomycota", "Microsporidia",
      tax$Classification[which(grepl("Microsporidiomycota", tax$Classification))])

# Fix Glugeida
tax$Classification[which(grepl("Glugeidae$", tax$Classification))] <-
  sub("Glugeidae", "Glugeida, Glugeidae",
      tax$Classification[which(grepl("Glugeidae$", tax$Classification))])

tax$Classification[which(grepl("Glugeidae,", tax$Classification))] <-
  sub("Glugeidae,", "Glugeida,",
      tax$Classification[which(grepl("Glugeidae,", tax$Classification))])

# Fix Malasseziomycetes
tax$Classification[which(grepl("Malasseziales", tax$Classification))] <-
  sub("Exobasidiomycetes,", "Malasseziomycetes,",
      tax$Classification[which(grepl("Malasseziales,", tax$Classification))])

# Fix 'Ascomycetes'
tax$Classification[
  which(grepl("Fungi, Dikarya, Ascomycota, Ascomycetes, Roesleriaceae",
              tax$Classification))
] <-
  "Fungi, Dikarya, Ascomycota, Pezizomycotina, Leotiomycetes, Leotiomycetidae,
Helotiales, Helotiaceae"

tax$Classification[
  which(grepl("Fungi, Dikarya, Ascomycota, Ascomycetes, Pseudeurotiaceae",
              tax$Classification))
] <-
  "Fungi, Dikarya, Ascomycota, Pezizomycotina, Leotiomycetes, Leotiomycetidae,
Thelebolales, Pseudeurotiaceae"

tax$Classification[
  which(grepl("Fungi, Dikarya, Ascomycota, Ascomycetes, Myxotrichaceae",
              tax$Classification))
] <-
  "Fungi, Dikarya, Ascomycota, Pezizomycotina, Leotiomycetes, Leotiomycetidae,
Myxotrichaceae"

tax$Classification[
  which(grepl("Fungi, Dikarya, Ascomycota, Ascomycetes", tax$Classification))
] <-
  "Fungi, Dikarya, Ascomycota"

# Update incertae sedis or ambiguous 'Ascomycota' genera
tax$Classification[tax$Current_genus == "Microcyclospora"] <-
  "Fungi, Dikarya, Ascomycota, Pezizomycotina, Dothideomycetes,
Dothideomycetidae"
tax$Classification[tax$Current_genus == "Microcyclosporella"] <-
  "Fungi, Dikarya, Ascomycota, Pezizomycotina, Dothideomycetes,
Dothideomycetidae"
tax$Classification[tax$Current_genus == "Emarellia"] <-
  "Fungi, Dikarya, Ascomycota, Pezizomycotina, Dothideomycetes,
Pleosporomycetidae, Pleosporales, Trematosphaeriaceae"
tax$Classification[tax$Current_genus == "Basidioascus"] <-
  "Fungi, Dikarya, Basidiomycota, Wallemiomycotina, Wallemiomycetes,
Geminibasidiales, Geminibasidiaceae"
tax$Classification[tax$Current_genus == "Arachnomyces"] <-
  "Fungi, Dikarya, Ascomycota, Pezizomycotina, Eurotiomycetes,
Eurotiomycetidae, Arachnomycetales, Arachnomycetaceae"
tax$Classification[tax$Current_genus == "Pyrenophora"] <-
  "Fungi, Dikarya, Ascomycota, Pezizomycotina, Dothideomycetes,
Pleosporomycetidae, Pleosporales, Pleosporaceae"

# Add missing genus
tax <- tax %>%
  ungroup() %>%
  add_row(Current_genus="Cochliobolus",
          Classification="Fungi, Dikarya, Ascomycota, Pezizomycotina,
          Dothideomycetes, Pleosporomycetidae, Pleosporales, Pleosporaceae")

# List remaining genera with no classification
tax %>%
  filter(Classification %in% c("-", "?")) %>%
  print(n=50)

## Split taxonomic ranks ##

subkingdoms <- c("Basidiobolomyceta", "Chytridiomyceta", "Dikarya",
                 "Mucoromyceta", "Olpidiomyceta", "Rozellomyceta",
                 "Zoopagomyceta")

suffixes <- c("myceta", "mycota", "mycotina", "mycetes", "mycetidae", "les",
              "ceae", "idae", "ida")

# Separate ranks into own columns based on suffixes
tax.2 <- tax %>%
  separate_rows(Classification, sep=",") %>%
  # Trim whitespace
  mutate(Classification=str_trim(Classification)) %>%
  # Identify suffix
  mutate(suffix=str_extract(
    Classification, paste0("(", paste0(suffixes, collapse="|"), ")$"))
  ) %>%
  # Assign column names based on suffix
  mutate(col=case_when(
    Classification == "Fungi" ~ "Current_kingdom",
    Classification %in% subkingdoms ~ "Current_subkingdom",
    suffix == "mycota" | Classification == "Microsporidia" ~ "Current_phylum",
    suffix == "mycotina" ~ "Current_subphylum",
    suffix == "mycetes" | Classification %in%
      c("Microsporea", "Rudimicrosporea") ~ "Current_class",
    suffix == "mycetidae" ~ "Current_subclass",
    suffix %in% c("les", "ida") ~ "Current_order",
    suffix %in% c("ceae", "idae") ~ "Current_family"
  )) %>%
  select(-suffix) %>%
  filter(!is.na(col)) %>%
  group_by(Current_genus) %>%
  pivot_wider(names_from=col, values_from=Classification) %>%
  ungroup() %>%
  # Add rank to incertae sedis to make unique
  mutate(across(everything(), ~replace_na(.x, "Incertae sedis")),
         Current_class=ifelse(Current_class == "Incertae sedis",
                              paste(Current_class, Current_phylum),
                              Current_class),
         Current_order=ifelse(Current_order == "Incertae sedis",
                              paste(Current_order, Current_class),
                              Current_order))

# Repeat for unprocessed MycoBank dataframe
all.tax <- mb.check %>%
  select(`Taxon name`, Current_genus) %>%
  left_join(tax, by="Current_genus") %>%
  distinct(.keep_all=TRUE) %>%
  separate_rows(Classification, sep=",") %>%
  # Trim whitespace
  mutate(Classification=str_trim(Classification)) %>%
  # Identify suffix category
  mutate(suffix=
           str_extract(Classification,
                       paste0("(", paste0(suffixes, collapse="|"), ")$"))) %>%
  # Assign column names based on suffix
  mutate(col=case_when(
    Classification == "Fungi" ~ "Current_kingdom",
    Classification %in% subkingdoms ~ "Current_subkingdom",
    suffix == "mycota" | Classification == "Microsporidia" ~ "Current_phylum",
    suffix == "mycotina" ~ "Current_subphylum",
    suffix == "mycetes" | Classification %in%
      c("Microsporea", "Rudimicrosporea") ~ "Current_class",
    suffix == "mycetidae" ~ "Current_subclass",
    suffix %in% c("les", "ida") ~ "Current_order",
    suffix %in% c("ceae", "idae") ~ "Current_family"
  )) %>%
  select(-suffix) %>%
  filter(!is.na(col)) %>%
  group_by(Current_genus, `Taxon name`) %>%
  pivot_wider(names_from=col, values_from=Classification) %>%
  ungroup() %>%
  # Add rank to incertae sedis to make unique
  mutate(across(everything(), ~replace_na(.x, "Incertae sedis")),
         Current_class=ifelse(Current_class == "Incertae sedis",
                              paste(Current_class, Current_phylum),
                              Current_class),
         Current_order=ifelse(Current_order == "Incertae sedis",
                              paste(Current_order, Current_class),
                              Current_order))

# Write taxonomy to file
write.csv(tax.2,
          paste0("MycoBank_processed_taxonomy-", Sys.Date(), ".csv"),
          row.names=FALSE, quote=FALSE)

# ── Map taxonomy to genome data ───────────────────────────────────────────────

df <- read.csv("NCBI_and_JGI_genomes_shared_metadata_260203.csv") %>%
  # Clear any existing taxonomy/mycobank columns
  select(-c(contains("Current"), "mycobank"))

# Read in flagged accessions to be excluded from plots
toremove <- read.csv("Assembly_accession_to_remove_N=72.csv")

df.2 <- df %>%
  # Exclude accessions
  filter(!Assembly_Accession %in% toremove$Assembly_Accession) %>%
  # Fix names/formatting
  mutate(
    Organism_Name=sub("\\[Candida\\]", "Candida", Organism_Name),
    Organism_Name=sub("$uncultured ", "", Organism_Name),
    Organism_Name=gsub(r"{\s*\([^\)]+\)}","", Organism_Name),
    Organism_Name=sub("Rhodosporidium", "Rhodotorula", Organism_Name),
    Organism_Name=sub("Porosterum", "Porostereum", Organism_Name),
    Organism_Name=sub("Setosphaeria turcica", "Exserohilum turcicum", Organism_Name),
    Database_name=str_extract(
      Organism_Name,
      "([A-Z][a-z]+\\s[a-z]+(?:\\s(?:subsp\\.|var\\.|aff\\.|f\\.|forma|sp\\.)\\s[a-z]+)*)"
    ),
    Database_genus=word(Database_name, 1)
  ) %>%
  # Add MycoBank taxonomy
  left_join(tax.2, by=c("Database_genus"="Current_genus"))

# For unmatched genera, check for synonyms to recover taxonomy
for (i in 1:nrow(df.2)) {
  
  for (level in c("family", "order", "subclass", "class",
                  "subphylum", "phylum", "subkingdom", "kingdom")) {
    
    if (is.na(df.2[i, paste0("Current_", level)]) &
        df.2$Genus[i] %in% all.tax$`Taxon name`) {
      
      df.2[i, paste0("Current_", level)] <- 
        unique(
          all.tax[intersect(grep(df.2$Genus[i], all.tax$`Taxon name`),
                            grep(df.2$Phylum[i], all.tax$Current_phylum)),
                  paste0("Current_", level)]
        )
      
    }
    
  }
  
}

# Use existing NCBI taxonomy if unable to match to MycoBank
df.3 <- df.2 %>%
  mutate(mycobank=ifelse(is.na(Current_kingdom),
                         FALSE,
                         TRUE),
         Current_kingdom=ifelse(is.na(Current_kingdom),
                                Kingdom,
                                Current_kingdom),
         Current_phylum=ifelse(is.na(Current_phylum),
                               Phylum,
                               Current_phylum),
         Current_class=ifelse(is.na(Current_class),
                              Class,
                              Current_class),
         Current_order=ifelse(is.na(Current_order),
                              Order,
                              Current_order),
         Current_family=ifelse(is.na(Current_family),
                               Family,
                               Current_family))

# Manually check cases where the old and updated phyla differ
df.check <- df.2 %>% 
  filter(Phylum != Current_phylum & Phylum != "") %>%
  select(-c(C, S, D, F, M, n, E, contains("Assembly")))

# Write data to file
write.csv(df.3,
          paste0("NCBI_and_JGI_genomes_updated_taxonomy-", Sys.Date(), ".csv"),
          row.names=FALSE, quote=FALSE)

# ── Plot taxonomy tree ────────────────────────────────────────────────────────

# Format to class-level
tax.3 <- tax.2 %>%
  mutate(across(everything(), as_factor)) %>%
  select(Current_phylum, Current_class) %>%
  unique()

# Generate tree object
tree <- as.phylo(~Current_phylum/Current_class, data=tax.3)

# Generate dataframe for tip labels
labels.df <- data.frame(class=tree$tip.label, Label=tree$tip.label)
labels.df$Label[grep("Incertae", labels.df$Label)] <- "Incertae sedis"

# Plot base tree
gg.tree <- ggtree(tree, layout="circular", linewidth=NA) %<+% labels.df

# Generate dataframe for nodes
phyla.nodes <- data.frame(
  phylum=unique(
    tax.3$Current_phylum[match(
      # Extract nodes from tree object
      gg.tree$data %>%
        filter(isTip) %>%
        arrange(y) %>%
        pull(label),
      tax.3$Current_class
    )]
  ),
  node=NA,
  # Generate field for alternating highlights
  box=rep(c(0, 1),
          # Total number of nodes
          length.out=length(
            unique(
              tax.3$Current_phylum[match(
                gg.tree$data %>%
                  filter(isTip) %>%
                  arrange(y) %>%
                  pull(label),
                tax.3$Current_class)]
            )
          )
  )
)

# Find most recent common ancestor for each phylum
for (i in 1:length(phyla.nodes$phylum)) {
  
  phyla.nodes$node[i] <-
    MRCA(tree, tax.3$Current_class[tax.3$Current_phylum == phyla.nodes$phylum[i]])
  
}

# Add to tree plot
gg.tree.2 <- gg.tree +
  xlim(-5, 20) +
  # Alternating highlights for phyla
  geom_highlight(data=phyla.nodes,
                 aes(node=node, fill=as.factor(box)),
                 extend=17,
                 alpha=0.2,
                 show.legend=FALSE) +
  scale_fill_manual(values=c("#A3A3A3", "white")) +
  # Lines for phyla
  geom_cladelab(data=phyla.nodes,
                mapping=aes(node=node, label=phylum),
                angle="auto",
                fontsize=0,
                fontface="italic",
                barsize=0.3,
                offset=17,
                offset.text=0.05) +
  # Tree branches
  geom_tree(linewidth=0.3) +
  # Tip labels
  geom_tiplab(aes(label=Label),
              size=2.8,
              offset=4,
              fontface="italic",
              align=TRUE,
              linesize=0.2) +
  ggpreview(width=7, height=7)

# Summarise number of long- and short-read assemblies per class
df.counts <- df.3 %>%
  group_by(Current_class, Seq_type) %>%
  summarise(num=n())

# Dummy plot to extract legend
gg.legend <- gg.tree.2 +
  new_scale_fill() +
  # Add assembly counts
  geom_fruit(data=df.counts,
             geom=geom_point,
             mapping=aes(y=Current_class, x=Seq_type, fill=Seq_type, size=num),
             shape=21,
             colour="black",
             offset=0.05,
             pwidth=0.15,
             alpha=0.5) +
  scale_fill_manual(values=c("#755cc5ff", "#da57a1ff"),
                    labels=c("long-read", "short-read")) +
  scale_size_continuous(trans="sqrt",
                        range=c(0.05, 5),
                        breaks=c(1, 10, 100, 1000),
                        labels=comma) +
  guides(size=guide_legend(title=NULL,
                           label.position="bottom",
                           direction="horizontal"),
         fill=guide_legend(title="# genome\nassemblies",
                           title.position="top",
                           title.hjust=0.5,
                           direction="vertical",
                           order=1)) +
  theme(legend.direction="horizontal",
        legend.margin=margin(0, 0, 0, 0),
        legend.text=element_text(size=8, margin=margin(0, 0, 0, 0)),
        legend.title=element_text(size=8, face="bold"),
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
             mapping=aes(y=Current_class, x=Seq_type, fill=Seq_type, size=num),
             shape=21,
             colour="black",
             offset=0.1,
             pwidth=0.3,
             alpha=0.5,
             show.legend=FALSE) +
  scale_fill_manual(values=c("#755cc5ff", "#da57a1ff"),
                    labels=c("long-read", "short-read")) +
  # Add boxplot ring with BUSCO completeness
  geom_fruit(data=df.2,
             geom=geom_boxplot,
             mapping=aes(y=Current_class, x=C),
             outliers=FALSE,
             linewidth=0.2,
             offset=0.5,
             pwidth=0.8,
             grid.params=list(linewidth=0.1, color="grey")) +
  scale_size_continuous(trans="sqrt",
                        range=c(0.05, 5),
                        breaks=c(1, 10, 100, 1000),
                        labels=comma) +
  theme(panel.background=element_blank()) +
  # Insert legend in center
  inset_element(size.legend.grob,
                l=0.44, b=0.445, r=0.56, t=0.555, on_top=FALSE) +
  ggpreview(width=7, height=7)

# Write to file
pdf(paste0("fig_5-", Sys.Date(), ".pdf"), width=7, height=7)
gg.tree.3
dev.off()

# ── Assemblies per rank inset barplot ─────────────────────────────────────────

# Generate dataframe with percentage of each rank with a genome
rank.df <- data.frame(
  rank=c("Class", "Order", "Family", "Genus"),
  prop=c(
    length(unique(df.3$Current_class)) / length(unique(tax.2$Current_class)),
    length(unique(df.3$Current_order)) / length(unique(tax.2$Current_order)),
    length(unique(df.3$Current_family)) / length(unique(tax.2$Current_family)),
    length(unique(df.3$Genus)) / length(unique(tax.2$Current_genus))
  )
)

rank.df$rank <- factor(rank.df$rank,
                       levels=c("Class", "Order", "Family", "Genus"))

# Plot
gg.rank <- ggplot(rank.df, aes(x=rank, y=prop)) +
  geom_col(width=0.7) +
  geom_text(aes(label=round(prop*100)),
            size=2.8,
            vjust=-1) +
  scale_y_continuous(labels=percent, limits=c(0, 1)) +
  labs(y="Percentage with\n≥1 assembly", x=NULL) +
  theme_minimal() +
  theme(axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=8, angle=45, hjust=1, vjust=1.2),
        axis.title=element_text(size=8)) +
  ggpreview(width=1.5, height=1.5)

# Write to file
cairo_pdf(paste0("fig_5_inset-", Sys.Date(), ".pdf"), width=1.5, height=1.5)
gg.rank
dev.off()

# ──────────────────────────────────────────────────────────────────────────────