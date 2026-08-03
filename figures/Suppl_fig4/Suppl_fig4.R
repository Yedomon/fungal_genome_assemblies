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
library(packcircles)  # v0.3.7
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


# ── Plot genus bubble plot ────────────────────────────────────────────────────

# Calculate total number of genera per phylum
tmp <- tax.2 %>%
  group_by(Current_phylum) %>%
  summarise(num.gen=n_distinct(Current_genus))

# Calculate number of genera with a genome assembly
df.gen <- df.3 %>%
  filter(Genus != "") %>%
  group_by(Current_phylum) %>%
  summarise(num.gen.seq=n_distinct(Genus)) %>%
  # Combine with total number of genera
  left_join(tmp) %>%
  mutate(perc=round(num.gen.seq/num.gen*100), # calculate percentage
         # add colours for phyla
         colour=c("#0072b2", "#c8cdda", "#de6765", "#bf987c", "#f8e6d6",
                  "#f4f1b3", "#755cc5ff", "#ba7bb4", "#5699d1", "#da57a1ff",
                  "#4abcbd", "#ebc74f", "#878365", "#cc9f28", "#60b345",
                  "#eda4a6", "#006779", "#657719", "#a0ca78", "#f59943"))

# Get radius and x and y coordinates for centre of larger circles
circle.layout <- circleProgressiveLayout(df.gen$num.gen,
                                         sizetype="area")

# Add a small gap between circles
circle.layout$radius <- circle.layout$radius * 0.95

# Create a dataframe of vertices to draw each 'circle'
circle.vertices <- circleLayoutVertices(circle.layout, npoints=50)

# Get radius and x and y coordinates for centre of nested circles
circle.layout.pub <- 
  circleProgressiveLayout(df.gen$num.gen.seq,
                          sizetype="area")

# Again add small gap between circles
circle.layout.pub$radius <- circle.layout.pub$radius * 0.95

# Replace x and y with that of the larger circles, but keep the same radius
circle.layout.pub <- data.frame(x=circle.layout$x,
                                y=circle.layout$y,
                                radius=circle.layout.pub$radius)

# Create a dataframe of vertices to draw each nested 'circle'
circle.vertices.pub <- circleLayoutVertices(circle.layout.pub,
                                            npoints=50)

# Combine original dataframe with the layout dataframe
circle.labels <- cbind(df.gen, circle.layout)

#Plot bubble plot
gg.circles.nested <- ggplot() +
  # Add circles for total genera
  geom_polygon(data=circle.vertices,
               aes(x, y, group=id, fill=as.factor(id)),
               colour=NA,
               alpha=0.3) +
  # Add circles for genomes
  geom_polygon(data=circle.vertices.pub,
               aes(x, y, group=id, fill=as.factor(id)),
               colour=NA) +
  # Phylum name slightly above center
  geom_text(data=circle.labels,
            aes(x, y + 1.5, size=num.gen.seq, label=Current_phylum),
            fontface="bold",
            show.legend=FALSE) +
  # Percentage slightly below center
  geom_text(data=circle.labels,
            aes(x, y - 1.5, size=num.gen.seq, label=paste0(perc, "%")),
            fontface="bold",
            show.legend=FALSE) +
  scale_size_continuous(range=c(1.5, 3.5)) +
  scale_fill_manual(values=circle.labels$colour,
                    labels=circle.labels$Current_phylum) +
  guides(fill=guide_legend(
    nrow=3,
    direction="horizontal",
    title=NULL,
    label.theme=element_text(size=7, margin=margin(l=-3)),
    keywidth=unit(7, "pt"),
    keyheight=unit(7, "pt"))
  ) +
  coord_equal() +
  theme_void() + 
  theme(legend.position="none") +
  ggpreview(width=10, height=8, unit="in")

# Write to file
pdf(paste0("fig_5c-", Sys.Date(), ".pdf"), width=10, height=8)
gg.circles.nested
dev.off()

# ──────────────────────────────────────────────────────────────────────────────