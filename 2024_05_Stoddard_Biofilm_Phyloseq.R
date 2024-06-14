### PHYLOSEQ ANALYSIS

## Load Libraries
library(phyloseq)
library(tidyverse)
library(csv)
library(ggplot2)
library(extrafont)
library(dunn.test)
library(ape)
library(vegan)


### Section 1: Loading Data

## Import Sample Names and Data from First Sequencing Run
sample_names <- read.csv("/data/home/htstodda/Biofilm_Fastq/workpath/Biofilm_sample_names.csv")
sdata <- as.csv("/data/home/htstodda/Biofilm_Fastq/Biofilm_Code/Biofilm_Scripts/Final Files/2024_05_Stoddard_Biofilm_Metadata/Biofilm_Metadata.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)

## Import Sample Names and Data from Second Sequencing Run
sample_namesRR <- read.csv("/data/home/htstodda/Biofilm_Fastq/workpath/Biofilm_sample_namesRR.csv")
sdataRR <- as.csv("/data/home/htstodda/Biofilm_Fastq/Biofilm_Code/Biofilm_Scripts/Final Files/2024_05_Stoddard_Biofilm_Metadata/Biofilm_MetadataRR.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)

## Combine Sample Names and Data from the Two Runs
sample_namesFINAL <- rbind(sample_names, sample_namesRR)
sample.names <- as.character(sample_namesFINAL)
sdataFINAL <- rbind(sdata, sdataRR)

## Mutate sdata to Have Factors
sdata2FINAL <- sdataFINAL %>%
  mutate(sample.names = Sample_ID,
         Sample_Type = factor(Sample_Type, levels = c("Planktonic", "Film")),
         Sample_Conditions = factor(Sample_Conditions, levels = c("Seed Culture", "Control", "Increased Agitation", "Decreased Temperature", "Increased Airflow", "Decreased Temperature and Increased Airflow", "pH 8")),
         Sample_Temperature = factor(Sample_Temperature, levels = c("30C", "40C")),
         Sample_Agitation_Speed = factor(Sample_Agitation_Speed, levels = c("100rpm", "750rpm")),
         Sample_Air_Flowrate = factor(Sample_Air_Flowrate, levels = c("0sL/h", "10sL/h"))) %>%
  full_join(sample_namesFINAL, by = "sample.names") %>%
  unique() %>%
  rownames_to_column(var = "row_number") %>%
  column_to_rownames(var = "sample.names") %>%
  filter(!is.na(Sample_Type))
#View(sdata2Final)

## Check for Duplicated Sample_IDs
duplicated(sdata2FINAL$Sample_ID)

## Load Sequence Table
seqtabFINAL <- readRDS("/data/home/htstodda/Biofilm_Fastq/workpath/seqtabFINAL.rds")
colnames(seqtabFINAL) <- NULL

## Load Taxa Table
taxaFINAL <- readRDS("/data/home/htstodda/Biofilm_Fastq/workpath/taxtabFINAL.rds")    
taxaFINAL[2:5,2:5]


### Section 2: Making Phyloseq Objects

## Make the First Phyloseq Object
samdataFINAL = sample_data(sdata2FINAL) # Sample Data
seqtabFINAL = otu_table(seqtabFINAL, taxa_are_rows = FALSE) # Sequence Data
taxtabFINAL = tax_table(taxaFINAL) # Taxa Data

#View which samdatFINAL to include
view(samdataFINAL)
samdataFINAL = samdataFINAL[c(1:12,17:44,55), ] # Remove Sequencing Controls Sequencing Controls and pH 8 Data From Sample Data
levels(samdataFINAL$Sample_Conditions)[match("Seed Culture", levels(samdataFINAL$Sample_Conditions))] <- "Inoculation Culture"
seqtabFINAL = seqtabFINAL[c(1:12,17:44,55), ] # Remove Sequencing Controls and pH 8 Data From Sample Data

psFINAL = phyloseq(otu_table(seqtabFINAL), tax_table(taxtabFINAL), sample_data(samdataFINAL)) # Combine the Data into one Object
#View(psFINAL)

write_rds(psFINAL, "/data/home/htstodda/Biofilm_Fastq/workpath/psFINAL_unrarefied_Biofilm.rds") #Save Unrarefied Phyloseq Data

## Filter our Eukaryotes, Mitochondria, and Chloroplast
ps_filtFINAL <- psFINAL %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
#View(ps_filtFINAL)


### Section 3: Normalize the Data

samplesover1000_all <- subset_samples(ps_filtFINAL, sample_sums(ps_filtFINAL) > 1000)
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)

## For Reproducible Data
set.seed(81)

rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))
ppsFINAL <- rarefy_samplesover1000_all 
#View(ppsFINAL)


### Section 4: Plotting

## Box Plot Comparing Alpha Diveristy Metrics "Observed" and "Shannon" Between Biofilm and Planktonic Samples
ppsFINAL %>% # Phyloseq Object
  plot_richness(
    x = "Sample_Type", # Compare Diveristy of this Data Type
    measures = c("Observed", "Shannon")) + # Specify Diversity Measures
  geom_boxplot(aes(x = Sample_Type, fill = Sample_Type), show.legend = FALSE, outlier.shape = NA)+ # Make a Boxplot, Set aes(fill to Data Type)
  theme_linedraw()+ # Change Theme to Classic
  xlab(NULL)+ # No x-axis Label
  theme(axis.text.y.left = element_text(size = 8), # Adjust y-axis Text
        axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 1, angle = 0), # Ajust x-label Position
        axis.title.y = element_text(size = 8))+ # Adjust y-axis Title
  theme(strip.text = element_text(face = "bold", size = 8))+ # Adjust Headings
  scale_fill_manual(values = c("#d3d3d3", "#d3d3d3")) +
  theme(plot.title=element_text(size = 8, face = "bold", hjust = 0.5)) + #Change Title Size, Face, and Position
  geom_point(aes(color = Sample_Conditions),
             position = position_jitterdodge()) + # Layer on Points for Other Data Type
  scale_color_manual(name  = "Sample Conditions", values = viridis(6), labels = c("Inoculation\nCulture", "Control", "Increased\nAgitation", "Decreased\nTemperature", "Increased\nAirflow", "Decreased\nTemperature\nand\nIncreased\nAirflow")) + # Set Colors for Data Points
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "bottom") + # Change Legend Title and Text Size and Position
  guides(color = guide_legend(nrow = 1, byrow = FALSE, title.position = "top", title.hjust = 0.5)) # Change Legend Rows
  
## Make Data Frames for Alpha Diversity Metrics
alphadivFINAL <- estimate_richness(ppsFINAL) %>% 
  rownames_to_column(var = "Sample_ID") %>% 
  left_join(sdata2FINAL, by = "Sample_ID")
#View(alphadivFINAL)

alphadivFINAL <- write_csv(alphadivFINAL, "/data/home/htstodda/Biofilm_Fastq/workpath/alpha_diversityFINAL.csv") # Save Alpha Diversity Data

## Making a PCoA Plot
## Adding a Phylogenetic Tree to Phyloseq Object Using ape Library
random_treeFINAL = rtree(ntaxa(ppsFINAL), rooted=TRUE, tip.label=taxa_names(ppsFINAL))
plot(random_treeFINAL)

justbacteriaFINAL = merge_phyloseq(ppsFINAL, samdataFINAL, random_treeFINAL)
justbacteriaFINAL

## Ordination
uni_distanceFINAL <- ordinate(
  physeq = justbacteriaFINAL, 
  method = "PCoA", 
  distance = "unifrac"
)

plot_ordination(
  physeq = justbacteriaFINAL, # Phyloseq Object
  ordination = uni_distanceFINAL)+ # Ordination
  geom_point(aes(fill = Sample_Conditions, shape = Sample_Type), size = 3) + # Sets the Fill Color to the Data Type
  scale_fill_manual(values = viridis(6), labels = c("Inoculation\nCulture", "Control", "Increased\nAgitation", "Decreased\nTemperature", "Increased\nAirflow", "Decreased\nTemperature\nand\nIncreased\nAirflow")) +
  scale_shape_manual(values = c(22,21))+
  theme_linedraw() + # Changes the Theme and Removes the Background
  theme(                             
    legend.title = element_blank(), #Removes the Legend Title
    legend.position = "right",
    legend.text = element_text(size = 8),                                 
    axis.text.y.left = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 8))+
  guides(fill = guide_legend(override.aes = list(shape = c(21)))) #Fills in the Legend Points

unifrac_distanceFINAL <- as.data.frame(uni_distanceFINAL$vectors) %>%
  rownames_to_column(var = "Sample_ID") %>%
  left_join(sdata2FINAL, by = "Sample_ID") %>%
  select(Sample_ID, Sample_Conditions, Sample_Type, row_number, everything())

write_csv(unifrac_distanceFINAL, "/data/home/htstodda/Biofilm_Fastq/workpath/unifrac_distanceFINAL.csv") #Save UniFrac Distance Data to a File

## Exploring Taxa
## Summarize the Abundance of Each Class
genusabundanceFINAL <- ppsFINAL %>%
  tax_glom(taxrank = "Genus") %>%                      # Agglomerate at Class Level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to Rel. Abundance
  psmelt() %>%                                         # Melt to Long Format
  arrange(Genus) 
#View(genusabundanceFINAL)

## Select and Summarize the Neccesary Variable
allFINAL <- genusabundanceFINAL %>%
  select(Phylum, Class, Family, Genus, Sample, row_number, Sample_Type, Sample_Conditions, Sample_Temperature, Sample_Agitation_Speed, Sample_Air_Flowrate, Abundance) %>%
  filter(Abundance > 0) %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus))
#View(allFINAL)

genusFINAL <- allFINAL %>%
  filter(Abundance > 0) %>%
  group_by(Sample_Conditions, Sample_Type) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  dplyr::group_by(Sample_Conditions, Sample_Type, Phylum, Class, Family, Genus, totalSum) %>%
  unite(Genus_Species, Genus, sep = " ", remove = FALSE) %>%
  summarise(
    Abundance = sum(Abundance),
    Genus = ifelse(Abundance < 0.05, "< 5 %", Genus),
    Genus_Species = ifelse(Abundance < 0.05, "< 5 %", Genus_Species)) %>% # Change Genus Label to Group Low Abundance Taxa Together
  group_by(Sample_Conditions, Sample_Type, Genus, Genus_Species, totalSum) %>%  # Now Group and Summarize Again to Group Newly Labeled Low Abundance Taxa Together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()
#View(genusFINAL)

## Color Function for Figure
colFuncFINAL <- colorRampPalette(c('#000000', '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#fffac8', '#800000', '#aaffc3', '#a9a9a9'))
colorsFINAL <- colFuncFINAL(length(unique(genusFINAL$Genus_Species)))
names(colorsFINAL) <- unique(genusFINAL$Genus_Species)[order(unique(genusFINAL$Genus_Species))]
length(unique(genusFINAL$Genus_Species))
length(colorsFINAL)

ggplot(genusFINAL)+
  geom_col(mapping = aes(x = Sample_Type, y = RelAb, fill = Genus_Species), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Sample_Conditions), labeller = labeller(Sample_Conditions = label_wrap_gen(width = 15)))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colorsFINAL)+                        
  xlab(NULL)+
  theme_linedraw()+
  labs(fill = "Genus")+
  theme(axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        legend.spacing.x = unit(1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 8, face = "bold", angle = 0),
        strip.text.y = element_text(size = 8, face = "bold", angle = 270),
        legend.title = element_text(size = 8, face = "bold"))+
  guides(fill=guide_legend(nrow=3,byrow=FALSE))

write_csv(genusFINAL, "/data/home/htstodda/Biofilm_Fastq/workpath/genus_relative_abundances_w_speciesFINAL.csv") # Save Rel. Abundance Data to a File


### Section 5: Statistics

## Alpha Diversity Stats
## Make Alpha Diversity Metrics into Data Frames
observed <- alphadivFINAL[2]
shannon <- alphadivFINAL[7]

## Append Those Objects to the Alpha Diversity Data
sdataPERMANOVA <- alphadivFINAL[12:13]
sdata2Div <- cbind(sdataPERMANOVA, observed, shannon)

## Krukall-Wallis and Dunn Tests
dunn.test(sdata2Div$Observed, sdata2Div$Sample_Type, method="bh", kw=TRUE, 
          table=TRUE, list=TRUE, altp = TRUE) # Kruskall-Wallis and Dunn Test Comparing Sample Type Using Observed Metric
dunn.test(sdata2Div$Shannon, sdata2Div$Sample_Type, method="bh", kw=TRUE, 
          table=TRUE, list=TRUE, altp = TRUE) # Kruskall-Wallis and Dunn Test Comparing Sample Type Using Shannon Metric

dunn.test(sdata2Div$Observed, sdata2Div$Sample_Conditions, method="bh", kw=TRUE, 
          table=TRUE, list=TRUE, altp = TRUE) # Kruskall-Wallis and Dunn Test Comparing Sample Conditions Using Observed Metric
dunn.test(sdata2Div$Shannon, sdata2Div$Sample_Conditions, method="bh", kw=TRUE, 
          table=TRUE, list=TRUE, altp = TRUE) # Kruskall-Wallis and Dunn Test Comparing Sample Conditions Using Shannon Metric

## Beta Diversity Stats
## PEMANOVA Test on UniFrac Distance
## Calculate distance matrix
dist=phyloseq::distance(justbacteriaFINAL, method="uunifrac")

## Calculate PERMANOVA
permanova=adonis(dist ~ Sample_Type + Sample_Conditions, data=sdataPERMANOVA)
permanova
