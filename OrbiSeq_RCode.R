# OrbiSeq: Optimized library preparation, sequencing, and data analysis for the generation of Orbivirus consensus sequences
# BTV Extraction/Sequencing Optimization - Combined Results
## Made by TD 05.15.24/08.03.25
#######################################################

#Load Required Libraries
library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)
library(forcats)
library(stats)
library(rstatix)
library(ggpubr)

packageVersion("rstatix")

#######################################################


############################################################################### RNA Extraction & Treatments

# read in data in spreadsheet
df <- read_excel("OrbiSeq_Raw_Data.xlsx", sheet = "sample_prep_op_final")

# add line break
df$ID <- gsub("\\(\\+\\) Secondary Purification", "(+) Secondary\nPurification", df$ID)
df$ID <- gsub("\\(\\-\\) Secondary Purification", "(-) Secondary\nPurification", df$ID)

# reorder samples
df$Group <- factor(df$Group)
df$Group <- factor(df$Group, levels = c('Primary RNA Purification', 
            'Carrier RNA', 'DNase Treatment', 
            'LiCl Treatment', 'Secondary RNA Purification'), ordered = TRUE) 

# subset data
dfqubit<- subset(df, Assay == "Qubit")
dfaln<- subset(df, Assay == "Percent Aligned")
dfcov<- subset(df, Assay == "Coverage")


# Perform Statistics


# QUBIT
# calculate mean, standard deviation, and 95% CI
dfqubit_DS <- dfqubit %>%
  group_by(Assay, Group, ID) %>%
  summarise(
    mean = mean(Output, na.rm = TRUE),
    sd = sd(Output, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se)

# PERCENT ALIGNED
# calculate mean, standard deviation, and 95% CI
dfaln_DS <- dfaln %>%
  group_by(Assay, Group, ID) %>%
  summarise(
    mean = mean(Output, na.rm = TRUE),
    sd = sd(Output, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se)

# COVERAGE
# calculate mean, standard deviation, and 95% CI
dfcov_DS <- dfcov %>%
  group_by(Assay, Group, ID) %>%
  summarise(
    mean = mean(Output, na.rm = TRUE),
    sd = sd(Output, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se)

# perform Mann-Whitney U test on all samples
df_MWU <- df %>% group_by(Assay, Group) %>% wilcox_test(Output ~ ID)
# *** All comparisons deemed not significant.


# Make combined (bar/point + error bar) plots with descriptive statistics


# QUBIT
# create plot
dfqubit_plot <- ggplot(dfqubit_DS, aes(x = ID, y = mean)) +
  geom_jitter(width = 0.05, height = NULL, data = dfqubit, aes(x = ID, y = Output), 
              stat = "identity", size = 4, stroke = 0.25, alpha = 0.5) + 
  geom_col(aes(fill = Group), colour = "black", alpha = 0.3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  xlab("") +
  ylab("Concentration (ng/uL)") +
  facet_grid(vars(Assay), vars(Group), scales = "free") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(strip.text.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(plot.title = element_text(size = 14)) +
  theme(legend.position = "none") + 
  theme(strip.text = element_text(size = 10))

# view plot
dfqubit_plot

# PERCENT ALIGNED
# create plot
dfaln_plot <- ggplot(dfaln_DS, aes(x = ID, y = mean)) +
  geom_jitter(width = 0.05, height = NULL, data = dfaln, aes(x = ID, y = Output), 
              stat = "identity", size = 4, stroke = 0.25, alpha = 0.5) + 
  geom_col(aes(fill = Group), colour = "black", alpha = 0.3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  xlab("") +
  ylab("Percent Aligned (%)") + 
  facet_grid(vars(Assay), vars (Group), scales = "free") + 
  scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
  theme_bw() +
  theme(strip.text.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(plot.title = element_text(size = 10)) +
  theme(legend.position = "none") + 
  theme( strip.text.x = element_blank() )

# view plot
dfaln_plot

# COVERAGE
# create plot
dfcov_plot <- ggplot(dfcov_DS, aes(x = ID, y = mean)) +
  geom_jitter(width = 0.05, height = NULL, data = dfcov, aes(x = ID, y = Output), 
              stat = "identity", size = 4, stroke = 0.25, alpha = 0.5) + 
  geom_col(aes(fill = Group), colour = "black", alpha = 0.3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  xlab("") +
  ylab("Coverage Depth (X)") + 
  facet_grid(vars(Assay), vars (Group), scales = "free") + 
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(strip.text.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + 
  theme(plot.title = element_text(size = 10)) +
  theme(legend.position = "none") + 
  theme( strip.text.x = element_blank() )

# view plot
dfcov_plot

# use patchwork to string plots together
big_DS_p <- ((dfqubit_plot)/(dfaln_plot)/(dfcov_plot)) + plot_annotation(tag_levels = 'A')

# view big plot
big_DS_p


# save plot to PDF
# ggsave("Illumina_sample_prep_op_combined_colored_082825.pdf", height=10, width=13, units="in")


############################################################################### END RNA Extraction & Treatments










########################################################################## Library Preparation Optimization
# read in data in spreadsheet
df2 <- read_excel("OrbiSeq_Raw_Data.xlsx", sheet = "lib_prep_op_final")

# add line break
df2$ID <- gsub("65ºC for 1 min.", "65ºC for\n1 min.", df2$ID)
df2$ID <- gsub("85ºC for 3 min.", "85ºC for\n3 min.", df2$ID)

# reorder facets
df2$Group <- factor(df2$Group, levels = c('Library Preparation Kit', 'KAPA Fragmentation')) 

# convert the 'ID' column back to factor
df2$ID <- factor(df2$ID, levels = unique(df2$ID))

# subset data
df2libprepkit<- subset(df2, Group == "Library Preparation Kit")
df2frag<- subset(df2, Group == "KAPA Fragmentation")


# Perform Statistics


# LIBPREP KIT
# calculate mean, standard deviation, and 95% CI
df2libprepkit_DS <- df2libprepkit %>%
  group_by(Group, ID) %>%
  summarise(
    mean = mean(Output, na.rm = TRUE),
    sd = sd(Output, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se)

# FRAGMENTATION
# calculate mean, standard deviation, and 95% CI
df2frag_DS <- df2frag %>%
  group_by(Group, ID) %>%
  summarise(
    mean = mean(Output, na.rm = TRUE),
    sd = sd(Output, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se)

# perform Mann-Whitney U test on all samples
df2_MWU <- df2 %>% group_by(Group) %>% wilcox_test(Output ~ ID)


# Make combined (bar/point + error bar) plots with descriptive statistics


# LIBPREP KIT
# create plot
df2libprepkit_plot <- ggplot(df2libprepkit_DS, aes(x = ID, y = mean)) +
  geom_jitter(width = 0.1, height = NULL, data = df2libprepkit, aes(x = ID, y = Output), 
              stat = "identity", size = 4, stroke = 0.25, alpha = 0.5) + 
  geom_col(fill = "#440555", colour = "black", alpha = 0.3, width = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  xlab("") +
  ylab("Concentration (ng/uL)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(plot.title = element_text(size = 14)) +
  facet_wrap(~Group, scales = "free_x") +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none") + 
  theme( strip.text.x = element_blank())

# view plot
df2libprepkit_plot

# FRAGMENTATION
# create plot
df2frag_plot <- ggplot(df2frag_DS, aes(x = ID, y = mean)) +
  geom_jitter(width = 0.1, height = NULL, data = df2frag, aes(x = ID, y = Output), 
              stat = "identity", size = 4, stroke = 0.25, alpha = 0.5) + 
  geom_col(fill = "#3A558A" ,colour = "black", alpha = 0.3, width = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  xlab("") +
  ylab("Concentration (ng/uL)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(plot.title = element_text(size = 14)) +
  facet_wrap(~Group, scales = "free_x") +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none") +
  theme( strip.text.x = element_blank())

# view plot
df2frag_plot

# use patchwork to string plots together
big_DS_p3 <- ((df2libprepkit_plot)/(df2frag_plot)) + plot_annotation(tag_levels = 'A')

# view big plot
big_DS_p3

# save plot to PDF
# ggsave("Illumina_lib_prep_op_combined_colored_101325.pdf", height=7, width=5, units="in")

########################################################################## END Library Preparation Optimization






################################################################################## Long vs Short Sequencing
# read in data in spreadsheet
df3 <- read_excel("OrbiSeq_Raw_Data.xlsx", sheet = "long_v_short_final")

# subset data
df3aln<- subset(df3, Assay == "Percent Aligned")
df3normcov<- subset(df3, Assay == "Normalized Coverage")
df3gp<- subset(df3, Assay == "Genome Completeness")


# Perform Statistics


# GENOME COMPLETENESS
# calculate mean, standard deviation, and 95% CI
df3gp_DS <- df3gp %>%
  group_by(Assay, Group, ID) %>%
  summarise(
    mean = mean(Output, na.rm = TRUE),
    sd = sd(Output, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se)

# perform Mann-Whitney U test on all samples
# df3gp_MWU <- df3gp %>% wilcox_test(Output ~ Group)
# ***STATS ON GENOME COMPLETENESS IS NOT POSSIBLE AS ALL VALUES ARE 100%

# PERCENT ALIGNED
# calculate mean, standard deviation, and 95% CI
df3aln_DS <- df3aln %>%
  group_by(Assay, Group, ID) %>%
  summarise(
    mean = mean(Output, na.rm = TRUE),
    sd = sd(Output, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se)

# perform Mann-Whitney U test on all samples
df3aln_MWU <- df3aln %>% wilcox_test(Output ~ Group)

# NORMALIZED COVERAGE
# calculate mean, standard deviation, and 95% CI
df3normcov_DS <- df3normcov %>%
  group_by(Assay, Group, ID) %>%
  summarise(
    mean = mean(Output, na.rm = TRUE),
    sd = sd(Output, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se)

# perform Mann-Whitney U test on all samples
df3normcov_MWU <- df3normcov %>% wilcox_test(Output ~ Group)


# Make combined (bar/point + error bar) plots with descriptive statistics


# set group colors 
group_colors <- c("Illumina" = "#440555", "Nanopore" = "#3A558A")

# GENOME COMPLETENESS
# create plot
df3gp_DS_barplot <- ggplot(df3gp_DS, aes(x = ID, y = mean)) +
  geom_point(data = df3gp, aes(x = ID, y = Output), stat = "identity", size = 4, stroke = 0.25, alpha = 0.5) + 
  geom_col(fill = group_colors , colour = "black", alpha = 0.3, width = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  xlab("") +
  ylab("Percent Genome\nCompletness") + 
  ggtitle("") +
  facet_grid(vars(Assay), vars (Group), scales = "free") + 
  scale_y_continuous(limits = c(0, NA),labels = scales::percent) +
  theme_bw() +
  theme(strip.text.y = element_blank()) +
  theme( strip.text.x = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none") 

# view plot
df3gp_DS_barplot

# PERCENT ALIGNED
# create plot
df3aln_DS_barplot <- ggplot(df3aln_DS, aes(x = ID, y = mean)) +
  geom_jitter(width = 0.05, height = NULL, data = df3aln, aes(x = ID, y = Output), 
              stat = "identity", size = 4, stroke = 0.25, alpha = 0.5) + 
  geom_col(fill = group_colors , colour = "black", alpha = 0.3, width = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  xlab("") +
  ylab("Percent\nAligned (%)")+ 
  facet_grid(vars(Assay), vars (Group), scales = "free") + 
  scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
  theme_bw() +
  theme(strip.text.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none") + 
  theme( strip.text.x = element_blank() )

# view plot
df3aln_DS_barplot

# NORMALIZED COVERAGE
# create plot
df3normcov_DS_barplot <- ggplot(df3normcov_DS, aes(x = ID, y = mean)) +
  geom_jitter(width = 0.05, height = NULL, data = df3normcov, aes(x = ID, y = Output), 
              stat = "identity", size = 4, stroke = 0.25, alpha = 0.5) + 
  geom_col(fill = group_colors , colour = "black", alpha = 0.3, width = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  xlab("") +
  ylab("Coverage Depth(X)\n/Reads(M)") + 
  facet_grid(vars(Assay), vars (Group), scales = "free") + 
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(strip.text.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + 
  theme(legend.position = "none") + 
  theme( strip.text.x = element_blank() )

# view plot
df3normcov_DS_barplot

# use patchwork to string plots together
big_DS_p2 <- ((df3gp_DS_barplot)/(df3aln_DS_barplot)/(df3normcov_DS_barplot)) + plot_annotation(tag_levels = 'A')

# view big plot
big_DS_p2

# save plot to PDF
# ggsave("BTV_op_long_v_short_combined_colored_101325.pdf", height=6, width=8, units="in")

################################################################################## END Long vs Short Sequencing