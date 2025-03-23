library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(RColorBrewer)

# Load data
df <- read_tsv("samples_metadata_and_essentials.tsv")
df <- df %>% filter(MotherBaby != "NegCtrl")
all_data <- read_tsv("species_mpa.txt")
all_data <- all_data %>%
  rename(microbe = MOB.001.BS.2)

# Convert 'Timepoint' to a factor
df_filtered <- df %>% drop_na(shannon) %>% mutate(Timepoint = as.factor(Timepoint))

# Create the box plot 
ggplot(df_filtered, aes(x = Timepoint, y = shannon, fill = Birthmode)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, color = "black") +  # Adding black borders
  stat_summary(fun = mean, geom = "point", shape = 17, size = 3, aes(group = Birthmode), 
               position = position_dodge(0.8), color = "#27AE60") +  # Elegant dark blue for mean points
  labs(title = "Shannon Diversity by Timepoint and Birthmode",
       x = "Timepoint",
       y = "Shannon Diversity",
       fill = "Birthmode") +
  scale_fill_manual(values = c("C-Section" = "#AEDFF7", "Vaginal" = "#F7B785")) +  # Elegant blue and red colors
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"), text = element_text(size = 12))

# Parameters
TIMELINE <- 3 # Time point to display
PERCENT_CUTOFF_OTHER <- 10 # Rate to classify as 'other'

# Extract the first three taxonomic levels
all_data <- all_data %>%
  mutate(first_3 = str_extract(microbe, "^[^|]+\\|[^|]+\\|[^|]+"))

# Group by first_3 and sum the counts
group <- all_data %>%
  select(-c(microbe, Taxa)) %>%
  group_by(first_3) %>%
  summarize(across(matches(paste0("BS\\.", TIMELINE, "$")), sum))

# Identify taxa that should be labeled as "Other"
group <- group %>%
  rowwise() %>%
  mutate(sum_across_m = sum(c_across(starts_with("M")))) %>%
  mutate(class = if_else(sum_across_m <= PERCENT_CUTOFF_OTHER, "Other", first_3)) %>%
  ungroup() %>%
  group_by(class) %>%
  summarize(across(starts_with("M"), sum))  # Sum the rows within each group (including "Other")

# Reshape the data from wide to long format and filter rows with non-zero sums
pivoted_group <- group %>%
  pivot_longer(cols = -class, names_to = "Subject", values_to = "Value") %>%
  pivot_wider(names_from = class, values_from = Value) %>%
  filter(rowSums(select(., -Subject)) != 0)

# Clean up the Subject column
pivoted_group <- pivoted_group %>%
  mutate(Subject = str_remove(Subject, paste0("\\.BS\\.", TIMELINE)))

# Reorder columns by mean abundance
column_order <- pivoted_group %>%
  select(-Subject) %>%
  summarize(across(everything(), mean)) %>%
  unlist() %>%
  sort(decreasing = TRUE) %>%
  names()

pivoted_group <- pivoted_group %>%
  select(c("Subject", all_of(column_order))) %>%
  left_join(df %>% select(Subject, Birthmode) %>% distinct(), by = "Subject")

# Pivot the data for ggplot
pivoted_group_long <- pivoted_group %>%
  pivot_longer(cols = -c(Subject, Birthmode), names_to = "Class", values_to = "Value") %>%
  filter(Value != 0)

# Create a dynamic color palette
unique_classes <- unique(pivoted_group_long$Class)
non_other_classes <- setdiff(unique_classes, "Other")
color_palette <- setNames(brewer.pal(min(length(non_other_classes), 12), "Paired"), non_other_classes)

# Add grey for the "Other" class
color_palette <- c(color_palette, Other = "grey")

# Plot with the legend placed at the bottom
ggplot(pivoted_group_long, aes(x = Subject, y = Value, fill = Class)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~Birthmode, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.ticks.x = element_line(color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.text.x = element_text(size = 12, face = "bold"),  # Increase facet label size
    panel.spacing = unit(2, "lines")  # Increase space between facets 
  ) +
  labs(
    y = "Percentage",
    x = "Subject",
    title = paste("Microbial Composition by Birth Mode - Time point: ", as.character(TIMELINE))
  ) +
  guides(fill = guide_legend(nrow = 3,byrow = TRUE))
  
