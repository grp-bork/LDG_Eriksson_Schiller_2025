### This script calculates the richness of European timeseries studies (obtained from Microbeatlas). It compares the December to June difference in richness to corresponding latitude ranges from the habitat modeling
### Author: Jonas Schiller & Enzo Faucher

# load dependencies
library(ggplot2)
library(data.table)
library(lubridate)
library(tibble)
library(vegan)

# Data
dir.create("../../Data")
dir.create("../../Data/Internal")
dir.create("../../Data/External") #should already contain Global_percentageIncrease_ensemble_mean_sd.csv.zip & microbeatlas_taxonomic_profile.csv.zip & microbeatlas_studies_meta.csv
dir.create("../../Data/Internal/FigS7_Europe_seasonality_JS")
# Output
dir.create("../Output")
dir.create("../Output/FigS7_Europe_seasonality_JS")

### read data
## model
modeled = fread(unzip("../../Data/External/Global_percentageIncrease_ensemble_mean_sd.csv.zip", 
                   files = "Global_percentageIncrease_ensemble_mean_sd.csv", 
                   exdir = tempdir()))
## european timeseries studies
# read otu
all_otu_df = fread(unzip("../../Data/External/microbeatlas_taxonomic_profile.csv.zip", 
                      files = "microbeatlas_taxonomic_profile.csv", 
                      exdir = tempdir()))
# read meta
all_metadata_df = fread("../../Data/External/microbeatlas_studies_meta.csv")

## Calculate richness of european timeseries studies
# filter
filtered_metadata = all_metadata_df %>%
  filter(sample_id %in% colnames(all_otu_df)) %>%
  arrange(match(sample_id, rownames(all_otu_df)[-1]))

## preprocess filtered meta
# Define a function to get the month name from the date
get_month = function(date) {
  tryCatch(
    {
      month_num = month(date)
      month_name = month.name[month_num]
    },
    error = function(e) {
      month_name = "Missing"
    }
  )
  return(month_name)
}
# Apply the function to create 'sampling_month' column
filtered_metadata = filtered_metadata %>%
  mutate(collection_date = ymd(collection_date),  # Convert to date format if not already
         sampling_month = sapply(collection_date, get_month))
#Add sampling season
filtered_metadata = filtered_metadata %>%
  mutate(sampling_season = case_when(
    sampling_month %in% c("December", "January", "February") ~ "Winter",
    sampling_month %in% c("March", "April", "May") ~ "Spring",
    sampling_month %in% c("June", "July", "August") ~ "Summer",
    sampling_month %in% c("September", "October", "November") ~ "Autumn",
    TRUE ~ "Unknown"
  ))

# Function to perform rarefaction and calculate richness
rarefy_and_richness = function(df, sample_size, iterations, seed = 42) {
  richness_results = matrix(NA, nrow = nrow(df), ncol = iterations)
  rownames(richness_results) = rownames(df)
  
  for (i in 1:nrow(df)) {
    if (sum(df[i, ]) < sample_size) {
      # Not enough reads, leave as NA
      next
    } else {
      for (j in 1:iterations) {
        set.seed(seed + j + i * 1000)  # unique seed per sample/iteration
        rarefied_row = rrarefy(df[i, , drop = FALSE], sample_size)
        richness_results[i, j] = sum(rarefied_row > 0)
      }
    }
  }
  
  median_richness = apply(richness_results, 1, median, na.rm = TRUE)
  return(median_richness)
}

otu_matrix = all_otu_df %>% column_to_rownames("OTU_id") %>%
  t()

# Define sample size for rarefaction and number of iterations
sample_size = 3000
iterations = 20

# Perform rarefaction and calculate median alpha diversity
median_alpha_diversity = rarefy_and_richness(otu_matrix, sample_size, iterations)

# join with meta
median_alpha_diversity_df = median_alpha_diversity %>%
  data.frame() %>%
  rename( "median_alpha_diversity"=".") %>%
  rownames_to_column("sample_id") %>%
  inner_join(filtered_metadata, by = "sample_id")

# Ensure sampling_month is a factor with levels in the correct order
median_alpha_diversity_df$sampling_month = factor(median_alpha_diversity_df$sampling_month, levels = month.name)

# plot
ggplot(median_alpha_diversity_df, aes(x = sampling_month, y = median_alpha_diversity, color = study_name, group = study_name)) +
  #geom_point()+
  geom_smooth(method = "loess", span=1, se = FALSE) +
  labs(title = "",
       x = "Sampling Month",
       y = "Alpha Diversity (Richness)",
       color = "Study name") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12))

fwrite(median_alpha_diversity_df, "../../Data/Internal/FigS7_Europe_seasonality_JS/median_alpha_diversity_df.csv")

## Calculate modeled December-June % difference in richness at "European" latitudes
# subseting to prokaryotic modeled richness
modeled_prokaryotes = modeled %>% filter(clade == "prokaryotes")

# copy
timeseries = median_alpha_diversity_df

# Compute timeseries % change
timeseries$sampling_month = factor(timeseries$sampling_month, 
                                    levels = c("January", "February", "March", "April", "May", "June",
                                               "July", "August", "September", "October", "November", "December"), 
                                    ordered = TRUE)

# percent diff: December to June (per study)
diversity_diff = timeseries %>%
  group_by(study_name, sampling_month) %>%
  summarise(median_alpha_div = median(median_alpha_diversity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sampling_month, values_from = median_alpha_div) %>%
  mutate(
    diff_Dec_Jun = December - June,
    pct_diff_Dec_Jun = ifelse(
      (December + June) == 0, NA, 
      ((December - June) / ((December + June) / 2)) * 100
    )
  ) %>%
  dplyr::select(study_name, pct_diff_Dec_Jun)

# Get lat/lon for each study
unique_location = timeseries %>%
  group_by(study_name) %>%
  summarise(
    latitude = if (n_distinct(latitude) == 1) unique(latitude) else mean(latitude, na.rm = TRUE),
    longitude = if (n_distinct(longitude) == 1) unique(longitude) else mean(longitude, na.rm = TRUE),
    .groups = "drop"
  )

final_timeseries = diversity_diff %>%
  left_join(unique_location, by = "study_name")

# Compute latitudes of interest from timeseries
latitudes_of_interest = unique(round(final_timeseries$latitude, 0))

# Subtract 0.5 to match modeled latitudes like 41.5, 42.5, etc.
modeled_prokaryotes_lat = modeled_prokaryotes %>%
  filter(y %in% (latitudes_of_interest - 0.5))

# Aggregate
modeled_prokaryotes_lat_avg = modeled_prokaryotes_lat %>%
  group_by(y) %>%
  summarise(
    median_pct_increase = median(as.numeric(mean_pct_increase), na.rm = TRUE),
    .groups = "drop"
  )

## Combine both datasets
# Timeseries with clean study name labels
timeseries_labeled = final_timeseries %>%
  mutate(
    study_label = gsub("_", " ", study_name),
    source = "Timeseries",
    value = pct_diff_Dec_Jun
  ) %>%
  dplyr::select(source, value, study_label)

# Modeled data
modeled_labeled = modeled_prokaryotes_lat_avg %>%
  mutate(source = "Modeled", value = median_pct_increase) %>%
  dplyr::select(source, value)

# Merge
plot_df = bind_rows(timeseries_labeled, modeled_labeled)

# Plot
ggplot(plot_df, aes(x = source, y = value)) +
  geom_boxplot(outlier.shape = NA, fill = "gray90", color = "black") +
  geom_jitter(
    data = filter(plot_df, source == "Timeseries"),
    aes(color = study_label),
    width = 0.2, size = 3
  ) +
  geom_jitter(
    data = filter(plot_df, source == "Modeled"),
    color = "dodgerblue", width = 0.2, size = 3
  ) +
  labs(
    x = NULL,
    y = "% Change (December vs June)",
    title = "Comparison of Timeseries and Modeled % Change",
    color = "Study"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+
  ylim(c(0,100))+
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
ggsave("../Output/FigS7_Europe_seasonality_JS/boxplot.png")
ggsave("../Output/FigS7_Europe_seasonality_JS/boxplot.svg")

# print medians shown in plot
plot_df %>%
  group_by(source) %>%
  summarise(median_value = median(value, na.rm = TRUE))
