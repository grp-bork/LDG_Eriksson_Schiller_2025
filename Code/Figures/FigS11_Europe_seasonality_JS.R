### This script calculates the richness of European timeseries studies (obtained from Microbeatlas). It compares the December to June difference in richness to corresponding latitude ranges from the habitat modeling
### Author: Jonas Schiller & Enzo Faucher

# load dependencies
library(ggplot2)
library(data.table)
library(lubridate)
library(tibble)
library(vegan)

# Data
dir.create("../Data")
dir.create("../Data/Internal")
dir.create("../Data/External") #should already contain Global_percentageIncrease_ensemble_mean_sd.csv & microbeatlas_taxonomic_profile.csv & microbeatlas_studies_meta.csv
dir.create("../Data/Internal/FigS11_Europe_seasonality_JS")
# Output
dir.create("./Output")
dir.create("./Output/FigS11_Europe_seasonality_JS")

### read data
## model
modeled = fread("../Data/External/Global_percentageIncrease_ensemble__mean_sd.csv")
## european timeseries studies
# read otu (microbeatlas)
all_otu_df = fread("../Data/External/microbeatlas_taxonomic_profile.csv")
# read meta (microbeatlas)
all_metadata_df = fread("../Data/External/microbeatlas_studies_meta.csv")

# metaG
## read data
# meta data
sheet2 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 4))
meta = merge(sheet2, sheet4, by = "biosample", all = TRUE)
# alpha diversity
file = "../Data/External/species_diversity_taxonomic_rank.xlsx"
sheets = excel_sheets(file)
alpha_diversity = NULL
for (s in sheets) {
  dt = as.data.table(read_excel(file, sheet = s))
  if (is.null(alpha_diversity)) {
    alpha_diversity = dt
  } else {
    # keep only columns not already present
    new_cols = setdiff(names(dt), names(alpha_diversity))
    if (length(new_cols) > 0) {
      alpha_diversity = merge(
        alpha_diversity,
        dt[, c("biosample", new_cols), with = FALSE],
        by = "biosample",
        all = TRUE
      )
    }
  }
}
# filter for Richness
richness_cols = c("biosample",colnames(alpha_diversity)[grep("Richness", colnames(alpha_diversity))])
alpha_diversity = alpha_diversity %>% dplyr::select(all_of(richness_cols))

# join contextual and alpha diversity data
meta_alpha_div = left_join(meta,
                           alpha_diversity)
# filter for water layer (MLD or mesopelagic)
meta_alpha_div = meta_alpha_div %>%
  filter(analysis_type == "MLD") %>%
  filter(rarefaction_included == "rarefaction_5000 & rarefaction_1000")

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

fwrite(median_alpha_diversity_df, "../Data/Internal/FigS11_Europe_seasonality_JS/median_alpha_diversity_df.csv")

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
ggsave("./Output/FigS11_Europe_seasonality_JS/boxplot.png")
ggsave("./Output/FigS11_Europe_seasonality_JS/boxplot.svg")

# print medians shown in plot
plot_df %>%
  group_by(source) %>%
  summarise(median_value = median(value, na.rm = TRUE))



### December/June change in diversity ALOHA
# biosamples belonging to the HOT ALOHA study
biosamples = c(
  "SAMN05991650-S001","SAMN05991651-S001","SAMN05991658-S001","SAMN05991659-S001",
  "SAMN05991666-S001","SAMN05991667-S001","SAMN05991673-S001","SAMN05991674-S001",
  "SAMN05991680-S002","SAMN05991686-S002","SAMN05991693-S001","SAMN05991700-S002",
  "SAMN05991707-S002","SAMN05991714-S001","SAMN05991721-S002","SAMN05991728-S001",
  "SAMN05991729-S001","SAMN11874501-S002","SAMN11874502-S001","SAMN11874503-S001",
  "SAMN11874513-S001","SAMN11874514-S001","SAMN11874515-S001","SAMN11874516-S001",
  "SAMN11874525-S002","SAMN11874526-S002","SAMN11874527-S001","SAMN11874537-S002",
  "SAMN11874538-S002","SAMN11874539-S001","SAMN11874546-S001","SAMN11874547-S001",
  "SAMN11874548-S001","SAMN11874558-S001","SAMN11874559-S001","SAMN11874570-S001",
  "SAMN11874571-S001","SAMN11874582-S001","SAMN11874583-S001","SAMN11874584-S001",
  "SAMN11874594-S001","SAMN11874595-S001","SAMN11874596-S001","SAMN11874606-S001",
  "SAMN11874607-S001","SAMN11874608-S001","SAMN11874618-S001","SAMN11874619-S001",
  "SAMN11874620-S001","SAMN11874630-S001","SAMN11874631-S001","SAMN11874632-S001",
  "SAMN11874633-S001","SAMN11874642-S001","SAMN11874643-S001","SAMN11874644-S001",
  "SAMN11874645-S001","SAMN11874654-S001","SAMN11874655-S001","SAMN11874656-S001",
  "SAMN11874666-S001","SAMN11874667-S001","SAMN11874668-S001","SAMN11874678-S001",
  "SAMN11874679-S001","SAMN11874680-S001"
)
meta_alpha_div_aloha = meta_alpha_div %>% filter(biosample %in% biosamples)

# ALOHA map
world = map_data("world")
ggplot() +
  geom_map(
    data = world,
    map = world,
    aes(long, lat, map_id = region),
    color = "white",
    fill = "darkgrey",
    size = 0.1
  ) +
  geom_point(
    data = meta_alpha_div_aloha,
    aes(x = lon, y = lat),
    size = 2,
    color = "black"
  ) +
  geom_text(
    data = meta_alpha_div_aloha,
    aes(x = lon, y = lat, label = "HOT ALOHA"),
    color = "black",
    hjust = -0.1,
    vjust = 0.9,
    size = 3
  ) +
  labs(x = "Longitude", y = "Latitude") +
  coord_equal(xlim = c(-180, -70)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("./Output/FigS11_Europe_seasonality_JS/aloha_map.png", width = 5, height = 7, units = "cm")
ggsave("./Output/FigS11_Europe_seasonality_JS/aloha_map.svg", width = 5, height = 7, units = "cm")


#### % Change in Richness
# response column
ycol = "prokaryotes_Richness"

# LOESS span
loess_span = 1

# loess fit
loess_mod = loess(
  formula = as.formula(paste0(ycol, " ~ month")),
  data = meta_alpha_div_aloha,
  span = loess_span,
  degree = 2,
  family = "gaussian",
  na.action = na.exclude,
  control = loess.control(surface = "direct")
)

# prediction of smooth + CI on a dense grid
grid_df = data.frame(month = seq(1, 12, by = 0.05))
grid_pred = predict(loess_mod, newdata = grid_df, se = TRUE)

grid_df = grid_df %>%
  mutate(
    fit = as.numeric(grid_pred$fit),
    se  = as.numeric(grid_pred$se.fit),
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se
  )

# Medians (raw data) for June + December
median_df = meta_alpha_div_aloha %>%
  filter(month %in% c(6, 12)) %>%
  group_by(month) %>%
  summarise(
    med = median(.data[[ycol]], na.rm = TRUE),
    n   = sum(!is.na(.data[[ycol]])),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      month.name[month], " median: ",
      round(med, 0),
      "\n(n=", n, ")"
    )
  )

# symmetric % difference based on medians
pct_diff_median_Dec_Jun = median_df %>%
  mutate(month_name = month.name[month]) %>%
  select(month_name, med) %>%
  pivot_wider(names_from = month_name, values_from = med) %>%
  mutate(
    diff_Dec_Jun = December - June,
    pct_diff_Dec_Jun = ifelse(
      (December + June) == 0, NA_real_,
      ((December - June) / ((December + June) / 2)) * 100
    )
  ) %>%
  select(pct_diff_Dec_Jun)

pct_diff_median_Dec_Jun


# plot: raw points + LOESS smooth + CI, no red points
ggplot(meta_alpha_div_aloha, aes(x = month, y = .data[[ycol]])) +
  geom_point() +
  geom_ribbon(
    data = grid_df,
    aes(x = month, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = grid_df,
    aes(x = month, y = fit),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = median_df,
    aes(
      x = month,
      y = med,
      label = label,
      hjust = ifelse(month == 12, 1.05, 0.5)
    ),
    color = "red",
    vjust = -1,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  coord_cartesian(ylim = c(0, 1400)) +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.name,
    expand = expansion(mult = c(0.02, 0.19))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = "", y = "Prokaryotic Richness") +
  annotate(
    "text",
    x = 1,
    y = 100,
    label = paste0("LOESS: span = ", loess_span),
    hjust = 0,
    color = "grey40",
    size = 3.5
  )
ggsave("./Output/FigS11_Europe_seasonality_JS/aloha_timeseries.png", width = 10, height = 9, units = "cm")
ggsave("./Output/FigS11_Europe_seasonality_JS/aloha_timeseries.svg", width = 10, height = 9, units = "cm")


# modeled pct increase at that latitude (unchanged)
latitudes_of_interest = unique(round(meta_alpha_div_aloha$lat, 0))

modeled_prokaryotes_lat = modeled_prokaryotes %>%
  filter(y %in% (latitudes_of_interest - 0.5))

modeled_prokaryotes_lat_avg = modeled_prokaryotes_lat %>%
  group_by(y) %>%
  summarise(
    median_pct_increase = median(as.numeric(mean_pct_increase), na.rm = TRUE),
    .groups = "drop"
  )

# barplot: use median-based % difference for timeseries
median_pct = pct_diff_median_Dec_Jun %>%
  summarise(value = as.numeric(pct_diff_Dec_Jun)) %>%
  mutate(type = "Timeseries")

modeled_pct = modeled_prokaryotes_lat_avg %>%
  summarise(value = median(median_pct_increase, na.rm = TRUE)) %>%
  mutate(type = "Modeled")

plot_df = bind_rows(modeled_pct, median_pct) %>%
  mutate(
    type = factor(
      type,
      levels = rev(c("Timeseries", "Modeled"))
    )
  )

ggplot(plot_df, aes(x = type, y = value)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f%%", value)),
    angle = 90,
    vjust = 0.5,
    hjust = -0.3,
    size = 3.5
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "", y = "% Change (December vs June)") +
  theme_minimal()

ggsave("./Output/FigS11_Europe_seasonality_JS/aloha_modeled_timeseries.png", width = 5, height = 7, units = "cm")
ggsave("./Output/FigS11_Europe_seasonality_JS/aloha_modeled_timeseries.svg", width = 5, height = 7, units = "cm")