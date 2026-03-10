Title: Code Repository for Global Marine Microbiome Diversity Analyses
Date: March 2026
Authors: Dominic Eriksson, Jonas Schiller

Description:
This directory contains all R scripts used to compute diversity metrics, extract model results, compute ensemble means and uncertainties, identify diversity hotspots, perform statistical tests, and generate all figures included in the associated manuscript.  
Scripts are grouped by analysis step.

--------------------------------------------------------------------------------
1_Alpha_diversity_computation
--------------------------------------------------------------------------------

01_alpha_diversity_clade_level_JS.R
Author: Jonas Schiller
Summary:
This script calculates species-level alpha diversity metrics (Richness, Shannon, Chao1) after rarefaction (i.e., 1,000 & 5,000 mOTU counts) for all prokaryotes and per taxonomic rank (domain–genus).  
Output files are used as target variables for habitat modeling and subsequent figure generation.

--------------------------------------------------------------------------------
2_Extract_model_output
--------------------------------------------------------------------------------

1_Extract_predictor_ranking_DE.R
Author: Dominic Eriksson
Summary:
Extracts and compiles predictor importance rankings from CEPHALOPOD species distribution model (SDM) outputs for prokaryotic taxa. Aggregates predictor rankings by variable and clade, generating feature importance tables for downstream analyses.

2_Extract_ensembleMembers_DE.R
Author: Dominic Eriksson
Summary:
Builds ensemble raster projections for CEPHALOPOD SDMs. Processes model output folders, constructs ensemble raster stacks for each algorithm and bootstrap replicate, and saves ensemble results for further use.

3_Extract_ensemble_means_with_predictivePower_DE.R
Author: Dominic Eriksson
Summary:
Calculates monthly and annual mean and standard deviation of species richness projections from CEPHALOPOD SDMs using parallel computation. Outputs ensemble raster stacks with predictive power metrics.

4_Get_successful_modeled_taxa_DE.R
Author: Dominic Eriksson
Summary:
Extracts successfully modeled taxa from CEPHALOPOD SDM outputs. Identifies failed modeling runs using log files and compiles a list of valid taxa for downstream analyses.

4b_Table_successfully_modeled_clades_DE.R
Author: Dominic Eriksson
Summary:
The script identifies which taxa have successfully generated diversity metric models by parsing directory names from the directories_with_maps.txt file. It filters for directories with a rarefaction depth of 5000, extracts clade and diversity metric information, creates presence–absence matrices, and generates tabular summaries and heatmap visualizations showing modeling success across taxonomic groups and diversity indices.

5 and 5b__Extract_ensemble_means_with_predictive_Power_genusLevel_DE
Author: Dominic Eriksson
Summary: 
The workflow aggregates ensemble model outputs across algorithms, time, and bootstraps to estimate species richness patterns per genus, quantify predictive performance and uncertainty, and produce standardized summaries for large-scale diversity analyses and visualization.

--------------------------------------------------------------------------------
3_Compute_ensembles_and_uncertainties
--------------------------------------------------------------------------------

1_AllEnsembleMembers_monthly_oneTable_DE.R
Author: Dominic Eriksson
Summary:
Compiles ensemble member outputs for diversity indices (Richness, Shannon, Chao1) into long-format tables. Filters, reshapes, and saves monthly ensemble data across all successfully modeled taxa.

2_Global_annual_richness_ensemble_mean_sd_DE.R
Author: Dominic Eriksson
Summary:
Computes annual ensemble mean and standard deviation for each diversity index (Richness, Shannon, Chao1) across all clades.

2b_Global_annual_richness_ensemble_mean_sd_genusLevel_DE.R
Author: Dominic Eriksson
Summary:
This R script computes global annual richness and Shannon diversity ensemble means and standard deviations for genus-level prokaryotic taxa.

2c_Remove_marginalSeas_using_shapefile_genusLevel_DE.R
Author: Dominic Eriksson
Summary: 
This R script removes marginal seas from prokaryotic diversity ensemble datasets using IHO ocean region shapefiles.

2d_Global_maps_annual_ensemble_mnean_richness_genusLevel_DE.R
Author: Dominic Eriksson
Summary: 
This script creates comprehensive multi-plot maps for prokaryotic genera showing annual mean richness and Shannon diversity patterns as global projections.

2e_LDG_ensemble_means_genusLevel_DE.R
Author: Dominic Eriksson
Summary: 
This R script creates comprehensive Latitudinal Diversity Gradients (LDG) for prokaryotic genera by processing ensemble datasets with marginal seas filtering and parallel computing optimization.

3_Difference_December_Junse_ensemble_mean_sd_DE.R
Author: Dominic Eriksson
Summary:
Computes the December–June difference in modeled diversity indices (Richness, Shannon, Chao1) across global grids. Aggregates across clades and algorithms to calculate mean, SD, coefficient of variation, and agreement in directional change.

4_LDG_annual_mean_sd_DE.R
Author: Dominic Eriksson
Summary:
Computes latitudinal diversity gradients (LDGs) using annual ensemble richness data. Aggregates by latitude and clade, computing mean, standard deviation, and confidence intervals.

5_LDG_monthly_mean_sd_DE.R
Author: Dominic Eriksson
Summary:
Computes LDGs for each month using ensemble richness data. Aggregates across models and months to derive mean and standard deviation per latitude, focusing on June and December for seasonal comparisons.

--------------------------------------------------------------------------------
4_Hotspots_diversity
--------------------------------------------------------------------------------

1_Top_richness_25Percentile_EnvConditions_DE.R
Author: Dominic Eriksson
Summary:
Identifies global diversity hotspots for each microbial clade (top 25th percentile of richness) and extracts corresponding environmental conditions at those locations.

--------------------------------------------------------------------------------
5_Wilcox_test
--------------------------------------------------------------------------------

1_LDG_wilcoxTest_across_bins_DE.R
Author: Dominic Eriksson
Summary:
Performs Wilcoxon tests comparing modeled richness between June and December within each 1° latitudinal bin to assess significance of seasonal LDG differences.

--------------------------------------------------------------------------------
6_Coastal_influence
--------------------------------------------------------------------------------

1_coastal_to_ocean_format_input_data__DE.R
Author: Dominic Eriksson
Summary:
This script prepares prokaryotic diversity data for coastal influence analysis by integrating richness metrics with environmental metadata and organizing samples into distance-based groups to assess how coastal proximity shapes marine microbial diversity patterns.
1b_Create_chunked_files_DE.R
Author: Dominic Eriksson
Summary:
This script partitions large prokaryotic diversity datasets by taxonomic group into smaller files to enable efficient parallel processing in the CEPHALOPOD species distribution modeling pipeline.
2_Extract_ensemble_members_DE.R
Author: Dominic Eriksson
Summary:
This script extracts individual ensemble members from CEPHALOPOD habitat model outputs for prokaryotic taxa, assembling bootstrap and temporal predictions into raster stacks for downstream uncertainty and variability analyses.
3_All_ensemble_members_monthly_oneTable_DE.R
Author: Dominic Eriksson
Summary:
This script compiles ensemble member outputs from CEPHALOPOD habitat models into a unified long-format dataset, integrating raster-based predictions across bootstrap replicates and time to enable efficient calculation of ensemble statistics and comparative analyses.
4_ensemble_mean_sd_richness_DE.R
Author: Dominic Eriksson
Summary:
This script computes ensemble summary statistics from prokaryotic richness predictions across coastal distance groups to quantify model uncertainty and generate spatial representations of prediction confidence.
4b_Remove_marginal_seas_DE.R
Author: Dominic Eriksson
Summary:
This script filters prokaryotic richness model outputs using IHO ocean region classifications to exclude marginal seas and retain only open ocean regions for analysis.
5_Visualization_DE.R
Author: Dominic Eriksson
Summary:
This script visualizes model uncertainty in prokaryotic richness predictions across coastal–open ocean gradients by mapping spatial patterns in the coefficient of variation.

--------------------------------------------------------------------------------
Figures
--------------------------------------------------------------------------------

Fig1_Sampling_coverage.R
Author: Dominic Eriksson, Jonas Schiller
Summary:
Visualizes the spatial, temporal, and vertical distribution of metagenomic samples, illustrating global sampling coverage.

Fig2_samplewise_alphadiv_latitude_JS.R
Author: Jonas Schiller
Summary:
Plots sample-wise alpha diversity versus latitude for surface (MLD) and mesopelagic samples. Includes statistical tests comparing tropical (0–40°) and temperate/polar (>40°) regions.

Fig3_GlobalMaps_absoluteDifference_December_June_Bacteria_Archaea_DE.R
Author: Dominic Eriksson
Summary:
Generates global maps of December–June differences in bacterial and archaeal richness, visualizing spatial trends and model uncertainty.

Fig3_GlobalMaps_annualBacrerial_archaeal_richness_DE.R
Author: Dominic Eriksson
Summary:
Generates annual global maps of bacterial and archaeal richness, showing mean and uncertainty distributions.

Fig3_LDG_absolute_difference_DE.R
Author: Dominic Eriksson
Summary:
Plots LDGs comparing June and December richness for Bacteria and Archaea, highlighting statistically significant latitudinal bins.

Fig3_LDG_annual_richness_BacteriaArchaea_DE.R
Author: Dominic Eriksson
Summary:
Plots global LDGs of annual richness for Bacteria and Archaea, showing mean and standard deviation along latitudinal gradients.

Fig4_annual_seasonal_LDG_JS.R
Author: Jonas Schiller
Summary:
Compares class-level and domain-level LDGs using pairwise Pearson correlations between modeled 1° latitude richness bins.

Fig5_Bubbleplot_predictor_ranking_DE.R
Author: Dominic Eriksson
Summary:
Generates bubble plots of top environmental predictors per microbial clade, showing median importance and IQR for each predictor variable.

Fig5_Global_maps_consensus_hotspots_DE.R
Author: Dominic Eriksson
Summary:
Generates global consensus hotspot maps, showing the number of clades exceeding the 75th percentile of richness per grid cell.

Fig5_Spider_plots_envParameters_DE.R
Author: Dominic Eriksson
Summary:
Generates spider plots of environmental parameters associated with diversity hotspots for microbial clades.

FigS1_Model_sensitivity_DE.R
Author: Dominic Eriksson
Summary:
FigS2_Sankey_JS.R
Author: Jonas Schiller
Summary:
Mean species richness across taxonomic ranks, for the surface mixed layer and the mesopelagic layer.
FigS3_LDG_surface_mesopelagic_stat_JS.R
Author: Jonas Schiller
Summary:
Tests robustness of observed surface LDGs by randomly subsampling surface samples to match mesopelagic sample counts.
FigS4_alpha_diversity_clade_level_GENUS_analysis_JS.R
Author: Jonas Schiller
Summary:
Analyzes genus-level richness results and compares them to species-level richness, focusing on contributions of Alphaproteobacteria and Cyanobacteriia.
FigS4_alpha_diversity_clade_level_GENUS_JS.R
Author: Jonas Schiller
Summary:
Calculates genus-level richness (after rarefaction) for overall prokaryotes and selected classes (Alphaproteobacteria, Cyanobacteriia).
FigS5_taxonomic_representation_JS.R
Author: Jonas Schiller
Summary:
Summarizes the fraction of the microbiome that is represented by marker gene based operational taxonomic untis (mOTUs) and the completeness of the community at the chosen rarefaction level.
FigS6_Global_maps_environmental_predictors_DE.R
Author: Dominic Eriksson
Summary:
Global maps of environmental conditions, including physiochemical parameters.
FigS7_Correlation_climatologies_vs_inSitu_measurements_DE.R
Author: Dominic Eriksson
Summary:
Correlation of monthly climatologies with in situ measurements from the Tara Ocean expeditions.
FigS8_Correlation_rarefaction_cutOffs_DE.R
Author: Dominic Eriksson
Summary:
Correlation of diversity metrics at two different rarefaction cutoffs.
FigS9_Correlation_bacterial_archaeal_richness_DE.R
Author: Dominic Eriksson
Summary:
Correlation fo archaeal and bacterial richness (global and within (sub)tropical area).
FigS10_LDGs_varying_ocean_basins_DE.R
Author: Dominic Eriksson
Summary:
Prokaryotic LDGs at latitude bands from different ocean basins.
FigS11_Europe_seasonality_JS.R
Author: Jonas Schiller
Summary:
Comparing the predicted seasonal richness change (prokaryotic) to stationary time series (Europe, HOT AlOHA).
FigS12_GlobalMaps_annual_Clade_richness_DE.R
Author: Dominic Eriksson
Summary:
Annual richness patterns (based on habitat models).
FigS13_LDG_per_genus_JS.R
Author: Jonas Schiller
Summary:
The similarity of species-richness based LDGs per genus.
FigS13_panel_C_DE.R
Author: Dominic Eriksson
Summary:
Latitudinal diversity gradients of the most species rich genera, Pelagibacter and Prochlorococcus_A.
FigS14_tanglegram_LDG_phylogeny_JS.R
Author: Jonas Schiller
Summary:
Compares LDG similarity among classes via pairwise Pearson correlations, visualized using phylogenetic tanglegrams.
FigS15_LDG_within_classes_JS.R
Author: Jonas Schiller
Summary:
Species-richness based LDGs of genera within classes.
FigS16_Environmental_predictor_multicollinearity_DE.R
Author: Dominic Eriksson
Summary:
Multicollinearity of environmental conditions.
FigS17_Environmental_predictor_usage_DE.R
Author: Dominic Eriksson
Summary:
Environmental predictor usage and median importance (across model ensemble members).

--------------------------------------------------------------------------------
Functions
--------------------------------------------------------------------------------

subset_ocean_regions.R
Summary:
Helper function used to subset spatial data (e.g., raster or shapefile objects) to specific ocean regions based on polygon masks. Utilized in figure generation scripts.
