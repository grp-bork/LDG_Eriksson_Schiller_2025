Title: Code Repository for Global Marine Microbiome Diversity Analyses
Date: October 2025
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
This script calculates alpha diversity metrics (Richness, Shannon, Chao1) after rarefaction (i.e., 1,000 & 5,000 mOTU counts) at the species and higher taxonomic levels (domain, phylum, class).  
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
Figures
--------------------------------------------------------------------------------

Fig1_Sampling_coverage_DE.R
Author: Dominic Eriksson
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

FigS1_alpha_diversity_clade_level_GENUS_JS.R
Author: Jonas Schiller
Summary:
Calculates genus-level richness (after rarefaction) for overall prokaryotes and selected classes (Alphaproteobacteria, Cyanobacteriia).

FigS1_alpha_diversity_clade_level_GENUS_analysis_JS.R
Author: Jonas Schiller
Summary:
Analyzes genus-level richness results and compares them to species-level richness, focusing on contributions of Alphaproteobacteria and Cyanobacteriia.

FigS2_LDG_surface_mesopelagic_stat_JS.R
Author: Jonas Schiller
Summary:
Tests robustness of observed surface LDGs by randomly subsampling surface samples to match mesopelagic sample counts.

FigS7_Europe_seasonality_JS.R
Author: Jonas Schiller, Enzo Faucher
Summary:
Calculates richness for European time-series (Microbeatlas) samples and compares December–June differences to modeled habitat data at corresponding latitudes.

FigS8_tanglegram_LDG_phylogeny_JS.R
Author: Jonas Schiller
Summary:
Compares LDG similarity among classes via pairwise Pearson correlations, visualized using phylogenetic tanglegrams.

--------------------------------------------------------------------------------
Functions
--------------------------------------------------------------------------------

subset_ocean_regions.R
Summary:
Helper function used to subset spatial data (e.g., raster or shapefile objects) to specific ocean regions based on polygon masks. Utilized in figure generation scripts.
