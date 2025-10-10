Description:
This directory contains all data files and contextual information used to reproduce the figures and analyses presented in the associated manuscript. Below is a description of each file, its origin, purpose, and associated analysis or figure.

--------------------------------------------------------------------------------
Root Directory
--------------------------------------------------------------------------------

samples_contextual.xlsx
Summary:
Contains contextual information for all metagenomic samples analyzed in this study.

Contents:
- Samples: Biosample IDs and references to associated studies
- In situ information: Filter sizes, coordinates, sampling depth, sampling time, ocean basin
- Climatologies: Matched climatology data (see “Environmental predictors” in the Materials and Methods section)
- Rarefaction analysis: Information specifying which analysis each sample was included in (rarefaction cutoff, depth layer)

Usage:
Input for environmental matching, habitat modeling, and figure-specific sample selection.

--------------------------------------------------------------------------------

diversity_metrics.xlsx
Summary:
Diversity metrics calculated from mOTUs4_OMDB.csv.zip.

Contents:
- Diversity metrics: Richness, Shannon index, Chao1
- Two rarefaction cutoffs: 1,000 and 5,000 mOTU observations
- Includes all taxonomic groups with model reliability (R² > 0.25, see Materials and Methods)

Generated in: 01_alpha_diversity_clade_level_JS.R
Used in: Habitat modeling and Figure 2, Supplementary Figures 1–2

--------------------------------------------------------------------------------
Data_generated Directory
--------------------------------------------------------------------------------

alpha_div_genus.csv
Summary:
Diversity metrics (Richness, Shannon, Chao1) for Alphaproteobacteria and Cyanobacteriia at the genus level, under two rarefaction cutoffs (1,000 & 5,000 mOTU observations).

Generated in: FigS1_alpha_diversity_clade_level_GENUS_JS.R
Used in: FigS1_alpha_diversity_clade_level_GENUS_analysis_JS.R
Associated Figure: Supplementary Fig. 1

--------------------------------------------------------------------------------

df_latGradients_alphaDiv_annual.csv
Summary:
Modeled annual diversity (surface mixed layer) summarized at 1° latitude bins (mean and standard deviation) for all modeled taxonomic groups.

Generated in: 4_LDG_annual_mean_sd_DE.R
Used in: Fig4_annual_seasonal_LDG_JS.R
Associated Figure: Figure 4

--------------------------------------------------------------------------------

df_latGradients_alphaDiv_month.csv
Summary:
Modeled seasonal (monthly resolution) diversity (surface mixed layer) summarized at 1° latitude bins (mean and standard deviation) for all modeled taxonomic groups.

Generated in: 5_LDG_monthly_mean_sd_DE.R
Used in: Fig4_annual_seasonal_LDG_JS.R
Associated Figure: Figure 4

--------------------------------------------------------------------------------

Wilcox_pairwise_test_june_December.csv
Summary:
Results of Wilcoxon tests assessing differences in modeled diversity between December and June (averaged per 1° latitude) at European latitudes.

Generated in: 1_LDG_wilcoxTest_across_bins_DE.R
Used in: Fig4_annual_seasonal_LDG_JS.R
Associated Figure: Figure 4

--------------------------------------------------------------------------------

Global_percentageIncrease_ensemble_mean_sd.csv.zip
Summary:
Mean modeled percentage change in richness, comparing December to June at European time-series latitudes.

Generated in: 3_Difference_December_Junse_ensemble_mean_sd_DE.R
Used in: FigS7_Europe_seasonality_JS.R
Associated Figure: Supplementary Fig. S7

--------------------------------------------------------------------------------

dist_bac.csv / dist_arc.csv
Summary:
Patristic distances of modeled bacterial (dist_bac.csv) and archaeal (dist_arc.csv) taxonomic groups, obtained and subsetted from GTDB v220.

Used in: FigS8_tanglegram_LDG_phylogeny_JS.R
Associated Figure: Supplementary Fig. S8

--------------------------------------------------------------------------------

microbeatlas_studies_meta.csv
Summary:
Contextual information about samples in Microbeatlas_taxonomic_profile.csv.zip, obtained and curated from the Microbeatlas database.

Used in: FigS7_Europe_seasonality_JS.R
Associated Figure: Supplementary Fig. S7

--------------------------------------------------------------------------------

External Datasets (Not Included)
--------------------------------------------------------------------------------

mOTUs4 taxonomic profile
Summary:
Species-level mOTU taxonomic profiles (v4) for all samples, obtained from the OMDB database.

Accessibility: https://omdb.microbiomics.io/repository/ocean/ (Downloaded 2025)

--------------------------------------------------------------------------------

mOTUs4 to GDTB mapping
Summary:
Mapped GTDB taxonomy for each metagenomic operational taxonomic unit (mOTU).

Accessibility: https://omdb.microbiomics.io/repository/ocean/ (Downloaded 2025)

--------------------------------------------------------------------------------

Microbeatlas taxonomic profile
Summary:
Taxonomic profile (97% OTU identity) for European time-series studies, obtained from the Microbeatlas database.

Used in: FigS7_Europe_seasonality_JS.R
Associated Figure: Supplementary Fig. S7


Ocean_regions_shapefile
Summary:
Shapefiles of global oceans and seas.
Accessibility: https://www.marineregions.org/downloads.php#goas (v01; Downloaded 09/2022)

--------------------------------------------------------------------------------

Surface_predictors
Summary:
Global surface ocean climatologies (temperature, salinity, nutrients, etc.).
Accessibility: https://data.up.ethz.ch/shared/bluecloud/OBS/CLIMATOLOGIES/ (Downloaded 2024)

Climatology files used:
climatology_A_0_50.nc, climatology_A_CHLA_regridded.nc, climatology_A_PAR_regridded.nc, climatology_eke_aviso.nc,
climatology_fsle_aviso_2001_2020.nc, climatology_i_0_50.nc, climatology_M_0_0.nc, climatology_n_0_50.nc,
climatology_Nstar_0_50.nc, climatology_O_0_50.nc, climatology_p_0_50.nc, climatology_s_0_50.nc, climatology_t_0_50.nc

Log10-transformed climatologies for normalization:
log10_climatology_eke_aviso.nc, log10_climatology_i_0_50.nc, log10_climatology_n_0_50.nc, log10_climatology_p_0_50.nc,
log10_climatology_S_PP_regridded.nc, log10_climatology_TOT_POC_CMEMS.nc
