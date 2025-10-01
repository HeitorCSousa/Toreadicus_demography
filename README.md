# Data and code for: Precipitation and temperature-dependent locomotor performance drive the abundance of a Neotropical lizard

**Authors:** Raquel da Silva Acácio, Heitor Campos de Sousa, Thiago Costa Gonçalves Portelinha, Guarino Rinaldi Colli, and Adriana Malvasio.

**Publication:** *Biotropica* (2025)

**DOI:** [INSERT DOI HERE ONCE AVAILABLE]

**Contact:** Heitor Campos de Sousa (heitorsousa.bio@gmail.com or heitor.sousa@uft.edu.br)

---

## Repository Content

This repository contains the data and R code necessary to reproduce the analyses and figures for the manuscript "Precipitation and temperature-dependent locomotor performance drive the abundance of a Neotropical lizard." The data and output with slow-running models are in separate subfolders.

## System Requirements

The analyses were conducted in R (version 4.2.0 or later). The following R packages are required to run the script. You will also need to have JAGS (Just Another Gibbs Sampler) installed on your system to run the Bayesian models.

* `jagsUI`
* `R2ucare`
* `INLA` (requires special installation, see https://www.r-inla.org/download-install)
* `RJDemetra`
* `GAMM4`
* `dplyr`
* `ggplot2`

---

## Metadata: File Descriptions

### Data Files

#### 1. `Toreadicus_PEL_imp.rds`
This file contains the complete capture data for all individuals of *Tropidurus oreadicus* from February 2018 to November 2021. Each row represents a single capture event.

* **local**: Locality where the study was conducted – _Parque Estadual do Lajeado_ (PEL).
* **campaign**: Number of the sampling campaign. One for every month (1-46)
* **date**: The date of capture in YYYY-MM-DD format.
* **monht.orig**: Original month of the sampling campaign.
* **year.orig**: Original month of the sampling campaign.
* **species**: Species: _Tropidurus oreadicus_.
* **fieldtag**: A unique numeric identifier for each individual lizard.
* **trap**: The numeric identifier (1-25) of the pitfall trap where the capture occurred.
* **toes**: Toes clipped to generate the unique ID.
* **id**: A unique numeric identifier for each individual lizard, based on the toes clipped.
* **recapture**: Whether the animal was already marked (Y or N).
* **mass**: The mass of the individual, measured in grams (g).
* **svl**: The snout-vent length of the individual, measured in millimeters (mm).
* **tl**: Tail length of the individual, measured in millimeters (mm).
* **tb** Tail base, _i.e.,_ tail length until the tail breakage or regenerated portion, measured in millimeters (mm).
* **tail** Whether the tail is broken (Y or N).
* **sex**: The sex of the individual, recorded as `F` (Female), `M` (Male), or `I` (Immature or unidentified).
* **eggs**: Whether the individual contained eggs (Y or N).
* **dead**: Whether the individual was found dead in the trap (Y or N).
* **colectors**: Names of the researchers who collected the data in the field.
* **obs**: Other observations.
* **ecophysio**: Whether the animal was taken to the laboratory to perform ecophysiological experiments (Y or N). This field is incomplete and was not filled precisely.
* **id_duplicated**: Whether the individual had a duplicated id with another individual (Y or N). The fieldtag rectified these situations.
* **year**: Year according to the date.
* **month**: Month according to the date.
* **day**: Day according to the date.
* **month.n**: Month since the first sampling campaign according to the date.


#### 2. `Toreadicus_data.rds`
This file contains the data organized to run the Bayesian Pradel Jolly-Seber Model.
Please, refer to the R code `Script_Dempgraphy_Toreadicus.R` for more information.

#### 3. `Points_Traps.txt`
This file contains the geographic coordinates of the pitfall traps
* **plot**: Plot. We only have one plot in this study.
* **trap**: The numeric identifier (1-25) of the pitfall trap where the capture occurred.
* **lat**: Latitude (in degrees)
* **long**: Longitude (in degrees).
* **local**: Locality where the study was conducted – _Parque Estadual do Lajeado_ (PEL).

#### 3. `perf_loc_Pel`, `data_Tpref`, and `CT_PEL`
These files contain the raw data from the laboratory ecophysiological experiments used to generate the thermal performance curve and hours of activity based on the thermal preferences.

* **pel**: A unique numeric identifier for the individual tested.
* **species or sp**: Species: _Tropidurus oreadicus_.
* **sex**: The sex of the individual, recorded as `F` (Female), `M` (Male), or `I` (Immature or unidentified).
* **CTmin and CTmax**: Minimum and maximum critical temperatures (°C).
* **date**: date when the experiment was conducted.
* **temp**: Body temperature of the individual.
* **veloc and acel**: Velocity and acceleration achieved by the lizard during the trial (cm/s and cm/s^2).
* **SVL**: The snout-vent length of the individual, measured in millimeters (mm).
* **TL**: Tail length of the individual, measured in millimeters (mm).
* **TB** Tail base, _i.e.,_ the tail length until the tail breakage or regenerated portion, measured in millimeters (mm).

#### 4. `macro_micro_clima_PEL_imp.rds`
This file contains the summarized monthly environmental, microclimatic, and ecophysiological predictor variables for each trap location.

* **trap**: The numeric identifier (1-25) of the pitfall trap.
* **year.GMT_3**: The year of the record (GMT-3).
* **month.GMT_3**: The month of the record (numeric, 1-12).
* **hour.GMT_3**: The hour of the record (numeric, 0-23).
* **temp_med**: Hourly mean air temperature recorded at 1m height at the trap location (°C).
* **temp_sd**: Hourly air temperature standard deviation recorded at 1m height at the trap location (°C).
* **rh_med**: Hourly mean of the relative humidity at 1m height (%).
* **rh_sd**: Hourly standard deviation of the daily relative humidity at 1m height (%).
* **dewpoint_temp.y**: Hourly Dewpoint temperature (mm) from ERA5 database (°C).
* **temp_2m.y**: Hourly air temperature at 2m height from ERA5 database (°C).
* **temp_skin.y**: Hourly skin surface temperature from ERA5 database (°C).
* **temp_soil.y**: Hourly soil temperature (10 cm) from ERA5 database (°C).
* **surf_pressure.y**: Hourly surface pressure from ERA5 database (°C).
* **solrad.y**: Hourly solar radiation from ERA5 database (W/m²).
* **precip**: Total hourly precipitation (mm).
* **Toreadicus_perf**: Mean estimated locomotor performance based on microclimate data (m/s).
* **Toreadicus_ha90**: Estimated monthly hours of activity considering the 90th percentile thermal preference range (total hours where microclimatic temperature was within the species' preferred thermal range).
* **Toreadicus_ha50**: Estimated monthly hours of activity considering the 50th percentile thermal preference range (total hours where microclimatic temperature was within the species' preferred thermal range).

### Code File

#### `Script_Demography_Toreadicus.R`
This R script contains all the code to perform the analyses presented in the manuscript. The script will:
1.  Load and process the data files.
2.  Generate the thermal performance curve using a GAMM.
3.  Run the Bayesian Pradel-Jolly-Seber (PJS) model using JAGS to estimate survival, recruitment, and capture probability.
4.  Perform time-series decomposition on the demographic parameters.
5.  Run the spatio-temporal N-mixture models using INLA to estimate abundance and its drivers.
6.  Generate all figures and tables included in the final manuscript and its supplementary materials.

---

## License
The code in this repository is released under the [MIT License](https://opensource.org/licenses/MIT). The data are available under the [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/) license.
