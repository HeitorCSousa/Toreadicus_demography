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

* **ID**: A unique numeric identifier for each individual lizard.
* **date**: The date of capture in YYYY-MM-DD format.
* **trap**: The numeric identifier (1-25) of the pitfall trap where the capture occurred.
* **svl**: The snout-vent length of the individual, measured in millimeters (mm).
* **sex**: The sex of the individual, recorded as `F` (Female), `M` (Male), or `I` (Immature or unidentified).

#### 2. `Toreadicus_data.rds`
This file contains the data organized to run the Bayesian Pradel Jolly-Seber Model.

* **ID**: A unique numeric identifier for each individual lizard.
* **date**: The date of capture in YYYY-MM-DD format.
* **trap**: The numeric identifier (1-25) of the pitfall trap where the capture occurred.
* **svl**: The snout-vent length of the individual, measured in millimeters (mm).
* **sex**: The sex of the individual, recorded as `F` (Female), `M` (Male), or `I` (Immature or unidentified).

#### 3. `Pontos_Armadilhas.txt`
This file contains the geographic coordinates of the pitfall traps
* **parcela**: Plot. We only have one plot in this study.
* **armadilha**: The numeric identifier (1-25) of the pitfall trap where the capture occurred.
* **trap**: The numeric identifier (1-25) of the pitfall trap where the capture occurred.
* **lat**: Latitude (in degrees)
* **long**: Longitude (in degrees).
* **local**: Locality where the study was conducted – _Parque Estadual do Lajeado_ (PEL).

#### 3. `desemp_loc_Pel`, `Dados_Tpref`, and `CT_PEL`
These files contain the raw data from the laboratory ecophysiological experiments used to generate the thermal performance curve and hours of activity based on the thermal preferences.

* **pel**: A unique numeric identifier for the individual tested.
* **especie or sp**: Species: _Tropidurus oreadicus_.
* **sexo**: The sex of the individual, recorded as `F` (Female), `M` (Male), or `I` (Immature or unidentified).
* **CTmin and CTmax**: Minimum and maximum critical temperatures (°C).
* **data**: date when the experiment was conducted.
* **temp**: Body temperature of the individual.
* **veloc and acel**: Velocity and acceleration achieved by the lizard during the trial (cm/s and cm/s^2).
* **CRC**: _Comprimento rostro-cloacal_ – The snout-vent length of the individual, measured in millimeters (mm).
* **CC**: _Comprimento da cauda_ – tail length of the individual, measured in millimeters (mm).
* **BC** _Base da cauda_ – tail base, _i.e.,_ the until the tail breakage or regenerated portion, measured in millimeters (mm).

#### 4. `macro_micro_clima_PEL_imp.rds`
This file contains the summarized monthly environmental, microclimatic, and ecophysiological predictor variables for each trap location.

* **trap**: The numeric identifier (1-25) of the pitfall trap.
* **month**: The month of the record (numeric, 1-12).
* **year**: The year of the record.
* **precip**: Total monthly precipitation (mm).
* **tmin_2m**: Monthly mean of the minimum daily air temperature at 2m height (°C).
* **tmed_sup**: Monthly mean of the surface (skin) temperature (°C).
* **insol**: Monthly mean of the maximum daily solar insolation (W/m²).
* **tmed_1m**: Monthly mean air temperature recorded at 1m height at the trap location (°C).
* **rh_1m**: Monthly mean of the maximum daily relative humidity at 1m height (%).
* **perf_med**: Mean estimated locomotor performance based on microclimate data (m/s).
* **ha**: Estimated monthly hours of activity (total hours where microclimatic temperature was within the species' preferred thermal range).
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
