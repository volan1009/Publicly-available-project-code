# Spatial Heatmap Generation: Global Asthma Burden Analysis (1990–2019)

This subdirectory contains source code and data required to generate spatial heatmaps visualizing Global Burden of Disease (GBD) data for asthma outcomes across countries and time periods.

## 📊 Data Source

The primary dataset utilized for spatial heatmaps was downloaded from the **Global Burden of Disease (GBD)** data visualization hub:

- **Link:** [vizhub.healthdata.org/gbd-results/](https://vizhub.healthdata.org/gbd-results/)
- **Metrics:** 
  - ASMR (Age-Standardized Mortality Rate)
  - ASDR (Age-Standardized Disability Rate / DALY Rate)
- **Time Period:** 1990, 2019
- **Cause:** Asthma
- **Age:** Age-standardized estimates

## ⚠️ Important Notes on Data Cleaning

Please carefully review the data processing steps if you intend to reuse this code, as several country naming inconsistencies exist between GBD and standard geographic databases.

### Known Data Quality Issues

* **Naming Discrepancies:** There are known inconsistencies between the country names provided in the GBD dataset and the standardized country names used within the `rnaturalearth` package.
* **Custom Logic Application:** To resolve these discrepancies, specific data cleaning logic was applied to the country names in this code. **This logic may not be directly applicable to your own analyses** without careful review.

### Country Name Mapping

The script includes a comprehensive `case_when()` mapping that harmonizes GBD country names to Natural Earth geographic names. Key transformations include:

- `Vietnam` → `Viet Nam`
- `South Korea` → `Republic of Korea`
- `Russia` → `Russian Federation`
- `Iran` → `Iran (Islamic Republic of)`
- And 60+ additional country-specific mappings

## 📂 Repository Structure

* `heatmap.csv` - GBD asthma data (DALY rates and mortality rates, 1990 & 2019)
* `heatmap.R` - Main R script for map generation
* `map_asmr.pdf` - Output: Asthma-specific mortality rates (two-panel: 1990 vs 2019)
* `map_asdr.pdf` - Output: Asthma disability-adjusted life years (two-panel: 1990 vs 2019)

## 🚀 Getting Started

### Dependencies

Install the required R packages:

```r
install.packages(c(
  "dplyr",
  "ggplot2",
  "sf",
  "readr",
  "stringr",
  "rnaturalearth",
  "cowplot",
  "ggpubr",
  "scales"
))
```

### Execution

1. Ensure `heatmap.csv` is in your working directory.
2. Run the main script:

```r
source("heatmap.R")
```

### Script Output

The script generates two PDF files:
- **`map_asmr.pdf`**: Age-standardized asthma mortality rates
  - Side-by-side comparison of 1990 and 2019
  - Includes top 5 countries table for each year
  
- **`map_asdr.pdf`**: Age-standardized asthma DALY rates
  - Side-by-side comparison of 1990 and 2019
  - Includes top 5 countries table for each year

## 📝 Visualization Details

- **Map Projection:** Robinson projection for global context
- **Color Scale:** Log10-transformed diverging scale (blues for low burden, reds for high burden)
- **Missing Data:** Countries without data are rendered in light grey
- **Data Quality:** Top 5 countries by metric are displayed in a table overlay

## ⚖️ License

MIT License - see LICENSE file for details.

## 📚 References

Global Burden of Disease Collaborative Network. Global Burden of Disease Study 2019 (GBD 2019) Results. Seattle, United States: Institute for Health Metrics and Evaluation (IHME), 2020. Available from http://vizhub.healthdata.org/gbd-results.

