# Network Meta-Analysis (NMA) Project

## Overview

This project implements a comprehensive **Network Meta-Analysis (NMA)** workflow using R and the `gemtc` package. Network Meta-Analysis is an advanced statistical technique that synthesizes evidence from multiple treatment comparison studies, allowing researchers to compare multiple treatments simultaneously, even when they have not been directly compared in clinical trials.

This implementation covers both **continuous outcomes** (e.g., ADAS-Cog scores) and **binary outcomes** (e.g., adverse events), with automated model fitting, diagnostics, and result generation.

## Project Structure

```
Network-Meta/
├── README.md                          # This file
├── Network Meta.pdf                   # Reference documentation
├── appendix_nma_workflow.R            # Main R analysis script
├── appendix_network_workflow.do       # Stata workflow (reference)
├── ADAS_Cog.csv                       # Continuous outcome data
├── AEs.csv                            # Binary outcome data (adverse events)
├── figures/                           # Output directory for plots
├── results/                           # Output directory for results tables
```

## Required Dependencies

### R Packages

The analysis requires the following R packages:

- **gemtc** - Network Meta-Analysis package implementing Bayesian hierarchical models
- **coda** - Analysis tools for Markov Chain Monte Carlo (MCMC) samples
- **ggplot2** - Data visualization
- **tidyr** - Data manipulation and reshaping
- **rjags** - Interface to JAGS (Just Another Gibbs Sampler)

### External Software

- **JAGS 4.x** - Required for Bayesian MCMC sampling
  - Download from: https://sourceforge.net/projects/mcmc-jags/files/
  - Windows users: Set `JAGS_HOME` environment variable or update the path in the R script

### Installation

```r
# Install required R packages
install.packages(c("gemtc", "coda", "ggplot2", "tidyr", "rjags"))

# For Linux users (example):
# sudo apt-get install jags

# For macOS users:
# brew install mcmc-jags
```

## Data Format

### Continuous Outcomes (ADAS_Cog.csv)

Required columns (case-insensitive, underscores/spaces ignored):

| Column | Type | Description |
|--------|------|-------------|
| **study** | text | Study identifier (unique per row + treatment combination) |
| **treatment** | text | Treatment name (lowercase alphanumeric + underscores) |
| **mean** | numeric | Mean outcome value for the arm |
| **stddev** | numeric | Standard deviation for the arm |
| **samplesize** | integer | Number of participants in the arm |

**Example:**
```
study,treatment,mean,stddev,samplesize
STUDY001,placebo,22.5,5.2,45
STUDY001,drug_a,18.3,6.1,43
STUDY002,placebo,23.1,4.9,50
STUDY002,drug_b,15.2,5.8,48
```

### Binary Outcomes (AEs.csv)

Required columns (case-insensitive, underscores/spaces ignored):

| Column | Type | Description |
|--------|------|-------------|
| **study** | text | Study identifier (unique per row + treatment combination) |
| **treatment** | text | Treatment name (lowercase alphanumeric + underscores) |
| **responders** | integer | Number of events (adverse events, responders, etc.) |
| **samplesize** | integer | Total number of participants in the arm |

**Example:**
```
study,treatment,responders,samplesize
TRIAL_A,placebo,12,100
TRIAL_A,drug_x,8,102
TRIAL_B,placebo,15,95
TRIAL_B,drug_x,5,98
```

## Validation Rules

The script enforces the following data quality checks:

1. **No blank values** - All study and treatment fields must be non-empty
2. **No duplicate rows** - Each study-treatment combination must appear only once
3. **Valid treatment IDs** - Treatment names must be lowercase alphanumeric or underscores
4. **Reference treatment** - "placebo" must be present in the data
5. **Multi-arm studies** - Each study must have at least 2 treatment arms
6. **Numeric validity** - Outcome values must convert successfully to numeric type
7. **Outcome range** - For binary outcomes, responder counts must be between 0 and sample size

## Analysis Workflow

### Step 1: Model Fitting

The script fits **two competing models**:

#### Consistency Model
- Assumes consistency between direct and indirect evidence
- Uses `mtc.model()` with `type = "consistency"`
- Recommended for main analyses

#### Unrelated Mean Effects (UME) Model
- Assumes no consistency between different comparisons
- Uses `mtc.model()` with `type = "ume"`
- Useful for sensitivity analysis
- **Note:** Not fit for networks with multi-arm trials due to gemtc limitations

### Step 2: MCMC Sampling

Each model is fit using Markov Chain Monte Carlo (MCMC):

- **n.adapt** = 20,000 (adaptation iterations, discarded)
- **n.iter** = 50,000 (posterior sampling iterations per chain)
- **n.chain** = 4 (number of independent chains)
- **thin** = 1 (keep every iteration)
- **seed** = 20240511 (reproducibility)

### Step 3: Diagnostics & Output

For each outcome, the script generates:

#### Model Summaries
- `summary_[outcome].txt` - Full model summary (consistency model)
- `summary_ume_[outcome].txt` - Full model summary (UME model)

#### Trace Plots
- `trace_[outcome].png` - MCMC chain traces (visual inspection for convergence)

#### Gelman-Rubin Convergence Diagnostics
- `gelman_[outcome].png` - Potential Scale Reduction Factor (PSRF) plots
  - PSRF < 1.05 indicates good convergence

#### Forest Plot
- `forest_[outcome].png` - Treatment effects vs. reference (placebo)

#### Ranking Analysis
- `rank_probability_[outcome].png` - Probability of each rank for each treatment
- `sucra_[outcome].png` - Surface Under the Cumulative RAnking curve (SUCRA)

#### League Tables
- `league_table_[outcome].csv` - All pairwise relative effects (log scale for binary outcomes)
- `league_table_or_[outcome].csv` - All pairwise odds ratios (binary outcomes only)
- `sucra_[outcome].csv` - SUCRA values for treatment ranking

#### Node-Splitting Analysis
- `nodesplit_[outcome].txt` - Inconsistency assessment results
- `nodesplit_[outcome].png` - Node-splitting diagnostics
- **Note:** Skipped if network structure doesn't permit

#### Heterogeneity Assessment
- `heterogeneity_[outcome].txt` - Tau parameter estimates (between-study heterogeneity)
- `heterogeneity_[outcome].png` - Heterogeneity visualization

### Step 4: Network Visualization

- `network_[outcome].png` - Network graph showing treatment connections and comparisons

## Running the Analysis

### Basic Execution

```r
# Source the script from R or RStudio
source("appendix_nma_workflow.R")
```

### Command Line Execution

```bash
# On Windows or Linux
Rscript appendix_nma_workflow.R
```

### Configuration

Modify the following lines in `appendix_nma_workflow.R` before running:

```r
# Line 278: Set JAGS installation path (Windows users especially)
Sys.setenv(JAGS_HOME = "D:/Program Files/JAGS/JAGS-4.3.2")

# Line 565-566: Add or modify outcome analyses
run_continuous_nma(file.path(project_dir, "ADAS_Cog.csv"), "ADAS_Cog")
run_binary_nma(file.path(project_dir, "AEs.csv"), "AEs")
```

## Statistical Methods

### Continuous Outcomes (ADAS-Cog Example)

- **Likelihood:** Normal distribution
- **Link function:** Identity
- **Model type:** Random-effects hierarchical Bayesian model
- **Prior:** Non-informative (gemtc defaults)
- **Interpretation:** Mean differences between treatments

### Binary Outcomes (Adverse Events Example)

- **Likelihood:** Binomial distribution
- **Link function:** Logit
- **Model type:** Random-effects hierarchical Bayesian model
- **Interpretation:** Log odds ratios (and odds ratios) between treatments

## Key Outputs Interpretation

### League Table

A symmetric matrix showing all pairwise comparisons:
- **Diagonal:** Reference (1.0 for odds ratios, 0 for log odds ratios)
- **Upper/lower triangle:** Relative effect of row treatment vs. column treatment
- **95% CI:** Included in results; CI not crossing 1 (or 0 for log scale) indicates statistical significance

### SUCRA Values

Ranking metric based on cumulative rank probabilities:
- **100%** = Best performing treatment
- **0%** = Worst performing treatment
- Higher values indicate better efficacy/safety

### Forest Plot

Visual representation of relative effects:
- **Point estimate:** Treatment effect size
- **Horizontal line:** 95% credible interval
- **Vertical line at 0:** No effect

## Troubleshooting

### JAGS Installation Issues

**Error:** `"rjags is installed but JAGS 4.x is not available"`

**Solution:**
1. Download JAGS 4.x from https://sourceforge.net/projects/mcmc-jags/files/
2. Install JAGS
3. Set environment variable or update the path in line 278:
   ```r
   Sys.setenv(JAGS_HOME = "path/to/JAGS-4.x.x")
   ```

### Data Validation Failures

**Error:** `"Blank study values detected"`

**Solution:** Check that all study and treatment fields are non-empty; trim whitespace using Excel or R

### Model Convergence Issues

**Solution:**
- Increase `n.adapt` and `n.iter` parameters
- Check data quality and consistency
- Review trace plots for signs of mixing problems
- Ensure sufficient trials in the network

## References

- Salanti G. (2012). "Indirect and mixed-treatment comparison, network, or multiple-treatments meta-analysis: several names, several benefits." *Journal of Clinical Epidemiology*, 65(5), 475-483.
- Dias S, Sutton AJ, Ades AE, et al. (2013). "Evidence synthesis for decision making." *Research Synthesis Methods*, 4(3), 230-246.
- Dias S, Welton NJ, Caldwell DM, et al. (2010). "Checking consistency in mixed treatment comparison meta-analysis." *Statistics in Medicine*, 29(7-8), 932-944.

## Notes

- All analyses use non-informative priors
- Results are fully Bayesian; credible intervals represent 95% probability ranges
- Random-effects models account for between-study heterogeneity
- Results include full MCMC diagnostics for transparency

## Contact & Attribution

This NMA workflow implements standard Bayesian hierarchical meta-analysis methods using the R `gemtc` package. Adaptations and documentation for this specific project structure.

---

*Last updated: 2026-05-11*
