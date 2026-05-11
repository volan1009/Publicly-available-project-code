# Methodological Simulation: Stress-Testing GBTM for Recurrent Count Outcomes

> **Author:** *Wandering Data Ghost #2*
> **Context:** Independent methodological framework developed to address analytical bottlenecks in psychiatric epidemiology.

## 📌 Overview

Group-Based Trajectory Modeling (GBTM) is widely used in psychiatric epidemiology to identify latent trajectory classes in recurrent count outcomes. However, as a fundamentally data-driven method, the robustness of estimated trajectories depends critically on data completeness.

This subdirectory contains the R codebase and methodological protocol for a **Monte Carlo stress-test simulation**. It is designed to quantify the magnitude of trajectory misclassification bias when data are degraded by:

1. **Observation Loss** (e.g., treatment success, lack of formal diagnosis, competing risks).
2. **Administrative Censoring** (e.g., study end date).

## 💡 Key Methodological Highlights

* **Empirically-Anchored DGM:** A robust Data Generating Mechanism (DGM) utilizing Poisson outcomes with quadratic mean trajectories (`y ~ t + t^2`), anchored to plausible clinical parameters.
* **Sequential Degradation Logic:** A chronologically realistic degradation pipeline simulating clinical "under-ascertainment" followed by study truncation.
* **Optimal Label Alignment:** An automated greedy matching algorithm to align un-ordered latent class predictions from `flexmix` with the underlying Ground Truth.
* **Quantification of Boundary Shift:** The simulation visually and numerically demonstrates how intermediate-burden patients are bidirectionally misclassified due to trajectory "flattening".

## 📊 Visualizing the "Attrition Cliff"

The output visualization demonstrates the measurement bias:
- **Black/blue lines:** True underlying mechanism
- **Orange line:** "Flattened" trajectory estimated by GBTM when facing observation loss and censoring
- **Key finding:** Intermediate-burden subjects are systematically misallocated to adjacent risk classes

## 📂 Repository Structure

* `Illustrative code.R` - Main R script for simulation (pilot study with 5,000 subjects)
* `PMD_trajectory_simulation.pdf` - Methodological research proposal and protocol
* `Rplot.pdf` - High-resolution output figure (trajectory recovery comparison)
* `LICENSE` - MIT License

## 🚀 Getting Started

### Dependencies

Install required R packages:

```r
install.packages(c("dplyr", "tidyr", "flexmix", "ggplot2"))
```

### Execution

1. Ensure all files are in your working directory.
2. Run the simulation:

```r
source("Illustrative code.R")
```

### What the Script Does

1. **Generates synthetic data:**
   - 5,000 subjects followed over 10 time periods
   - Three latent trajectory classes with Poisson-distributed count outcomes
   - Quadratic trajectory form with class-specific coefficients

2. **Applies data degradation:**
   - Observation loss (cumulative event threshold = 3)
   - Administrative censoring (random follow-up truncation)

3. **Fits GBTM models:**
   - Model 1: Complete (ground truth) data
   - Model 2: Degraded (observed) data

4. **Evaluates recovery:**
   - Greedy label alignment algorithm
   - Class-specific attrition rates
   - Visual comparison of trajectory curves

5. **Outputs:**
   - `Rplot.pdf`: Three-panel figure comparing theoretical vs. fitted trajectories
   - Console output: Attrition rates by class

## 📝 Model Specification

### Data Generation Process

For subject $i$ in class $k$ at time $t$:

$$\log(\lambda_{ikt}) = \beta_{k0} + \beta_{k1}t + \beta_{k2}t^2$$

$$Y_{ikt} \sim \text{Poisson}(\lambda_{ikt})$$

**Class-specific parameters:**
- Class 1 (Low Burden): $\beta_0 = -1.5$, $\beta_1 = 0.1$, $\beta_2 = -0.01$
- Class 2 (Intermediate): $\beta_0 = -0.5$, $\beta_1 = 0.3$, $\beta_2 = -0.05$
- Class 3 (High Burden): $\beta_0 = 0.5$, $\beta_1 = 0.2$, $\beta_2 = -0.02$

### Degradation Mechanisms

1. **Observation Loss:** Subject exits observation if cumulative events ≥ 3
2. **Administrative Censoring:** Follow-up duration ~ Negative Binomial(size=2, μ=2), capped at $T_{\max} + 2$

## 📊 Expected Results

The pilot simulation typically demonstrates:
- **Intermediate-burden class attrition:** 10–15% misclassification
- **Low/high-burden stability:** 2–5% misclassification
- **Trajectory flattening effect:** Orange curves lie below theoretical trajectories

## ⚖️ License

MIT License - see LICENSE file for details.

## 📚 Citation

```
Author: Wandering Data Ghost #2
Title: Methodological Simulation - Stress-Testing GBTM for Recurrent Count Outcomes
Repository: https://github.com/volan1009/Publicly-available-project-code
Date: 2026-05-11
```

