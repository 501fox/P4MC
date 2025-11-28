Markdown

# P4MC: A Monte Carlo Classification Method Based on Partial Variables for Postoperative Complication Prediction

<div align="center">
  <img src="framework/framework.png" width="100%" alt="P4MC Framework Flowchart">
</div>

<br>

## 📝 Abstract

Predicting postoperative complications is a critical task for improving patient care and optimizing resource management. We frame this challenge as a binary classification problem utilizing two-stage clinical data: preoperative and intraoperative variables. A major obstacle for this task is the unavailability of intraoperative data for preoperative prediction. Existing approaches, which typically rely solely on preoperative data or single intraoperative indicator imputation, often fail to capture critical predictive information.

To provide a more robust solution, we propose a novel method, **P4MC** (**P**reoperative, **P**artial intrao**P**erative, **P**ostoperative **M**onte **C**arlo classification). This approach explicitly models the dependency between the two stages to recover latent intraoperative information. Specifically, P4MC leverages Partial Dimension Reduction to compress high-dimensional intraoperative variables into a low-dimensional sufficient predictor. By employing Monte Carlo sampling on the conditional distribution of two-stage data, we generate latent features to serve as inputs for a classification model. The effectiveness of P4MC is validated through extensive simulation experiments. Furthermore, we demonstrate its practical utility on two real-world cancer cohorts: laparoscopic pancreaticoduodenectomy and hepatocellular carcinoma, showing that our method outperforms traditional machine learning models and yields clinically meaningful insights.


## 📦 Installation

You can install the development version of P4MC directly from GitHub:

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install P4MC
devtools::install_github("501fox/P4MC")
```

## 🚀 Quick Start
Here is how to load the data and run the prediction model:

```r

library(P4MC)

# 1. Load the example datasets
data("data_simulation")
data("data_LPD")

# Check the data structure
head(data_simulation)
head(data_LPD)

# 2. Run the P4MC evaluation
# Arguments:
#   datasets: The input dataframe
#   k: Number of folds for cross-validation (default: 10)
#   S: Number of MCMC sampling iterations (default: 200)
#   grid: Number of grid points for KDE (default: 25)
#   seed: Random seed for reproducibility (default: 123)
#   calculate_auc: Whether to compute AUC (default: FALSE)

results_sim <- run_P4MC(
  datasets = data_simulation,
  k = 10,             # Using 5-fold CV for this example
  csv_path = "P4MC_results_sim.csv",
  calculate_auc = TRUE
)

results_LPD <- run_P4MC(
  datasets = data_LPD,
  k = 5,              # Using 5-fold CV for this example
  csv_path = "P4MC_results_LPD.csv",
  calculate_auc = TRUE
)

# 3. View the final summary (Mean ± SD)
print(results_sim$final_table)
print(results_LPD$final_table)

# 4. Access detailed confusion matrices for each fold
# print(results_sim$confusion_tables)
# print(results_LPD$confusion_tables)
```
