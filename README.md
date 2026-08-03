# Evaluating Statistical Methods for Modeling and forecasting Hourly Precipitation

This repository contains the R code and selected outputs associated with the paper:

> **Evaluating statistical methods for modeling and forecasting hourly precipitation**  
> Rohan Hemant Chhatre, James O'Donnell, Nalini Ravishanker
> Japanese Journal of Statistics and Data Science

The project develops and evaluates probabilistic precipitation forecasting models based on a zero-adjusted gamma (ZAGA) distribution fitted using component-wise gradient boosting. Forecast performance is compared with benchmark and machine-learning approaches across weather stations in Connecticut, USA.

## Repository structure

```text
.
├── code/
│   ├── revisions/              # Code prepared during manuscript revisions.
│       ├── run_gamboost_meriden_all_final_rolling.R # Code for Model 3 from the paper.
│       ├── run_gamboost_meriden_rainonly_final_rolling.R # Code for Model 1 from the paper.
│       ├── run_gamboost_meriden_exogenous_final_rolling.R # Code for Model 2 from the paper.
│       ├── run_arima_meriden_rolling.R # Code for ARIMA modelfrom the paper.
│       ├── run_persistence_meriden_rolling.R # Code for Persistence model from the paper.
│       ├── run_random_forest_meriden_rolling.R # Code for Random Forest model from the paper.
│       ├── run_lstm_meriden_rolling.R # Code for LSTM from the paper.
│   ├── old code/               # Earlier development scripts retained for reference
│   ├── ct_map.R                # Map of the Connecticut weather stations
│   ├── events.R                # Event-based forecast comparisons
│   └── lead_times_rmse.R       # Forecast-error analysis by lead time
├── merged_imputed_plots/       # Diagnostic plots for imputed station data
├── output_github/output/       # Selected model outputs and figures
├── residuals/                  # Residual diagnostic outputs
└── README.md
```

## Forecasting models

The repository includes code for comparing the proposed ZAGA boosting models with forecasting benchmarks such as:

- persistence forecasting;
- ARIMA/ARIMAX;
- random forest;
- LSTM; and
- Three Models proposed under the ZAGA-GAMLSS framework.

The event-level evaluation includes selected high-impact precipitation events, including Tropical Storm Elsa, Hurricane Henri, and the remnants of Hurricane Ida.

## Software requirements

The analysis was conducted in **R**. Packages used across the scripts include, among others:

```r
install.packages(c(
  "data.table",
  "ggplot2",
  "lubridate",
  "maps"
))
```

Additional packages required by individual model-fitting scripts should be installed from CRAN or the appropriate source. In particular, the forecasting scripts may require packages for GAMLSS boosting and model-specific benchmark methods.

## Data

The original weather data are not publicly available due to data-use restrictions and file size limitations. Researchers interested in accessing the data may request it by contacting the corresponding author.
Corresponding author:
Rohan Chhatre
Email: `fbs24003@uconn.edu`

## Running the analysis

1. Clone or download this repository.
2. Open the repository as an R project or set the repository root as the working directory.
3. Install the required packages.
4. Obtain and prepare the input data.
5. Replace any remaining machine-specific file paths with project-relative paths.
6. Run the relevant scripts under `revisions`.
7. Run `code/lead_times_rmse.R` and `code/events.R` to reproduce the comparative evaluation and event-based figures.

Note: Avoid paths such as `/Users/name/Documents/...`, because they only work on the original author's computer.

## Results

Selected figures, predictions, residual diagnostics, and imputation plots are included in the output folders. These files are provided to facilitate comparison with the results reported in the paper.

## Citation

When using this code, please cite the accompanying paper:
A machine-readable `CITATION.cff` file is included in the repository to cite the paper and the code.

## Contact

For questions about the code or paper, please contact:

**Rohan Chhatre**  
Department of Statistics, University of Connecticut  
GitHub: [Rohan-Chhatre](https://github.com/Rohan-Chhatre)
Email: [fbs24003@uconn.edu](mailto:fbs24003@uconn.edu)

## Acknowledgements

The authors are grateful to the anonymous reviewers and the editors for their valuable suggestions, which helped us improve our paper. We also thank Dr. Zheng Ren and Dr. Marc De Vos for helpful suggestions and for assistance in obtaining the data. The authors would like to express thanks to the UCONN HPC support staff for their efforts.
