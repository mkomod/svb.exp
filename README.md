# Variational Bayes for high-dimensional proportional hazards models

## Install and setup

```
git clone https://github.com/mkomod/svb.exp.git
cd svb.exp
guzip data/*
```

## Directory structure

```
svb.exp/
├── data		# datasets
├── figures		# saved figure (initially empty)
├── R			# R code for experiments
└── RData		# saved results
```

## Simulation study

All simulations can be replicated by running

**Highly recommended that these scripts are ran on a cluster**.

```
Rscript R/simulations/01-mcmc_comparison.R
Rscript R/simulations/02-simulations.R
Rscript R/simulations/03-mcmc_models.R
Rscript R/simulations/04-sensitivity_design.R
Rscript R/simulations/05-sensitivity_params.R
```

All output is saved to `RData/simulations/`, `RData/comparison/` and `RData/sensitivity/`.

Figures and tables can be reproduced by running

```
Rscript R/simulations/06-figures.R
Rscript R/simulations/07-tables.R
```

**Note**: you may need to run `Rscript R/simulations/03-mcmc_models.R` beforehand.

## Application to real data

To reproduce our models run

**Highly recommended that these scripts are ran on a cluster**

```
Rscript R/application/tcga_data.R
Rscript R/application/yau_models.R
```

Prefit models can be found under `RData/models/`

To reproduce tables and figures run

```
Rscript R/application/tcga_tables_figures.R
Rscript R/application/yau_tables_figures.R
```

## Results

All our results are saved in `RData/`.

```
RData/
├── comparison          # comparison to MCMC
├── data                # datasets
├── models              # fit models (real datasets)
├── sensitivity         # sensitivity analysis
└── simulations         # simulation study
```

