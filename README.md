# Variational Bayes for high-dimensional proportional hazards models

## Install and setup

```
git clone https://github.com/mkomod/svb.exp.git
```

```
cd svb.exp
guzip data/*
```

## Directory structure

```
svb.exp/
├── data		# compressed datasets
├── figures		# saved figure (initially empty)
├── logs		# logs for scripts
├── R			# R code for experiments
├── RData		# saved results
└── scripts		# cluster scripts
```

## Results

We have saved all our results under `RData/`. These can be used to generate the figures and table from our paper. The following sections describe how our results can be reproduced

```
RData/
├── comparison		# comparison to MCMC
├── data                # real datasets
├── models              # models fit to real datasets
├── sensitivity		# sensitivity analysis
└── simulations		# simulation study, comparison to other methods
```

## Simulation study

All simulations can be replicated by running

```
Rscript R/simulations/01-mcmc_comparison.R
Rscript R/simulations/02-simulations.R
Rscript R/simulations/03-mcmc_models.R
Rscript R/simulations/04-sensitivity_design.R
Rscript R/simulations/05-sensitivity_params.R
```

**Highly recommended that these scripts on a cluster**. Accompanying PBS scripts can be found under `scripts/`

Figures and tables can be reproduced by running

```
Rscript R/simulations/06-figures.R
Rscript R/simulations/07-tables.R
```

**Note**: you may need to run `Rscript R/simulations/03-mcmc_models.R` beforehand.

## Application to real data

To reproduce our models run

```
Rscript R/application/tcga_data.R
Rscript R/application/yau_models.R
```

**Highly recommended that scripts on a cluster**

To reproduce tables and figures run

```
Rscript R/application/tcga_tables_figures.R
Rscript R/application/yau_tables_figures.R
```


