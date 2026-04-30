
<!-- README.md is generated from README.Rmd. Please edit that file -->

# intCPUE <img src="man/figures/logo.png" align="right" height="136" alt="intCPUE logo" />

> Integrated CPUE standardization with spatiotemporal models in TMB

<!-- badges: start -->

[![R-CMD-check](https://github.com/RujiaBi/intCPUE/workflows/R-CMD-check/badge.svg)](https://github.com/RujiaBi/intCPUE/actions)
<!-- badges: end -->

**intCPUE** is a TMB-based framework for integrated CPUE standardization
across multiple fisheries or surveys, with optional preferential
sampling correction (under development).

## Table of Contents

- [Overview](#overview)
- [Contact](#contact)
- [Citation](#citation)
- [Installation](#installation)
- [Data structure](#data-structure)
- [Coordinate projection](#coordinate-projection)
- [Build spatial mesh](#build-spatial-mesh)
- [Fit the model](#fit-the-model)
- [Getting index with bias
  correction](#getting-index-with-bias-correction)
- [Next steps](#next-steps)

## Overview

The model supports:

- spatiotemporal random fields via an SPDE mesh (using `fmesher`)

- multiple fleets / surveys through catchability components

- mgcv-style `s()` smooth terms parsed into mixed-effects from inside
  TMB

The package supports three observation-data structures through
`data_type = c("auto", "delta", "positive", "mixed")`:

- `"auto"` (default): detect the data structure by `flagid`

- `"delta"`: every flag/type has both zero and positive catches

- `"positive"`: every flag/type is positive-only

- `"mixed"`: some flags/types are delta and some are positive-only

For `data_type = "delta"`, `delta_link = c("Poisson", "traditional")`
chooses between the Poisson-link delta model and a traditional
delta-lognormal model. If any positive-only flag/type is present, the
model uses the traditional positive-lognormal interpretation for those
observations. For positive-only data, `delta_link` is ignored because
there is no encounter submodel.

The baseline population density can include:

- temporal effect (`tid`)

- spatial random field (`omega`)

- spatiotemporal random field (`epsilon`)

- population factor effects

- population smooth effects

In the current package design:

- `tid` is always present

- the default `tid` effect is a fixed effect

- if `formula_population` contains `f(tid)`, the default fixed `tid`
  effect is replaced by `intercept + random tid effect`

- `omega` and `epsilon` are controlled by `pop_spatial` and
  `pop_spatiotemporal`; wrapper templates remove disabled SPDE modules
  from the compiled objective

- vessel effects are part of the core catchability model and are always
  included; if `n_v = 1`, the vessel effect is fixed to 0

- systematic fishery differences (`flag_f`) are part of the core model
  and can be `q_diffs_system_type = "random"` or `"fixed"`

Smooth terms can be assigned either to catchability only or to the
population layer. Population-layer smooths and factor effects are
included both in the observation model and in projection, so they affect
the standardized index. Catchability-layer smooths and factor effects
affect the observation model only.

Here provides a minimal workflow: install → data preparation →
coordinate projection → mesh → model fitting → index extraction.

## Contact

For questions, suggestions, or collaboration, please contact:

- **Rujia Bi** — <rbi@iattc.org>

## Citation

If you use **intCPUE** in your work, please cite it as software:

> Bi, R. (2026). *intCPUE: Integrated CPUE standardization with TMB*. R
> package (v0.1.0). <https://github.com/RujiaBi/intCPUE>.

Once a paper or DOI is available, this section will be updated.

## Installation (development version)

``` r
# install.packages("remotes")
remotes::install_github("RujiaBi/intCPUE")
```

``` r
library(intCPUE)
```

## Data structure

An intCPUE model requires a data frame containing the following columns:

- `cpue` — positive catch rate

- `encounter` — encounter indicator (0 = zero catch, 1 = positive catch)

- `lon`, `lat` — geographic coordinates

- `vesid` — vessel ID (0-based)

- `tid` — time index (0-based)

- `flagid` — fishery / survey ID (0-based; 0 = reference fishery)

These columns are required. Additional columns may be included as
needed.

``` r
# Example structure based on a two-fishery YFT CPUE data set.
# Replace column names on the right-hand side with the names in your data.
combined_data[combined_data$flag == "JPN", "flagid"] <- 0
combined_data[combined_data$flag == "PS", "flagid"] <- 1
combined_data$vesid <- as.numeric(as.factor(as.character(combined_data$vesno))) - 1
combined_data$tid <- combined_data$yq - min(combined_data$yq)

data_input <- data.frame(
  cpue = combined_data$cpue,
  encounter = combined_data$encounter,
  lon = combined_data$lonc1,
  lat = combined_data$latc1,
  vesid = combined_data$vesid,
  tid = combined_data$tid,
  flagid = combined_data$flagid,
  hbf = combined_data$hbf
)
```

### Important

- `vesid`, `tid` and `flagid` must be 0-based contiguous integers

- `flagid` must use 0 as the reference fishery / survey

## Coordinate projection (lon/lat → UTM)

Longitude may be in either -180..180 or 0..360.

`make_utm()` automatically:

- detects longitude convention

- selects an appropriate UTM zone

- scales coordinates for numerical stability

``` r
utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
data_utm <- utm$data_utm
```

## Build spatial mesh

The mesh must be constructed using the scaled projected coordinates:
`utm_x_scale` and `utm_y_scale`

### K-means mesh

``` r
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), type = "kmeans", n_knots = 50)
plot(mesh)
```

### Cutoff mesh

``` r
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), type = "cutoff", cutoff = 0.1)
plot(mesh)
```

### Tailor mesh

``` r
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), type = "tailored",
    convex = -0.1,         # for a finer boundary
    max.edge = c(0.5, 2),   # max triangle edge length; inner and outer meshes
    offset = c(0.1, 0.5),  # inner and outer border widths
    cutoff = 0.05)
plot(mesh)
```

### Custom mesh

``` r
bnd <- INLA::inla.nonconvex.hull(cbind(data_utm$utm_x_scale, data_utm$utm_y_scale), convex = -0.1)
mesh_inla <- INLA::inla.mesh.2d(
  boundary = bnd,
  max.edge = c(0.5, 2)
)
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), mesh = mesh_inla)
plot(mesh)
```

## Fit the model

``` r
formula_catchability <- ~ s(hbf)
formula_population <- ~ s(temp) + s(chl)
```

`intCPUE()` uses separate one-sided formulas for catchability and
population covariates. The response is not written in either formula;
`cpue` and `encounter` are read directly from `data_utm`.

Full interface:

``` r
fit <- intCPUE(
  data_utm,
  mesh,
  formula_catchability = NULL,
  formula_population = NULL,
  projection_data = NULL,
  data_type = c("auto", "delta", "positive", "mixed"),
  delta_link = c("Poisson", "traditional"),
  pop_spatial = c("on", "off"),
  pop_spatiotemporal = c("on", "off"),
  pop_spatiotemporal_type = c("rw", "ar1"),
  q_diffs_system_type = c("random", "fixed"),
  q_diffs_time = c("on", "off"),
  q_diffs_spatial = c("on", "off"),
  obs_sd = c("shared", "flag"),
  control = list(eval.max = 1e5, iter.max = 1e5),
  ncores = NULL,
  ...,
  silent = FALSE,
  restart_max = 1L,
  newton_max = 2L,
  coord_max = 5L
)
```

Parameter guide:

- `data_utm`: output from `make_utm()`; must contain `cpue`,
  `encounter`, `lon`, `lat`, `vesid`, `tid`, `flagid`, `utm_x_scale`,
  and `utm_y_scale`.
- `mesh`: an `intCPUEmesh` from `make_mesh()` or a compatible `fmesher`
  mesh.
- `formula_catchability`: one-sided formula for covariates that affect
  catchability/observation process only, for example `~ s(hbf)` or
  `~ f(gear_group) + s(depth)`.
- `formula_population`: one-sided formula for covariates that affect the
  latent population density surface, for example `~ s(temp) + s(chl)` or
  `~ f(tid) + s(temp)`.
- `projection_data`: optional covariate table on the extrapolation grid
  for variables used in `formula_population`.
- `data_type`: observation-data structure. `"auto"` detects it by
  `flagid`; `"delta"` means every flag has zeros and positives;
  `"positive"` means all observations are positive-only; `"mixed"` means
  some flags are delta and some are positive-only.
- `delta_link`: delta-model link. `"Poisson"` is only allowed for pure
  `data_type = "delta"`. `"traditional"` is required for
  `data_type = "mixed"`. For `data_type = "positive"`, `delta_link` is
  ignored because only the positive lognormal component is fitted.
- `pop_spatial`: turns the time-constant population spatial field
  (`omega`) `"on"` or `"off"`.
- `pop_spatiotemporal`: turns the population spatiotemporal field
  (`epsilon`) `"on"` or `"off"`.
- `pop_spatiotemporal_type`: temporal structure for `epsilon`; `"rw"` is
  a random walk and `"ar1"` is autoregressive.
- `q_diffs_system_type`: controls only the main flag/system effect
  (`flag_f`); `"random"` treats it as a random effect and `"fixed"`
  treats it as an ordinary fixed effect.
- `q_diffs_time`: turns flag-specific temporal catchability differences
  (`flag_t`) `"on"` or `"off"`.
- `q_diffs_spatial`: turns flag-specific spatial catchability
  differences (`flag_s`) `"on"` or `"off"`.
- `obs_sd`: `"shared"` uses one lognormal SD for positive catches;
  `"flag"` estimates one lognormal SD per flag.
- `control`: passed to `stats::nlminb()`.
- `ncores`: optional number of OpenMP threads for TMB.
- `...`: passed to `make_data()`, for example `area_scale`.
- `silent`: passed to `TMB::MakeADFun()`.
- `restart_max`, `newton_max`, `coord_max`: post-fit polishing controls.
  Set all three to `0` to run only the initial `nlminb()` optimization.

Formula rules:

- formulas must be one-sided, for example
  `formula_catchability = ~ s(hbf)`;
- use `s(x)` for smooth terms, `x` for a fixed factor effect, and `f(x)`
  for a random factor effect;
- do not write `cpue ~ ...` because the response is already supplied
  through `data_utm$cpue` and `data_utm$encounter`;
- vessel effects do not need to be written as `f(vesid)` because they
  are always included in the core catchability model;
- if `n_v = 1`, the vessel effect is fixed to 0;
- `formula_population = ~ f(tid) + ...` replaces the default fixed `tid`
  effect with `intercept + random tid effect`.

Factor terms:

- `formula_catchability = ~ season` adds `season` as a fixed factor
  effect on catchability;
- `formula_catchability = ~ f(season)` adds `season` as a random factor
  effect on catchability, with an estimated standard deviation;
- the same syntax can be used in `formula_population`, where factor
  effects enter both the observation model and projection;
- fixed and random factor terms must be simple variable names. Inline
  interactions such as `season:flag_factor` or `f(season:flag_factor)`
  are not supported by the parser;
- if a factor effect should differ by flag, create the interaction as a
  new factor column first. For example:
  `data_input$season_flag <- interaction(data_input$season, data_input$flagid, drop = TRUE)`,
  then use `formula_catchability = ~ season_flag` for a fixed
  interaction or `formula_catchability = ~ f(season_flag)` for a random
  interaction;
- do not include `vesid` in either formula because vessel effects are
  included by default in the catchability model;
- do not include `tid` as a fixed population factor because fixed `tid`
  is already included by default. Use `formula_population = ~ f(tid)`
  only when you want to replace the default fixed `tid` effect with
  `intercept + random tid effect`.

Smooth terms with missing covariates:

- `intCPUE` evaluates each smooth term on the rows where the variables
  required by that smooth are non-missing;
- rows with `NA` for a smooth covariate are expanded back into the
  design matrix with zeros, so the contribution of that smooth is 0 for
  those rows;
- this is useful when a covariate is available for only some flags. For
  example, if `hbf` is not defined for one flag, set it to `NA` for that
  flag and use `formula_catchability = ~ s(hbf)`. For that flag, the
  model uses `s(hbf) = 0`;
- this does not create a different smooth curve for each flag. It
  creates one shared smooth curve that is active only for rows with
  non-missing `hbf`;
- if the smooth relationship should differ by flag, use a factor `by`
  smooth. Prefer creating the factor explicitly, for example
  `data_input$flag_factor <- factor(data_input$flagid)` and then
  `formula_catchability = ~ s(hbf, by = flag_factor)`. In that case,
  each flag-specific smooth needs enough non-missing observations to
  estimate its curve. Avoid using numeric `by = flagid`, which would
  create a numeric varying-coefficient smooth rather than separate
  flag-specific smooths.

Mixed data handling:

- `data_type = "mixed"` must use `delta_link = "traditional"`;
  `delta_link = "Poisson"` is only valid for pure `data_type = "delta"`;
- the reference flag (`flagid = 0`) must be a delta flag in mixed
  models;
- delta flags contribute both encounter and positive-catch likelihoods;
- positive-only flags have encounter probability fixed at 1 and
  contribute only the positive-catch likelihood;
- `ves_v_1` is compressed to vessels that appear in delta observations;
- `flag_f_1`, `flag_t_1`, and `flag_s_1` are compressed to delta
  non-reference flags only;
- component-1 smooth and factor design matrices are built separately in
  `make_data()`;
- a smooth/factor term supported only by positive-only observations is
  dropped from component 1 but can remain in component 2;
- component 2 keeps the full positive-catch smooth/factor design, so if
  component 1 has no `s(x)` but component 2 has `s(x)`,
  `ln_smooth_sigma_*` in the shared smooth block refers to component 2.

Interpretation of the baseline model:

- `data_type = "delta"` means every flag/type has zero and positive
  catches; it can use `delta_link = "Poisson"` or
  `delta_link = "traditional"`

- `data_type = "positive"` fits only the positive-lognormal component,
  ignores `delta_link`, and projection uses `exp(eta2_proj)`

- `data_type = "mixed"` uses the traditional delta-lognormal model;
  delta flags contribute encounter + positive likelihood, while
  positive-only flags contribute only positive likelihood

- the population-density linear predictor is:
  `temporal effect (tid) + omega + epsilon + population factor effects + population smooth effects`

- `tid` is always included; `omega` and `epsilon` are included when
  `pop_spatial = "on"` and `pop_spatiotemporal = "on"`

- by default, `tid` is a fixed effect

- `formula_population = ~ f(tid) + ...` switches the model to
  `intercept + random tid effect` instead of the default fixed `tid`

- vessel effects are always included in the catchability model; if
  `n_v = 1`, the vessel effect is fixed to 0

- `q_diffs_system_type = "random"` treats the main flag/system effect
  (`flag_f`) as a random effect; `"fixed"` estimates `flag_f` as
  ordinary fixed effects

- `formula_catchability` adds catchability-side smooth and factor
  effects only

- `formula_population` adds population-side smooth and factor effects,
  and these enter both the observation model and projection

When population covariates are not uniquely defined within each
extrapolation cell, provide them explicitly through `projection_data`.
If they also vary over time, include a `tid` column in `projection_data`
so the projection covariates are matched by grid cell and time.

### `projection_data` format

If `formula_population` is used, population covariates must be available
on the extrapolation grid:

- Static population covariates: `projection_data` should contain one row
  per extrapolation grid cell, with columns `utm_x_scale`,
  `utm_y_scale`, and the covariates used in `formula_population`.

- Time-varying population covariates: `projection_data` should contain
  one row per grid cell-time combination, with columns `utm_x_scale`,
  `utm_y_scale`, `tid`, and the covariates used in `formula_population`.

- If `formula_population` mixes static and time-varying covariates, use
  the grid cell-time format and repeat the static covariate values
  across `tid`.

- Optional `region` can be supplied in `data_utm` (for example via
  `data_input` before `make_utm()`). If no `region` column is supplied,
  `get_index()` returns only the total index over the full extrapolation
  grid. If `region` is supplied, each extrapolation grid cell must
  belong to exactly one region, and `get_index(..., regions = "all")`
  returns each region-specific index plus the total index.
  Region-specific indices and the total index are handled as one joint
  index vector for bias correction.

`tid` must use the same 0-based coding as `data_utm$tid`.

Region example:

``` r
data_input$region <- combined_data$region
utm <- make_utm(data_input, utm_zone = NULL, hemisphere = "auto", coord_scale = "auto")
data_utm <- utm$data_utm
```

``` r
ncores <- 4
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), type = "cutoff", cutoff = 0.1)
fit <- intCPUE(
  data_utm = data_utm,
  mesh = mesh,
  formula_catchability = formula_catchability,
  formula_population = formula_population,
  data_type = "auto",
  delta_link = "Poisson",
  pop_spatial = "on",
  pop_spatiotemporal = "on",
  pop_spatiotemporal_type = "rw",
  q_diffs_system_type = "random",
  q_diffs_time = "off",  
  q_diffs_spatial = "off",
  obs_sd = "shared",
  ncores = ncores
)
```

After fitting:

``` r
# Without a region request, this returns the total index.
index_total <- get_index(fit)

# With regions = "all", this returns region-specific indices plus total.
index_region <- get_index(fit, regions = "all")
```

### Core model components

- `data_type` — `"auto"`, `"delta"`, `"positive"`, or `"mixed"`

- `delta_link` — `"Poisson"` or `"traditional"`; `"Poisson"` is only
  available for pure `data_type = "delta"`

- `pop_spatiotemporal_type` — temporal dependence for the population
  spatiotemporal field (`"rw"` or `"ar1"`)

The current package version always includes:

- temporal effect (`tid`)

- vessel effect in catchability

- systematic fishery catchability differences

Optional core components:

- population spatial random field (`omega`), controlled by `pop_spatial`

- population spatiotemporal random field (`epsilon`), controlled by
  `pop_spatiotemporal`

By default:

- `tid` is a fixed effect

- if `formula_population` contains `f(tid)`, the model uses
  `intercept + random tid effect`

- if `n_v = 1`, the vessel effect is fixed to 0

### User-facing switches

- `q_diffs_system_type` controls whether the main flag/system effect
  (`flag_f`) is `"random"` or `"fixed"`

- `q_diffs_time` — time-varying catchability difference

- `q_diffs_spatial` — spatial catchability difference

For the reference fishery (`flagid = 0`), the `q_diffs_*` terms are
constrained to 0.

### Observation error

- `obs_sd = "shared"` uses one lognormal observation SD across all flags

- `obs_sd = "flag"` estimates one lognormal observation SD for each flag

## Getting index with bias correction

``` r
index <- get_index(fit)
index_by_region <- get_index(fit, regions = "all")
plot_index(index)
```

## Diagnostics

``` r
check_convergence(fit)
calc_marginal_aic(fit)

pred <- get_predicted(fit, data = data_utm)
plots <- plot_residuals(pred, observed_col = "cpue")

plots$observed_predicted
plots$spatial_residual

plot_anisotropy(fit)
```

## Next steps

- preferential sampling correction

- length frequency part
