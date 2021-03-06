Replication code for Gneiting and Resin (2021)
================

### This repository provides replication code for the article “Regression Diagnostics meets Forecast Evaluation: Conditional Calibration, Reliability Diagrams, and Coefficient of Determination” by Tilmann Gneiting and Johannes Resin ([arXiv:2108.03210](https://arxiv.org/abs/2108.03210v2)).

## Required packages

``` r
library(openxlsx)
library(RColorBrewer)
library(GoFKernel)
library(isotone)
library(quantreg)
```

## (Down)load the data

``` r
source("load_data.R")
```

Download and load the data. This will save some files in the working
directory.

The following data is download:

  - Data for the butterfly population size example from Tredennick et
    al. (2021) available at <https://doi.org/10.5281/zenodo.4311358>
  - The Bank of England forecasts of CPI inflation rates available at
    <https://www.bankofengland.co.uk/-/media/boe/files/monetary-policy-report/2021/february/mpr-february-2021-chart-slides-and-data.zip>
  - Historical quarterly CPI inflation rates published by the UK Office
    for National Statistics (ONS) available at
    <https://www.ons.gov.uk/economy/inflationandpriceindices/timeseries/d7g7>

## Load necessary functions and variables

``` r
source("functions.R")
```

## Plot figures from main text and supplement

``` r
source("figures.R")
```

A new folder ‘figs’ is created in the working directory and all figures
are saved there as png files.

## Figure 2

<p align="middle">

<img src="figs/fig2a.png" width="400" />
<img src="figs/fig2b.png" width="400" />

</p>

## Figure 3

<p align="middle">

<img src="figs/fig3.png" width="525" />

</p>

## Figure 4

<p align="middle">

<img src="figs/fig4a.png" width="315" />
<img src="figs/fig4b.png" width="270" />
<img src="figs/fig4c.png" width="270" />

</p>

<p align="middle">

<img src="figs/fig4d.png" width="315" />
<img src="figs/fig4e.png" width="270" />
<img src="figs/fig4f.png" width="270" />

</p>

<p align="middle">

<img src="figs/fig4g.png" width="315" />
<img src="figs/fig4h.png" width="270" />
<img src="figs/fig4i.png" width="270" />

</p>

## Figure 5

<p align="middle">

<img src="figs/fig5.png" width="525" />

</p>

## Figure 6

<p align="middle">

<img src="figs/fig6a.png" width="315" />
<img src="figs/fig6b.png" width="270" />
<img src="figs/fig6c.png" width="270" />

</p>

<p align="middle">

<img src="figs/fig6d.png" width="315" />
<img src="figs/fig6e.png" width="270" />
<img src="figs/fig6f.png" width="270" />

</p>

<p align="middle">

<img src="figs/fig6g.png" width="315" />
<img src="figs/fig6h.png" width="270" />
<img src="figs/fig6i.png" width="270" />

</p>

## Figure 7

<p align="middle">

<img src="figs/fig7.png" width="500" />

</p>

## Figure 8

<p align="middle">

<img src="figs/fig8a.png" width="270" />
<img src="figs/fig8b.png" width="270" />
<img src="figs/fig8c.png" width="270" />

</p>

<p align="middle">

<img src="figs/fig8d.png" width="270" />
<img src="figs/fig8e.png" width="270" />
<img src="figs/fig8f.png" width="270" />

</p>

## Figure 9

<p align="middle">

<img src="figs/fig9a.png" width="400" />
<img src="figs/fig9b.png" width="400" />

</p>

## Figure S1

<p align="middle">

<img src="figs/figS1a.png" width="315" />
<img src="figs/figS1b.png" width="270" />
<img src="figs/figS1c.png" width="270" />

</p>

## Figure S2

<p align="middle">

<img src="figs/figS2a.png" width="315" />
<img src="figs/figS2b.png" width="270" />
<img src="figs/figS2c.png" width="270" />

</p>

<p align="middle">

<img src="figs/figS2d.png" width="315" />
<img src="figs/figS2e.png" width="270" />
<img src="figs/figS2f.png" width="270" />

</p>

## Figure S3

<p align="middle">

<img src="figs/figS3a.png" width="270" />
<img src="figs/figS3b.png" width="270" />
<img src="figs/figS3c.png" width="270" />

</p>

<p align="middle">

<img src="figs/figS3d.png" width="270" />
<img src="figs/figS3e.png" width="270" />
<img src="figs/figS3f.png" width="270" />

</p>

## Figure S4

<p align="middle">

<img src="figs/figS4a.png" width="270" />
<img src="figs/figS4b.png" width="270" />
<img src="figs/figS4c.png" width="270" />

</p>

<p align="middle">

<img src="figs/figS4d.png" width="270" />
<img src="figs/figS4e.png" width="270" />
<img src="figs/figS4f.png" width="270" />

</p>

## Figure S5

<p align="middle">

<img src="figs/figS5a.png" width="270" />
<img src="figs/figS5b.png" width="270" />
<img src="figs/figS5c.png" width="270" />

</p>

<p align="middle">

<img src="figs/figS5d.png" width="270" />
<img src="figs/figS5e.png" width="270" />
<img src="figs/figS5f.png" width="270" />

</p>

## Figure S6

<p align="middle">

<img src="figs/figS6a.png" width="270" />
<img src="figs/figS6b.png" width="270" />
<img src="figs/figS6c.png" width="270" />

</p>

<p align="middle">

<img src="figs/figS6d.png" width="270" />
<img src="figs/figS6e.png" width="270" />
<img src="figs/figS6f.png" width="270" />

</p>

## Figure S7

<p align="middle">

<img src="figs/figS7a.png" width="270" />
<img src="figs/figS7b.png" width="270" />
<img src="figs/figS7c.png" width="270" />

</p>

<p align="middle">

<img src="figs/figS7d.png" width="270" />
<img src="figs/figS7e.png" width="270" />
<img src="figs/figS7f.png" width="270" />

</p>

## Figure S8

<p align="middle">

<img src="figs/figS8a.png" width="270" />
<img src="figs/figS8b.png" width="270" />
<img src="figs/figS8c.png" width="270" />

</p>

<p align="middle">

<img src="figs/figS8d.png" width="270" />
<img src="figs/figS8e.png" width="270" />
<img src="figs/figS8f.png" width="270" />

</p>
