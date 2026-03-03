# SimuRg

## Overview

**SimuRg** provides a comprehensive workflow for non-linear mixed-effects model
development in pharmacometrics. The package provides the entire modeling pipeline:
from model calibration with Monolix fitter(2023) and output processing to goodness-of-fit 
visualization, simulation, and sensitivity analysis. To use Monolix, it should 
be installed. 

Key features:
- **Model calibration** via the Monolix fitter
- **Automated output conversion** into a generalized fit output
- **Diagnostic visualization** for model assessment
- **Simulation capabilities** for model predictions
- **Sensitivity analysis** tools for parameter exploration

## Installation

### From CRAN
```r
install.packages("SimuRg")
```

## Illustrated example

First of all, the model should be calibrated with Monolix fitter. For this goal,
Monolix should be installed on the computer. As this software have commercial 
license, we start our example with the conversion from the Monolix output files 
into the generalized fit output.

```r
library("SimuRg")
library(stringr)
# Convert Monolix project results
test_folder <- system.file("extdata", "Monolix_objects", package = "SimuRg")
if (substr(test_folder, nchar(test_folder), nchar(test_folder)) != "/")
   test_folder <- str_c(test_folder, "/")
pro_name <- "proj-solo"
result <- sg_converter(folder_path = test_folder, proj_name = pro_name)

# Running goodness-of-fit objects
sg_gof_obpr(result)
sg_gof_res(
  fpath_i = result,
  res_type = "IWRES",
  vs_time = TRUE
)
sg_gof_par_cov(
  fpath_i = result,
  ptype = "IndParvsCov",
  cont_cov = cont_cov,
  cat_cov = cat_cov
)
sg_gof_par_dist(fpath_i = result)
sg_gof_res_dist(fpath_i = result, res_type = "IWRES")
sg_gof_res(
  fpath_i = result,
  res_type = "IWRES",
  vs_time = TRUE
)
```

