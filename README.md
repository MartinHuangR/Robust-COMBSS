# Robust Best Subset Selection via Fast Approximate MM-Estimation

This repository contains code to reproduce results as given in the paper "Robust Best Subset Selection via Fast Approximate MM-Estimation". For any questions please email martin.huang@sydney.edu.au.

# Simulations
The repository is split into three folders, "RMS", "COMBSS", and "NBA", corresponding to sections 4.1, 4.2 and 4.3 respectively. The row-wise contamination models results are under "MSSim.R" and "COMBSS.R". The response contamination model results are under "MSSimC-NS.R" and "COMBSSC-NS.R". Their respective code to produce figures fall under "MSPlots.R", "COMBSSPlots.R", "MSPlotsC-NS.R", and "COMBSSC-NSPlots.R" respectively.
\\
The Moka style benchmark is found under the file COMBSSBenchmark.R and COMBSSBenchmarkPlots.R. 
\\
The files "Functions.R" and "CppFuncs.cpp" contain functions for all three folders. We have also included the simulation results in .RData files. The file paths assumes "Robust COMBSS" is the home working directory.

# Data Availability
The raw NBA data csv files is given in the "NBA" folder, and the code to process the data into the used dataset is provided in NBA/NBA.R. It can also be extracted as a csv file given as NBAOutlier.csv.

# Package References
This analysis could not have been done without the following packages.

* Eddelbuettel D, Francois R, Allaire J, Ushey K, Kou Q, Russell N, Ucar I, Bates D, Chambers J (2025). Rcpp: Seamless R and C++ Integration. R package version 1.0.14, https://CRAN.R-project.org/package=Rcpp.
* Kepplinger D, Salibián-Barrera M, Cohen Freue G (2024). pense: Penalized Elastic Net S/MM-Estimator of Regression. R package version 2.2.2 https://CRAN.R-project.org/package=pense.
* Martin Maechler, Peter Rousseeuw, Christophe Croux, Valentin Todorov, Andreas Ruckstuhl, Matias Salibian-Barrera, Tobias Verbeke, Manuel Koller, c("Eduardo", "L. T.") Conceicao and Maria Anna di Palma (2024). robustbase: Basic Robust Statistics R package version 0.99-4-1. URL http://CRAN.R-project.org/package=robustbase
* Pedersen T (2024). patchwork: The Composer of Plots. R package version 1.3.0, https://CRAN.R-project.org/package=patchwork.
* Salibian-Barrera M (2025). FRBmodelselection: Model Selection with the Fast and Robust Bootstrap. R package version 0.0.0.9000, commit 5b1e977ed6c7bd5e8cecbb6b8ca9206a8314f009, https://github.com/msalibian/FRBmodelselection.
* van den Brand T (2025). ggh4x: Hacks for 'ggplot2'. R package version 0.3.1, https://CRAN.R-project.org/package=ggh4x.
* Van Aelst S, Willems G (2013). "Fast and Robust Bootstrap for Multivariate Inference: The R Package FRB." Journal of Statistical Software, 53(3), 1-32. doi:10.18637/jss.v053.i03 https://doi.org/10.18637/jss.v053.i03, https://www.jstatsoft.org/v53/i03/.
* Vaughan D, Dancho M (2022). furrr: Apply Mapping Functions in Parallel using Futures. R package version 0.3.1, https://CRAN.R-project.org/package=furrr.
* Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0
* Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). "Welcome to the tidyverse." Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686 https://doi.org/10.21105/joss.01686.
* Wickham H, Henry L (2025). purrr: Functional Programming Tools. R package version 1.0.4, https://CRAN.R-project.org/package=purrr.