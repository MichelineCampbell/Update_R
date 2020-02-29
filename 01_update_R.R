## Script Name: 01_update_R.R
## Purpose: Gets names of all installed packages, concatenates them, installs them into new version of R
## Initial Author: Dr. Micheline Campbell
## Contributing Authors: Authors
## Date Created: 20200229
## Last Edited: 20200229
##  By: MC
## Email: michelineleecampbell@gmail.com
## Notes: Follow these steps to get a shiny new version of R with shiny updated packages. I did write this after a couple of beers, so it may not be bullet proof.

# 1: In Rstudio
ip <- as.data.frame(installed.packages()[,c(1)]) # get installed packages
rownames(ip) <- NULL # no rownames thankyou
dput(as.character(ip[[1]])) # this is awesome and puts ip into concatenated format. COPY THIS OUTPUT

#2: Copy output into: 
packages <- c() # see example at bottom

#3 Now, in base R gui (not RStudio) run: 
if(!require(installr)) {
  install.packages("installr"); require(installr)} #load / install+load installr
updateR()

# 4. Reopen R Studio. Check you are using latest installed version of R
# Now copy your packages object you created in step 2.
#

# packages <-

# 5 run install command below, sit back and watch everything become beautiful and up to date. Note that you're fresh out of luck for anything you haven't sourced from CRAN. 

install.packages(packages)


## example packages object e.g. as at 20200229
# # packages <- c("abc", "abc.data", "abind", "abnormality", "accelerometry", 
#               "acepack", "acss.data", "ada", "additivityTests", "ade4", "adegenet", 
#               "adegraphics", "adehabitatHR", "adehabitatLT", "adehabitatMA", 
#               "adephylo", "ahaz", "airGR", "alabama", "AlgDesign", "animation", 
#               "Anthropometry", "anytime", "ape", "archetypes", "archivist", 
#               "aroma.affymetrix", "aroma.apd", "aroma.core", "assertthat", 
#               "AUC", "audio", "automap", "awsMethods", "backports", "base64enc", 
#               "BB", "Benchmarking", "BH", "bibtex", "biclust", "bigmemory", 
#               "bigmemory.sri", "binaryLogic", "bindr", "bindrcpp", "bio3d", 
#               "Bios2cor", "bit", "bit64", "bitops", "biwavelet", "blob", "bmp", 
#               "bold", "brew", "broom", "BsMD", "Cairo", "callr", "capushe", 
#               "car", "carData", "caret", "caTools", "cellranger", "changepoint", 
#               "checkmate", "CircStats", "circular", "classInt", "cli", "clipr", 
#               "clue", "clues", "clusterGeneration", "clv", "cmprsk", "coda", 
#               "colorspace", "colourpicker", "combinat", "ComICS", "compare", 
#               "coneproj", "conf.design", "corrplot", "cowplot", "cranlogs", 
#               "crayon", "crosstalk", "crul", "cubature", "curl", "cyclocomp", 
#               "daewr", "data.table", "data.tree", "DBI", "dbscan", "deldir", 
#               "DEoptimR", "depth", "desc", "DescTools", "devtools", "DiagrammeR", 
#               "dichromat", "digest", "directlabels", "DoE.base", "doParallel", 
#               "dotCall64", "downloader", "dplR", "dplyr", "drgee", "DT", "dtw", 
#               "dummies", "dummy", "dvmisc", "dygraphs", "dynamicTreeCut", "e1071", 
#               "EbayesThresh", "effects", "elasticnet", "ellipse", "ellipsis", 
#               "emmeans", "emojifont", "emplik", "equivalence", "estimability", 
#               "etm", "evaluate", "exactRankTests", "expint", "expm", "extraDistr", 
#               "extrafont", "extrafontdb", "FactoMineR", "fansi", "fastcluster", 
#               "fastmap", "fastmatch", "fda", "fftwtools", "fields", "filehash", 
#               "filematrix", "fit.models", "flashClust", "flexclust", "flextable", 
#               "flock", "FNN", "forcats", "foreach", "forecast", "formattable", 
#               "Formula", "fracdiff", "fractaldim", "FrF2", "furrr", "future", 
#               "GA", "gamlss", "gamlss.data", "gamlss.dist", "gbRd", "gdata", 
#               "GDINA", "gdtools", "geepack", "genalg", "GenBinomApps", "GENEAread", 
#               "generics", "geometry", "geosphere", "ggdendro", "ggformula", 
#               "GGIR", "ggmap", "ggplot2", "ggpubr", "ggrepel", "ggsci", "ggsignif", 
#               "ggsn", "ggstance", "ggthemes", "ggvis", "git2r", "gld", "glm2", 
#               "globals", "glue", "gmodels", "goftest", "gower", "gplots", "Grid2Polygons", 
#               "gridExtra", "gsl", "gstat", "gtable", "gtools", "HaploSim", 
#               "haven", "hexbin", "highlight", "highr", "Hmisc", "hms", "htmlTable", 
#               "htmltools", "htmlwidgets", "httpcode", "httpuv", "httr", "hydroGOF", 
#               "hydroTSM", "ICGE", "igraph", "imager", "imputeTS", "imputeYn", 
#               "ineq", "influenceR", "inline", "inlmisc", "insol", "installr", 
#               "intervals", "ipred", "irlba", "Iso", "iterators", "itertools", 
#               "ivtools", "jomo", "jpeg", "jsonlite", "Kendall", "kernelFactory", 
#               "klaR", "km.ci", "kml", "KMsurv", "knitr", "labeling", "labelled", 
#               "lars", "later", "latex2exp", "latticeExtra", "lava", "lawstat", 
#               "lazyeval", "ldbounds", "leaflet", "leaps", "LearnBayes", "lhs", 
#               "lifecycle", "linprog", "lintr", "listenv", "lme4", "lmerTest", 
#               "lmom", "lmtest", "locfit", "lomb", "longitudinalData", "lpSolve", 
#               "lpSolveAPI", "lubridate", "magic", "magrittr", "MALDIquant", 
#               "manipulate", "manipulateWidget", "mapdata", "mapproj", "maps", 
#               "maptools", "maptree", "markdown", "matchingMarkets", "matlab", 
#               "matrixcalc", "MatrixModels", "matrixStats", "maxstat", "mcmc", 
#               "MCMCpack", "mda", "memoise", "metap", "MHadaptive", "mice", 
#               "mime", "miniUI", "minpack.lm", "minqa", "misc3d", "miscF", "missForest", 
#               "mitml", "mitools", "mixexp", "mixtools", "mmap", "mnormt", "ModelMetrics", 
#               "modeltools", "moments", "mosaic", "mosaicCore", "mosaicData", 
#               "munsell", "mvtnorm", "natserv", "ncdf4", "neldermead", "networkD3", 
#               "nleqslv", "nloptr", "nnls", "nortest", "numbers", "numDeriv", 
#               "officer", "openssl", "openxlsx", "optimbase", "optimsimplex", 
#               "ordinal", "oz", "padr", "PairedData", "palr", "pan", "PAutilities", 
#               "pbapply", "pbkrtest", "pcaPP", "pdist", "pedigree", "pegas", 
#               "permute", "phangorn", "pheatmap", "phylobase", "phylosignal", 
#               "phytools", "pillar", "pkgbuild", "pkgconfig", "pkgload", "plgp", 
#               "plogr", "plotly", "plotrix", "pls", "plyr", "png", "polyclip", 
#               "polynom", "poweRlaw", "pracma", "praise", "prettyunits", "processx", 
#               "prodlim", "progress", "promises", "prophet", "protiq", "proto", 
#               "proxy", "ps", "PSCBS", "pscl", "psych", "purrr", "quadprog", 
#               "quantmod", "quantreg", "questionr", "QUIC", "R.cache", "R.devices", 
#               "R.filesets", "R.huge", "R.methodsS3", "R.oo", "R.rsp", "R.utils", 
#               "R2admb", "R2jags", "R2WinBUGS", "R6", "RandomFields", "RandomFieldsUtils", 
#               "randomForest", "randomNames", "randtoolbox", "ranger", "rappdirs", 
#               "rARPACK", "raster", "rasterImage", "rbenchmark", "RColorBrewer", 
#               "Rcpp", "RcppArmadillo", "RcppEigen", "RcppProgress", "RCurl", 
#               "Rdpack", "readbitmap", "readODS", "readr", "readxl", "recipes", 
#               "RefManageR", "rematch", "remotes", "rentrez", "repr", "reshape", 
#               "reshape2", "reticulate", "rex", "rgdal", "rgeos", "rgexf", "rgl", 
#               "Rglpk", "RgoogleMaps", "rintrojs", "rio", "ritis", "rjags", 
#               "rJava", "rjson", "RJSONIO", "rlang", "RMariaDB", "rmarkdown", 
#               "rncl", "RNeXML", "rngWELL", "robust", "robustbase", "ROCR", 
#               "Rook", "rootSolve", "rotl", "rpanel", "rprojroot", "rrcov", 
#               "rredlist", "Rsolnp", "RSpectra", "RSQLite", "rstan", "rstudioapi", 
#               "Rttf2pt1", "rvest", "rworldmap", "sandwich", "scales", "scatterplot3d", 
#               "scico", "SDMTools", "seewave", "segmented", "selectr", "SemiPar", 
#               "seqinr", "sf", "sfsmisc", "shape", "shapes", "shiny", "shinyalert", 
#               "shinyBS", "shinydashboard", "shinyjs", "shinythemes", "shinyWidgets", 
#               "showtext", "showtextdb", "signal", "slam", "sn", "solrium", 
#               "sourcetools", "sp", "spacetime", "spam", "sparseLDA", "SparseM", 
#               "spatstat", "spatstat.data", "spatstat.utils", "spData", "spdep", 
#               "spectral", "SQUAREM", "stabledist", "StanHeaders", "stdReg", 
#               "stinepack", "stringdist", "stringi", "stringr", "SuppDists", 
#               "survey", "survminer", "survMisc", "sysfonts", "systemfonts", 
#               "taxize", "tcltk2", "tensor", "testthat", "tgp", "tibble", "tidyr", 
#               "tidyselect", "tiff", "timeDate", "timeline", "timelineS", "timetk", 
#               "tinytex", "tkrplot", "tolerance", "toOrdinal", "treeClust", 
#               "triebeard", "truncnorm", "TSA", "tseries", "TTR", "tuneR", "turboEM", 
#               "ucminf", "units", "univOutl", "unmarked", "urca", "urltools", 
#               "uroot", "utf8", "uuid", "V8", "vcd", "vctrs", "vegan", "vegawidget", 
#               "verification", "versions", "VGAM", "viridis", "viridisLite", 
#               "visNetwork", "wavethresh", "webshot", "wellknown", "wesanderson", 
#               "WGCNA", "whisker", "WikidataR", "WikipediR", "wikitaxa", "withr", 
#               "worrms", "writexl", "xfun", "xlsx", "xlsxjars", "XML", "xml2", 
#               "xmlparsedata", "xtable", "xts", "yaImpute", "yaml", "zeallot", 
#               "zip", "zoo", "base", "boot", "class", "cluster", "codetools", 
#               "compiler", "datasets", "foreign", "graphics", "grDevices", "grid", 
#               "KernSmooth", "lattice", "MASS", "Matrix", "methods", "mgcv", 
#               "nlme", "nnet", "parallel", "rpart", "spatial", "splines", "stats", 
#               "stats4", "survival", "tcltk", "tools", "translations", "utils"
# )
