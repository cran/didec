
<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# didec 1.1.0

## Main changes

-   Add a new method for estimating the directed dependence coefficient (didec, mfoci, VarClustPartition)
-   Add the data standardization procedure for predictor variables (didec, mfoci, VarClustPartition)
-   Add a new method for selecting features (mfoci)
-   Add new methods for permuting the response variables (mfoci, VarClustPartition)

## Minor changes

-   Change the return of the function into a list (mfoci)
-   Add automatic removing missing values (NA) (didec, mfoci, VarClustPartition)
-   Update the bioclimatic dataset with 19 variables 
-   Add warning when variable is constant (didec, mfoci, VarClustPartition)
-   Use function cor.fk in R package pcaPP to calculate bivariate kendall tau (VarClustPartition)
-   Change terms of the return of the function (mfoci)
-   Fix typos in the manual