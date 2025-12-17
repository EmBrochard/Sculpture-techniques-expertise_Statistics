# Sculpture-techniques-expertise_Statistics
This code was used to carry out the study titled <b> "Carved in stone: Experimental criteria for identifying Paleolithic <i>bas-relief<i> production techniques and sculptors’ expertise" </b>, by Émilie Brochard, Luc Doyon, Lloyd A. Courtenay, Gilles Tosello, Lila Geis and Francesco d'Errico, currently under review at <i> PLoS ONE </i>.

-----------------------------------------------------------------------------------------------------------------

## <b> Author Details </b>

Authors listed in alphabetical order

<b> Author </b>: Émilie Brochard

<b> Email </b>: emilie.brochard@u-bordeaux.fr

<b> ORCID </b>: https://orcid.org/0009-0003-5864-0844

<b> Current Afiliation </b>: University of Bordeaux [CNRS, PACEA UMR5199]

<b> Contributions </b>: Wrote code for "Roughness".

<br>

<b> Author </b>: Lloyd A. Courtenay

<b> Email </b>: ladc1995@gmail.com

<b> ORCID </b>: https://orcid.org/0000-0002-4810-2001

<b> Current Afiliation </b>: University of Bordeaux [CNRS, PACEA UMR5199]

<b> Contributions </b>: Wrote code for "Roughness", "Engraving_Analysis", "Engraving_Extract_Data", "EFA Functions", "Profile_Code".

<br>

<b> Author </b>: Luc Doyon

<b> Email </b>: luc.doyon@u-bordeaux.fr

<b> ORCID </b>: https://orcid.org/0000-0001-7163-6186

<b> Current Afiliation </b>: University of Bordeaux [CNRS, PACEA UMR5199]

<b> Contributions </b>: Wrote code for "Roughness"

---------------------------------------------------------------------------------------------------

This code has been designed for the open-source free R programming languages.

---------------------------------------------------------------------------------------------------

## <b> Repository Details </b>

The present repository contains:

* <b> Codes </b>
   * <b>Source</b>
    * EFA Functions
      * This code implements Elliptic Fourier Analysis (EFA) for outline-based shape analysis, including normalization, harmonic selection, PCA-based visualization of shape variation, and statistical diagnostics.
    * Profile_Code.R
      * This R script provides a set of functions for the extraction, normalization, visualization, and quantitative analysis of 2D profile data, including geometric measurements, asymmetry indices, landmark-based shape characterization, and circular–linear statistical analyses.
      <br>
  * <b>Engraving_Extract_Data.R</b>
    * This R script provides functions for generating profile images, extracting geometric measurements, computing landmarks, and exporting data for morphometric analyses, including writing Morphologika-compatible files and CSV tables.
      <br>
  * <b>Engraving_Analysis.R</b>
    * This main R script implements a complete analytical pipeline for 2D profile data, combining geometric measurements, circular and linear statistics, multivariate analyses, and elliptic Fourier–based shape and form analyses to investigate morphological variability and group differences.
      <br>
  * <b>Roughness.R</b>
    * Code used to .......

--------------------------------------------------------

## <b> System Requirements for Deep Learning </b>

<i> Note that here we specify the versions of the libraries that we used for the present study. We are unaware if earlier or later versions of the same libraries will work or present problems, because we have not tested this, the objective here is simply to state how we ran each of the codes presented </i>

* R - <i> v.4.4.1 </i>
* The following R libraries
  * GraphGMM - <i>v.1.0.0</i>
  * pValueRobust - <i>v.0.1.0</i>
  * geomorph - <i>v.4.0.9</i>
  * shapes - <i>v.1.2.7</i>
  * ggplot2 - <i>v.3.5.1</i>
  * circular - <i>v.0.5.1</i>
  * RVAideMemoire - <i>v.0.9-83-7</i>
  * dplyr - <i>v.1.1.4</i>
  * ggpubr - <i>v.0.6.0</i>
  * writexl - <i>v.1.5.1</i>
  * corrplot - <i>v.0.95</i>
  * FactoMineR - <i>v.2.11</i>
  * factoextra - <i>v.1.0.7</i>
  * caret - <i>v.7.0-1</i>
  * MASS  - <i>v.7.3-60.2</i>
  * Morpho - <i>v.2.12</i>
  
--------------------------------------------------------

## <b> Repository Citation </b>

Please cite this repository as:

 <b>  NAME (2024) Code for ...    . https://github.com/EmBrochard/Sculpture-techniques-expertise_Statistics </b>

--------------------------------------------------------

Comments, questions, doubts, suggestions and corrections can all be directed to É. Brochard at the email provided above.

---------------------------------------------------------------------------------------------------

## License

This project is licensed under the GNU Affero General Public License v3.0.
See the LICENSE file for details.

<b> Copyright (C) 2025 Emilie Brochard, Lloyd Courtenay, Luc Doyon </b>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, version 3.
