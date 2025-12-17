# Code written by Lloyd Austin Courtenay

# Copyright (C) 2025 Lloyd Courtenay
# SPDX-License-Identifier: AGPL-3.0

# load libraries and functions ------------

library(GraphGMM) # Lloyd Courtenay Library
library(pValueRobust) # Lloyd Courtenay Library
library(geomorph) # for geometric morphometrics
library(shapes) # for geometric morphometrics
library(ggplot2) # for plotting
library(circular) # for processing angles
library(RVAideMemoire) # for MANOVA

source(".\\Source\\EFA Functions.R")
source(".\\Source\\Profile_Code.R")

#

# load data ------------------

# write here what the separator is
# "" is a tab
# "," is a comma
# ";" is a semi-colon
# " " is a space

separator = "\t"

results <- extract_measurements("Gravures_fines", sep = separator)
measurements <- results$metric_results
landmarks <- results$landmarks

measurements$D <- as.numeric(measurements$D)
measurements$WIS <- as.numeric(measurements$WIS)
measurements$A <- as.numeric(measurements$AOA)
measurements$Sample <- as.factor(measurements$Sample)
measurements$OA_lin <- cos(measurements$OA) + sin(measurements$OA)

measurements <- measurements[,-c(3,4)]

landmarks <- superimpose_outlines(landmarks)

calculate_harmonic_power(landmarks)

n_harmonics <- 6

shape_variables <- elliptic_fourier_analysis(landmarks, n_harmonics, scale = TRUE)
form_variables <- elliptic_fourier_analysis(landmarks, n_harmonics, scale = FALSE)

#

# metric statistics linear ----------------

summary_stats(measurements$D)
summary_stats(measurements$WIS)
summary_stats(measurements$A)

# Depth

visualize_boxplot(measurements$D, measurements$Sample, "Depth (\u03BCm)")

shapiro.test(measurements$D)

# if inhomogeneous (shapiro p < 0.003)

kruskal.test(measurements$D, measurements$Sample)

# if homogeneous (shapiro p > 0.003)

anova(lm(measurements$D ~ measurements$Sample))

# WIS

visualize_boxplot(measurements$WIS, measurements$Sample, "WIS (\u03BCm)")

shapiro.test(measurements$WIS)

# if inhomogeneous (shapiro p < 0.003)

kruskal.test(measurements$WIS, measurements$Sample)

# if homogeneous (shapiro p > 0.003)

anova(lm(measurements$WIS ~ measurements$Sample))

# A

visualize_boxplot(measurements$A, measurements$Sample, "Asymmetry")

shapiro.test(measurements$A)

# if inhomogeneous (shapiro p < 0.003)

kruskal.test(measurements$A, measurements$Sample)

# if homogeneous (shapiro p > 0.003)

anova(lm(measurements$A ~ measurements$Sample))

#

# descriptive statistics for separate variables

target_variables <- "SetG"

summary_stats(measurements$D[measurements$Sample == target_variables])
summary_stats(measurements$WIS[measurements$Sample == target_variables])
summary_stats(measurements$A[measurements$Sample == target_variables])


x <- measurements$D[measurements$Sample == target_variables]
summary(x)
quantile(x, probs = c(0.25, 0.5, 0.75), type = 7)

x <- measurements$WIS[measurements$Sample == target_variables]
summary(x)
quantile(x, probs = c(0.25, 0.5, 0.75), type = 7)

x <- measurements$A[measurements$Sample == target_variables]
summary(x)
quantile(x, probs = c(0.25, 0.5, 0.75), type = 7)


#

# metric statistics circular ---------------

# for all notches 

angles <- circular::circular(measurements$OA, type = "angles", units = "radians")
circular::rose.diag(angles, col = "grey", bins = 15)

# preferential orientation (p < 0.003 means there is preferential orientation)

rayleigh.test(
  circular(angles,
           type = "angles", units = "radians")
)

descriptive_circular_analysis(angles)
circular::median.circular(angles)[[1]] * (180/pi)

# specific sets

target_sample <- "SetI" # CHANGE THIS!

angles <- circular::circular(measurements$OA[measurements$Sample == target_sample], type = "angles", units = "radians")
circular::rose.diag(angles, col = "grey", bins = 15)

# preferential orientation (p < 0.003 means there is preferential orientation)

rayleigh.test(
  circular(angles,
           type = "angles", units = "radians")
)

descriptive_circular_analysis(angles)
circular::median.circular(angles)[[1]] * (180/pi)

# compare two sets
# if p < 0.003 then the two sets of angles are different

set1 <- "SetH"
set2 <- "SetI"

sample1 <- circular(
  measurements$OA[measurements$Sample == set1],
  type = "angles",
  units = "radians"
)
sample2 <- circular(
  measurements$OA[measurements$Sample == set2],
  type = "angles",
  units = "radians"
)

target <- as.numeric(c(sample1, sample2))
ndat <- c(length(sample1), length(sample2))
g <- 2

PgVal(target = target, ndat = ndat, g = g) # test statistic
pchisq(PgVal(target = target, ndat = ndat, g = g), g-1, lower.tail = FALSE) # p value

#

# metric statistics pca ------------

pca_plot(measurements[,c(1,2,4,6)], measurements$Sample,
         label_colours = c('red1', 'palevioletred1','palegreen3', 'saddlebrown', 'purple', 'skyblue2',
                           'yellow4', 'orange', '#caaa84', '#d4b277',
                           '#ddbb6a', '#e6c35b', '#efcc49', '#f7d532', '#ffde00'),
         Chull = TRUE)

biplot_data <- prcomp(measurements[,c(1,2,4,6)], scale = TRUE)
X11(); biplot(biplot_data)

cumsum(pca_plot(measurements[,c(1,2,4,6)])$variance * 100)
pc_scores <- pca_plot(measurements[,c(1,2,4,6)])$pc_scores

# Univariate analysis of PC1

shapiro.test(pc_scores[,1])

# if inhomogeneous (shapiro p < 0.003)

kruskal.test(pc_scores[,1], measurements$Sample)

# if homogeneous (shapiro p > 0.003)

anova(lm(pc_scores[,1] ~ measurements$Sample))

#

# MANOVA

pairwise.perm.manova(pc_scores[,1:2],
                     measurements$Sample,
                     test = c("Wilks"),
                     nperm = 999,
                     progress = TRUE)

#

# EFA (elliptic fourier analysis) Allometry ------------

shape_size_data <- list(
  coordinates = as.matrix(shape_variables),
  CSize = measurements$D,
  Sample = as.factor(measurements$Sample)
)

allometric_fit <- lm.rrpp(coordinates ~ CSize + Sample, data = shape_size_data)

plot_residuals(allometric_fit, plot_labels = FALSE)
plot_residuals_qq(allometric_fit, plot_labels = FALSE)

anova.lm.rrpp(allometric_fit)$table

#

# EFA PCA Shape -----------------------------

GraphGMM::pca_plot(shape_variables, measurements$Sample, Chull = TRUE,
                   label_colours =c('red1', 'palevioletred1','palegreen3', 'saddlebrown', 'purple', 'skyblue2',
                                    'yellow4', 'orange', '#caaa84', '#d4b277',
                                    '#ddbb6a', '#e6c35b', '#efcc49', '#f7d532', '#ffde00'))
visualise_shape_change(landmarks, GraphGMM::pca_plot(shape_variables)$pc_scores, TRUE)

pc_scores <- GraphGMM::pca_plot(shape_variables)$pc_scores
visualize_boxplot(pc_scores[,1], measurements$Sample, "PC1")

cumsum(pca_plot(shape_variables)$variance * 100)
pc_scores <- pca_plot(shape_variables)$pc_scores

# Univariate analysis of PC1

shapiro.test(pc_scores[,1])

# if inhomogeneous (shapiro p < 0.003)

kruskal.test(pc_scores[,1], measurements$Sample)

# if homogeneous (shapiro p > 0.003)

anova(lm(pc_scores[,1] ~ measurements$Sample))

#

# MANOVA

pairwise.perm.manova(pc_scores[,1:5],
                     measurements$Sample,
                     test = c("Wilks"),
                     nperm = 999,
                     progress = TRUE)

#

# EFA PCA Form -----------------------------


GraphGMM::pca_plot(form_variables, measurements$Sample, Chull = TRUE,
                   label_colours = c('red1', 'palevioletred1','palegreen3', 'saddlebrown', 'purple', 'skyblue2',
                                     'yellow4', 'orange', '#caaa84', '#d4b277',
                                     '#ddbb6a', '#e6c35b', '#efcc49', '#f7d532', '#ffde00'))
visualise_shape_change(landmarks, GraphGMM::pca_plot(form_variables)$pc_scores, TRUE)

pc_scores <- GraphGMM::pca_plot(form_variables)$pc_scores
visualize_boxplot(pc_scores[,1], measurements$Sample, "PC1")

cumsum(pca_plot(form_variables)$variance * 100)
pc_scores <- pca_plot(form_variables)$pc_scores

# Univariate analysis of PC1

shapiro.test(pc_scores[,1])

# if inhomogeneous (shapiro p < 0.003)

kruskal.test(pc_scores[,1], measurements$Sample)

# if homogeneous (shapiro p > 0.003)

anova(lm(pc_scores[,1] ~ measurements$Sample))

#

# MANOVA

pairwise.perm.manova(pc_scores[,1:5],
                     measurements$Sample,
                     test = c("Wilks"),
                     nperm = 999,
                     progress = TRUE)

#

# pick a specific coordinate for tps ---------------------

target_space <- form_variables
pc_scores <- pca_plot(target_space)$pc_scores

target_coordinate <- c(1,0)

target_pc <- morphological_predictor(
  landmarks,
  pc_scores[,1:2],
  target_coordinate
)
plot_outline(target_pc)

#

# evaluate a p value ----------------

p_value <- 0.0028
FPR(p_value) * 100

# FPR is the False Positive Risk - the probability of being wrong
# look Courtenay_2024 in Quaternary Environment and Humans
# only calculate FPR for p values under 0.368

#


