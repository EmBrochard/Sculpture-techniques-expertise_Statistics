# Code written by Lloyd Austin Courtenay

# Copyright (C) 2025 Lloyd Courtenay
# SPDX-License-Identifier: AGPL-3.0

# libraries and functions -------------------------------

library(GraphGMM)

source(".\\Source\\Profile_Code.R")

#

# your folders -----------------------------------

# write here what the separator is
# "\t" is a tab
# "," is a comma
# ";" is a semi-colon
# " " is a space

separator = "\t"


# write here the name of the folder with the profiles

folder_name = "Gravures_fines"

#

# run code -----------------

create_images(folder_name, sep = separator)

results <- extract_measurements("Dataset", sep = separator)
metric_results <- results$metric_results
landmarks <- results$landmarks

write_morphologika_file("landmark_data", landmarks, metric_results$Sample)
write.table(metric_results, "metric_results.csv", col.names = TRUE, row.names = FALSE,
            sep = ";")

#



