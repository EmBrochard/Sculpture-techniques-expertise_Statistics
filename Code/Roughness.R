# Code written by Emilie Brochard, Lloyd Austin Courtenay, Luc Doyon (alphabetical order)

# Copyright (C) 2025 Émilie Brochard, Lloyd Courtenay, Luc Doyon
# SPDX-License-Identifier: AGPL-3.0

# External functions -------------------------------

circular_variance <- function(x) {
  
  return(suppressWarnings(1 - circular::trigonometric.moment(x, p = 1)$rho))
  
}

circular_std_skewness <- function(x) {
  
  return(
    suppressWarnings((
      circular::trigonometric.moment(x, p = 2, center = TRUE)$sin
    )/(
      (1 - circular::trigonometric.moment(x, p = 1)$rho)^(3/2)
    ))
  )
  
}

circular_std_kurtosis <- function(x) {
  
  return(
    suppressWarnings((
      circular::trigonometric.moment(
        x, p = 2, center = TRUE
      )$cos - circular::trigonometric.moment(
        x, p = 1
      )$rho^4)/(
        (1-circular::trigonometric.moment(x, p = 1)$rho)^2
      ))
  )
  
}

label_imbalance <- function(labels) {
  
  n <- length(labels)
  g <- length(levels(labels))
  c <- c()
  entropy_index <- c()
  for (ci in 1:g) {
    c <- c(
      c, table(labels)[ci][[1]]
    )
    entropy_index <- c(
      entropy_index, (c[ci] / n) * log(c[ci] / n)
    )
  }
  H <- -sum(entropy_index)
  balance_index <- H / log(g)
  return(balance_index)
  
}


mean_theta <- function(x) {
  
  n = length(x)
  numerator_x = 0
  numerator_y = 0
  
  for (i in 1:n) {
    
    numerator_x <- numerator_x + cos(x[i])
    numerator_y <- numerator_y + sin(x[i])
    
  }
  
  Y = numerator_y / n
  X = numerator_x / n
  
  r = sqrt(X^2 + Y^2)
  
  cosa = X / r
  sina = Y / r
  
  return(atan2(sina, cosa))
  
}

add_alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

`%!in%` = Negate(`%in%`)

split_data = function(data,p = 0.7, s = 666) {
  
  set.seed(s)
  
  index = sample (1:dim(data)[1])
  train = data [index[1: floor(dim(data)[1]*p)],]
  test = data [index[((ceiling(dim(data)[1]*p)) + 1) :dim(data)[1]],]
  
  return(list(train = train, test = test))
  
}

# Jammalamadaka Sarma Correlation Coefficient
# With randomised calculation of p-values
# an asymptotic test best suited for small sample sizes
# to compare two sets of circular variables

# circular-circular association correlation coefficients
# have received limited attention so far and could be further
# explored, but at present this test appears to be the
# most robust across a wide range of simulations

JammalamadakaSarma_Rand <- function(x, y, n = 1000) {
  
  observed_correlation <- suppressWarnings(circular::cor.circular(x, y))
  statistic = 1
  
  for (iteration in 1:n) {
    
    x_prima <- sample(x)
    
    randomised_correlation <- suppressWarnings(circular::cor.circular(x_prima, y))
    
    if (abs(randomised_correlation) >= abs(observed_correlation)) {
      
      statistic = statistic + 1
      
    }
    
  }
  
  p <- statistic / (n + 1)
  
  return(list(
    cor = observed_correlation,
    p = p
  ))
  
  
}

# Johnson Wehrly Mardia Correlation Coefficient
# for the comparison of circular variables with linear variables
# This is a parametric test
# p-values are calculated using randomised permutations for 
# the computation of p

JohnsonWehrlyMardia_Correlation <- function(linear, circular) {
  
  r_xc = cor(linear, cos(circular))
  r_xs = cor(linear, sin(circular))
  r_cs = cor(cos(circular), sin(circular))
  coefficient <- (r_xc^2 + r_xs^2 - (2 * r_xc * r_xs * r_cs)) / (1 - r_cs^2)
  return(coefficient)
  
}

JohnsonWehrlyMardia_p <- function(linear, circular, n = 1000) {
  
  coefficient <- JohnsonWehrlyMardia_Correlation(linear, circular)
  
  statistic = 1
  
  for (iteration in 1:n) {
    
    x_prima <- sample(linear)
    
    randomised_correlation <- JohnsonWehrlyMardia_Correlation(x_prima, circular)
    
    if (abs(randomised_correlation) >= abs(coefficient)) {
      
      statistic = statistic + 1
      
    }
    
  }
  
  p <- statistic / (n + 1)
  
  return(list(
    cor = coefficient,
    p = p
  ))
  
}

# the Mardia rank based correlation coefficient with randomised
# permutations for the calculation of p-values
# this is a non parametric test to compare circular variables with
# linear variables

MardiaRank_Correlation <- function(linear, circular, n = 1000) {
  
  U_prima_calc <- function(uniform_scores) {
    
    T_c = 0
    T_s = 0
    
    for (individual in 1:N) {
      
      T_c = T_c + individual * cos(uniform_scores[individual])
      T_s = T_s + individual * sin(uniform_scores[individual])
      
    }
    
    U_prima = (T_s * T_s) + (T_c * T_c)
    
    return(U_prima)
    
  }
  
  N = length(linear)
  circular_2 = 0
  
  rank_linear <- rank(linear, ties.method = "random")
  
  for (individual in 1:N) {
    
    circular_2[rank_linear[individual]] <- circular[individual]
    
  }
  
  rank_circular = rank(circular_2, ties.method = "random")
  
  uniform_scores <- rank_circular * 2 * pi / N
  
  observed_U <-  U_prima_calc(uniform_scores)
  statistic = 1
  
  for (iteration in 1:n) {
    
    randomised_ranks <- sample(uniform_scores)
    randomised_U <- U_prima_calc(randomised_ranks)
    
    if (randomised_U > observed_U) {
      
      statistic = statistic + 1
      
    }
    
  }
  
  p = statistic / (n + 1)
  
  return(list(
    U = observed_U,
    p = p
  ))
  
}

# Mardia Watson Wheeler test
# non parametric test to compare two sets of angles and
# see if they are different

CosSinUniScores <- function(target) {
  
  N = length(target)
  ranks = rank(target, ties.method= "random")
  CosUniScores = cos((ranks*2*pi)/N)
  SinUniScores = sin((ranks*2*pi)/N)
  
  return(list(
    CosUniScores,
    SinUniScores
  ))
  
}

WgVal <- function(CSUScores, ndat, g) {
  
  CosUScores = CSUScores[[1]]
  SinUScores = CSUScores[[2]]
  N = length(CosUScores)
  ndatcsum = cumsum(ndat)
  Wg = 0
  
  for (k in 1:g) {
    
    CosUScoresk = 0
    SinUScoresk = 0
    
    if (k==1) {
      
      low = 0
      
    } else if (k > 1) {
      
      low = ndatcsum[k - 1]
        
    }
    
    for (j in 1:ndat[k]) {
      
      CosUScoresk[j] = CosUScores[j+low] ; SinUScoresk[j] = SinUScores[j+low]
      
    }
    
    sumCkSq = (sum(CosUScoresk))**2
    sumSkSq = (sum(SinUScoresk))**2
    Wg = Wg + (sumCkSq + sumSkSq) / ndat[k]
    
  }
  
  Wg = 2 * Wg
  
  return(Wg)
  
}

# Fisher's non-parametric test for comparing the medians of two
# sets of angles.

MinusPiPi<-function(sample) {
  
  n = length(sample)
  
  for (j in 1:n) {
    
    if (sample[j] < -pi) {
      
      sample[j] = sample[j]+(2*pi)
      
    } else
      
      if (sample[j] > pi) {
        
        sample[j] = sample[j]-(2*pi)
        
      }
    
  }
  
  return(sample)
  
}

PgVal<-function(target, ndat, g) {
  
  N = length(target)
  sumterms = 0 ; M = 0
  ndatcsum = cumsum(ndat)
  gmedian = suppressWarnings(circular::median.circular(target)[[1]])
  
  for (k in 1:g) {
    
    if (k==1) {
      
      low = 0
      
    } else if (k > 1) {
      
      low = ndatcsum[k-1]
      
    }
    
    sample = suppressWarnings(circular::circular(0))
    
    for (j in 1:ndat[k]) {
      
      sample[j] = target[j+low]
      
    }
    
    shiftdat = MinusPiPi(sample - gmedian)
    m = length(shiftdat[shiftdat<0])
    M = M+m
    sumterms = sumterms + m*m/ndat[k]
    
  }
  
  term1 = ((N*N)/(M*(N-M)))
  term2 = (N*M)/(N-M)
  
  Pg = term1*sumterms - term2
  
  return(Pg)
  
}

#

# Step 1 -------------------------------------------

# A. Read the CSV file
data <- read.csv("BDD_V16.csv",sep=";")

# B. Create an output directories
output_dir = "results"
# Ensure the output directory exists
if (!dir.exists(output_dir)) dir.create(output_dir)
output_dir = "results/Step1"
# Ensure the output directory exists
if (!dir.exists(output_dir)) dir.create(output_dir)


# C. Eliminate rows with "****"
data_char <- data.frame(lapply(data, as.character), stringsAsFactors = FALSE)
rows_to_remove <- apply(data_char, 1, function(x) any(x == "****"))
data_clean <- data[!rows_to_remove, ]
write.csv(data_clean, "results/Step1/data_clean.csv", row.names = FALSE)


#

# Step 2 ------------------------------

# Step 2.A ------------------------------
  # Stats results (Shapiro, Kruskal, ANOVA), Summary stats (mean, median, CI, etc), 
  # Pairwise results (Wilcoxon, Tukey), Boxplots

# Before running the function, ensure to indicate the right start and end column number

generate_boxplots_statistics <- function(data = data_clean, group_col = "Technique", 
                                         start_col = 10, end_col = 38, output_dir = "results/Step2") {
  
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
  if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(writexl)
  
  # Input validation
  if (!group_col %in% colnames(data)) {
    stop("The group column (", group_col, ")is not present in the dataset.")
  }
  if (ncol(data) < end_col) {
    stop("The specified column range exceeds the number of columns in the dataset.")
  }
  
  # Verify that the output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Prepare dataframes to store test results
  statistical_results <- data.frame(
    variable = character(),
    test = character(),
    shapiro_p = numeric(),
    shapiro_w = numeric(),
    Df = numeric(),
    test_statistic = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
 
  
  pairwise_results <- data.frame(
    Variable = character(),
    Test = character(),
    Group1 = character(),
    Group2 = character(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Prepare the table for boxplot descriptive statistics
  boxplot_summary_results <- data.frame(
    Variable = character(),
    Technique = character(),
    Min = numeric(),
    CI_lower = numeric(),
    Q1 = numeric(),
    Median = numeric(),
    Mean = numeric(),
    SD = numeric(), 
    NMAD = numeric(),
    sqrtBWMV = numeric(),
    Q3 = numeric(),
    CI_upper = numeric(),
    Max = numeric(),
    stringsAsFactors = FALSE
  )
  
  all_plots <- c()
  plot_index <- 1
  
  # Loop through the specified columns
  for (col in colnames(data)[start_col:end_col]) {
    
    if (col %!in% c("Std", "First_Direction", "Second_direction", "Third_direction")) {
      
      # 1: Ensure that the column is numeric
      data[[col]] <- as.numeric(as.character(data[[col]]))
      
      # 2: Check for NA, Inf, or non-numeric values
      if (any(is.na(data[[col]]) | is.infinite(data[[col]]) | !is.numeric(data[[col]]))) {
        print(paste("⚠️ The column", col, "- contains NA, Inf, or non-numeric values"))
        next
      }
      
      # 3: Replace NA values with the median
      data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
      
      # 4: Shapiro test
      shapiro_results <- shapiro.test(data[[col]])
      
      if (shapiro_results$p.value < 0.003) {
        
        kruskal_result <- tryCatch({
          kruskal.test(as.formula(paste(col, "~", group_col)), data = data)
        }, error = function(e) {
          print(paste("⚠️ The column", col, "- - Kruskal-Wallis test failed:", e$message))
          return(NULL)
        })
        
        statistical_results <- rbind(
          statistical_results,
          data.frame(
            variable = col,
            test = "kruskal",
            shapiro_p = shapiro_results$p.value,
            shapiro_w = shapiro_results$statistic,
            Df = kruskal_result$parameter[[1]],
            test_statistic = kruskal_result$statistic[[1]],
            p_value = kruskal_result$p.value
          )
        )
        
        test_statistic <- kruskal_result$statistic[[1]]
        used_test <- "kruskal"
        df <- kruskal_result$parameter[[1]]
        p_value <- kruskal_result$p.value
        
        # 5: Pairwise Wilcoxon test with Bonferroni correction
        pairwise_table <- tryCatch({
          pairwise.wilcox.test(data[[col]], data[[group_col]], p.adjust.method = "bonferroni")
        }, error = function(e) {
          print(paste("⚠️ The column", col, "- Wilcoxon test failed:", e$message))
          return(NULL)
        })
        
        comparisons <- list()
        wilcoxon_matrix <- matrix(NA, ncol = length(unique(data[[group_col]])), nrow = length(unique(data[[group_col]])))
        rownames(wilcoxon_matrix) <- unique(data[[group_col]])
        colnames(wilcoxon_matrix) <- unique(data[[group_col]])
        
        if (!is.null(pairwise_table)) {
          
          # Store Wilcoxon test results in a matrix for each variable
          all_pairs <- which(!is.na(pairwise_table$p.value), arr.ind = TRUE)
          for (i in seq_len(nrow(all_pairs))) {
            group1 <- rownames(pairwise_table$p.value)[all_pairs[i, 1]]
            group2 <- colnames(pairwise_table$p.value)[all_pairs[i, 2]]
            p_val <- pairwise_table$p.value[all_pairs[i, 1], all_pairs[i, 2]]
            
            pairwise_results <- rbind(
              pairwise_results,
              data.frame(
                Variable = col,
                Test = "Wilcoxon",
                Group1 = group1,
                Group2 = group2,
                P_value = p_val
              )
            )
            
            wilcoxon_matrix[group1, group2] <- p_val
            wilcoxon_matrix[group2, group1] <- p_val
            
            if (!is.na(p_val) && p_val < 0.003) {
              comparisons <- append(comparisons, list(c(group1, group2)))
            }
          }
          
          # Export the Wilcoxon matrix for this column
          write.table(wilcoxon_matrix, file.path(output_dir, paste0("wilcoxon_matrix_", col, ".csv")),
                      sep = ",", row.names = TRUE)
          
        }
        
      } else {
        
        anova_results <- anova(lm(data[[col]] ~ data[[group_col]]))
        
        statistical_results <- rbind(
          statistical_results,
          data.frame(
            variable = col,
            test = "ANOVA",
            shapiro_p = shapiro_results$p.value,
            shapiro_w = shapiro_results$statistic,
            Df = anova_results$Df[1],
            test_statistic = anova_results$`F value`[1],
            p_value = anova_results$`Pr(>F)`[1]
          )
        )
        
        test_statistic <- anova_results$`F value`[1]
        used_test <- "ANOVA"
        df <- anova_results$Df[1]
        p_value <- anova_results$`Pr(>F)`[1]
        
        # 6: Test pairwise Tukey
        
        tukey_results <- TukeyHSD(aov(lm(data[[col]] ~ data[[group_col]])), conf.level = 0.95)$`data[[group_col]]`
        comparisons <- list()
        
        if (!is.null(tukey_results)) {
          
          all_pairs <- which(!is.na(tukey_results[,4]), arr.ind = TRUE)
          n_comparisons <- nrow(tukey_results)
          
          for (pair_set in 1:n_comparisons) {
            
            comparison_string <- rownames(tukey_results)[pair_set]
            
            separation_symbol <- gregexpr("-", comparison_string)[[1]][1]
            group1 <- substr(comparison_string, 1, separation_symbol - 1)
            group2 <- substr(comparison_string, separation_symbol + 1, nchar(comparison_string))
            p_val <- tukey_results[pair_set,4]
            
            pairwise_results <- rbind(
              pairwise_results,
              data.frame(
                Variable = col,
                Test = "Tukey HSD",
                Group1 = group1,
                Group2 = group2,
                P_value = p_val
              )
            )
            
            if (!is.na(p_val) && p_val < 0.003) {
              
              comparisons <- append(comparisons, list(c(group1, group2)))
              
            }
            
          }
          
          # Export the Tukey matrix for this column 
          write.table(tukey_results, file.path(output_dir, paste0("tukey_matrix_", col, ".csv")),
                      row.names = TRUE, col.names = TRUE,
                      sep = ",")
          
        }
        
      }
      
      if (is.null(kruskal_result)) next
      
      
      # Calculate parameter statistics (mean, median, etc)
      summary_stats <- data %>%
        group_by(!!sym(group_col)) %>%
        summarise(
          Median = median(.data[[col]], na.rm = TRUE),
          Mean = mean(.data[[col]], na.rm = TRUE),
          Q1 = quantile(.data[[col]], 0.25, na.rm = TRUE),
          Q3 = quantile(.data[[col]], 0.75, na.rm = TRUE),
          Min = min(.data[[col]], na.rm = TRUE),
          Max = max(.data[[col]], na.rm = TRUE),
          SD = sd(.data[[col]], na.rm = TRUE),
          NMAD = pValueRobust::median_absolute_deviation(.data[[col]]),
          sqrtBWMV = pValueRobust::biweight_midvariance(.data[[col]]),
          CI_lower = quantile(.data[[col]], 0.025, na.rm = TRUE), 
          CI_upper = quantile(.data[[col]], 0.975, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(Variable = col) %>%
        relocate(Variable, .before = 1) %>%
        rename(Technique = !!group_col)
      
      boxplot_summary_results <- bind_rows(boxplot_summary_results, summary_stats)
      
  # Generate the boxplot
      data[[group_col]] <- factor(data[[group_col]], levels = c("CS","PIPB","PIPP", "PDPP","PDPC","PSc",
                                                                "ScS","ScB","PScPo","PoHS","PoS","PoC"))
      
      if (used_test == "kruskal") {
        
        if (p_value < 0.0001) {
          
          p_value_for_string <- formatC(p_value, format = "e", digits = 2)
          
        } else {
          
          p_value_for_string <- round(p_value, 4)
          
        }
        
        plot_title <- paste0("Kruskal-Wallis: \U03C7\U00B2(", df, ") = ", round(test_statistic, 2), "; p = ", p_value_for_string)
        
      } else {
        
        if (p_value < 0.0001) {
          
          p_value_for_string <- formatC(p_value, format = "e", digits = 2)
          
        } else {
          
          p_value_for_string <- round(p_value, 4)
          
        }
        
        plot_title <- paste0("ANOVA: F(", df, ") = ", round(test_statistic, 2), "; p = ", p_value_for_string)
        
      }
      
      p <- ggplot(data, aes_string(x = group_col, y = col)) +
        geom_boxplot(aes_string(fill = group_col), outlier.shape = NA) +
        geom_jitter(width = 0.1, size = 1, color = "black") +
        scale_fill_manual(values = c("grey90", "navyblue","royalblue3","skyblue1","paleturquoise",
                                     "lightseagreen", "darkgreen","olivedrab2","yellow3","gold",
                                     "yellow1","lightgoldenrod1")) +
        labs(x = group_col, y = col, title = plot_title) +
        theme_classic() +
        theme(
          plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
          axis.text = element_text(size = 10),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1,face = "bold", size = 12, color = "black",
                                     margin = margin(t = 7)),
          axis.title.x = element_text(size = 15, face = "bold",
                                      margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
          axis.title.y = element_text(size = 15, face = "bold",
                                      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
          plot.title = element_text(hjust = 0, size = 12, face = "bold"),
          legend.position = "none"
        )
      
      # Add significant comparisons with asterisks — DO NOT RUN if you want clean boxplots
      #if (length(comparisons) > 0) {
        
       # p <- p + stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", size = 2)
        
        
     # }
      
      
      
      
      
      all_plots[[plot_index]] <- p
      
      plot_index <- plot_index + 1
      
      # Save the boxplot 
      ggsave(
        filename = file.path(output_dir, paste0("boxplot_", col, ".png")),
        plot = p,
        width = 7, height = 5
      )
    
    }
    
  }
  
  # Export the results to CSV file
  write.csv(statistical_results, file = file.path(output_dir, "stat_results.csv"), row.names = FALSE)
  write.csv(pairwise_results, file = file.path(output_dir, "pairwise_results.csv"), row.names = FALSE)
  write_xlsx(pairwise_results, file.path(output_dir, "pairwise_results.xlsx"))
  write.csv(boxplot_summary_results, file = file.path(output_dir, "summary_stats.csv"), row.names = FALSE)
  
  
  
  message("✅ Analysis complete! The results have been saved in the directory", output_dir)
  
  return(all_plots)  
  
}

results <- generate_boxplots_statistics()

#Print results in the console if you want to display the boxplots
print (results)


# Step 2.B ---------------------
 # Study of circular data : Circular stat results (Mardia-Watson-Wheeler), Summary stats (mean, median, CI, etc), 
 # Rose diagram

data[data$ID != "CS",]

generate_boxplots_circular_statistics <- function(
    data = data_clean, group_col = "Technique", start_col = 10, end_col = 38,
    output_dir = "results/Step2"
)
{
  
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
  if (!requireNamespace("circular", quietly = TRUE)) install.packages("circular")
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(circular)
  
  # Input validation
  if (!group_col %in% colnames(data)) {
    stop("The group column (", group_col, ") is not present in the dataset.")
  }
  if (ncol(data) < end_col) {
    stop("The specified column range exceeds the number of columns in the dataset.")
  }
  
  # Check if the output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Prepare dataframes to store test results
  
  statistical_results <- data.frame(
    variable = character(),
    test = character(),
    Df = numeric(),
    Wg = numeric(),
    p_value = numeric(),
    test2 = character(),
    Pg = numeric(),
    p_value2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  circular_summary_results <- data.frame(
    Variable = character(),
    Technique = character(),
    Min = numeric(),
    CI_lower = numeric(),
    Q1 = numeric(),
    Median = numeric(),
    Mean = numeric(),
    v = numeric(), 
    shat = numeric(),
    khat = numeric(),
    Q3 = numeric(),
    CI_upper = numeric(),
    Max = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  
  for (col in colnames(data[start_col:end_col])) {
    
    if (col %in% c("Std", "First_Direction", "Second_direction", "Third_direction")) {
      
      data[[col]] <- as.numeric(as.character(data[[col]]))
      
      if (any(is.na(data[[col]]) | is.infinite(data[[col]]) | !is.numeric(data[[col]]))) {
        print(paste("⚠️ The colum", col, "- contains NA, Inf, or non-numeric values"))
        next
      }
      
      data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
      
      collected_p_values <- c()
      collected_Wg_values <- c()
      collected_df_values <- c()
      collected_Pg_values <- c()
      collected_p_Pg_values <- c()
      
      for (individual_sample_i in levels(as.factor(data[[group_col]]))) {
        
        for (individual_sample_j in levels(as.factor(data[[group_col]]))) {
          
          if (individual_sample_i != individual_sample_j) {
            
            sample_1 <- data[data[[group_col]] == individual_sample_i,][[col]]
            sample_2 <- data[data[[group_col]] == individual_sample_j,][[col]]
            
            sample_1 <- sample_1 * (pi / 180)
            sample_2 <- sample_2 * (pi / 180)
            
            target_sample <- as.numeric(c(sample_1, sample_2))
            ndat <- c(length(sample_1), length(sample_2))
            g <- 2
            
            pairwise_test_statistic <- WgVal(CosSinUniScores(target_sample), ndat = ndat, g = g)
            pairwise_p_value <- pchisq(WgVal(CosSinUniScores(target_sample), ndat = ndat, g = g),
                                       2 * (g - 1), lower.tail = FALSE)
            
            collected_p_values <- c(collected_p_values, pairwise_p_value)
            collected_Wg_values <- c(collected_Wg_values, pairwise_test_statistic)
            collected_df_values <- c(collected_df_values, 2 * (g - 1))
            
            pairwise_test_statistic_pg <- PgVal(target = target_sample, ndat = ndat, g = g)
            pairwise_p_value_pg <- pchisq(
              PgVal(target = target_sample, ndat = ndat, g = g),
              g - 1, lower.tail = FALSE
            )
            
            collected_Pg_values <- c(collected_Pg_values, pairwise_test_statistic_pg)
            collected_p_Pg_values <- c(collected_p_Pg_values, pairwise_p_value_pg)
            
          }
          
        }
        
      }
      
      statistical_results <- rbind(
        statistical_results,
        data.frame(
          variable = col,
          test = "Mardia-Watson-Wheeler",
          Df = median(collected_df_values),
          Wg = median(collected_Wg_values),
          p_value = median(collected_p_values),
          test2 = "Fisher's Non-Parametric",
          Pg = median(collected_Pg_values),
          p_value2 = median(collected_p_Pg_values),
          stringsAsFactors = FALSE
        )
      )
      
      plot_colours <-  c(
        "grey90", "navyblue","royalblue3","skyblue1","paleturquoise",
        "lightseagreen", "darkgreen","olivedrab2","yellow3","gold",
        "yellow1","lightgoldenrod1"
      )
     # Uncomment the following line to save as PNG
      # png(file.path(output_dir, paste0("rose_diagram_", col, ".png")), width = 2000, height = 2000, res = 300)
      
      rose.diag(suppressWarnings(circular(data[[col]], units = "degrees")), col = "grey", bins = 15)
      index <- 1
      
      n_samples <- length(levels(as.factor(data[[group_col]])))
      offset <- seq(-0.05, -0.3, length = n_samples)
      
      for (individual_sample in levels(as.factor(data[[group_col]]))) {
        
        individual_values <- data[data[[group_col]] == individual_sample,][[col]]
        
        points(suppressWarnings(circular(individual_values, units = "degrees")), col = plot_colours[index],
               next.points = offset[index])
        
        index <- index + 1
        
      }
      # Uncomment the following line to save as PNG
      # dev.off()
      
      summary_stats <- data %>%
        group_by(!!sym(group_col)) %>%
        summarise(
          Median = suppressWarnings(median.circular(.data[[col]] * (pi / 180), na.rm = TRUE)[[1]]) * (180 / pi),
          Mean = suppressWarnings(mean_theta(.data[[col]] * (pi / 180))[[1]]) * (180 / pi),
          Q1 = pValueRobust::quantile_CI(.data[[col]] * (pi / 180),
                                         q = 0.25)[[1]] * (180 / pi),
          Q3 = pValueRobust::quantile_CI(.data[[col]] * (pi / 180),
                                         q = 0.75)[[1]] * (180 / pi),
          Min = min(.data[[col]], na.rm = TRUE)[[1]],
          Max = max(.data[[col]], na.rm = TRUE)[[1]],
          v = suppressWarnings(circular_variance(.data[[col]] * (pi / 180))[[1]]),
          shat = suppressWarnings(circular_std_skewness(.data[[col]] * (pi / 180))[[1]]),
          khat = suppressWarnings(circular_std_kurtosis(.data[[col]] * (pi / 180))[[1]]),
          CI_lower = pValueRobust::quantile_CI(.data[[col]] * (pi / 180),
                                               0.025)[[1]] * (180 / pi), 
          CI_upper = pValueRobust::quantile_CI(.data[[col]] * (pi / 180),
                                               0.975)[[1]] * (180 / pi),
          .groups = "drop"
        ) %>%
        mutate(Variable = col) %>%
        relocate(Variable, .before = 1) %>%
        rename(Technique = !!group_col)
      
      circular_summary_results <- bind_rows(circular_summary_results, summary_stats)
      
    }
    
  }
  

  
  write.csv(statistical_results,
            file = file.path(output_dir, "circular_stat_results.csv"), row.names = FALSE)
  write.csv(circular_summary_results,
            file = file.path(output_dir, "circular_summary_stat_results.csv"), row.names = FALSE)
  
}

generate_boxplots_circular_statistics()




#

# Step 3 -------------------------------

# Step 3.A. Filter statistical results -------------------------------

filter_significant_variables <- function(stat_results_file = "results/Step2/stat_results.csv", 
                                         circular_results_file = "results/Step2/circular_stat_results.csv",
                                         original_data_file = "results/Step1/data_clean.csv", 
                                         threshold = 0.003, output_dir = "results/Step3") {
  # Load required library
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Read the Kruskal-Wallis results
  stat_results <- read.csv(stat_results_file)
  circular_stat_results <- read.csv(circular_results_file)
  
  # Filter rows where the p-value is less than the threshold
  significant_results <- stat_results %>%
    filter(p_value < threshold)
  significant_circ_results <- circular_stat_results %>%
    filter(p_value < threshold)
  
  if (nrow(significant_circ_results) != 0) {
    
    significant_circ_results <- data.frame(
      variable = significant_circ_results$variable,
      test = significant_circ_results$test,
      shapiro_p = 0,
      shapiro_w = 0,
      Df = significant_circ_results$Df,
      test_statistic = significant_circ_results$Wg,
      p_value = significant_circ_results$p_value
    )
    
  }
  
  significant_results <- rbind(significant_results, significant_circ_results)
  
  # Ensure Variable is a character vector
  significant_results$Variable <- as.character(significant_results$variable)
  
  # Print significant variables
  message("Significant variables (p < ", threshold, "):")
  print(significant_results$variable)
  
  # Save significant variables to a CSV file
  output_file <- file.path(output_dir, "significant_variables.csv")
  write.csv(significant_results$variable, file = output_file, row.names = FALSE)
  
  # Read the original data
  original_data <- read.csv(original_data_file)
  
  # Check for missing columns
  missing_columns <- setdiff(significant_results$variable, names(original_data))
  if (length(missing_columns) > 0) {
    stop("The following columns are missing in the original data: ", paste(missing_columns, collapse = ", "))
  }
  
  # Create a new dataframe with only the significant variables
  significant_data <- original_data %>%
    dplyr::select(dplyr::all_of(c("Technique", significant_results$variable))) #Modify according to the parameter you want to analyze
  
  # Save the new dataframe to a CSV file
  output_file <- file.path(output_dir, "significant_data.csv")
  write.csv(significant_data, file = output_file, row.names = FALSE)
  
  message("Filtered significant data saved to: ", output_file)
  return(significant_data)
}

filter_significant_variables()

#

# Step 3.B. Filter out highly correlated variables ------------------------------------

filter_highly_correlated <- function(sigdata = "results/Step3/significant_data.csv", 
                                     threshold = 0.7, 
                                     output_dir = "results/Step3") {
  
  # Ensure that the output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Read the significant data file
  sigdata <- read.csv(sigdata)
  
  # Display content for diagnostic purposes
  message("Contenu de significant_data.csv :")
  print(head(sigdata))  # Show the first rows
  
  # Verify the data structure
  message("Structure de significant_data.csv :")
  str(sigdata)
  
  # If the dataframe is empty, stop execution
  if (nrow(sigdata) == 0) {
    stop("The file significant_data.csv is empty. Check its contents before continuing.")
  }
  
  # Check if the "Technique" column exists and set it aside
  if ("Technique" %in% colnames(sigdata)) {
    technique_col <- sigdata["Technique"]
    sigdata <- sigdata[, colnames(sigdata) != "Technique"]  # Temporarily remove "Technique"
  } else {
    technique_col <- NULL
  }
  
  # Calculate the correlation matrix (numeric variables only)
  numeric_data <- sigdata[, sapply(sigdata, is.numeric)]
  
  circular_stats <- numeric_data
  col_indices <- c()
  for (column in 1:ncol(circular_stats)) {
    
    if (colnames(circular_stats)[column] %in% c("Std",
                                                "First_Direction",
                                                "Second_direction",
                                                "Third_direction")) {
      
      col_indices <- c(col_indices, column)
      
    }
    
  }
  
  if(!is.null(col_indices)) {
    
    circular_stats <- circular_stats[,col_indices] * (pi / 180)
    
    numeric_data <- numeric_data[,-col_indices]
    
  } else {
    
    circular_stats <- NULL
    
  }
  
  corr_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
  corr_matrix_names <- rownames(corr_matrix)
  
  # Circular statistics part
  # step a - calculate correlation between linear and circular variables
  
  linear_circular_cor_results <- c()
  
  if (!is.null(circular_stats)) {
    
    if (ncol(matrix(circular_stats)) != 0) {
      
      if (ncol(matrix(circular_stats)) == 1) {
        
        circular_stats <- data.frame(circular_stats)
        colnames(circular_stats) <- colnames(sigdata[col_indices])
        
      }
      
      circular_circular_cor_matrix <- matrix(NA, ncol = ncol(circular_stats), nrow = ncol(circular_stats))
      
      # Compute correlations and store in upper triangle
      for (i in 1:ncol(circular_stats)) {
        for (j in i:ncol(circular_stats)) {  # Start from i to avoid redundant comparisons
          if (i != j) {
            circular_circular_cor_matrix[i, j] <- JammalamadakaSarma_Rand(
              circular_stats[, i], circular_stats[, j]
            )$cor
          }
        }
      }
      
      # Add column and row names
      rownames(circular_circular_cor_matrix) <- colnames(circular_stats)
      colnames(circular_circular_cor_matrix) <- colnames(circular_stats)
      
      linear_circular_cor_matrix <- matrix(NA, nrow = ncol(numeric_data), ncol = ncol(circular_stats))
      
      # Compute linear-circular correlations
      for (i in 1:ncol(numeric_data)) {
        for (j in 1:ncol(circular_stats)) {
          linear_circular_cor_matrix[i, j] <- JohnsonWehrlyMardia_Correlation(
            numeric_data[, i], circular_stats[, j]
          )
        }
      }
      
      # Add row and column names
      rownames(linear_circular_cor_matrix) <- colnames(numeric_data)
      colnames(linear_circular_cor_matrix) <- colnames(circular_stats)
      
      circular_circular_cor_matrix[
        lower.tri(circular_circular_cor_matrix)
      ] <- t(circular_circular_cor_matrix)[lower.tri(circular_circular_cor_matrix)]
      
      num_linear <- ncol(numeric_data)
      num_circular <- ncol(circular_stats)
      
      full_correlation_matrix <- matrix(NA, nrow = num_linear + num_circular, ncol = num_linear + num_circular)
      
      full_correlation_matrix[1:num_linear, 1:num_linear] <- corr_matrix
      
      full_correlation_matrix[
        (num_linear + 1):(num_linear + num_circular),
        (num_linear + 1):(num_linear + num_circular)
      ] <- circular_circular_cor_matrix
      
      full_correlation_matrix[1:num_linear, (num_linear + 1):(num_linear + num_circular)] <- linear_circular_cor_matrix
      
      full_correlation_matrix[(num_linear + 1):(num_linear + num_circular), 1:num_linear] <- t(linear_circular_cor_matrix)
      
      diag(full_correlation_matrix) <- 1
      
      full_variable_names <- c(colnames(numeric_data), colnames(circular_stats))
      rownames(full_correlation_matrix) <- full_variable_names
      colnames(full_correlation_matrix) <- full_variable_names
      
    } else {
      
      full_correlation_matrix <- corr_matrix
      
    }
    
  } else {
    
    full_correlation_matrix <- corr_matrix
    
  }
  
  # step b - calculation correlation between circular and circular variables
  
  # Display the correlation matrix
  message("Correlation matrix:")
  print(full_correlation_matrix)
  
  # Save the correlation matrix
  write.csv(full_correlation_matrix, file = file.path(output_dir, "correlation_matrix.csv"))
  
  # visualization of the correlation matrix
  if (requireNamespace("corrplot", quietly = TRUE)) {
    
    library(corrplot)
    #Uncomment the following line to save as PNG
    #png(file.path(output_dir, "corr_plot.png"), width = 2000, height = 2000, res = 300)
    corrplot(full_correlation_matrix, method = "circle", title = "Correlation matrix", mar = c(0,0,1,0),
             tl.col = "#003366") 
    #Uncomment the following line to save as PNG
    #dev.off()
    
  } else {
    
    message("The 'corrplot' package is not installed. Install it using install.packages('corrplot').")
    
  }
  
  # Identify highly correlated variables
  high_corr <- caret::findCorrelation(full_correlation_matrix, cutoff = threshold, verbose = FALSE)
  
  # Convert indices to column names
  high_corr_vars <- colnames(full_correlation_matrix)[high_corr]
  
  # Check which variables will be removed
  message("Variables fortement corrélées détectées (seuil = ", threshold, "):")
  print(high_corr_vars)
  
  # Save the list of removed variables
  write.csv(high_corr_vars, file = file.path(output_dir, "highly_correlated_variables.csv"), row.names = FALSE)
  
  # Check if epLsar is detected
  message("epLsar removed?", "epLsar" %in% high_corr_vars)
  message("NewEpLsar removed?", "NewEpLsar" %in% high_corr_vars)
  
  # Remove the detected variables
  if (length(high_corr_vars) > 0) {
    filtered_data <- sigdata[, !(colnames(sigdata) %in% high_corr_vars), drop = FALSE]
  } else {
    filtered_data <- sigdata
  }
  
  # Ensure there are still variables left
  if (ncol(filtered_data) == 0) {
    stop("All variables have been removed. Try lowering the threshold.")
  }
  
  # Reintegrate the "Technique" column only if it existed originally
  if (!is.null(technique_col)) {
    filtered_data <- cbind(technique_col, filtered_data)
  }
  
  # Save the filtered file
  write.csv(filtered_data, file = file.path(output_dir, "filtered_data.csv"), row.names = FALSE)
  
  message("Number of variables removed due to high correlation: ", length(high_corr_vars))
  return(filtered_data)
}

# Run the function
filter_highly_correlated()

#

# Step 4 ------------------------
 # PCA 

# Additional categories may be added in the additional_cols = c("ID", "Technique", "Expertise")
# Confidence threshold may be changed in the confidence_level = 0.95
# PCA with custom dimensions (dimensions may be changed manually in the dims = c(x, y))

pca_with_custom_dimensions <- function(filtered_data_file = "results/Step3/filtered_data.csv", 
                                       original_data_file = "results/Step1/data_clean.csv", 
                                       group_col = "Technique", 
                                       additional_cols = c("ID", "Expertise","Technique"), 
                                       confidence_level = 0.95, 
                                       dims = c(1, 2), 
                                       output_dir = "results/Step4") {
 
   # Load required libraries
  if (!requireNamespace("FactoMineR", quietly = TRUE)) install.packages("FactoMineR")
  if (!requireNamespace("factoextra", quietly = TRUE)) install.packages("factoextra")
  library(FactoMineR)
  library(factoextra)

  # Ensure the output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir)

  # Read the filtered and original data file
  filtered_data <- read.csv(filtered_data_file)
  original_data <- read.csv(original_data_file)
  
  # Ensure the grouping column exists
  if (!group_col %in% colnames(filtered_data)) {
    stop("The grouping column (", group_col, ") is not found in the data.")
  }
  
  # Separate the grouping column and numeric data for PCA
  group <- filtered_data[[group_col]]
  
  for (column in 1:ncol(filtered_data)) {
    
    if (colnames(filtered_data)[column] %in% c("Std",
                                               "First_Direction",
                                               "Second_direction",
                                               "Third_direction")) {
      
      filtered_data[column] <- sin(filtered_data[column] * (pi / 180)) + cos(filtered_data[column] * (pi / 180))
      
    }
    
  }
  
  numeric_data <- filtered_data[, !(colnames(filtered_data) %in% group_col)]
  
  # Perform PCA on normalized data
  pca_result <- PCA(numeric_data, scale.unit = TRUE, graph = FALSE)
  
  # Extract the individual coordinates for all dimensions
  individuals_coord <- as.data.frame(pca_result$ind$coord)
  
  # Identify the index for the additional columns in the original data file
  col_ind <- which(colnames(original_data) %in% additional_cols)
  
  # Extract additional columns
  additional_data <- as.data.frame(original_data[, col_ind])
  
  # Merge PCA results with additional data using the common key
  individuals_coordinates <- cbind(individuals_coord, additional_data)

    # Save the individual coordinates to a CSV file
  individuals_file <- file.path(output_dir, "pca_individuals_coordinates.csv")
  write.csv(individuals_coordinates, individuals_file, row.names = FALSE)
  message("Individual PCA coordinates saved to: ", individuals_file)
  
  # Extract the percentage of variance explained by the selected PCs
  variance_explained <- pca_result$eig[dims, "percentage of variance"]
  
  # Axis labels with variation explained
  x_label <- paste0("Principal Component ", dims[1], " (", round(variance_explained[1], 2), "%)")
  y_label <- paste0("Principal Component ", dims[2], " (", round(variance_explained[2], 2), "%)")
  
  
  # Plot 1: Individuals with confidence ellipses
  
  

  custom_colors <- c("CS"="#D0C4DF", "PIPB"="#5D675B", "PIPP"="#ffa8cb", "PDPP"="skyblue1", 
                     "PDPC"="#991932", "PSc"="lightseagreen", "ScS"="darkgreen", 
                     "ScB"="olivedrab2", "PScPo"="#e31a1c", "PoHS"="yellow", 
                     "PoS"= "orange", "PoC"="darkmagenta")

  
  
  custom_shape <- c(0:6,8,15:18)
  
  #hulls <- individuals_coordinates %>%
  #  dplyr::group_by(!!sym(group_col)) %>%
  #  dplyr::slice(chull(!!sym(names(individuals_coord)[dims[1]]), 
  #                     !!sym(names(individuals_coord)[dims[2]])))
  #
  # Plot points with uniform shape for now and custom colors
  individuals_plot <- fviz_pca_ind(
    pca_result,
    axes = dims,        # Specify the axes to visualize
    geom.ind = "point", # Use points for individuals
    col.ind = group,    # Color individuals by group
    addEllipses = TRUE, # Add confidence ellipses
    ellipse.level = confidence_level, # Confidence level for ellipses
    legend.title = group_col
  ) +
    # Customize point
    scale_color_manual(values = custom_colors) +  # Set custom colors for points
    scale_shape_manual(values = custom_shape) +  # Set custom shapes for points
    scale_fill_manual(values = custom_colors) +   # Set custom colors for confidence ellipse
    labs(title = paste("PCA Individuals: Dimensions PC", dims[1], "and PC", dims[2]),
         subtitle = paste0("Confidence Level: ", confidence_level * 100, "%"),
         x = x_label,
         y = y_label) +
    theme_classic()
  
  # Display the plot
  print(individuals_plot)
  
  

  
  
  # Save the individuals plot
  individuals_file <- file.path(output_dir, paste0("pca_individuals_PC", dims[1], "_PC", dims[2], ".png"))
  ggsave(filename = individuals_file, plot = individuals_plot, width = 8, height = 6)
  
  # Plot 2: Variables (eigenvectors)
  variables_plot <- fviz_pca_var(
    pca_result,
    axes = dims,         # Specify the PCs to visualize
    col.var = "contrib", # Color variables by contribution
    #gradient.cols = c("blue", "green", "gold"), # Gradient for contributions
    gradient.cols = c("blue", "green", "orange"), # Gradient for contributions
    repel = TRUE         # Avoid overlapping labels
  ) +
    labs(title = paste("PCA Variables: Dimensions PC", dims[1], "and PC", dims[2]),
         x = x_label,
         y = y_label) +
    theme_classic()
  
  # Save the variables plot
  variables_file <- file.path(output_dir, paste0("pca_variables_PC", dims[1], "_PC", dims[2], ".png"))
  ggsave(filename = variables_file, plot = variables_plot, width = 8, height = 6)
  
  message("PCA individuals plot saved to: ", individuals_file)
  message("PCA variables plot saved to: ", variables_file)
  
  return(list(individuals_plot = individuals_plot, variables_plot = variables_plot))
}

pca_with_custom_dimensions()

#


# Step 5 ---------------------------------------
# Step 5.A. Linear discriminant analysis---------------------------------------

library(caret)
library(MASS)
library(Morpho) # The only reliable library known for performing cross-validated CVA

data <- read.csv("BDD_V16.csv",sep=";")

study_data <- read.csv("results/Step3/filtered_data.csv")

# Find optimal number of pcs according to bishop 1999 (plus we checked and for this data is best)

initial_pca <- prcomp(study_data[, -1], scale = TRUE)
initial_perc_variance <- initial_pca$sdev^2 / sum(initial_pca$sdev^2)

sum <- 0; for (pc in 1:length(initial_perc_variance)) {
  
  sum <- sum + initial_perc_variance[pc]
  
  if (sum > 0.95) {
    
    n_pcs <- pc
    break
      
  }
  
}


output_dir <- "results/Step5"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Calculate if your data is balanced or not

label_imbalance(as.factor(study_data$Technique))

#Study_data$Technique <- data$Expertise

scaling_params <- preProcess(study_data[, -1], method = c("center", "scale"))

train_test <- split_data(study_data, p = 0.8, s = 666)
train_data <- train_test$train
test_data <- train_test$test

train_scaled <- predict(scaling_params, train_data[, -1])
test_scaled <- predict(scaling_params, test_data[, -1])

train_data_pca <- prcomp(train_scaled)



train_data_pcs <- train_data_pca$x[,1:n_pcs]
test_data_pcs <- predict(train_data_pca, as.matrix(test_scaled))[,1:n_pcs]

train_data <- data.frame(train_data_pcs, Sample = as.factor(train_test$train$Technique))
test_data <- data.frame(test_data_pcs, Sample = as.factor(train_test$test$Technique))

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     classProbs = TRUE)

lda_model <- train(
  Sample ~ .,
  data = train_data,
  method = "lda", # svmRadial = support vector machine, lda = linear discriminant analysis
  trControl = ctrl
)

confusionMatrix(predict(lda_model, test_data), test_data$Sample)

# LDA without PCA

study_data <- read.csv("results/Step3/filtered_data.csv")

train_test <- split_data(study_data, p = 0.8, s = 666)
train_data <- train_test$train
test_data <- train_test$test

train_data$Technique <- as.factor(train_data$Technique)
test_data$Technique <- as.factor(test_data$Technique)

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     classProbs = TRUE)

lda_model <- train(
  Technique ~ .,
  data = train_data,
  method = "lda",
  trControl = ctrl
)

confusionMatrix(predict(lda_model, test_data), test_data$Technique)

conf_matrix <-confusionMatrix(predict(lda_model, test_data), test_data$Technique)


write.csv(conf_matrix$table, file.path(output_dir, "lda_no_pca_confusion_matrix.csv"))
write.csv(conf_matrix$overall, file.path(output_dir, "lda_no_pca_overall_statistics.csv"))
write.csv(conf_matrix$byClass, file.path(output_dir, "lda_no_pca_class_statistics.csv"))

#

# Step 5.B.Cross Validated CVA -----------------------

# The same as an LDA plot, but cross-validated which means
# it is corrected to be as reliable as possible without forcing separation

study_data <- read.csv("results/Step3/filtered_data.csv")

pca <- prcomp(scale(study_data[,2:ncol(study_data)]))  # Full PCA
pc_scores <- pca$x[,1:n_pcs]                               # Keep components that explain 95% of variance
sum(summary(pca)$importance[2,1:n_pcs])                   # Cumulative % variance of the selected PCs

non_cv_CVA <- CVA(pc_scores, study_data$Technique, cv = FALSE)
cv_CVA <- CVA(pc_scores, study_data$Technique, cv = TRUE)

custom_colors <- c("CS"="#D0C4DF", "PIPB"="#5D675B", "PIPP"="#ffa8cb", "PDPP"="skyblue1", 
                   "PDPC"="#991932", "PSc"="lightseagreen", "ScS"="darkgreen", 
                   "ScB"="olivedrab2", "PScPo"="#e31a1c", "PoHS"="yellow", 
                   "PoS"= "orange", "PoC"="darkmagenta")

custom_shape <- c(0:6,8,15:18)

output_dir <- "results/Step5"
dir.create(output_dir, showWarnings = FALSE)

conf_matrix <- confusionMatrix(cv_CVA$groups, cv_CVA$class)
write.csv(conf_matrix$table, file.path(output_dir, "cvaconfusion_matrix.csv"))  # Save the confusion matrix
write.csv(conf_matrix$overall, file.path(output_dir, "cvaoverall_statistics.csv"))  # Save overall statistics
write.csv(conf_matrix$byClass, file.path(output_dir, "cvaclass_statistics.csv"))  # Save per-class statistics


cva_data <- data.frame(
  cv1 = non_cv_CVA$CVscores[,1],
  cv2 = non_cv_CVA$CVscores[,2],
  Sample = study_data$Technique
)

ggplot2::ggplot(data = cva_data, ggplot2::aes(x = cv1, y = cv2, colour = Sample)) +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::geom_point(stat = "identity", size = 4) +
  ggplot2::stat_ellipse(size = 1) +
  #conf_interval +
  ggplot2::xlab(paste0("CV1 (", round(non_cv_CVA$Var[1,2][[1]], 2), "%)")) +
  ggplot2::ylab(paste0("CV2 (", round(non_cv_CVA$Var[2,2][[1]], 2), "%)")) +
  ggplot2::ggtitle("Non Cross Validated CVA") +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
    plot.title = ggplot2::element_text(face = "bold", size = 20),
    plot.subtitle = ggplot2::element_text(size = 15),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
    axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 13),
    legend.title = ggplot2::element_text(size = 18, face = "bold"),
    legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
    legend.box.background = ggplot2::element_rect(colour = "black"),
    legend.position = "bottom"
  ) +
  ggplot2::geom_vline(xintercept = 0,
                      colour = "black",
                      size = 0.5,
                      linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0,
                      colour = "black",
                      linetype = "dashed",
                      size = 0.5)




ggsave(file.path(output_dir, "non_cv_CVA_plot.png"), 
       plot = last_plot(), 
       width = 10, height = 8)




cvcva_data <- data.frame(
  cv1 = cv_CVA$CVcv[,1],
  cv2 = cv_CVA$CVcv[,2],
  Sample = study_data$Technique
)


ggplot2::ggplot(data = cvcva_data, ggplot2::aes(x = cv1, y = cv2, colour = Sample, shape = Sample)) +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::scale_shape_manual(values = custom_shape) +
  ggplot2::geom_point(stat = "identity", size = 2) +
  ggplot2::stat_ellipse(size = 1) +
  #conf_interval +
  ggplot2::xlab(paste0("cv-CV1 (", round(cv_CVA$Var[1,2][[1]], 2), "%)")) +
  ggplot2::ylab(paste0("cv-CV2 (", round(cv_CVA$Var[2,2][[1]], 2), "%)")) +
  ggplot2::ggtitle("Cross Validated CVA") +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
    plot.title = ggplot2::element_text(face = "bold", size = 20),
    plot.subtitle = ggplot2::element_text(size = 15),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
    axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 13),
    legend.title = ggplot2::element_text(size = 18, face = "bold"),
    legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
    legend.box.background = ggplot2::element_rect(colour = "black"),
    legend.position = "bottom"
  ) +
  ggplot2::geom_vline(xintercept = 0,
                      colour = "black",
                      size = 0.5,
                      linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0,
                      colour = "black",
                      linetype = "dashed",
                      size = 0.5)


ggsave(file.path(output_dir, "cv_CVA_plot.png"), 
       plot = last_plot(), 
       width = 10, height = 8)


confusionMatrix(cv_CVA$groups,
                cv_CVA$class)




# Step 5.C. cross-validation CVA : technique by level of expertise -----------------------

study_data <- read.csv("results/Step3/filtered_data.csv")
study_data$Expertise <- read.csv("results/Step1/data_clean.csv")$Expertise
study_data <- study_data[study_data$Technique != "CS",]
output_dir <- "results/Step5"
dir.create(output_dir, showWarnings = FALSE)
all_plots <- c()
plot_index <- 1
for (expertise in levels(as.factor(study_data$Expertise))) {
  
  technique_data <- study_data[study_data$Expertise == expertise,]
  numeric_data <- technique_data[,-c(1, ncol(technique_data))]
  pc_scores <- prcomp(scale(numeric_data))$x[,1:n_pcs]
  cv_CVA <- CVA(pc_scores, as.factor(technique_data$Technique), cv = TRUE)

  custom_colors <- c("CS"="#D0C4DF", "PIPB"="#5D675B", "PIPP"="#ffa8cb", "PDPP"="skyblue1", 
                     "PDPC"="#991932", "PSc"="lightseagreen", "ScS"="darkgreen", 
                     "ScB"="olivedrab2", "PScPo"="#e31a1c", "PoHS"="yellow", 
                     "PoS"= "orange", "PoC"="darkmagenta")
  
  custom_shape <- c(0:6,8,15:18)
  
  conf_matrix <- confusionMatrix(cv_CVA$groups, cv_CVA$class)
  write.csv(conf_matrix$table, file.path(output_dir, paste0(expertise, "_cvaconfusion_matrix.csv")))  # Sauvegarde de la confusion matrix
  write.csv(conf_matrix$overall, file.path(output_dir, paste0(expertise, "_cvaoverall_statistics.csv")))  # Sauvegarde des statistiques globales
  write.csv(conf_matrix$byClass, file.path(output_dir, paste0(expertise, "_cvaclass_statistics.csv")))  # Sauvegarde des statistiques par classe
  
  cvcva_data <- data.frame(
    cv1 = cv_CVA$CVcv[,1],
    cv2 = cv_CVA$CVcv[,2],
    Sample = as.factor(technique_data$Technique)
  )
  
  current_plot <- ggplot2::ggplot(data = cvcva_data, ggplot2::aes(x = cv1, y = cv2, colour = Sample, shape = Sample)) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_shape_manual(values = custom_shape) +
    ggplot2::geom_point(stat = "identity", size = 2) +
    ggplot2::stat_ellipse(size = 1) +
    #conf_interval +
    ggplot2::xlab(paste0("cv-CV1 (", round(cv_CVA$Var[1,2][[1]], 2), "%)")) +
    ggplot2::ylab(paste0("cv-CV2 (", round(cv_CVA$Var[2,2][[1]], 2), "%)")) +
    ggplot2::ggtitle(paste0("Cross Validated CVA ", expertise)) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 18, face = "bold"),
      legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
      legend.box.background = ggplot2::element_rect(colour = "black"),
      legend.position = "bottom"
    ) +
    ggplot2::geom_vline(xintercept = 0,
                        colour = "black",
                        size = 0.5,
                        linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0,
                        colour = "black",
                        linetype = "dashed",
                        size = 0.5)
  
  all_plots[[plot_index]] <- current_plot
  plot_index <- 1 + plot_index
  print (all_plots)
  
  #Uncomment to display in the console
  #ggsave(file.path(output_dir, paste0(expertise, "_cv_CVA_plot.png")), 
  #plot = last_plot(), 
  # width = 10, height = 8)
  
  
  
}


pc_scores <- prcomp(scale(study_data[,2:ncol(study_data)]))$x[,n_pcs]

non_cv_CVA <- CVA(pc_scores, study_data$Technique, cv = FALSE)
cv_CVA <- CVA(pc_scores, study_data$Technique, cv = TRUE)

custom_colors <- c("CS"="#D0C4DF", "PIPB"="#5D675B", "PIPP"="#ffa8cb", "PDPP"="skyblue1", 
                   "PDPC"="#991932", "PSc"="lightseagreen", "ScS"="darkgreen", 
                   "ScB"="olivedrab2", "PScPo"="#e31a1c", "PoHS"="yellow", 
                   "PoS"= "orange", "PoC"="darkmagenta")

custom_shape <- c(0:6,8,15:18)

output_dir <- "results/Step5"
dir.create(output_dir, showWarnings = FALSE)

conf_matrix <- confusionMatrix(cv_CVA$groups, cv_CVA$class)
write.csv(conf_matrix$table, file.path(output_dir, "cvaconfusion_matrix.csv"))  # Save the confusion matrix
write.csv(conf_matrix$overall, file.path(output_dir, "cvaoverall_statistics.csv"))  # Save overall statistics
write.csv(conf_matrix$byClass, file.path(output_dir, "cvaclass_statistics.csv"))  # Save per-class statistics


cva_data <- data.frame(
  cv1 = non_cv_CVA$CVscores[,1],
  cv2 = non_cv_CVA$CVscores[,2],
  Sample = study_data$Technique
)

ggplot2::ggplot(data = cva_data, ggplot2::aes(x = cv1, y = cv2, colour = Sample)) +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::geom_point(stat = "identity", size = 4) +
  ggplot2::stat_ellipse(size = 1) +
  #conf_interval +
  ggplot2::xlab(paste0("CV1 (", round(non_cv_CVA$Var[1,2][[1]], 2), "%)")) +
  ggplot2::ylab(paste0("CV2 (", round(non_cv_CVA$Var[2,2][[1]], 2), "%)")) +
  ggplot2::ggtitle("Non Cross Validated CVA") +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
    plot.title = ggplot2::element_text(face = "bold", size = 20),
    plot.subtitle = ggplot2::element_text(size = 15),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
    axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 13),
    legend.title = ggplot2::element_text(size = 18, face = "bold"),
    legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
    legend.box.background = ggplot2::element_rect(colour = "black"),
    legend.position = "bottom"
  ) +
  ggplot2::geom_vline(xintercept = 0,
                      colour = "black",
                      size = 0.5,
                      linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0,
                      colour = "black",
                      linetype = "dashed",
                      size = 0.5)




ggsave(file.path(output_dir, "non_cv_CVA_plot.png"), 
       plot = last_plot(), 
       width = 10, height = 8)




cvcva_data <- data.frame(
  cv1 = cv_CVA$CVcv[,1],
  cv2 = cv_CVA$CVcv[,2],
  Sample = study_data$Technique
)


ggplot2::ggplot(data = cvcva_data, ggplot2::aes(x = cv1, y = cv2, colour = Sample, shape = Sample)) +
  ggplot2::scale_color_manual(values = custom_colors) +
  ggplot2::scale_shape_manual(values = custom_shape) +
  ggplot2::geom_point(stat = "identity", size = 2) +
  ggplot2::stat_ellipse(size = 1) +
  #conf_interval +
  ggplot2::xlab(paste0("cv-CV1 (", round(cv_CVA$Var[1,2][[1]], 2), "%)")) +
  ggplot2::ylab(paste0("cv-CV2 (", round(cv_CVA$Var[2,2][[1]], 2), "%)")) +
  ggplot2::ggtitle("Cross Validated CVA") +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
    plot.title = ggplot2::element_text(face = "bold", size = 20),
    plot.subtitle = ggplot2::element_text(size = 15),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
    axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 13),
    legend.title = ggplot2::element_text(size = 18, face = "bold"),
    legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
    legend.box.background = ggplot2::element_rect(colour = "black"),
    legend.position = "bottom"
  ) +
  ggplot2::geom_vline(xintercept = 0,
                      colour = "black",
                      size = 0.5,
                      linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0,
                      colour = "black",
                      linetype = "dashed",
                      size = 0.5)


ggsave(file.path(output_dir, "cv_CVA_plot.png"), 
       plot = last_plot(), 
       width = 10, height = 8)


confusionMatrix(cv_CVA$groups,
                cv_CVA$class)





