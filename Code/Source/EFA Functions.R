
# Code written by Lloyd Austin Courtenay

# this code requires the libraries ggplot2

`%!in%` = Negate(`%in%`)

superimpose_outlines <- function(coordinate_array) {
  
  superimposed_coords <- coordinate_array
  for (coord in 1:dim(coordinate_array)[3]) {
    
    target_coords <- coordinate_array[,,coord]
    
    # center
    
    target_coords <- apply(target_coords, 2, function(x) x - mean(x))
    
    # align
    
    target_coords <- target_coords %*% svd(var(target_coords))$u
    
    superimposed_coords[,,coord] <- target_coords
    
    
  }
  
  return(superimposed_coords)
  
}

elliptic_fourier_coefficients <- function(coordinates, n_harmonics, scale = FALSE) {
  
  n_points <- nrow(coordinates)
  
  dx <- coordinates[,1] - coordinates[,1][c(n_points, (1:(n_points - 1)))]
  dy <- coordinates[,2] - coordinates[,2][c(n_points, (1:(n_points - 1)))]
  dt <- sqrt(dx^2 + dy^2)
  dt[dt < 1e-10] <- 1e-10  # to avoid Nan - recomendation Momocs
  
  tp <- cumsum(dt) # in momocs tp is called t1
  tp_1 <- c(0, tp[-n_points]) # in momocs tp_1 is called t1m1
  perim_T <- sum(dt) # in momocs this is called T (changed for compatability reasons)
  
  an <- bn <- cn <- dn <- numeric(n_harmonics)
  
  for (harmonic in 1:n_harmonics) {
    Ti <- (perim_T / (2 * pi^2 * harmonic^2))
    numerator <- 2 * harmonic * pi # in momocs this is called r
    an[harmonic] <- Ti * sum(
      (dx / dt) * (cos((numerator * tp) / perim_T) - cos((numerator * tp_1) / perim_T))
    )
    bn[harmonic] <- Ti * sum(
      (dx / dt) * (sin((numerator * tp) / perim_T) - sin((numerator * tp_1) / perim_T))
    )
    cn[harmonic] <- Ti * sum(
      (dy / dt) * (cos((numerator * tp) / perim_T) - cos((numerator * tp_1) / perim_T))
    )
    dn[harmonic] <- Ti * sum(
      (dy / dt) * (sin((numerator * tp) / perim_T) - sin((numerator * tp_1) / perim_T))
    )
  }
  
  if (scale == TRUE) {
    
    a1 <- an[1]
    b1 <- bn[1]
    c1 <- cn[1]
    d1 <- dn[1]
    
    # calculate rotation angle
    
    psi <- 0.5 * atan(
      2 * (a1 * b1 + c1 * d1) /
        (a1^2 + c1^2 - b1^2 - d1^2)
    ) %% pi # originally theta - not sure why the modular division with pi
    
    phaseshift_matrix <- matrix(
      c(cos(psi), sin(psi), -sin(psi), cos(psi)), 2, 2
    ) # originally phaseshift
    
    transposed_matrix <- matrix(c(a1, c1, b1, d1), 2, 2) %*% phaseshift_matrix # originally M2
    
    v <- apply(transposed_matrix^2, 2, sum) # in momocs but i'm not sure why
    if (v[1] < v[2]) {
      psi <- psi + pi/2
    }
    
    psi <- ((psi + pi/2) %% pi) - (pi/2)
    
    a_prima <- a1 * cos(psi) + b1 * sin(psi) # originally Aa
    c_prima <- c1 * cos(psi) + d1 * sin(psi) # originally Cc
    
    Lambda <- sqrt(a_prima^2 + c_prima^2) # originally scale
    
    theta <- atan(c_prima / a_prima) %% pi # originally psi
    
    if (a_prima < 0) {
      theta <- theta + pi
    }
    
    size_coefficient <- 1 / Lambda
    
    rotation_matrix <- matrix(
      c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2
    )
    
    norm_an <- norm_bn <- norm_cn <- norm_dn <- numeric(n_harmonics)
    
    for (harmonic in 1:n_harmonics) {
      normalised_matrix <- size_coefficient * rotation_matrix %*%
        matrix(c(an[harmonic], cn[harmonic], bn[harmonic], dn[harmonic]),
               2, 2) %*%
        matrix(c(
          cos(harmonic * psi), sin(harmonic * psi),
          -sin(harmonic * psi), cos(harmonic * psi)
        ), 2, 2)
      norm_an[harmonic] <- normalised_matrix[1,1]
      norm_bn[harmonic] <- normalised_matrix[1,2]
      norm_cn[harmonic] <- normalised_matrix[2,1]
      norm_dn[harmonic] <- normalised_matrix[2,2]
    }
    
    coe_ind <- c(norm_an, norm_bn, norm_cn, norm_dn)
    
    if (coe_ind[1] < 0) {
      coe_ind <- coe_ind * -1 # recomendation of momocs
    }
    coe_ind[abs(coe_ind) < 1e-12] = 0 # recomendation of momocs
    
  } else {
    
    coe_ind <- c(an, bn, cn, dn)
    
  }
  
  return(coe_ind)
  
}

elliptic_fourier_analysis <- function(coordinate_array, n_harmonics, scale = FALSE) {
  
  coefficients <- array(numeric(), dim = c(0, n_harmonics * 4))
  
  for (individual in 1:dim(coordinate_array)[3]) {
    
    #coordinates <- coo[individual][[1]]
    coordinates <- coordinate_array[,,individual]
    
    coe_ind <- elliptic_fourier_coefficients(coordinates,
                                             n_harmonics = n_harmonics,
                                             scale = scale)
    
    coefficients <- abind::abind(coefficients, coe_ind, along = 1)
    
  }
  
  coefficients <- as.data.frame(coefficients)
  coef_names <- c()
  for(coef_val in c("A", "B", "C", "D")) {
    for(harmonic in 1:n_harmonics) {
      coef_names <- c(coef_names, paste0(coef_val, harmonic))
    }
  }
  colnames(coefficients) <- c(coef_names)
  
  if (scale) {
    
    indices <- c()
    
    for (column in 1:ncol(coefficients)) {
      
      if (mean(coefficients[,column]) == median(coefficients[,column])) {
        
        indices <- c(indices, column)
        
      }
      
    }
    
    coefficients <- coefficients[, -indices]
    
  }
  
  return(coefficients)
  
}

harmonic_power <- function(an, bn, cn, dn) {
  power <- (an^2 + bn^2 + cn^2 + dn^2)/2
  return(power)
}

calculate_harmonic_power <- function(coordinate_array, lower_limits = 0.95, upper_limits = 0.99) {
  
  #possible_n_harmonics <- floor(nrow(coo_object[1]) / 2) - 1
  possible_n_harmonics <- floor(dim(coordinate_array)[1] / 2) - 1
  
  results <- array(numeric(), dim = c(0, (possible_n_harmonics - 1)))
  x <- 2:possible_n_harmonics
  
  for (individual in 1:dim(coordinate_array)[3]) {
    x_coefs <- elliptic_fourier_coefficients(coordinate_array[,,individual],
                                             n_harmonics = possible_n_harmonics,
                                             scale = TRUE)
    
    an <- x_coefs[1:possible_n_harmonics]
    bn <- x_coefs[(possible_n_harmonics + 1):(possible_n_harmonics * 2)]
    cn <- x_coefs[((possible_n_harmonics * 2) + 1):(possible_n_harmonics * 3)]
    dn <- x_coefs[((possible_n_harmonics * 3) + 1):(possible_n_harmonics * 4)]
    
    power <- harmonic_power(an, bn, cn, dn)[x]
    results <- abind::abind(results, power, along = 1)
    
  }
  
  results <- t(apply(results, 1, function(pow) {cumsum(pow) / sum(pow)})) * 100
  
  opt_rank <- ceiling(
    (find_optimal_harmonic(results, lower_limits) + find_optimal_harmonic(results, upper_limits)) / 2
  )
  
  decomp_c <- c(); decomp_samp <- c(); for(i in 1:ncol(results)) {
    decomp_c <- c(decomp_c, results[,i])
    decomp_samp <- c(decomp_samp, rep(i, nrow(results)))
  }
  
  power_plot_data <- data.frame(
    Harmonic_Power = decomp_c,
    Harmonic_Rank = as.factor(decomp_samp)
  )
  
  result_plot <- ggplot2::ggplot(data = power_plot_data, ggplot2::aes(x = Harmonic_Rank, y = Harmonic_Power)) +
    ggplot2::geom_boxplot(fill = "lightgrey", size = 1) +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(xintercept = opt_rank, col = "red", size = 1) +
    ggplot2::xlab("Harmonic") +
    ggplot2::ylab("Power") +
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"),
      axis.title = ggplot2::element_text(face = "bold", size = 20),
      axis.text = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, b = 5)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(l = 10))
    )
  
  power_value <- median(power_plot_data$Harmonic_Power[power_plot_data$Harmonic_Rank == opt_rank])
  
  cat(paste0(opt_rank, " harmonics has a power of ", power_value, " %"))
  
  return(result_plot)
  
}

find_optimal_harmonic <- function(harmonic_results, target_threshold = 0.95) {
  
  median_power <- apply(harmonic_results, 2, median)
  
  return(which(median_power > target_threshold * 100)[1])
  
}

plot_residuals <- function(lmmodel, label_values = NULL, plot_labels = TRUE){
  
  # modification of the plot.lm.rrpp() function from RRPP library
  
  residuals <- lmmodel$LM$residuals
  fitted <- lmmodel$LM$fitted
  
  residuals <- scale(residuals, scale = FALSE)
  fitted <- scale(fitted, scale = FALSE)
  residuals <- sqrt(diag(tcrossprod(residuals)))
  fitted <- sqrt(diag(tcrossprod(fitted)))
  
  if(length(unique(round(fitted, 7))) <= 2) {
    
    lfr <- list()
    lfr$x <- fitted
    lfr$y <- scale(residuals)
    fit <- lm(scale(residuals) ~ fitted)
    lfr$fitted <- fit$fitted.values
    lfr$residuals <- fit$residuals
    
  } else {
    
    options(warn = -1)
    lfr <- loess(scale(residuals)~fitted, span = 1)
    options(warn = 0)
    
  }
  
  lfr <- cbind(lfr$x, lfr$y, lfr$fitted)
  lfr <- lfr[order(lfr[,1]),]
  
  lfr <- as.data.frame(lfr)
  
  main_plot <- ggplot2::ggplot(
    data = lfr
  ) + 
    ggplot2::geom_line(ggplot2::aes(x = fitted, y = V3), color = "red", size = 1) +
    ggplot2::geom_hline(yintercept = 0, col = "grey", size = 1,
                        linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(x = fitted, y = V2), size = 3) +
    ggplot2::xlab("Standardized Fitted Values") +
    ggplot2::ylab("Standardized Residuals") +
    ggplot2::ggtitle("Residuals vs Fitted") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = ggplot2::element_text(size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
    )
  
  if (plot_labels) {
    
    if (is.null(label_values)) {
      stop("label_values have not been provided")
    } else {
      
      xlim <- range(lfr$fitted)
      ylim <- range(lfr$V2)
      xlim[1] <- xlim[1] - (sd(lfr$fitted) / 4)
      xlim[2] <- xlim[2] + (sd(lfr$fitted) / 4)
      ylim[1] <- ylim[1] - (sd(lfr$V2) / 4)
      ylim[2] <- ylim[2] + (sd(lfr$V2) / 4)
      
      final_plot <- main_plot + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
        ggplot2::geom_text(aes(label = label_values,
                               x = fitted, y = V2),
                           nudge_x = 0.005, nudge_y = 0.005)
      
    }
    
    return(final_plot)
    
  } else {
    
    return(main_plot)
    
  }
  
}

plot_residuals_qq <- function(lmmodel, label_values = NULL, plot_labels = TRUE){
  
  residuals <- lmmodel$LM$residuals
  
  residuals <- sqrt(diag(tcrossprod(scale(residuals, scale = FALSE))))
  residuals <- sort(residuals)
  n <- length(residuals)
  tq <- (seq(1, n) - 0.5) / n
  tq <- qnorm(tq)
  
  qq_data <- data.frame(x = tq, y = residuals)
  
  main_plot <- ggplot2::ggplot(
    data = qq_data
  ) + 
    ggplot2::geom_point(ggplot2::aes(x = x, y = y), size = 3) +
    ggplot2::xlab("Theoretical Quantiles") +
    ggplot2::ylab("Standardized Residuals") +
    ggplot2::ggtitle("Q-Q Plot") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = ggplot2::element_text(size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
    )
  
  if (plot_labels) {
    
    if (is.null(label_values)) {
      stop("label_values have not been provided")
    } else {
      
      xlim <- range(qq_data$x)
      ylim <- range(qq_data$y)
      xlim[1] <- xlim[1] - (sd(qq_data$x) / 4)
      xlim[2] <- xlim[2] + (sd(qq_data$x) / 4)
      ylim[1] <- ylim[1] - (sd(qq_data$y) / 4)
      ylim[2] <- ylim[2] + (sd(qq_data$y) / 4)
      
      final_plot <- main_plot + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
        ggplot2::geom_text(aes(label = label_values,
                               x = x, y = y),
                           nudge_x = 0.005, nudge_y = 0.005)
      
    }
    
    print(shapiro.test(residuals))
    
    return(final_plot)
    
  } else {
    
    print(shapiro.test(residuals))
    
    return(main_plot)
    
  }
  
}

plot_calibration_curves <- function(p, priors = 0.2, limits = c(1, 1e-05)) {
  
  par(mfrow = c(1,2))
  pValueRobust::calibration_curve(method = "FPR", priors = priors, lwd = 2, col = "black",
                                  p_values = limits)
  #abline(v = p, col = "orange", lwd = 2)
  points(p, pValueRobust::FPR(p, priors = priors), pch = 19, col = "red", cex = 2)
  abline(h = 0.1, col = "blue", lwd = 2)
  pValueRobust::calibration_curve(method = "p_BFB vs FPR", priors = priors, lwd = 2,
                                  p_values = limits)
  points(p, pValueRobust::FPR(p, priors = priors), pch = 19, col = "red", cex = 2)
  abline(h = 0.1, col = "blue", lwd = 2)
  par(mfrow = c(1,1))
  
}

visualize_boxplot <- function(variable, samples, variable_name,
                              ylim = NULL, col = "lightgrey") {
  
  plot_data <- data.frame(
    variable = variable,
    samples = samples
  )
  
  if (is.null(ylim)) {
    ggplot2::ggplot(data = plot_data, ggplot2::aes(x = samples, y = variable)) +
      ggplot2::geom_boxplot(fill = col, size = 1, outlier.shape = NA) +
      ggplot2::theme_bw() +
      ggplot2::xlab("Sample") +
      ggplot2::ylab(variable_name) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"),
        axis.title = ggplot2::element_text(face = "bold", size = 20),
        axis.text = ggplot2::element_text(size = 15),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, b = 5)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(l = 10, r = 10))
      )
  } else{
    ggplot2::ggplot(data = plot_data, ggplot2::aes(x = samples, y = variable)) +
      ggplot2::geom_boxplot(fill = "lightgrey", size = 1, outlier.shape = NA) +
      ggplot2::theme_bw() +
      ggplot2::xlab("Sample") +
      ggplot2::ylab(variable_name) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"),
        axis.title = ggplot2::element_text(face = "bold", size = 20),
        axis.text = ggplot2::element_text(size = 15),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, b = 5)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(l = 10, r = 10))
      ) +
      ggplot2::coord_cartesian(ylim = ylim)
  }
  
}

add_alpha <- function(col, alpha=1){ # Function for setting the alpha of colours
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}


visualise_shape_change <- function(coordinates, pc_scores, gmm = TRUE, size_correct = FALSE) {
  
  min_pc_1 <- min(pc_scores[,1])
  min_pc_2 <- min(pc_scores[,2])
  max_pc_1 <- max(pc_scores[,1])
  max_pc_2 <- max(pc_scores[,2])
  
  min_pc_1_individual <- pc_scores[pc_scores[,1] == min_pc_1, 1:2]
  min_pc_2_individual <- pc_scores[pc_scores[,2] == min_pc_2, 1:2]
  max_pc_1_individual <- pc_scores[pc_scores[,1] == max_pc_1, 1:2]
  max_pc_2_individual <- pc_scores[pc_scores[,2] == max_pc_2, 1:2]
  
  negative_pc_1 <- GraphGMM::morphological_predictor(
    coordinates,
    pc_scores[,1:2],
    min_pc_1_individual
  )
  
  positive_pc_1 <- GraphGMM::morphological_predictor(
    coordinates,
    pc_scores[,1:2],
    max_pc_1_individual
  )
  
  negative_pc_2 <- GraphGMM::morphological_predictor(
    coordinates,
    pc_scores[,1:2],
    min_pc_2_individual
  )
  
  positive_pc_2 <- GraphGMM::morphological_predictor(
    coordinates,
    pc_scores[,1:2],
    max_pc_2_individual
  )
  
  par(mfrow = c(2,2))
  
  if (!size_correct) {
    
    plot_outline(negative_pc_1)
    
    title(main = paste0("-PC", 1))
    
    plot_outline(positive_pc_1)
    
    title(main = paste0("PC", 1))
    
    plot_outline(negative_pc_2)
    
    title(main = paste0("-PC", 2))
    
    plot_outline(positive_pc_2)
    
    title(main = paste0("PC", 2))
    
  } else {
    
    xlim_pc1 = range(rbind(negative_pc_1, positive_pc_1)[,1])
    ylim_pc1 = range(rbind(negative_pc_1, positive_pc_1)[,2])
    xlim_pc2 = range(rbind(negative_pc_2, positive_pc_2)[,1])
    ylim_pc2 = range(rbind(negative_pc_2, positive_pc_2)[,2])
    
    plot_outline(negative_pc_1, xlim_pc1, ylim_pc1)
    
    title(main = paste0("-PC", 1))
    
    plot_outline(positive_pc_1, xlim_pc1, ylim_pc1)
    
    title(main = paste0("PC", 1))
    
    plot_outline(negative_pc_2, xlim_pc2, ylim_pc2)
    
    title(main = paste0("-PC", 2))
    
    plot_outline(positive_pc_2, xlim_pc2, ylim_pc2)
    
    title(main = paste0("PC", 2))
      
  }
  
  par(mfrow = c(1,1))
  
}

plot_outline <- function(coordinates, xlim = NULL, ylim = NULL) {
  
  if (is.null(xlim) | is.null(ylim)) {
    
    plot(coordinates[,1], coordinates[,2], asp = 1, col = NULL, bty = "n",
         xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    
  } else {
    
    plot(coordinates[,1], coordinates[,2], asp = 1, col = NULL, bty = "n",
         xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         xlim = xlim, ylim = ylim)
    
  }
  
  lines(coordinates[,1], coordinates[,2], lwd = 3)
  
}

#
