
# Code written by Lloyd Austin Courtenay

# Copyright (C) 2025 Lloyd Courtenay
# SPDX-License-Identifier: AGPL-3.0

plot_profile <- function(profile) {
  
  par(mfrow = c(1,1), mar = c(5.1, 5, 4.1, 2.))
  
  xlim = range(profile[,1])
  ylim = round(range(profile[,2]), 2)
  
  plot(profile[,1], profile[,2], asp = 1, type = "l", lwd = 0,
       xlab = "", ylab = "", xaxt = "none", yaxt = "none", col ="#0073e3")
  
  grid(col = "lightgray", lty = "dotted")
  abline(h = 0)
  
  lines(profile[,1], profile[,2], lwd = 2, col ="#0073e3")
  
  mtext(side = 1, line = 3, "X (\u03BCm)", cex = 1.25, font = 2)
  mtext(side = 2, line = 3, "Z (\u03BCm)", cex = 1.25, font = 2)
  
  axis(1, seq(xlim[1], xlim[2], by = 20), font = 1, cex.axis = 1.5)
  axis(2, seq(ylim[1], ylim[2], by = 20), font = 1, cex.axis = 1.5)
  
  par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
  
}

check_directory <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
    print("Directory Created")
  }
}

create_images <- function(target_sample, sep, level = FALSE) {
  
  folder_path <- paste0(".\\", target_sample, "\\Profiles")
  target_path <- paste0(".\\", target_sample, "\\Images")
  directory_contents <- list.dirs(folder_path, recursive = FALSE)
  
  for(subfolder in directory_contents){
    
    sample <- substr(subfolder, nchar(folder_path) + 2, nchar(subfolder))
    
    full_path <- list.files(subfolder, recursive = FALSE, full.names = TRUE)
    file_name <- list.files(subfolder, recursive = FALSE, full.names = FALSE)
    
    save_directory <- paste0(target_path, "\\", sample)
    
    check_directory(save_directory)
    
    for (file in 1:length(full_path)) {
      
      profile_target <- read.table(full_path[file], head = FALSE, sep = sep)
      
      if (is.character(profile_target[,2])) {
        profile_target <- profile_target[
          profile_target[,2] != "***",
        ]
        profile_target[,2] <- as.numeric(profile_target[,2])
      }
      
      if (level) {
        profile_target = helmert_rotation(profile_target)
      }
      
      save_file_name <- paste0(
        save_directory, "\\", file, "_",
        sub(" - .*", "", file_name[file]), ".png"
      )
      
      png(save_file_name, width = 1100, height = 1100 / 2)
      
      plot_profile(profile_target)
      
      dev.off()
      
    }
    
  }
  
}


normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

euclid <- function(x1, x2, y1, y2) {
  
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(distance)
  
}

herons_formula <- function(s1, s2, s3) {
  
  p <- (s1 + s2 + s3) / 2
  a <- sqrt(p * (p - s2) * (p - s3) * (p - s1))
  D <- (2 * a) / s1
  
  return(D)
  
}

calculate_angle <- function(s1, s2, s3) {
  
  upper = s2^2 + s3^2 - s1^2
  lower = 2 * s2 * s3
  return(acos(upper / lower))
  
}

helmert_rotation <- function(profile) {
  
  n = nrow(profile)
  
  transformed_profile = as.matrix(profile)
  
  z_translation = profile[1,2]
  transformed_profile[,2] = profile[,2] - z_translation
  
  end_of_mark = transformed_profile[n,]
  
  c = euclid(0, end_of_mark[1], 0, end_of_mark[2])
  b = end_of_mark[2]
  gamma = acos(b / c)
  
  rotation_matrix = matrix(
    c(cos(gamma), -sin(gamma), sin(gamma), cos(gamma)), ncol = 2, byrow = TRUE
  )
  
  rotated_profile = rotation_matrix %*% t(as.matrix(transformed_profile))
  
  final_profile = array(numeric(), dim = c(n, 2))
  final_profile[,1] = t(rotated_profile)[,2]
  final_profile[,2] = -(t(rotated_profile)[,1])
  
  return(final_profile)
  
}

circ_to_lin <- function(angle) {
  
  return(cos(angle) + sin(angle))
  
}

calculate_asymmetry <- function(profile) {
  
  profile <- normalize(profile)
  
  deepest_point <- profile[which(profile[,2] == min(profile[,2]))[1],]
  
  profile[,1] <- profile[,1] - deepest_point[1]
  profile[,2] <- profile[,2] - deepest_point[2]
  
  left_wall <- profile[profile[,1] < 0,]
  right_wall <- profile[profile[,1] > 0,]
  
  left_wall_geom <- sf::st_linestring(as.matrix(left_wall))
  right_wall_geom <- sf::st_linestring(as.matrix(right_wall))
  left_wall <- sf::st_coordinates(sf::st_line_sample(left_wall_geom, 30, type = "regular"))[,1:2]
  right_wall <- sf::st_coordinates(sf::st_line_sample(right_wall_geom, 30, type = "regular"))[,1:2]
  
  left_wall[,1] <- abs(left_wall[,1])
  
  d <- c()
  
  for (i in 1:30) {
    
    d <- c(d, euclid(left_wall[i,1], right_wall[i,1], left_wall[i,2], right_wall[i,2]))
    
  }
  
  sigma_d2 <- sum(d^2)
  A = sqrt(sigma_d2 / 30)
  
  return(A)
  
}

calculate_parameters <- function(profile) {
  
  # n = number of steps
  # D = Depth
  # WIS = Width at Surface
  # LWIS = Left hand WIS
  # RWIS = Right hand WIS
  # OAL = Opening Angle on the Left hand side
  # OAR = Opening Angle on the Right hand side
  # OA = Opening Angle
  # OA_lin = Linear transformation of OA
  # OAL_lin = Linear transformation of OAL
  # OAR_lin = Linear transformation of OAR
  # AWIS = Asymmetry Index based on WIS values
  # AOA = Asymmetry Index based on OA values
  # A = Asymmetry Index calculated by removing axis of symmetry
  
  n = nrow(profile)
  begining = profile[n,]
  end = profile[1,]
  
  WIS = euclid(begining[1], end[1], begining[2], end[2])
  
  deepest_point <- profile[which(profile[,2] == 0)[1],]
  
  corrected_D = herons_formula(
    euclid(begining[1], end[1], begining[2], end[2]),
    euclid(begining[1], deepest_point[1], begining[2], deepest_point[2]),
    euclid(end[1], deepest_point[1], end[2], deepest_point[2])
  )
  
  D_axis = profile[which(profile[,2] == 0)[1],1]
  
  LWIS = abs(begining[1])
  RWIS = end[1]
  
  OAL = atan(LWIS / corrected_D)
  OAR = atan(RWIS / corrected_D)
  OA = calculate_angle(
    euclid(begining[1], end[1], begining[2], end[2]),
    euclid(begining[1], deepest_point[1], begining[2], deepest_point[2]),
    euclid(end[1], deepest_point[1], end[2], deepest_point[2])
  )
  
  #OA_lin = circ_to_lin(OA)
  
  AWIS = abs(LWIS - RWIS) / (LWIS + RWIS)
  AOA = abs(circ_to_lin(OAL) - circ_to_lin(OAR)) / (circ_to_lin(OAL) + circ_to_lin(OAR))
  
  A = calculate_asymmetry(profile)
  
  individual_measurement <- c(abs(corrected_D), WIS, AOA, AWIS, OA, A)
  
  return(individual_measurement)
  
}


descriptive_circular_analysis <- function(circular_object) {
  
  
  calc_mean<-function(target){
    x = numeric()
    y = numeric()
    for(i in 1:length(target)){
      y = c(y, sin(target[i][[1]]))
      x = c(x, cos(target[i][[1]]))
      Y = sum(y)/length(target) # B1
      X = sum(x)/length(target) # A1
      r = sqrt((X^2) + (Y^2)) # equivalent to rbar
      cosa = X/r
      sina = Y/r
      mr = atan2(sina,cosa) # equivalent to tbar
    }
    return(mr)
  }
  
  `%!in%` = Negate(`%in%`)
  
  target <- circular::circular(circular_object,
                               type = "directions",
                               units = "radians")
  
  # shat according to Mardia 1972
  # datapoints will be concentrated about the mean if shat is close to 0
  
  s_hat <- (
    circular::trigonometric.moment(target,
                                   p = 2,
                                   center = TRUE)$sin)/((1 - circular::trigonometric.moment(
                                     target,
                                     p = 1
                                   )$rho)**(3/2)
                                   )
  
  # khat according to Mardia 1972
  # datapoints close to 0 are flat while peaks are seen the further away you get from 0
  k_hat <- (
    circular::trigonometric.moment(target,
                                   p = 2,
                                   center = TRUE)$cos - circular::trigonometric.moment(
                                     target,
                                     p = 1)$rho**4
  )/((1 - circular::trigonometric.moment(
    target,
    p = 1)$rho)**2
    
  )
  
  # Sample circular variance (V)
  # value between 0 and 1. the closer to 0 the value the more concentrated the datapoints
  V <- (
    1 - circular::trigonometric.moment(target, p = 1)$rho
  )
  
  # delta hat is the sample circular dispersion
  
  delta_hat <- (
    1-circular::trigonometric.moment(
      target, p = 2
    )$rho)/(2*circular::trigonometric.moment(
      target, p = 1
    )$rho**2
    )
  
  Opposite <- circular::circular(target - pi, type = "directions", units = "radians")
  
  return(list(
    Standardised_Skewness = s_hat,
    Standardised_Kurtosis = k_hat,
    Sample_Circular_Variance = V,
    Sample_Circular_Dispersion = delta_hat,
    Central_Orientation_Radians = c(calc_mean(target), calc_mean(Opposite)),
    Central_Orientation_Degrees = c(calc_mean(target) * (180/pi), calc_mean(Opposite) * (180/pi))
  ))
  
}

check_nas <- function(theo, emp) {
  
  data_set <- data.frame(theo = theo, emp = emp)
  data_set <- na.omit(data_set)
  return(list(
    theo = data_set[,1],
    emp = data_set[,2]
  ))
  
}

extract_measurements <- function(target_sample, sep) {
  
  folder_path <- paste0(".\\", target_sample, "\\Profiles")
  directory_contents <- list.dirs(folder_path, recursive = FALSE)
  
  metric_results <- array(dim = c(0, 7))
  landmarks <- array(numeric(), dim = c(30, 2, 0))
  
  for(subfolder in directory_contents){
    
    sample <- substr(subfolder, nchar(folder_path) + 2, nchar(subfolder))
    
    full_path <- list.files(subfolder, recursive = FALSE, full.names = TRUE)
    file_name <- list.files(subfolder, recursive = FALSE, full.names = FALSE)
    
    for (file in 1:length(full_path)) {
      
      profile_target <- read.table(full_path[file], head = FALSE, sep = sep)
      
      if (is.character(profile_target[,2])) {
        profile_target <- profile_target[
          profile_target[,2] != "***",
        ]
        profile_target[,2] <- as.numeric(profile_target[,2])
      }
      
      profile_target = helmert_rotation(profile_target)
      
      profile_target <- remove_orientation(profile_target)
      
      profile_target <- profile_target[order(profile_target[,1]), ]
      
      individual_results <- calculate_parameters(profile_target)
      
      individual_results <- c(individual_results, sample)
      
      metric_results <- abind::abind(metric_results, individual_results, along = 1)
      
    }
    
    landmarks <- abind::abind(
      landmarks,
      compute_landmarks(subfolder, 30, sep = sep),
      along = 3
    )
    
  }
  
  metric_results <- data.frame(metric_results)
  
  colnames(metric_results) <- c("D", "WIS", "AOA", "AWIS", "OA", "A", "Sample")
  
  for (column in 1:6) {
    metric_results[,column] <- as.numeric(metric_results[,column])
  }
  metric_results$Sample <- as.factor(metric_results$Sample)
  
  return(list(
    metric_results = metric_results,
    landmarks = landmarks
  ))
  
}

bootstrap_CIs <- function(distribution, n = 1000, robust = TRUE) {
  
  values <- c()
  
  if (robust) {
    
    for (i in 1:n) {
      
      target <- sample(distribution, replace = TRUE)
      values <- c(values, median(target))
      
    }
    
    return(pValueRobust::quantile_CI(values, q = c(0.025, 0.975)))
    
  } else {
    
    for (i in 1:n) {
      
      target <- sample(distribution, replace = TRUE)
      values <- c(values, mean(target))
      
    }
    
    return(pValueRobust::quantile_CI(values, q = c(0.025, 0.975)))
    
  }
  
}

summary_stats <- function(variable) {
  
  p_val <- shapiro.test(variable)$p.value
  w <- shapiro.test(variable)$statistic[[1]]
  fpr_val <- FPR(p_val) * 100
  
  if (p_val < 0.003) {
    
    cat(
      paste0(
        
        "\n",
        "Shaprio W:\t\t\t", w, "\n",
        "Shaprio p:\t\t\t", p_val, "\n",
        "Shaprio FPR (%):\t\t", fpr_val, "\n",
        "Minimum:\t\t\t", min(variable), "\n",
        "Lower Q:\t\t\t", pValueRobust::quantile_CI(variable, q = c(0.025, 0.975))[1], "\n",
        "Lower CI:\t\t\t", bootstrap_CIs(variable)[1], "\n",
        "Median:\t\t\t\t", median(variable), "\n",
        "Biweight Midvariance:\t\t", pValueRobust::biweight_midvariance(variable), "\n",
        "Upper CI:\t\t\t", bootstrap_CIs(variable)[2], "\n",
        "Upper Q:\t\t\t", pValueRobust::quantile_CI(variable, q = c(0.025, 0.975))[2], "\n",
        "Maximum:\t\t\t", max(variable), "\n\n"
        
      )
    )
    
  } else {
    
    cat(
      paste0(
        
        "\n",
        "Shaprio W:\t\t\t", w, "\n",
        "Shaprio p:\t\t\t", p_val, "\n",
        "Shaprio FPR (%):\t\t", fpr_val, "\n",
        "Minimum:\t\t\t", min(variable), "\n",
        "Lower Q:\t\t\t", pValueRobust::quantile_CI(variable, q = c(0.025, 0.975))[1], "\n",
        "Lower CI:\t\t\t", bootstrap_CIs(variable, robust = FALSE)[1], "\n",
        "Mean:\t\t\t\t", mean(variable), "\n",
        "Standard Deviation:\t\t", sd(variable), "\n",
        "Upper CI:\t\t\t", bootstrap_CIs(variable, robust = FALSE)[2], "\n",
        "Upper Q:\t\t\t", pValueRobust::quantile_CI(variable, q = c(0.025, 0.975))[2], "\n",
        "Maximum:\t\t\t", max(variable), "\n\n"
        
      )
    )
    
  }
  
}


MinusPiPi<-function(sample) {
  n = length(sample)
  for (j in 1:n) {
    if (sample[j] < -pi) {sample[j] = sample[j]+(2*pi)} else
      if (sample[j] > pi) {sample[j] = sample[j]-(2*pi)} }
  return(sample)
}

PgVal<-function(target, ndat, g) {
  N = length(target)
  sumterms = 0 ; M = 0
  ndatcsum = cumsum(ndat)
  gmedian = circular::median.circular(target)[[1]]
  for (k in 1:g) {
    if (k==1)
    {low = 0} else
      if (k > 1)
      {low = ndatcsum[k-1]}
    sample = circular::circular(0)
    for (j in 1:ndat[k])
    {sample[j] = target[j+low]}
    shiftdat = MinusPiPi(sample - gmedian)
    m = length(shiftdat[shiftdat<0]) ; M = M+m
    sumterms = sumterms + m*m/ndat[k]
  }
  term1 = ((N*N)/(M*(N-M)))
  term2 = (N*M)/(N-M)
  Pg = term1*sumterms - term2
  return(Pg)
}

compute_landmarks <- function(file_path,
                              n_landmarks, sep) {
  
  landmark_tensor <- array(numeric(),
                           dim = c(n_landmarks, 2, 0))
  
  for(profile in (list.files(file_path))) {
    
    profile_target <- read.table(paste0(
      file_path, "\\", profile
    ), head = FALSE, sep = sep)
    
    if (is.character(profile_target[,2])) {
      profile_target <- profile_target[
        profile_target[,2] != "***",
      ]
      profile_target[,2] <- as.numeric(profile_target[,2])
    }
    
    profile_target = helmert_rotation(profile_target)
    
    profile_target <- remove_orientation(profile_target)
    
    profile_target <- profile_target[order(profile_target[,1]), ]
    
    profile_geom <- sf::st_linestring(
      as.matrix(profile_target)
    )
    
    landmarks = sf::st_line_sample(
      profile_geom, n_landmarks, type = "regular"
    )
    
    landmark_tensor <- abind::abind(
      landmark_tensor,
      sf::st_coordinates(landmarks)[,1:2], along = 3
    )
    
  }
  
  return(landmark_tensor)
  
}

remove_orientation <- function(profile) {
  
  deepest_point <- profile[which(profile[,2] == min(profile[,2]))[1],]
  
  profile[,1] <- profile[,1] - deepest_point[1]
  profile[,2] <- profile[,2] - deepest_point[2]
  
  left_side = abs(range(profile[,1]))[1]
  right_side = abs(range(profile[,1]))[2]
  
  if (right_side < left_side) {
    
    profile[,1] = -profile[,1]
    
  }
  
  return(profile)
  
}

JammalamadakaSarma_Rand <- function(x, y, n = 1000) {
  
  observed_correlation <- circular::cor.circular(x, y)
  statistic = 1
  
  for (iteration in 1:n) {
    
    x_prima <- sample(x)
    
    randomised_correlation <- circular::cor.circular(x_prima, y)
    
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
