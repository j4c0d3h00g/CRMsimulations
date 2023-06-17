# Cellwise Robust M-regression (CRM) - Simulation Study Efficiency Levels -------------------------
rm(list = ls())

# Setup -------------------------------------------------------------------------------------------

# simulation setting - - - - - - - - - - - - - - -
n <- 400                                          # number of cases
p <- 50                                           # number of predictor variables
pct_case_out <- 0.05                              # percentage of casewise outliers
pct_cell_out <- 0.10                              # percentage of cellwise outliers for each casewise outlier
n_sims <- 50                                      # number of simulations
eff_levels <- c(0.85, 0.90, 0.95, 0.975, 0.99)    # efficiency levels weight functions
n_eff_levels <- length(eff_levels)                # number of evaluated efficiency levels
data_list_clean <- list()                         # list to store clean data of each iteration
data_list_contaminated <- list()                  # list to store contaminated data of each iteration
betas_list <- list()                              # list to store betas of each iteration
case_outliers_list <- list()                      # list to store casewise outliers of each iteration
outliers_list <- list()                           # list to store contaminated cells of each iteration


# CRM input parameters - - - - - - - - - - - - - -
maxiter             <- 100
tolerance           <- 0.01
outlyingness.factor <- 1.5
spadieta            <- seq(0.9, 0.1, -0.1)


# Load packages -----------------------------------------------------------------------------------
library(MASS)
library(cellWise)
library(lubridate)
library(robustbase)
library(RobStatTM)
library(rJava)
library(fsdaR)
library(CRMwf)


# Start simulation procedure ----------------------------------------------------------------------
resultmatrix <- data.frame(weightfunction = rep(c("Hampel", "Tukey", "Huber", "Gauss", "Quadratic"), each = 5),
                           efficiency = rep(c(0.85, 0.90, 0.95, 0.975, 0.99), times = 5))
resultmatrix$parameters = NA
resultmatrix$MSEP = NA
resultmatrix$MAE = NA
resultmatrix$precision = NA
resultmatrix$recall = NA

results_MAE         <- data.frame(matrix(nrow = n_sims, ncol = 5))
results_MSEP        <- data.frame(matrix(nrow = n_sims, ncol = 5))
results_precision   <- data.frame(matrix(nrow = n_sims, ncol = 5))
results_recall      <- data.frame(matrix(nrow = n_sims, ncol = 5))
results_time        <- data.frame(matrix(nrow = n_sims, ncol = 5))
results_avg_time    <- data.frame(matrix(nrow = n_eff_levels, ncol = 5))

names(results_MAE)          <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
names(results_MSEP)         <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
names(results_precision)    <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
names(results_recall)       <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
names(results_time)         <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
names(results_avg_time)     <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")


t_start <- proc.time()
set.seed(2019)
cat(paste("\n*", n_sims ,"simulations - started at", format(Sys.time(), "%X"),
          "============================================================\n\n"))


for (j_sim in 1:n_sims) {
  
  # Generate clean sample -------------------------------------------------------------------------
  mu <- rep(0, p)
  Sigma <- diag(p)
  Sigma[(row(Sigma) - col(Sigma)) == 1 | row(Sigma) - col(Sigma) == -1] <- 0.5
  X <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  
  slopes <- rnorm(n = p, mean = 0, sd = 1)
  slopes <- 10 * slopes / sqrt(sum(slopes^2))
  intercept <- 10
  noise <- rnorm(n = n, mean = 0, sd = 0.5)
  
  y <- intercept + X %*% slopes + noise
  
  betas <- c(intercept, slopes)
  betas_list[[j_sim]] <- betas
  
  
  # Add contamination in design matrix ------------------------------------------------------------
  Xc <- X
  contamination <- colMeans(X) + 6 * apply(X, 2, sd)
  
  case_outliers <- sample(n, size = n * pct_case_out)
  case_outliers_list[[j_sim]] <- case_outliers
  
  outliers_mat_flag <- matrix(FALSE, nrow = n, ncol = p)
  
  for (i_outliers in case_outliers) {
    cell_outliers <- sample(p, size = p * pct_cell_out)
    Xc[i_outliers, cell_outliers] <- contamination[cell_outliers] + rnorm(length(cell_outliers))
    outliers_mat_flag[i_outliers, cell_outliers] <- TRUE
  }
  outliers_list[[j_sim]] <- outliers_mat_flag
  
  
  # Collect data samples --------------------------------------------------------------------------
  data_clean        <- cbind.data.frame(y, X)
  data_contaminated <- cbind.data.frame(y, Xc)
  names(data_clean) <- names(data_contaminated) <- c("y", paste0("X", 1:p))
  
  # Save data of iteration such that the data are the same for each efficiency level --------------
  data_list_clean[[j_sim]] <- data_clean
  data_list_contaminated[[j_sim]] <- data_contaminated
}


i <- 1
for (eff_level in eff_levels) {
  cat(paste0("\n* Evaluated efficiency level ", eff_level, " (", i, "/", n_eff_levels, ")", " =======================================\n"))
  
  # Derive parameters for the weight functions corresponding to the efficiency level --------------
  parameters_hampel       <- psifun(eff = eff_level, fun = "hampel")$c1
  parameters_hampel       <- parameters_hampel[1] * parameters_hampel[2:4]
  parameters_tukey        <- bisquare(eff_level)
  parameters_huber        <- huber(eff_level)
  parameters_gauss        <- .psi.ggw.findc(ms = -0.5, b = 1.5, eff = eff_level)
  parameters_gauss[1]     <- parameters_gauss[4]
  parameters_gauss        <- parameters_gauss[1:3]
  parameters_quadratic    <- .psi.lqq.findc(ms = -0.5, b.c = 1.5, eff = eff_level)
  parameters_quadratic    <- parameters_quadratic[c(2, 1, 3)]
  
  # Save parameters in the resul tmatrix -----------------------------------------------------------
  resultmatrix$parameters[i]        <- list(c(parameters_hampel))
  resultmatrix$parameters[i+5]      <- parameters_tukey
  resultmatrix$parameters[i+10]     <- parameters_huber
  resultmatrix$parameters[i+15]     <- list(c(parameters_gauss))
  resultmatrix$parameters[i+20]     <- list(c(parameters_quadratic))
  
  
  for (j_sim in 1:n_sims) {
    cat(paste0("\n* Simulation ", j_sim, "/", n_sims, " =======================================\n"))
    
    # Fit cellwise robust M regression models ------------------------------------------------------
    crm_hampel <- suppressWarnings(crm_functional(formula   = y ~ .,
                                                  data      = data_list_contaminated[[j_sim]],
                                                  maxiter   = maxiter,
                                                  tolerance = tolerance,
                                                  outlyingness.factor = outlyingness.factor,
                                                  spadieta  = spadieta,
                                                  center    = "median",
                                                  scale     = "qn",
                                                  regtype   = "MM",
                                                  verbose   = FALSE,
                                                  weightfunction = "hampel",
                                                  weightthreshold = 1,
                                                  parameters = parameters_hampel))
    
    crm_tukey <- suppressWarnings(crm_functional(formula   = y ~ .,
                                                 data      = data_list_contaminated[[j_sim]],
                                                 maxiter   = maxiter,
                                                 tolerance = tolerance,
                                                 outlyingness.factor = outlyingness.factor,
                                                 spadieta  = spadieta,
                                                 center    = "median",
                                                 scale     = "qn",
                                                 regtype   = "MM",
                                                 verbose   = FALSE,
                                                 weightfunction = "tukey",
                                                 weightthreshold = 0.7,
                                                 parameters = parameters_tukey))
    
    crm_huber <- suppressWarnings(crm_functional(formula   = y ~ .,
                                                 data      = data_list_contaminated[[j_sim]],
                                                 maxiter   = maxiter,
                                                 tolerance = tolerance,
                                                 outlyingness.factor = outlyingness.factor,
                                                 spadieta  = spadieta,
                                                 center    = "median",
                                                 scale     = "qn",
                                                 regtype   = "MM",
                                                 verbose   = FALSE,
                                                 weightfunction = "huber",
                                                 weightthreshold = 1,
                                                 parameters = parameters_huber))
    
    crm_gauss <- suppressWarnings(crm_functional(formula   = y ~ .,
                                                 data      = data_list_contaminated[[j_sim]],
                                                 maxiter   = maxiter,
                                                 tolerance = tolerance,
                                                 outlyingness.factor = outlyingness.factor,
                                                 spadieta  = spadieta,
                                                 center    = "median",
                                                 scale     = "qn",
                                                 regtype   = "MM",
                                                 verbose   = FALSE,
                                                 weightfunction = "gauss",
                                                 weightthreshold = 1,
                                                 parameters = parameters_gauss))
    
    crm_quadratic <- suppressWarnings(crm_functional(formula   = y ~ .,
                                                     data      = data_list_contaminated[[j_sim]],
                                                     maxiter   = maxiter,
                                                     tolerance = tolerance,
                                                     outlyingness.factor = outlyingness.factor,
                                                     spadieta  = spadieta,
                                                     center    = "median",
                                                     scale     = "qn",
                                                     regtype   = "MM",
                                                     verbose   = FALSE,
                                                     weightfunction = "quadratic",
                                                     weightthreshold = 1,
                                                     parameters = parameters_quadratic))


  # Evaluate performance --------------------------------------------------------------------------
  
  # Mean Absolute Error - - - - - - - - - - - - -
  results_MAE$'CRM-Hampel'[j_sim]       <- mean(abs(crm_hampel$coefficients - betas_list[[j_sim]]))
  results_MAE$'CRM-Tukey'[j_sim]        <- mean(abs(crm_tukey$coefficients - betas_list[[j_sim]]))
  results_MAE$'CRM-Huber'[j_sim]        <- mean(abs(crm_huber$coefficients - betas_list[[j_sim]]))
  results_MAE$'CRM-Gauss'[j_sim]        <- mean(abs(crm_gauss$coefficients - betas_list[[j_sim]]))
  results_MAE$'CRM-Quadratic'[j_sim]    <- mean(abs(crm_quadratic$coefficients - betas_list[[j_sim]]))
  
  
  # Mean Squared Error of Prediction - - - - - - -
  results_MSEP$'CRM-Hampel'[j_sim]       <- mean((crm_hampel$fitted.values - data_list_clean[[j_sim]]$y)[-case_outliers_list[[j_sim]]]^2)
  results_MSEP$'CRM-Tukey'[j_sim]        <- mean((crm_tukey$fitted.values - data_list_clean[[j_sim]]$y)[-case_outliers_list[[j_sim]]]^2)
  results_MSEP$'CRM-Huber'[j_sim]        <- mean((crm_huber$fitted.values - data_list_clean[[j_sim]]$y)[-case_outliers_list[[j_sim]]]^2)
  results_MSEP$'CRM-Gauss'[j_sim]        <- mean((crm_gauss$fitted.values - data_list_clean[[j_sim]]$y)[-case_outliers_list[[j_sim]]]^2)
  results_MSEP$'CRM-Quadratic'[j_sim]    <- mean((crm_quadratic$fitted.values - data_list_clean[[j_sim]]$y)[-case_outliers_list[[j_sim]]]^2)
  
  
  # Precision & recall of detected cells - - - - -
  CM_hampel     <- table(Prediction = ifelse(c(crm_hampel$cellwiseoutliers) != 0, TRUE, FALSE),
                         Reference  = c(outliers_list[[j_sim]]))
  CM_tukey      <- table(Prediction = ifelse(c(crm_tukey$cellwiseoutliers) != 0, TRUE, FALSE),
                         Reference  = c(outliers_list[[j_sim]]))
  CM_huber      <- table(Prediction = ifelse(c(crm_huber$cellwiseoutliers) != 0, TRUE, FALSE),
                         Reference  = c(outliers_list[[j_sim]]))
  CM_gauss      <- table(Prediction = ifelse(c(crm_gauss$cellwiseoutliers) != 0, TRUE, FALSE),
                         Reference  = c(outliers_list[[j_sim]]))
  CM_quadratic  <- table(Prediction = ifelse(c(crm_quadratic$cellwiseoutliers) != 0, TRUE, FALSE),
                         Reference  = c(outliers_list[[j_sim]]))
  
  results_precision$'CRM-Hampel'[j_sim]     <- CM_hampel[2, 2] / (CM_hampel[2, 1] + CM_hampel[2, 2])
  results_precision$'CRM-Tukey'[j_sim]      <- CM_tukey[2, 2] / (CM_tukey[2, 1] + CM_tukey[2, 2])
  results_precision$'CRM-Huber'[j_sim]      <- CM_huber[2, 2] / (CM_huber[2, 1] + CM_huber[2, 2])
  results_precision$'CRM-Gauss'[j_sim]      <- CM_gauss[2, 2] / (CM_gauss[2, 1] + CM_gauss[2, 2])
  results_precision$'CRM-Quadratic'[j_sim]  <- CM_quadratic[2, 2] / (CM_quadratic[2, 1] + CM_quadratic[2, 2])
  
  results_recall$'CRM-Hampel'[j_sim]        <- CM_hampel[2, 2] / (CM_hampel[1, 2] + CM_hampel[2, 2])
  results_recall$'CRM-Tukey'[j_sim]         <- CM_tukey[2, 2] / (CM_tukey[1, 2] + CM_tukey[2, 2])
  results_recall$'CRM-Huber'[j_sim]         <- CM_huber[2, 2] / (CM_huber[1, 2] + CM_huber[2, 2])
  results_recall$'CRM-Gauss'[j_sim]         <- CM_gauss[2, 2] / (CM_gauss[1, 2] + CM_gauss[2, 2])
  results_recall$'CRM-Quadratic'[j_sim]     <- CM_quadratic[2, 2] / (CM_quadratic[1, 2] + CM_quadratic[2, 2])
  
  
  # Execution time - - - - - - - - - - - - - - - -
  results_time$'CRM-Hampel'[j_sim]      <- crm_hampel$time
  results_time$'CRM-Tukey'[j_sim]       <- crm_tukey$time
  results_time$'CRM-Huber'[j_sim]       <- crm_huber$time
  results_time$'CRM-Gauss'[j_sim]       <- crm_gauss$time
  results_time$'CRM-Quadratic'[j_sim]   <- crm_quadratic$time
  
  
  # Elapsed time - - - - - - - - - - - - - - - - -
  t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
  cat(paste0(" - ", round(100 * (j_sim + (i - 1) * n_sims)  / (n_sims * n_eff_levels), 2),
             "% - ", sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))
  }
  
  
  resultmatrix$MSEP[i]       <- round(mean(results_MSEP$'CRM-Hampel'), 4)
  resultmatrix$MSEP[i+5]     <- round(mean(results_MSEP$'CRM-Tukey'), 4)
  resultmatrix$MSEP[i+10]    <- round(mean(results_MSEP$'CRM-Huber'), 4)
  resultmatrix$MSEP[i+15]    <- round(mean(results_MSEP$'CRM-Gauss'), 4)
  resultmatrix$MSEP[i+20]    <- round(mean(results_MSEP$'CRM-Quadratic'), 4)
  
  resultmatrix$MAE[i]        <- round(mean(results_MAE$'CRM-Hampel'), 4)
  resultmatrix$MAE[i+5]      <- round(mean(results_MAE$'CRM-Tukey'), 4)
  resultmatrix$MAE[i+10]     <- round(mean(results_MAE$'CRM-Huber'), 4)
  resultmatrix$MAE[i+15]     <- round(mean(results_MAE$'CRM-Gauss'), 4)
  resultmatrix$MAE[i+20]     <- round(mean(results_MAE$'CRM-Quadratic'), 4)
  
  resultmatrix$precision[i]       <- round(mean(results_precision$'CRM-Hampel'), 4)
  resultmatrix$precision[i+5]     <- round(mean(results_precision$'CRM-Tukey'), 4)
  resultmatrix$precision[i+10]    <- round(mean(results_precision$'CRM-Huber'), 4)
  resultmatrix$precision[i+15]    <- round(mean(results_precision$'CRM-Gauss'), 4)
  resultmatrix$precision[i+20]    <- round(mean(results_precision$'CRM-Quadratic'), 4)
  
  resultmatrix$recall[i]          <- round(mean(results_recall$'CRM-Hampel'), 4)
  resultmatrix$recall[i+5]        <- round(mean(results_recall$'CRM-Tukey'), 4)
  resultmatrix$recall[i+10]       <- round(mean(results_recall$'CRM-Huber'), 4)
  resultmatrix$recall[i+15]       <- round(mean(results_recall$'CRM-Gauss'), 4)
  resultmatrix$recall[i+20]       <- round(mean(results_recall$'CRM-Quadratic'), 4)
  
  results_avg_time$'CRM-Hampel'[i]        <- round(mean(results_time$'CRM-Hampel'), 4)
  results_avg_time$'CRM-Tukey'[i]         <- round(mean(results_time$'CRM-Tukey'), 4)
  results_avg_time$'CRM-Huber'[i]         <- round(mean(results_time$'CRM-Huber'), 4)
  results_avg_time$'CRM-Gauss'[i]         <- round(mean(results_time$'CRM-Gauss'), 4)
  results_avg_time$'CRM-Quadratic'[i]     <- round(mean(results_time$'CRM-Quadratic'), 4)
  
  i <- i + 1
}


# Create matrix containing the performance evaluation of CRM with the Hampel weight function, 
# where the weights correspond to the quantiles of the standard normal distribution
results_hampel_norm       <- data.frame(matrix(nrow = n_sims, ncol = 5))
names(results_hampel_norm)   <- c("MSEP", "MAE", "precision", "recall", "time")

t_start <- proc.time()
cat(paste("\n*", n_sims ,"simulations - started at", format(Sys.time(), "%X"),
          "============================================================\n\n"))

# Evaluate CRM using the Hampel weight function with the standard normal quantiles as weights -----
for (j_sim in 1:n_sims) {
  
  # Fit regression models -------------------------------------------------------------------------
  crm_hampel_norm <- suppressWarnings(crm_functional(formula   = y ~ .,
                                                data      = data_list_contaminated[[j_sim]],
                                                maxiter   = maxiter,
                                                tolerance = tolerance,
                                                outlyingness.factor = outlyingness.factor,
                                                spadieta  = spadieta,
                                                center    = "median",
                                                scale     = "qn",
                                                regtype   = "MM",
                                                verbose   = FALSE,
                                                weightfunction = "hampel",
                                                weightthreshold = 1,
                                                parameters = c(qnorm(0.95), qnorm(0.975), qnorm(0.999))))
  
  
  # Mean Absolute Error - - - - - - - - - - - - -
  results_hampel_norm$'MSEP'[j_sim]  <- mean((crm_hampel_norm$fitted.values - data_list_clean[[j_sim]]$y)[-case_outliers_list[[j_sim]]]^2)
  
  # Mean Squared Error of Prediction - - - - - - -
  results_hampel_norm$'MAE'[j_sim]   <- mean(abs(crm_hampel_norm$coefficients - betas_list[[j_sim]]))
  
  # Precision & recall of detected cells - - - - -
  CM_hampel_norm                              <- table(Prediction = ifelse(c(crm_hampel_norm$cellwiseoutliers) != 0, TRUE, FALSE),
                                                                            Reference  = c(outliers_list[[j_sim]]))
  results_hampel_norm$'precision'[j_sim]      <- CM_hampel_norm[2, 2] / (CM_hampel_norm[2, 1] + CM_hampel_norm[2, 2])
  results_hampel_norm$'recall'[j_sim]         <- CM_hampel_norm[2, 2] / (CM_hampel_norm[1, 2] + CM_hampel_norm[2, 2])
  
  # Execution time - - - - - - - - - - - - - - - -
  results_hampel_norm$'time'[j_sim]   <- crm_hampel_norm$time 
  
  # Elapsed time - - - - - - - - - - - - - - - - -
  t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
  cat(paste0(" - ", round(100 * j_sim / n_sims, 2),
             "% - ", sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))
}

average_results_hampel_norm               <- data.frame(matrix(nrow = 1, ncol = 5))
names(average_results_hampel_norm)        <- c("MSEP", "MAE", "precision", "recall", "time")

average_results_hampel_norm$'MSEP'        <- round(mean(results_hampel_norm$'MSEP'), 4)
average_results_hampel_norm$'MAE'         <- round(mean(results_hampel_norm$'MAE'), 4)
average_results_hampel_norm$'precision'   <- round(mean(results_hampel_norm$'precision'), 4)
average_results_hampel_norm$'recall'      <- round(mean(results_hampel_norm$'recall'), 4)
average_results_hampel_norm$'time'        <- round(mean(results_hampel_norm$'time'), 4)


t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
cat(paste("\nTime elapsed:",
          sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))
