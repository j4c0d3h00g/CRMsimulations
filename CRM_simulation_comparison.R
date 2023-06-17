# Cellwise Robust M-regression (CRM) - Simulation Study Initial Comparison Weight Functions -------
rm(list = ls())

# Setup -------------------------------------------------------------------------------------------

# simulation setting - - - - - - - - - - - - - - -
n <- 400                            # number of cases
p <- 50                             # number of predictor variables
pct_case_out <- 0.05                # percentage of casewise outliers
pct_cell_out <- 0.10                # percentage of cellwise outliers for each casewise outlier
n_sims <- 100                       # number of simulations
data_list_clean <- list()           # list to store clean data of each iteration
data_list_contaminated <- list()    # list to store contaminated data of each iteration
betas_list <- list()                # list to store betas of each iteration
case_outliers_list <- list()        # list to store casewise outliers of each iteration
outliers_list <- list()             # list to store contaminated cells of each iteration


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
library(CRMwf)

# Start simulation procedure ----------------------------------------------------------------------
results_MAE         <- data.frame(matrix(nrow = n_sims, ncol = 5))
results_MSEP        <- data.frame(matrix(nrow = n_sims, ncol = 5))
results_precision   <- data.frame(matrix(nrow = n_sims, ncol = 5))
results_recall      <- data.frame(matrix(nrow = n_sims, ncol = 5))
results_time        <- data.frame(matrix(nrow = n_sims, ncol = 5))

names(results_MAE)          <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
names(results_MSEP)         <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
names(results_precision)    <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
names(results_recall)       <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
names(results_time)         <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")


t_start <- proc.time()
set.seed(2019)
cat(paste("\n* Simulations - started at", format(Sys.time(), "%X"),
          "============================================================\n"))


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
  
  for (i in case_outliers) {
    cell_outliers <- sample(p, size = p * pct_cell_out)
    Xc[i, cell_outliers] <- contamination[cell_outliers] + rnorm(length(cell_outliers))
    outliers_mat_flag[i, cell_outliers] <- TRUE
  }
  outliers_list[[j_sim]] <- outliers_mat_flag
  
  
  # Collect data samples --------------------------------------------------------------------------
  data_clean        <- cbind.data.frame(y, X)
  data_contaminated <- cbind.data.frame(y, Xc)
  names(data_clean) <- names(data_contaminated) <- c("y", paste0("X", 1:p))
  
  
  # Save data of iteration ------------------------------------------------------------------------
  data_list_clean[[j_sim]] <- data_clean
  data_list_contaminated[[j_sim]] <- data_contaminated
}

  
for (j_sim in 1:n_sims) {

  # Fit cellwise robust M regression models -------------------------------------------------------
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
                                                parameters = c(qnorm(0.95), qnorm(0.975), qnorm(0.999))))
  
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
                                               parameters = c(4.685)))
  
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
                                               parameters = c(1.345)))
  
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
                                               parameters = c(1.063, 1.387, 1.5)))
  
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
                                                   parameters = c(0.982, 1.473, 1.5)))
  
  
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
  cat(paste0(" - ", round(100 * j_sim / n_sims, 2),
             "% - ", sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))
}


t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
cat(paste("\nTime elapsed:",
          sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))



# Study results -----------------------------------------------------------------------------------

# Mean Absolute Error - - - - - - - - - - - - - -
boxplot(results_MAE, ylab = "MAE", col = "lightblue",
        ylim = c(0, max(results_MAE)), las = 1, cex.axis = 1.5, cex.lab = 1.6)
points(colMeans(results_MAE), col = "red", pch = 18, cex = 1.7)
text(rep(0, 8), labels = round(colMeans(results_MAE), 4), cex = 2,
     font = ifelse(1:5 == which.min(colMeans(results_MAE)), 2, 1))


# Mean Squared Error of Prediction - - - - - - - -
boxplot(results_MSEP, ylab = "MSEP", col = "lightblue",
        ylim = c(0, max(results_MSEP)), las = 0.99, cex.axis = 1.5, cex.lab = 1.6)
points(colMeans(results_MSEP), col = "red", pch = 18, cex = 1.7)
text(rep(0, 8), labels = round(colMeans(results_MSEP), 4), cex = 2,
     font = ifelse(1:5 == which.min(colMeans(results_MSEP)), 2, 1))


# Precision - - - - - - - - - - - - - - - - - - - 
boxplot(results_precision, ylab = "Precision", col = "lightblue",
        ylim = c(0, max(results_precision)), las = 1, cex.axis = 1.5, cex.lab = 1.6)
points(colMeans(results_precision), col = "red", pch = 18, cex = 1.7)
text(rep(0, 8), labels = round(colMeans(results_precision), 4), cex = 2,
     font = ifelse(1:5 == which.max(colMeans(results_precision)), 2, 1))


# Recall - - - - - - - - - - - - - - - - - - - - -
boxplot(results_recall, ylab = "Recall", col = "lightblue",
        ylim = c(0, max(results_recall)), las = 1, cex.axis = 1.5, cex.lab = 1.6)
points(colMeans(results_recall), col = "red", pch = 18, cex = 1.7)
text(rep(0, 8), labels = round(colMeans(results_recall), 4), cex = 2,
     font = ifelse(1:5 == which.max(colMeans(results_recall)), 2, 1))


# Execution time - - - - - - - - - - - - - - - - -
cat("CRM-Hampel average execution time:", round(mean(results_time$'CRM-Hampel'), 1), "seconds\n\n")
cat("CRM-Tukey average execution time:", round(mean(results_time$'CRM-Tukey'), 1), "seconds\n\n")
cat("CRM-Huber average execution time:", round(mean(results_time$'CRM-Huber'), 1), "seconds\n\n")
cat("CRM-Gauss average execution time:", round(mean(results_time$'CRM-Gauss'), 1), "seconds\n\n")
cat("CRM-Quadratic average execution time:", round(mean(results_time$'CRM-Quadratic'), 1), "seconds\n\n")
