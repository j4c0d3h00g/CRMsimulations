# Cellwise Robust M-regression (CRM) - Simulation Study Controlling Contamination -----------------
# Here, only the amount of cellwise contamination is increased. This corresponds to situation 2.
rm(list = ls())

# Setup -------------------------------------------------------------------------------------------

# simulation setting - - - - - - - - - - - - - - -
n <- 400                                # number of cases
p <- 60                                 # number of predictor variables
pct_case_out <- 0.10                    # percentage of casewise outliers
pct_cell_out_seq <- seq(0, 0.5, 0.05)   # sequence of percentages of cellwise outliers for each casewise outlier
n_sims <- 15                            # number of simulations for each value of pct_case_out
betas_list <- list()                    # list to store betas of each iteration
y_list <- list()                        # list to store dependent variable of each iteration
X_list <- list()                        # list to store data matrix of each iteration
Xc_list <- list()                       # list to store contaminated data matrix of each iteration


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
results_MAE       <- list()
results_MSEP      <- list()
results_RMSEI     <- list()
results_precision <- list()
results_recall    <- list()
results_time      <- list()

t_start <- proc.time()
set.seed(2019)
cat(paste("\n* Simulations started at", format(Sys.time(), "%X"),
          "============================================================\n"))

n_pct_cell_out <- length(pct_cell_out_seq)
for (j_sim in 1:n_sims) {
  
  # Generate clean sample -------------------------------------------------------------------------
  mu <- rep(0, p)
  Sigma <- diag(p)
  Sigma[(row(Sigma) - col(Sigma)) == 1 | row(Sigma) - col(Sigma) == -1] <- 0.5
  X <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  X_list[[j_sim]] <- X
  
  slopes <- rnorm(n = p, mean = 0, sd = 1)
  slopes <- 10 * slopes / sqrt(sum(slopes^2))
  intercept <- 10
  noise <- rnorm(n = n, mean = 0, sd = 0.5)
  
  y <- intercept + X %*% slopes + noise
  y_list[[j_sim]] <- y
  betas <- c(intercept, slopes)
  betas_list[[j_sim]] <- betas
  
  
  # Setup for contamination -----------------------------------------------------------------------
  Xc <- X
  Xc_list[[j_sim]] <- Xc
  contamination <- colMeans(X) + 6 * apply(X, 2, sd)
}  


for (j_sim in 1:n_sims) {
  cat(paste0("\n* Simulation ", j_sim, "/", n_sims, " =======================================\n"))
  
  results_MAE_j_sim       <- data.frame(matrix(nrow = n_pct_cell_out, ncol = 6))
  results_MSEP_j_sim      <- data.frame(matrix(nrow = n_pct_cell_out, ncol = 6))
  results_precision_j_sim <- data.frame(matrix(nrow = n_pct_cell_out, ncol = 6))
  results_recall_j_sim    <- data.frame(matrix(nrow = n_pct_cell_out, ncol = 6))
  results_time_j_sim      <- data.frame(matrix(nrow = n_pct_cell_out, ncol = 6))
  
  names(results_MAE_j_sim)   <- c("CRM-Hampel-norm", "CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
  names(results_MSEP_j_sim)  <- c("CRM-Hampel-norm", "CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
  names(results_precision_j_sim)  <- c("CRM-Hampel-norm", "CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
  names(results_recall_j_sim)  <- c("CRM-Hampel-norm", "CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
  names(results_time_j_sim)  <- c("CRM-Hampel-norm", "CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
  
  cell_outliers <- c()
  uncontaminated_rows <- 1:n
  uncontaminated_cells <- 1:p
  uncontaminated_cells_list <- list()
  outliers_mat_flag <- matrix(FALSE, nrow = n, ncol = p)
  
  case_outliers <- sample(uncontaminated_rows, size = round(n * pct_case_out))
  for (i in case_outliers) {
    uncontaminated_cells_list[[i]] <- uncontaminated_cells
  }
  
  
  for (cell_pct_index in 1:n_pct_cell_out) {
    new_cell_outliers <- list()
    
    # Add contamination in design matrix ----------------------------------------------------------
    
    # For each casewise outlier the additional cellwise contamination is added --------------------
    for (i in case_outliers) {
      extra_pct_cell_out <- ifelse(cell_pct_index == 1,
                                   pct_cell_out_seq[1],
                                   pct_cell_out_seq[cell_pct_index] - pct_cell_out_seq[cell_pct_index - 1])
      
      new_cell_outliers[[i]] <- sample(uncontaminated_cells_list[[i]], size = round(p * extra_pct_cell_out))
      
      cell_outliers       <- c(cell_outliers, new_cell_outliers[[i]])
      uncontaminated_cells_list[[i]] <- setdiff(uncontaminated_cells_list[[i]], new_cell_outliers[[i]])
      
      Xc_list[[j_sim]][i, new_cell_outliers[[i]]] <- contamination[new_cell_outliers[[i]]] + rnorm(length(new_cell_outliers[[i]]))
      outliers_mat_flag[i, new_cell_outliers[[i]]] <- TRUE
    }
    
    
    cat(paste0("\n # cellwise outliers = ", length(cell_outliers),
               " (",100 * pct_cell_out_seq[cell_pct_index], "%)\n"))
    
    
    # Collect data samples ------------------------------------------------------------------------
    data_clean        <- cbind.data.frame(y_list[[j_sim]], X_list[[j_sim]])
    data_contaminated <- cbind.data.frame(y_list[[j_sim]], Xc_list[[j_sim]])
    names(data_clean) <- names(data_contaminated) <- c("y", paste0("X", 1:p))
    
    
    # Fit cellwise robust M regression models -----------------------------------------------------
    crm_hampel_norm <- suppressWarnings(crm_functional(formula   = y ~ .,
                                                       data      = data_contaminated,
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
    
    crm_hampel <- suppressWarnings(crm_functional(formula   = y ~ .,
                                                  data      = data_contaminated,
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
                                                  parameters = c(1.382, 2.764, 5.528)))
    
    crm_tukey <- suppressWarnings(crm_functional(formula   = y ~ .,
                                                 data      = data_contaminated,
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
                                                 data      = data_contaminated,
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
                                                 data      = data_contaminated,
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
                                                     data      = data_contaminated,
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
    
    
    # Evaluate performance ------------------------------------------------------------------------
    
    # Mean Absolute Error - - - - - - - - - - - -
    results_MAE_j_sim$'CRM-Hampel-norm'[cell_pct_index]     <- mean(abs(crm_hampel_norm$coefficients - betas_list[[j_sim]]))
    results_MAE_j_sim$'CRM-Hampel'[cell_pct_index]          <- mean(abs(crm_hampel$coefficients - betas_list[[j_sim]]))
    results_MAE_j_sim$'CRM-Tukey'[cell_pct_index]           <- mean(abs(crm_tukey$coefficients - betas_list[[j_sim]]))
    results_MAE_j_sim$'CRM-Huber'[cell_pct_index]           <- mean(abs(crm_huber$coefficients - betas_list[[j_sim]]))
    results_MAE_j_sim$'CRM-Gauss'[cell_pct_index]           <- mean(abs(crm_gauss$coefficients - betas_list[[j_sim]]))
    results_MAE_j_sim$'CRM-Quadratic'[cell_pct_index]       <- mean(abs(crm_quadratic$coefficients - betas_list[[j_sim]]))
    
    
    # Mean Squared Error of Prediction - - - - - -
    if (length(case_outliers) == 0) {
      results_MSEP_j_sim$'CRM-Hampel-norm'[cell_pct_index]      <- mean((crm_hampel_norm$fitted.values - data_clean$y)^2)
      results_MSEP_j_sim$'CRM-Hampel'[cell_pct_index]           <- mean((crm_hampel$fitted.values - data_clean$y)^2)
      results_MSEP_j_sim$'CRM-Tukey'[cell_pct_index]            <- mean((crm_tukey$fitted.values - data_clean$y)^2)
      results_MSEP_j_sim$'CRM-Huber'[cell_pct_index]            <- mean((crm_huber$fitted.values - data_clean$y)^2)
      results_MSEP_j_sim$'CRM-Gauss'[cell_pct_index]            <- mean((crm_gauss$fitted.values - data_clean$y)^2)
      results_MSEP_j_sim$'CRM-Quadratic'[cell_pct_index]        <- mean((crm_quadratic$fitted.values - data_clean$y)^2)
    } else {
      results_MSEP_j_sim$'CRM-Hampel-norm'[cell_pct_index]      <- mean((crm_hampel_norm$fitted.values - data_clean$y)[-case_outliers]^2)
      results_MSEP_j_sim$'CRM-Hampel'[cell_pct_index]           <- mean((crm_hampel$fitted.values - data_clean$y)[-case_outliers]^2)
      results_MSEP_j_sim$'CRM-Tukey'[cell_pct_index]            <- mean((crm_tukey$fitted.values - data_clean$y)[-case_outliers]^2)
      results_MSEP_j_sim$'CRM-Huber'[cell_pct_index]            <- mean((crm_huber$fitted.values - data_clean$y)[-case_outliers]^2)
      results_MSEP_j_sim$'CRM-Gauss'[cell_pct_index]            <- mean((crm_gauss$fitted.values - data_clean$y)[-case_outliers]^2)
      results_MSEP_j_sim$'CRM-Quadratic'[cell_pct_index]        <- mean((crm_quadratic$fitted.values - data_clean$y)[-case_outliers]^2)
    }
    
    
    # Precision & recall of detected cells - - - -
    if (length(cell_outliers) == 0) {
      results_precision_j_sim$'CRM-Hampel-norm'[cell_pct_index]   <- NA
      results_precision_j_sim$'CRM-Hampel'[cell_pct_index]        <- NA
      results_precision_j_sim$'CRM-Tukey'[cell_pct_index]         <- NA
      results_precision_j_sim$'CRM-Huber'[cell_pct_index]         <- NA
      results_precision_j_sim$'CRM-Gauss'[cell_pct_index]         <- NA
      results_precision_j_sim$'CRM-Quadratic'[cell_pct_index]     <- NA
      
      results_recall_j_sim$'CRM-Hampel-norm'[cell_pct_index]      <- NA
      results_recall_j_sim$'CRM-Hampel'[cell_pct_index]           <- NA
      results_recall_j_sim$'CRM-Tukey'[cell_pct_index]            <- NA
      results_recall_j_sim$'CRM-Huber'[cell_pct_index]            <- NA
      results_recall_j_sim$'CRM-Gauss'[cell_pct_index]            <- NA
      results_recall_j_sim$'CRM-Quadratic'[cell_pct_index]        <- NA
    } else {
      CM_hampel_norm    <- table(Prediction = ifelse(c(crm_hampel_norm$cellwiseoutliers) != 0, TRUE, FALSE),
                                  Reference  = c(outliers_mat_flag))
      CM_hampel         <- table(Prediction = ifelse(c(crm_hampel$cellwiseoutliers) != 0, TRUE, FALSE),
                                  Reference  = c(outliers_mat_flag))
      CM_tukey          <- table(Prediction = ifelse(c(crm_tukey$cellwiseoutliers) != 0, TRUE, FALSE),
                                  Reference  = c(outliers_mat_flag))
      CM_huber          <- table(Prediction = ifelse(c(crm_huber$cellwiseoutliers) != 0, TRUE, FALSE),
                                  Reference  = c(outliers_mat_flag))
      CM_gauss          <- table(Prediction = ifelse(c(crm_gauss$cellwiseoutliers) != 0, TRUE, FALSE),
                                  Reference  = c(outliers_mat_flag))
      CM_quadratic      <- table(Prediction = ifelse(c(crm_quadratic$cellwiseoutliers) != 0, TRUE, FALSE),
                                  Reference  = c(outliers_mat_flag))
      
      results_precision_j_sim$'CRM-Hampel-norm'[cell_pct_index]     <- CM_hampel_norm[2, 2] / (CM_hampel_norm[2, 1] + CM_hampel_norm[2, 2])
      results_precision_j_sim$'CRM-Hampel'[cell_pct_index]          <- CM_hampel[2, 2] / (CM_hampel[2, 1] + CM_hampel[2, 2])
      results_precision_j_sim$'CRM-Tukey'[cell_pct_index]           <- CM_tukey[2, 2] / (CM_tukey[2, 1] + CM_tukey[2, 2])
      results_precision_j_sim$'CRM-Huber'[cell_pct_index]           <- CM_huber[2, 2] / (CM_huber[2, 1] + CM_huber[2, 2])
      results_precision_j_sim$'CRM-Gauss'[cell_pct_index]           <- CM_gauss[2, 2] / (CM_gauss[2, 1] + CM_gauss[2, 2])
      results_precision_j_sim$'CRM-Quadratic'[cell_pct_index]       <- CM_quadratic[2, 2] / (CM_quadratic[2, 1] + CM_quadratic[2, 2])
      
      results_recall_j_sim$'CRM-Hampel-norm'[cell_pct_index]        <- CM_hampel_norm[2, 2] / (CM_hampel_norm[1, 2] + CM_hampel_norm[2, 2])
      results_recall_j_sim$'CRM-Hampel'[cell_pct_index]             <- CM_hampel[2, 2] / (CM_hampel[1, 2] + CM_hampel[2, 2])
      results_recall_j_sim$'CRM-Tukey'[cell_pct_index]              <- CM_tukey[2, 2] / (CM_tukey[1, 2] + CM_tukey[2, 2])
      results_recall_j_sim$'CRM-Huber'[cell_pct_index]              <- CM_huber[2, 2] / (CM_huber[1, 2] + CM_huber[2, 2])
      results_recall_j_sim$'CRM-Gauss'[cell_pct_index]              <- CM_gauss[2, 2] / (CM_gauss[1, 2] + CM_gauss[2, 2])
      results_recall_j_sim$'CRM-Quadratic'[cell_pct_index]          <- CM_quadratic[2, 2] / (CM_quadratic[1, 2] + CM_quadratic[2, 2])
    }
    
    
    # Execution time - - - - - - - - - - - - - - -
    results_time_j_sim$'CRM-Hampel-norm'[cell_pct_index]      <- crm_hampel_norm$time
    results_time_j_sim$'CRM-Hampel'[cell_pct_index]           <- crm_hampel$time
    results_time_j_sim$'CRM-Tukey'[cell_pct_index]            <- crm_tukey$time
    results_time_j_sim$'CRM-Huber'[cell_pct_index]            <- crm_huber$time
    results_time_j_sim$'CRM-Gauss'[cell_pct_index]            <- crm_gauss$time
    results_time_j_sim$'CRM-Quadratic'[cell_pct_index]        <- crm_quadratic$time
    
    
    # Elapsed time - - - - - - - - - - - - - - - -
    t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
    cat(paste0(" - ", round(100 * (cell_pct_index + (j_sim - 1) * n_pct_cell_out) /
                              (n_sims * n_pct_cell_out), 2),
               "% - ", sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))
    
  } # end of forloop "cell_pct_index in 1:n_pct_cell_out"
  
  
  results_MAE[[j_sim]]       <- results_MAE_j_sim
  results_MSEP[[j_sim]]      <- results_MSEP_j_sim
  results_precision[[j_sim]] <- results_precision_j_sim
  results_recall[[j_sim]]    <- results_recall_j_sim
  results_time[[j_sim]]      <- results_time_j_sim
  
} # end of forloop "j_sim in 1:n_sims"


t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
cat(paste("\nTime elapsed:",
          sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))


# Study results -----------------------------------------------------------------------------------
results_MAE       <- cbind.data.frame(do.call(rbind, results_MAE),       pct_cell_out = rep(pct_cell_out_seq, n_sims))
results_MSEP      <- cbind.data.frame(do.call(rbind, results_MSEP),      pct_cell_out = rep(pct_cell_out_seq, n_sims))
results_precision <- cbind.data.frame(do.call(rbind, results_precision), pct_cell_out = rep(pct_cell_out_seq, n_sims))
results_recall    <- cbind.data.frame(do.call(rbind, results_recall),    pct_cell_out = rep(pct_cell_out_seq, n_sims))
results_time      <- cbind.data.frame(do.call(rbind, results_time),      pct_cell_out = rep(pct_cell_out_seq, n_sims))

results_MAE_mean       <- aggregate(. ~ pct_cell_out, data = results_MAE,       FUN = mean)
results_MSEP_mean      <- aggregate(. ~ pct_cell_out, data = results_MSEP,      FUN = mean)
results_precision_mean <- aggregate(. ~ pct_cell_out, data = results_precision, FUN = mean)
results_recall_mean    <- aggregate(. ~ pct_cell_out, data = results_recall,    FUN = mean)
results_time_mean      <- aggregate(. ~ pct_cell_out, data = results_time,      FUN = mean)


# The following plots are not used in our research - - -

# Mean Absolute Error - - - - - - - - - - - - - -
plot(100 * results_MAE_mean$pct_cell_out, results_MAE_mean$'CRM-Hampel-norm', type = "l", xlab = "% cellwise outliers",
     ylab = "Average MAE", ylim = range(results_MAE_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2)
lines(100 * results_MAE_mean$pct_cell_out, results_MAE_mean$'CRM-Hampel',      lwd = 2, lty = 2)
lines(100 * results_MAE_mean$pct_cell_out, results_MAE_mean$'CRM-Tukey',      lwd = 2, lty = 3)
lines(100 * results_MAE_mean$pct_cell_out, results_MAE_mean$'CRM-Huber',  lwd = 2, lty = 4)
lines(100 * results_MAE_mean$pct_cell_out, results_MAE_mean$'CRM-Gauss',     lwd = 2, lty = 5)
lines(100 * results_MAE_mean$pct_cell_out, results_MAE_mean$'CRM-Quadratic', lwd = 2, lty = 6)
legend("bottomright", legend = names(results_MAE_mean)[-1], lty = 1:6, lwd = 2, cex = 1, y.intersp = 0.5)


# Mean Squared Error of Prediction - - - - - - - -
plot(100 * results_MSEP_mean$pct_cell_out, results_MSEP_mean$'CRM-Hampel-norm', type = "l", xlab = "% cellwise outliers",
     ylab = "Average MSEP", ylim = range(results_MSEP_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2)
lines(100 * results_MSEP_mean$pct_cell_out, results_MSEP_mean$'CRM-Hampel',      lwd = 2, lty = 2)
lines(100 * results_MSEP_mean$pct_cell_out, results_MSEP_mean$'CRM-Tukey',      lwd = 2, lty = 3)
lines(100 * results_MSEP_mean$pct_cell_out, results_MSEP_mean$'CRM-Huber',  lwd = 2, lty = 4)
lines(100 * results_MSEP_mean$pct_cell_out, results_MSEP_mean$'CRM-Gauss',     lwd = 2, lty = 5)
lines(100 * results_MSEP_mean$pct_cell_out, results_MSEP_mean$'CRM-Quadratic', lwd = 2, lty = 6)
legend("topleft", legend = names(results_MSEP_mean)[-1], lty = 1:6, lwd = 2)


# Precision detected cells - - - - - - - - - - - -
par(mfrow = c(1, 2), pty = "s")
plot(100 * results_precision_mean$pct_cell_out, results_precision_mean$'CRM-Hampel-norm', type = "l", xlab = "% cellwise outliers",
     ylab = "Average precision", las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2, ylim = c(0, 1))
lines(100 * results_precision_mean$pct_cell_out, results_precision_mean$'CRM-Hampel',      lwd = 2, lty = 2)
lines(100 * results_precision_mean$pct_cell_out, results_precision_mean$'CRM-Tukey',      lwd = 2, lty = 3)
lines(100 * results_precision_mean$pct_cell_out, results_precision_mean$'CRM-Huber',  lwd = 2, lty = 4)
lines(100 * results_precision_mean$pct_cell_out, results_precision_mean$'CRM-Gauss',     lwd = 2, lty = 5)
lines(100 * results_precision_mean$pct_cell_out, results_precision_mean$'CRM-Quadratic', lwd = 2, lty = 6)
legend("topleft", legend = names(results_precision_mean)[-1], lty = 1:6, lwd = 2)


# Recall of detected cells - - - - - - - - - - - -
plot(100 * results_recall_mean$pct_cell_out, results_recall_mean$'CRM-Hampel-norm', type = "l", xlab = "% cellwise outliers",
     ylab = "Average recall", las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2, ylim = c(0, 1))
lines(100 * results_recall_mean$pct_cell_out, results_recall_mean$'CRM-Hampel',      lwd = 2, lty = 2)
lines(100 * results_recall_mean$pct_cell_out, results_recall_mean$'CRM-Tukey',      lwd = 2, lty = 3)
lines(100 * results_recall_mean$pct_cell_out, results_recall_mean$'CRM-Huber',  lwd = 2, lty = 4)
lines(100 * results_recall_mean$pct_cell_out, results_recall_mean$'CRM-Gauss',     lwd = 2, lty = 5)
lines(100 * results_recall_mean$pct_cell_out, results_recall_mean$'CRM-Quadratic', lwd = 2, lty = 6)
legend("topleft", legend = names(results_recall_mean)[-1], lty = 1:6, lwd = 2)


# Execution time - - - - - - - - - - - - - - - - -
cat("CRM-Hampel-norm average execution time:", round(mean(results_time$'CRM-Hampel-norm'), 1), "seconds")
cat("CRM-Hampel average execution time:", round(mean(results_time$'CRM-Hampel'), 1), "seconds")
cat("CRM-Tukey average execution time:", round(mean(results_time$'CRM-Tukey'), 1), "seconds")
cat("CRM-Huber average execution time:", round(mean(results_time$'CRM-Huber'), 1), "seconds")
cat("CRM-Gauss average execution time:", round(mean(results_time$'CRM-Gauss'), 1), "seconds")
cat("CRM-Quadratic average execution time:", round(mean(results_time$'CRM-Quadratic'), 1), "seconds")


# The following plot is used in our research - - - - - - -

# Mean Absolute Error for different weight functions - - - 
library(tidyverse)
library(hrbrthemes)
library(kableExtra)
options(knitr.table.format = "html")
library(streamgraph)
library(viridis)
library(DT)
library(plotly)

MAE_matrix <- data.frame(pct_cell_out = rep(100 * results_MAE_mean$pct_cell_out, times = 6),
                         crm_function = rep(c("CRM-Hampel-norm", "CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic"), each = 11),
                         average_MAE = c(results_MAE_mean$'CRM-Hampel-norm', results_MAE_mean$'CRM-Hampel', results_MAE_mean$'CRM-Tukey', results_MAE_mean$'CRM-Huber',
                                         results_MAE_mean$'CRM-Gauss', results_MAE_mean$'CRM-Quadratic'))

tmp <- MAE_matrix %>%
  mutate(crm_function_2 = crm_function)

tmp %>%
  ggplot(aes(x = pct_cell_out, y = average_MAE)) +
  geom_line(data = tmp %>% dplyr::select(-crm_function), aes(group = crm_function_2), color = "grey", size = 0.5, alpha = 0.9) +
  geom_line(aes(color = crm_function), color = "#69b3a2", size = 1.2) +
  scale_color_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=14),
    panel.grid = element_blank()
  ) +
  ylab("Average MAE") + xlab("% cellwise contamination") +
  facet_wrap(~factor(crm_function, levels = c("CRM-Hampel-norm", "CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")))
