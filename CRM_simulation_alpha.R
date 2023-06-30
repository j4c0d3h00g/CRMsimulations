# Cellwise Robust M-regression (CRM) - Simulation Study Increasing Robustness ---------------------
rm(list = ls())

# Setup -------------------------------------------------------------------------------------------

# simulation setting - - - - - - - - - - - - - - -
n <- 400                                    # number of cases
p <- 60                                     # number of predictor variables
pct_case_out_seq <- seq(0, 0.5, 0.05)       # sequence of percentages of casewise outliers
n_pct_case_out <- length(pct_case_out_seq)  # number of considered casewise outliers levels
pct_cell_out <- 0.10                        # percentage of cellwise outliers for each casewise outlier
n_sims <- 10                                # number of simulations for each value of pct_case_out
alphas <- c(0.6, 0.7, 0.8, 0.9, 1.0)        # values of alpha to adjust the parameters of the weight functions
number_of_alphas <- length(alphas)          # number of considered alpha values
betas_list <- list()                        # list to store betas of each iteration
y_list <- list()                            # list to store dependent variable of each iteration
X_list <- list()                            # list to store data matrix of each iteration
Xc_list <- list()                           # list to store contaminated data matrix of each iteration


# CRM input parameters - - - - - - - - - - - - - -
maxiter             <- 100
tolerance           <- 0.01
outlyingness.factor <- 1.5
spadieta            <- seq(0.9, 0.1, -0.1)


# Load packages -----------------------------------------------------------------------------------
library(CRMwf)
library(MASS)
library(cellWise)
library(lubridate)
library(robustbase)


# Start simulation procedure ----------------------------------------------------------------------
results_MAE_hampel            <- list()
results_MAE_tukey             <- list()
results_MAE_huber             <- list()
results_MAE_gauss             <- list()
results_MAE_quadratic         <- list()

names(results_MAE_hampel)     <- c("alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")
names(results_MAE_tukey)      <- c("alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")
names(results_MAE_huber)      <- c("alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")
names(results_MAE_gauss)      <- c("alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")
names(results_MAE_quadratic)  <- c("alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")

parameters_hampel     <- c(1.382, 2.764, 5.528)
parameters_tukey      <- c(4.685)
parameters_huber      <- c(1.345)
parameters_gauss      <- c(1.063, 1.387, 1.5)
parameters_quadratic  <- c(0.982, 1.473, 1.5)

t_start <- proc.time()
set.seed(2019)
cat(paste("\n* Simulations started at", format(Sys.time(), "%X"),
          "============================================================\n"))

n_pct_case_out <- length(pct_case_out_seq)
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
  
  results_MAE_j_sim          <- data.frame(matrix(nrow = n_pct_case_out, ncol = 5))
  names(results_MAE_j_sim)   <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")
  
  case_outliers <- c()
  uncontaminated_rows <- 1:n
  outliers_mat_flag <- matrix(FALSE, nrow = n, ncol = p)
  
  
  for (case_pct_index in 1:n_pct_case_out) {
    
    # Add contamination in design matrix ----------------------------------------------------------
    extra_pct_case_out  <- ifelse(case_pct_index == 1,
                                    pct_case_out_seq[1],
                                    pct_case_out_seq[case_pct_index] - pct_case_out_seq[case_pct_index - 1])
    
    new_case_outliers   <- sample(uncontaminated_rows, size = round(n * extra_pct_case_out))
    case_outliers       <- c(case_outliers, new_case_outliers)
    uncontaminated_rows <- setdiff(uncontaminated_rows, case_outliers)
    
    cell_outliers <- NULL
    for (i in new_case_outliers) {
      cell_outliers <- sample(p, size = p * pct_cell_out)
      Xc_list[[j_sim]][i, cell_outliers] <- contamination[cell_outliers] + rnorm(length(cell_outliers))
      outliers_mat_flag[i, cell_outliers] <- TRUE
    }
    
    cat(paste0("\n # casewise outliers = ", length(case_outliers),
               " (",100 * pct_case_out_seq[case_pct_index], "%)\n"))
    
    
    # Collect data samples ------------------------------------------------------------------------
    data_clean        <- cbind.data.frame(y_list[[j_sim]], X_list[[j_sim]])
    data_contaminated <- cbind.data.frame(y_list[[j_sim]], Xc_list[[j_sim]])
    names(data_clean) <- names(data_contaminated) <- c("y", paste0("X", 1:p))
    
    alpha_index <- 0
    for (alpha in alphas) {
      alpha_index <- alpha_index + 1
      cat(paste0("\n* Evaluated alpha ", alpha_index, "/", number_of_alphas, " =======================================\n"))
      
      # Lower the parameters with a factor alpha ---------------------------------------------------
      parameters_hampel_use           <- alpha * parameters_hampel
      parameters_tukey_use            <- alpha * parameters_tukey
      parameters_huber_use            <- alpha * parameters_huber
      parameters_gauss_use            <- parameters_gauss
      parameters_gauss_use[1:2]       <- alpha * parameters_gauss_use[1:2]
      parameters_quadratic_use        <- parameters_quadratic
      parameters_quadratic_use[1:2]   <- alpha * parameters_gauss_use[1:2]
    
      # Fit cellwise robust M regression models -----------------------------------------------------
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
                                                    parameters = parameters_hampel_use))
      
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
                                                   parameters = parameters_tukey_use))
      
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
                                                   parameters = parameters_huber_use))
      
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
                                                   parameters = parameters_gauss_use))
      
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
                                                       parameters = parameters_quadratic_use))
      
      
      # Evaluate performance ------------------------------------------------------------------------
      
      # Mean Absolute Error - - - - - - - - - - - -
      results_MAE_hampel_temp[case_pct_index, alpha_index]      <- mean(abs(crm_hampel$coefficients - betas_list[[j_sim]]))
      results_MAE_tukey_temp[case_pct_index, alpha_index]       <- mean(abs(crm_tukey$coefficients - betas_list[[j_sim]]))
      results_MAE_huber_temp[case_pct_index, alpha_index]       <- mean(abs(crm_huber$coefficients - betas_list[[j_sim]]))
      results_MAE_gauss_temp[case_pct_index, alpha_index]       <- mean(abs(crm_gauss$coefficients - betas_list[[j_sim]]))
      results_MAE_quadratic_temp[case_pct_index, alpha_index]   <- mean(abs(crm_quadratic$coefficients - betas_list[[j_sim]]))
    
      # Elapsed time - - - - - - - - - - - - - - - -
      t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
      cat(paste0(" - ", round(100 * (alpha_index + number_of_alphas * ((j_sim - 1) * n_pct_case_out + (case_pct_index - 1))) /
                                (n_sims * n_pct_case_out * number_of_alphas), 2),
                 "% - ", sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))
    
      } # end of forloop "alpha in alphas"
      
    } # end of forloop "case_pct_index in 1:n_pct_case_out"
  
  results_MAE_hampel[[j_sim]]       <- results_MAE_hampel_temp
  results_MAE_tukey[[j_sim]]        <- results_MAE_tukey_temp
  results_MAE_huber[[j_sim]]        <- results_MAE_huber_temp
  results_MAE_gauss[[j_sim]]        <- results_MAE_gauss_temp
  results_MAE_quadratic[[j_sim]]    <- results_MAE_quadratic_temp
  
} # end of forloop "j_sim in 1:n_sims"

results_MAE_hampel       <- cbind.data.frame(do.call(rbind, results_MAE_hampel),      pct_case_out = rep(pct_case_out_seq, n_sims))
results_MAE_tukey        <- cbind.data.frame(do.call(rbind, results_MAE_tukey),       pct_case_out = rep(pct_case_out_seq, n_sims))
results_MAE_huber        <- cbind.data.frame(do.call(rbind, results_MAE_huber),       pct_case_out = rep(pct_case_out_seq, n_sims))
results_MAE_gauss        <- cbind.data.frame(do.call(rbind, results_MAE_gauss),       pct_case_out = rep(pct_case_out_seq, n_sims))
results_MAE_quadratic    <- cbind.data.frame(do.call(rbind, results_MAE_quadratic),   pct_case_out = rep(pct_case_out_seq, n_sims))

results_MAE_hampel_mean           <- aggregate(. ~ pct_case_out, data = results_MAE_hampel,          FUN = mean)
results_MAE_tukey_mean            <- aggregate(. ~ pct_case_out, data = results_MAE_tukey,           FUN = mean)
results_MAE_huber_mean            <- aggregate(. ~ pct_case_out, data = results_MAE_huber,           FUN = mean)
results_MAE_gauss_mean            <- aggregate(. ~ pct_case_out, data = results_MAE_gauss,           FUN = mean)
results_MAE_quadratic_mean        <- aggregate(. ~ pct_case_out, data = results_MAE_quadratic,       FUN = mean)

names(results_MAE_hampel_mean)      <- c("pct_case_out", "alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")
names(results_MAE_tukey_mean)       <- c("pct_case_out", "alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")
names(results_MAE_huber_mean)       <- c("pct_case_out", "alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")
names(results_MAE_gauss_mean)       <- c("pct_case_out", "alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")
names(results_MAE_quadratic_mean)   <- c("pct_case_out", "alpha-0.6", "alpha-0.7", "alpha-0.8", "alpha-0.9", "alpha-1.0")


# The following plot is used in our research - - - - - - -

# Mean Absolute Error for different weight functions and alphas - - - 
library(tidyverse)
library(hrbrthemes)
library(kableExtra)
options(knitr.table.format = "html")
library(streamgraph)
library(viridis)
library(DT)
library(plotly)

# Transform data matrices of each weight function - - - - -

# The columns correspond to the percentage of casewise outliers,
# the weight function, the value of alpha, and the average MAE. - - - 
MAE_matrix_hampel <- data.frame(pct_case_out = rep(100 * results_MAE_hampel_mean$pct_case_out, times = 5),
                                crm_function = rep("CRM-Hampel", each = 55),
                                alpha_value = rep(c("0.6", "0.7", "0.8", "0.9", "1.0"), each = 11),
                                average_MAE = c(results_MAE_hampel_mean$'alpha-0.6', results_MAE_hampel_mean$'alpha-0.7', results_MAE_hampel_mean$'alpha-0.8',
                                                results_MAE_hampel_mean$'alpha-0.9', results_MAE_hampel_mean$'alpha-1.0'))
MAE_matrix_tukey <- data.frame(pct_case_out = rep(100 * results_MAE_tukey_mean$pct_case_out, times = 5),
                               crm_function = rep("CRM-Tukey", each = 55),
                               alpha_value = rep(c("0.6", "0.7", "0.8", "0.9", "1.0"), each = 11),
                               average_MAE = c(results_MAE_tukey_mean$'alpha-0.6', results_MAE_tukey_mean$'alpha-0.7', results_MAE_tukey_mean$'alpha-0.8',
                                               results_MAE_tukey_mean$'alpha-0.9', results_MAE_tukey_mean$'alpha-1.0'))
MAE_matrix_huber <- data.frame(pct_case_out = rep(100 * results_MAE_huber_mean$pct_case_out, times = 5),
                               crm_function = rep("CRM-Huber", each = 55),
                               alpha_value = rep(c("0.6", "0.7", "0.8", "0.9", "1.0"), each = 11),
                               average_MAE = c(results_MAE_huber_mean$'alpha-0.6', results_MAE_huber_mean$'alpha-0.7', results_MAE_huber_mean$'alpha-0.8',
                                               results_MAE_huber_mean$'alpha-0.9', results_MAE_huber_mean$'alpha-1.0'))
MAE_matrix_gauss <- data.frame(pct_case_out = rep(100 * results_MAE_gauss_mean$pct_case_out, times = 5),
                               crm_function = rep("CRM-Gauss", each = 55),
                               alpha_value = rep(c("0.6", "0.7", "0.8", "0.9", "1.0"), each = 11),
                               average_MAE = c(results_MAE_gauss_mean$'alpha-0.6', results_MAE_gauss_mean$'alpha-0.7', results_MAE_gauss_mean$'alpha-0.8',
                                               results_MAE_gauss_mean$'alpha-0.9', results_MAE_gauss_mean$'alpha-1.0'))
MAE_matrix_quadratic <- data.frame(pct_case_out = rep(100 * results_MAE_quadratic_mean$pct_case_out, times = 5),
                                   crm_function = rep("CRM-Quadratic", each = 55),
                                   alpha_value = rep(c("0.6", "0.7", "0.8", "0.9", "1.0"), each = 11),
                                   average_MAE = c(results_MAE_quadratic_mean$'alpha-0.6', results_MAE_quadratic_mean$'alpha-0.7', results_MAE_quadratic_mean$'alpha-0.8',
                                                   results_MAE_quadratic_mean$'alpha-0.9', results_MAE_quadratic_mean$'alpha-1.0'))

# Combine the matrices corresponding to the different weight functions in one large matrix - - - 
MAE_matrix <- rbind(MAE_matrix_hampel, MAE_matrix_tukey, MAE_matrix_huber, MAE_matrix_gauss, MAE_matrix_quadratic)

# Plot the average MAE for different alpha values for the different weight functions used in CRM - - - 
MAE_matrix %>%
  ggplot(aes(x = pct_case_out, y = average_MAE, group = alpha_value)) +
  geom_line(aes(color = alpha_value), size = 1) +
  scale_color_viridis(discrete = TRUE) +
  theme(
    legend.position = "none",
    plot.title = element_text(size=14)
  ) +
  theme_ipsum() +
  ylab("Average MAE") + xlab("% casewise outliers") + labs(colour = "Alpha") +
  facet_wrap(~factor(crm_function, levels = c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")))
