# Nutrient contents -------------------------------------------------------------------------------

# Load & preprocess data --------------------------------------------------------------------------
rm(list = ls())
data(nutrients, package = "robCompositions")
nutrients <- as.data.frame(nutrients)

# change rownames to the English name of the food product combined with the row numbers:
rownames(nutrients) <- paste(nutrients$name_E, 1:nrow(nutrients))

# select numerical variables:
nutrients <- nutrients[1:200, c(23, 14, 15, 17, 18, 20)]

# check for missing values:
print(sum(rowSums(is.na(nutrients)) > 0))
print(apply(nutrients, 2, function (x) sum(is.na(x))))
# 7 out of 200 food products have at least 1 missing value

nutrients <- na.omit(nutrients) # remove food products with at least 1 NA value
print(dim(nutrients))

# log-transform all variables:
transNutrients <- log(nutrients + 0.01)
colnames(transNutrients) <- paste0("log.", names(nutrients))
# transNutrients is the data set that we analyse with CRM.


# Cellwise Robust M-regression --------------------------------------------------------------------
library(CRMwf)

# CRM input parameters - - - - - - - - - - - - - -
maxiter             <- 100
tolerance           <- 0.01
outlyingness.factor <- 1
spadieta            <- seq(0.9, 0.1, -0.1)


crm_hampel <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                              data      = transNutrients,
                                              maxiter   = maxiter,
                                              tolerance = tolerance,
                                              outlyingness.factor = outlyingness.factor,
                                              spadieta  = spadieta,
                                              center    = "median",
                                              scale     = "sd",
                                              regtype   = "MM",
                                              seed      = 2019,
                                              verbose   = FALSE,
                                              weightfunction = "hampel",
                                              weightthreshold = 1,
                                              parameters = c(1.382, 2.764, 5.528)))

crm_tukey <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                             data      = transNutrients,
                                             maxiter   = maxiter,
                                             tolerance = tolerance,
                                             outlyingness.factor = outlyingness.factor,
                                             spadieta  = spadieta,
                                             center    = "median",
                                             scale     = "sd",
                                             regtype   = "MM",
                                             seed      = 2019,
                                             verbose   = FALSE,
                                             weightfunction = "tukey",
                                             weightthreshold = 0.7,
                                             parameters = c(4.685)))

crm_huber <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                             data      = transNutrients,
                                             maxiter   = maxiter,
                                             tolerance = tolerance,
                                             outlyingness.factor = outlyingness.factor,
                                             spadieta  = spadieta,
                                             center    = "median",
                                             scale     = "sd",
                                             regtype   = "MM",
                                             seed      = 2019,
                                             verbose   = FALSE,
                                             weightfunction = "huber",
                                             weightthreshold = 1,
                                             parameters = c(1.345)))

crm_gauss <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                             data      = transNutrients,
                                             maxiter   = maxiter,
                                             tolerance = tolerance,
                                             outlyingness.factor = outlyingness.factor,
                                             spadieta  = spadieta,
                                             center    = "median",
                                             scale     = "sd",
                                             regtype   = "MM",
                                             seed      = 2019,
                                             verbose   = FALSE,
                                             weightfunction = "gauss",
                                             weightthreshold = 1,
                                             parameters = c(1.063, 1.387, 1.5)))

crm_quadratic <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                                 data      = transNutrients,
                                                 maxiter   = maxiter,
                                                 tolerance = tolerance,
                                                 outlyingness.factor = outlyingness.factor,
                                                 spadieta  = spadieta,
                                                 center    = "median",
                                                 scale     = "sd",
                                                 regtype   = "MM",
                                                 seed      = 2019,
                                                 verbose   = FALSE,
                                                 weightfunction = "quadratic",
                                                 weightthreshold = 1,
                                                 parameters = c(0.982, 1.473, 1.5)))


# Analyse results ---------------------------------------------------------------------------------

# execution time - - - - - - - - - - - - - - - - -
print(crm_hampel$time)
print(crm_tukey$time)
print(crm_huber$time)
print(crm_gauss$time)
print(crm_quadratic$time)


# coefficients - - - - - - - - - - - - - - - - - -
print(round(crm_hampel$coefficients, 5))
print(round(crm_tukey$coefficients, 5))
print(round(crm_huber$coefficients, 5))
print(round(crm_gauss$coefficients, 5))
print(round(crm_quadratic$coefficients, 5))


# casewise outliers - - - - - - - - - - - - - - -
casesOfInterest_hampel <- which(crm_hampel$casewiseoutliers)
print(length(casesOfInterest_hampel))
print(rownames(nutrients)[casesOfInterest_hampel])
# 36 out of 193 food products are considered as casewise outliers by the CRM algorithm using the Hampel weight function

casesOfInterest_tukey <- which(crm_tukey$casewiseoutliers)
print(length(casesOfInterest_tukey))
print(rownames(nutrients)[casesOfInterest_tukey])
# 14 out of 193 food products are considered as casewise outliers by the CRM algorithm using the Tukeys's bisquare weight function

casesOfInterest_huber <- which(crm_huber$casewiseoutliers)
print(length(casesOfInterest_huber))
print(rownames(nutrients)[casesOfInterest_huber])
# 29 out of 193 food products are considered as casewise outliers by the CRM algorithm using the Huber weight function

casesOfInterest_gauss <- which(crm_gauss$casewiseoutliers)
print(length(casesOfInterest_gauss))
print(rownames(nutrients)[casesOfInterest_gauss])
# 0 out of 193 food products are considered as casewise outliers by the CRM algorithm using the Generalized Gauss weight function

casesOfInterest_quadratic <- which(crm_quadratic$casewiseoutliers)
print(length(casesOfInterest_quadratic))
print(rownames(nutrients)[casesOfInterest_quadratic])
# 103 out of 193 food products are considered as casewise outliers by the CRM algorithm using the linear quadratic quadratic weight function


# heatmap of cellwise outliers for CRM using the Hampel weight function - - - - -
library(crmReg)

# original nutrients data:
nutrients_casesOfInterest <- round(nutrients[casesOfInterest_hampel, -1], 1)
rownames(nutrients_casesOfInterest) <- substr(rownames(nutrients_casesOfInterest), 1,
                                              nchar(rownames(nutrients_casesOfInterest)) - c(2, rep(3, length(casesOfInterest_hampel)-1)))

binary_cellwise_outliers <- crm_hampel$cellwiseoutliers[casesOfInterest_hampel, ]
binary_cellwise_outliers[binary_cellwise_outliers < 0] <- -1
binary_cellwise_outliers[binary_cellwise_outliers > 0] <- 1

cellwiseheatmap(cellwiseoutliers = binary_cellwise_outliers,
                data = nutrients_casesOfInterest,
                margins = c(10, 22), notecex = 2)

cellwiseheatmap(cellwiseoutliers = crm_hampel$cellwiseoutliers[casesOfInterest_hampel, ],
                data = nutrients_casesOfInterest,
                col.scale.factor = 1/5,
                margins = c(10, 22), notecex = 2)


# imputed nutrients data by CRM:
nutrients.imputed <- exp(crm_hampel$data.imputed)
colnames(nutrients.imputed) <- colnames(nutrients)
nutrients.imputed_casesOfInterest <- round(nutrients.imputed[casesOfInterest_hampel, -1], 1)
rownames(nutrients.imputed_casesOfInterest) <- rownames(nutrients_casesOfInterest)

cellwiseheatmap(cellwiseoutliers = crm_hampel$cellwiseoutliers[casesOfInterest_hampel, ],
                col.scale.factor = 1/5,
                data = nutrients.imputed_casesOfInterest,
                margins = c(10, 22), notecex = 2)



# Cross-validation --------------------------------------------------------------------------------
library(CRMwf)
library(caret)
library(cellWise)
library(robustbase)

nfolds <- 10
set.seed(2019)
folds <- createFolds(y = transNutrients$log.cholesterol, k = nfolds)

RMSE <- data.frame("CRM-Hampel"     = rep(NA, nfolds),
                   "CRM-Tukey"      = rep(NA, nfolds),
                   "CRM-Huber"      = rep(NA, nfolds),
                   "CRM-Gauss"      = rep(NA, nfolds),
                   "CRM-Quadratic"  = rep(NA, nfolds))
names(RMSE) <- c("CRM-Hampel", "CRM-Tukey", "CRM-Huber", "CRM-Gauss", "CRM-Quadratic")

set.seed(2019)
for (k in 1:nfolds) {
  cat(paste("\n Fold", k, "out of", nfolds, "\n"))
  
  train <- transNutrients[-folds[[k]], ]
  test  <- transNutrients[ folds[[k]], ]
  
  
  crm_hampel_train <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                                      data      = train,
                                                      maxiter   = maxiter,
                                                      tolerance = tolerance,
                                                      outlyingness.factor = outlyingness.factor,
                                                      spadieta  = spadieta,
                                                      center    = "median",
                                                      scale     = "sd",
                                                      regtype   = "MM",
                                                      seed      = 2019,
                                                      verbose   = FALSE,
                                                      weightfunction = "hampel",
                                                      weightthreshold = 1,
                                                      parameters = c(1.382, 2.764, 5.528)))
  
  crm_tukey_train <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                                     data      = train,
                                                     maxiter   = maxiter,
                                                     tolerance = tolerance,
                                                     outlyingness.factor = outlyingness.factor,
                                                     spadieta  = spadieta,
                                                     center    = "median",
                                                     scale     = "sd",
                                                     regtype   = "MM",
                                                     seed      = 2019,
                                                     verbose   = FALSE,
                                                     weightfunction = "tukey",
                                                     weightthreshold = 0.7,
                                                     parameters = c(4.685)))
  
  crm_huber_train <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                                     data      = train,
                                                     maxiter   = maxiter,
                                                     tolerance = tolerance,
                                                     outlyingness.factor = outlyingness.factor,
                                                     spadieta  = spadieta,
                                                     center    = "median",
                                                     scale     = "sd",
                                                     regtype   = "MM",
                                                     seed      = 2019,
                                                     verbose   = FALSE,
                                                     weightfunction = "huber",
                                                     weightthreshold = 1,
                                                     parameters = c(1.345)))
  
  crm_gauss_train <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                                     data      = train,
                                                     maxiter   = maxiter,
                                                     tolerance = tolerance,
                                                     outlyingness.factor = outlyingness.factor,
                                                     spadieta  = spadieta,
                                                     center    = "median",
                                                     scale     = "sd",
                                                     regtype   = "MM",
                                                     seed      = 2019,
                                                     verbose   = FALSE,
                                                     weightfunction = "gauss",
                                                     weightthreshold = 1,
                                                     parameters = c(1.063, 1.387, 1.5)))
  
  crm_quadratic_train <- suppressWarnings(crm_functional(formula   = log.cholesterol ~ .,
                                                         data      = train,
                                                         maxiter   = maxiter,
                                                         tolerance = tolerance,
                                                         outlyingness.factor = outlyingness.factor,
                                                         spadieta  = spadieta,
                                                         center    = "median",
                                                         scale     = "sd",
                                                         regtype   = "MM",
                                                         seed      = 2019,
                                                         verbose   = FALSE,
                                                         weightfunction = "quadratic",
                                                         weightthreshold = 1,
                                                         parameters = c(0.982, 1.473, 1.5)))
  
  
  pred_test_crm_hampel      <- predict(crm_hampel_train,     newdata = test)
  pred_test_crm_tukey       <- predict(crm_tukey_train,      newdata = test)
  pred_test_crm_huber       <- predict(crm_huber_train,      newdata = test)
  pred_test_crm_gauss       <- predict(crm_gauss_train,      newdata = test)
  pred_test_crm_quadratic   <- predict(crm_quadratic_train,  newdata = test)
  
  RMSE$'CRM-Hampel'[k]      <- sqrt(mean((test$log.cholesterol - pred_test_crm_hampel)^2,     trim = 0.1))
  RMSE$'CRM-Tukey'[k]       <- sqrt(mean((test$log.cholesterol - pred_test_crm_tukey)^2,      trim = 0.1))
  RMSE$'CRM-Huber'[k]       <- sqrt(mean((test$log.cholesterol - pred_test_crm_huber)^2,      trim = 0.1))
  RMSE$'CRM-Gauss'[k]       <- sqrt(mean((test$log.cholesterol - pred_test_crm_gauss)^2,      trim = 0.1))
  RMSE$'CRM-Quadratic'[k]   <- sqrt(mean((test$log.cholesterol - pred_test_crm_quadratic)^2,  trim = 0.1))
}

boxplot(RMSE, ylab = "10% trimmed RMSEP", col = "olivedrab4",
        ylim = c(0, max(RMSE)), cex.lab = 1.4, cex.axis = 1.4, cex.main = 2)
points(colMeans(RMSE), pch = 18, cex = 1.5)
text(rep(0, ncol(RMSE)), labels = round(colMeans(RMSE), 4), cex = 2,
     font = ifelse(1:ncol(RMSE) == which.min(colMeans(RMSE)), 2, 1))

