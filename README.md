# CRMsimulations

This repository contains simulation studies where the performance of cellwise robust M (CRM) regression using different weight functions is evaluated. In these simulation studies the package `CRMwf` is used, which can be found on https://github.com/j4c0d3h00g/CRMwf. The repository `CRMsimulations`contains the following R files:
- `CRM_plot_weightfunctions.R`: creates plots of the different weight functions for varying parameters.
- `CRM_simulation_comparison.R`: compares the use of different weight functions in CRM based on the MSEP, MAE, precision, and recall. 
- `CRM_simulation_efficiency.R`: evaluates the use of different weight functions in CRM based on the MSEP, MAE, precision, and recall. Here the parameters are selected based on different efficiency levels. 
- `CRM_simulation_breakdown.R`: evaluates the use of different weight functions in CRM based on the MSEP, MAE, precision, and recall. Here the parameters are selected based on different breakdown points.
- `CRM_simulation_contamination_situation1.R`: evaluates the breakdown behavior of CRM using different weight functions. Here, only the amount of casewise contamination is increased.
- `CRM_simulation_contamination_situation2.R`: evaluates the breakdown behavior of CRM using different weight functions. Here, only the amount of cellwise contamination is increased.
- `CRM_simulation_contamination_situation3.R`: evaluates the breakdown behavior of CRM using different weight functions. Here, both the amount of casewise and cellwise contamination is increased.
- `CRM_simulation_alpha.R`: evaluates whether the robustness of the CRM regression estimator can increase when the parameters for the weight functions are lowered.
- `CRM_simulation_realdata.R`: evaluates the predictive performance of CRM using different weight functions in a real data application. To illustrate how CRM regression works in practice, a cellwise heatmap that shows which values are imputed for the outlying cells is included for CRM using the Hampel weight function. 
