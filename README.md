## Sharing the effort of the European Green Deal among countries
### Steininger, K.W., Williges, K., Meyer, L.H., Maczek, F., and Riahi, K.
This repository consists of a static version of all code used in creation of the results as reported and analysed in the paper, as well as the necessary datasets to reproduce all results.

### Code files
The calculations distributing emissions as described in the paper are programmed in python. After calculation of the range of possible combinations of equity interpretations via the python script and/or Jupyter notebook, assessment of negotiation convergence points, final analysis and figures is carried out in R.

All data and code required for the python and R portions of the work are contained in folders `python_code` and `r_code` respectively.

#### Python code
The distribution calculations rely on two modules. The `main.py` module contains functions that calculate carbon budget distributions based on the interpretation and higher level functions that combine several distributions and print results to csv files. The `util.py` module contains globally used variables and provides helper functions for main functions.
  * `EU_data_new.xlsx.xlsx` contains statistics for EU country emission and interpretation specific data
  * `EU_GDP_approach_distribution.xlsx` budget distribution data, based on the 2021 EU proposal (see C1 interpretation in paper)
  * `capacity_indicators_18_01_2022.csv` governance capacity indicators data used for "C3 - Government effectiveness" interpretation

##### Jupyter notebooks
The `ES_Tool_v4.ipynb` notebook has been used to calculate the results for the publication. It can be used to experiment with different interpretation combinations to see the effects on 2030 carbon budgets of EU member states. The notebook comes with a manual (in comments) on how to use the functions and parameters.

#### R code
Two R scripts are included. The first, `analysis.R`, contains all code needed to find the negotiation convergence points as discussed in the paper, as well as to produce all figures found in the main text and Supplementary Information.

The second script, `shiny_app.R`, consists of the code necessary to run the accompanying shiny app.

##### Data file index:
The following files are needed to fully run `analysis.R`
  * `./data/all_absolutes_with_RES.rds`
  * `./data/all_shares_with_RES.rds`
  * `./data/emi_pop_data.xlsx`
  * `./data/opt_results/max_optimals_new.rds`
  * `./data/opt_results/esd_optimals_new.rds`
  * `./data/single_interpretations.rds`
  * `./data/sensitivity/h_year_sens.csv`
  * `./data/sensitivity/reductions_2030_new_gov_indic_with_zero_res.csv`

##### Shiny app
The following files are required by the shiny app:
  * `all_shares_with_RES.rda`
  * `max_optimals_new.rda`
  * `esd_optimals_new.rda`
  * `fig4_boxData.rda`
  * `popData.rda`
