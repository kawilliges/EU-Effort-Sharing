# R script to calculate possible negotiation convergence points, and to 
# plot all figures included in the paper "Sharing the effort of the 
# European Green Deal among countries" in Nature Communications. 

wkDir <- 'C:/GitHub/python_effort_sharing_tool/r_analysis/'
setwd(wkDir)
getwd()

library(tidyverse)
library(ggthemes)
library(readxl)
library(stringr)
library(patchwork)
library(NMOF)
library(optimx)
library(viridis)
library(plyr)
library(Ternary)
library(readr)

# 1. Data import ---------------------------------------------------------------

# all scenarios run from python script, e.g. all combinations of equity 
# interpretations, saved into one dataset. Two datasets are included,
# "all_absolutes" and "all_shares", which hold emissions in quantities and
# shares of total EU ES emissions, respectively

# all_absolutes needed for calculation of negotiation points only
all_absolutes <- readRDS("./data/all_absolutes_with_RES.rds")

# all_shares used for plotting figures
all_shares <- readRDS("./data/all_shares_with_RES.rds")


# 1.3 Other needed data creation / import                                 -----
# Define the current Fit for 55 reduction targets
esdReduct <- data.frame('countries' = c('BEL', 'BGR', 'CZE', 'DNK', 'DEU', 
                                        'EST', 'IRL', 'GRC', 'ESP', 'FRA', 
                                        'HRV', 'ITA', 'CYP', 'LVA', 'LTU', 
                                        'LUX', 'HUN', 'MLT', 'NLD', 'AUT', 
                                        'POL', 'PRT', 'ROU', 'SVN', 'SVK', 
                                        'FIN', 'SWE'),
                        'value' = c(-0.47, -0.10, -0.26, -0.50, -0.50, 
                                    -0.24, -0.42, -0.227, -0.377, -0.475, 
                                    -0.165, -0.437, -0.32, -0.17, -0.21, 
                                    -0.50, -0.187, -0.19,  -0.48, -0.48, 
                                    -0.177, -0.287, -0.127, -0.27, -0.227,
                                    -0.50, -0.50)
)
esdReduct$type = 'esd'
esdReduct$fac <- '55'

# Define the previous (2018) agreed-upon reduction targets
esdReductOld <- data.frame('countries' = c('BEL', 'BGR', 'CZE', 'DNK', 'DEU', 
                                           'EST', 'IRL', 'GRC', 'ESP', 'FRA', 
                                           'HRV', 'ITA', 'CYP', 'LVA', 'LTU', 
                                           'LUX', 'HUN', 'MLT', 'NLD', 'AUT', 
                                           'POL', 'PRT', 'ROU', 'SVN', 'SVK', 
                                           'FIN', 'SWE'),
                           'value' = c(-0.35, 0, -0.14, -0.39, -0.38, -0.13, 
                                       -0.30, -0.16, -0.26, -0.37, -0.07, 
                                       -0.33, -0.24, -0.06, -0.09, -0.40, -0.07,
                                       -0.19,  -0.36, -0.36, -0.07, -0.17, 
                                       -0.02, -0.15, -0.12, -0.39, -0.40)
)
esdReductOld$type = 'esd'
esdReductOld$fac <- '2018'

# Import info on population data and 2005 emissions
popData <- read_xlsx('./data/emi_pop_data.xlsx')

## Define a list of country ISOs
countryISOs <- c('AUT','BEL','BGR','CYP','CZE','DEU','DNK','ESP','EST','FIN',
                 'FRA','GRC','HRV','HUN','IRL','ITA','LTU','LUX','LVA','MLT',
                 'NLD','POL','PRT','ROU','SVK','SVN','SWE')

# this function un-needed if the above datasets are imported. Only used
# if new scenarios are generated and should be added to main dataset.
dataProcessing <- function(){
  # 1.2. Importing all possible scenarios                                   -----
  # create a list with all datasets; first get vector of all filenames, then
  # import all into a dataframe
  
  # switch wkDir to /data to make it easier for the moment
  setwd('C:/GitHub/python_effort_sharing_tool/r_analysis/data/all_runs')
  
  input_share_csvs <- list.files(pattern = '\\-reduction-_zTrue.csv$')
  input_absolute_csvs <- list.files(pattern = '\\-absolute-_zTrue.csv$')
  all_shares <- lapply(input_share_csvs, read_csv)
  all_absolutes <- lapply(input_absolute_csvs, read_csv)
  
  names(all_shares) <- stringr::str_replace(input_share_csvs, pattern = ".csv", 
                                            replacement = "")
  
  names(all_absolutes) <- stringr::str_replace(input_absolute_csvs, 
                                               pattern = ".csv", replacement = "")
  
  setwd(wkDir)
  
  saveRDS(all_absolutes, file = "all_absolutes.rds")
  saveRDS(all_shares, file = 'all_shares.rds')
  
  all_absolutes <- readRDS(".data/all_absolutes.rds")
  all_shares <- readRDS(".data/all_shares.rds")
  
  # 1.2.1 Adding new scenarios if more are created
  # (read in and append to all_abs and all_shr)
  setwd('C:/GitHub/python_effort_sharing_tool/r_analysis/data/new_runs')
  
  input_share_csvs <- list.files(pattern = '\\-reduction-_zTrue.csv$')
  input_absolute_csvs <- list.files(pattern = '\\-absolute-_zTrue.csv$')
  all_new_shares <- lapply(input_share_csvs, read_csv)
  all_new_absolutes <- lapply(input_absolute_csvs, read_csv)
  
  names(all_new_shares) <- stringr::str_replace(input_share_csvs, pattern = ".csv", 
                                            replacement = "")
  
  names(all_new_absolutes) <- stringr::str_replace(input_absolute_csvs, 
                                               pattern = ".csv", replacement = "")
  
  all_absolutes_new <-c(all_absolutes, all_new_absolutes)
  all_shares_new <- c(all_shares, all_new_shares)
  
  setwd(wkDir)
  
  # now save the new dataset!!!
  saveRDS(all_absolutes_new, file = "all_absolutes_with_RES.rds")
  saveRDS(all_shares_new, file = 'all_shares_with_RES.rds')
}


# 1.1 Data functions -----------------------------------------------------------
# Function to find each country's maximum allocation - needed for some 
# data filtering and figure-making below
findTotQual <- function(inputData, capVar, needVar, respVar) {
  return(inputData %>% mutate(totQual = (as.numeric(inputData[[capVar]]) + 
                                           as.numeric(inputData[[needVar]]) + 
                                           as.numeric(inputData[[respVar]]))) %>% 
           subset(totQual < 1.04 & totQual > 0.96))
}

# function to find only the country preferred allocation, given inputs, 
# as above, needed for figure generation
getOnlyMaxes <- function(inputData, capVar, needVar, respVar) {
  return(findTotQual(inputData, capVar, needVar, respVar) %>% 
           group_by(countries) %>% slice(which.max(value)))
}


### Optimal calc - Min dist to preference --------------------------------------
### Find the optimal combination of qualifications which minimizes the average
### relative (%) loss of countries

# Function to return the average relative loss
# n: needs qualification level
# c: capability qualification level
# r: responsibility qualification level
# absRes: dataset of all possible results

# NOTE: This is included for documentation purposes only. Code execution
# takes considerable time, and for convenience the results are included
# as a dataset (imported below)

# get the max possible for each country
mainOpt <- function(scen_to_opt, maxOrFit, cVar, eVar, rVar){
  avgRelLoss <- function(x){

    if (maxOrFit == 'max'){
      # get the max possible for each country
      countryMax <- getOnlyMaxes(scen_to_opt, cVar, eVar, rVar)
    } else if (maxOrFit == 'fit'){
      # get the amount implied by the Fit for 55% 
      countryMax <- esdReductOld %>% left_join(popData, esdReductOld, 
                                            by = "countries") %>%
        mutate(value = emi_2005 * (1+value)) %>%
        subset(select = c(countries, value, type))
    }
    
    # select just the resulting distribution for countries given the inputs
    qualResult <- scen_to_opt %>% mutate(totQual = (scen_to_opt[[cVar]] +
                                                      scen_to_opt[[eVar]] + 
                                                      scen_to_opt[[rVar]])) #%>% 
    #  subset(totQual == 1.00) 
    
    qualResult <- qualResult %>% subset(totQual <  1.04 & totQual > 0.96)
    
    qualResult <- qualResult %>% filter((qualResult[[cVar]] > (x[[2]]-0.01) & 
                                           qualResult[[cVar]] < (x[[2]]+0.01)) &
                                          (qualResult[[eVar]] > (x[[1]]-0.01) &
                                             qualResult[[eVar]] < (x[[1]]+0.01)) &
                                          (qualResult[[rVar]] > (x[[3]]-0.01) &
                                             qualResult[[rVar]] < (x[[3]]+0.01))) %>% 
      subset(select = c(countries, value))
    
    names(qualResult)[names(qualResult) == "value"] <- "qual_result"
    
    # join the two dfs (maxes vs result of a specific weighting combo)
    join <- left_join(countryMax, qualResult, by = 'countries')
    
    # join the population data and calculate per capita values
    join <- left_join(join, popData, by = 'countries')
    
    join <- join %>% mutate(value = value / pop)
    join <- join %>% mutate(qual_result = qual_result / pop)
    
    # now calculate the difference bewtween the max value and the weight 
    # combo, giving more weight to higher differences
    join <- join %>% mutate(diff = ((value - qual_result) / value)^2)
    
    # find the average loss for the weighting combo
    avg_loss <- sum(join$diff) / nrow(join)
    return(avg_loss)
  }
  return(gridSearch(avgRelLoss, list(seq(0,1, 0.05), seq(0,1,0.05), 
                                     seq(0,1,0.05)))$minlevels)
}

# Grid search over params find that the optimal qualifications leading to 
# the smallest average relative loss in emissions budget is:

# 1. number: Equality
# 2. number: Capability
# 3. number: Responsibility

#  -----------------------------------------------------------------#
# create function to apply mainOpt() to all items in list, but to do so,
# first turn list entry (dataframe) name into character string and 
# retrieve the relevant parameters (e.g. what specifications used)

OptCall <- function(inputDF, mainOrFit = 'max'){
  # function takes input dataframe, returns a result of mainOpt function
  # using "max" allocation 
  dfToText <- deparse(names(inputDF))
  params <- dfToText %>% str_replace_all("\"", "") %>% 
    str_replace_all("-_zTrue", "") %>% str_split_fixed("-", 4)

  return(mainOpt(inputDF[[1]], mainOrFit, params[1], params[2], params[3]))
}

allOptMax <- function(){
  # Now run optimisation for all scenarios; first define the list which
  # will hold the results data, then run a loop to append all subsequent
  # results to the existing list (sloppy, but lapply doesn't work 
  # well here because of use of each scenario name - not passed by 
  # lapply)
  maxResultsTmp <- list(OptCall(all_absolutes[1]))

  for(x in 2:length(all_absolutes)) {
    maxResultsTmp <- c(maxResultsTmp, list(OptCall(all_absolutes[x])))    
  }
  
  # re-name the results to conform to scenarios
  names(maxResultsTmp) <- names(all_absolutes)
  return(maxResultsTmp)
}

# Now the same, for the 2018 ESD agreed targets
allOptEsd <- function(){
  esdOptResultsTmp <- list(OptCall(all_absolutes[1], 'fit'))
  
  for(x in 2:length(all_absolutes)) {
    esdOptResultsTmp <- c(esdOptResultsTmp, list(OptCall(all_absolutes[x], 'fit')))    
  }
  
  # re-name the results to conform to scenarios
  names(esdOptResultsTmp) <- names(all_absolutes)
  return(esdOptResultsTmp)
}

#maxResults <- allOptMax()
#esdOptResults <- allOptEsd()

# since we have run the functions already, load the saved data -----
maxResults <- readRDS(file = './data/opt_results/max_optimals_new.rds')
esdOptResults <- readRDS(file = './data/opt_results/esd_optimals_new.rds')

# FIGURE 2 ---------------------------------------------------------------------
# Plot the result of a 100% interpretation (for all interpretation), plus
# indicate the current Fit for 55% target
plot100 <- function(full_spec_df){
  
  # function to retrieve just the 100% specification data
  getAllocation <- function(cVar, cLvl, eVar, eLvl, rVar, rLvl, results)
  {
    return(results %>% filter((results[[cVar]] > (cLvl - 0.01) & 
                                 results[[cVar]] < (cLvl + 0.01)) & 
                                (results[[eVar]] > (eLvl - 0.01) & 
                                   results[[eVar]] < (eLvl + 0.01)) & 
                                (results[[rVar]] > (rLvl - 0.01) & 
                                   results[[rVar]] < (rLvl + 0.01))) %>%
             subset(select = c(countries, value)))
  }

  # now make point plots; divided into 3 panels based on if
  # majority of points are above/equal/below the Fit for 55% target
  
  # define list of breaks, labels, and colors for the charts
  breakList <- c("c_all", 'cv1_all', "cv2_all", 'cv3_all', "n_all", "nv1_all", 
                 'nv2_all', 'r_all', 'rv1_all', 'rv2_all', 'rv3_all', 
                 'rv4_all', 'esd')
  
  labelList <- c("C1-EU-capability", "C2-GDP/cap", "C3-Governance", 'C4-RES-cap',
                 "E1-Basic-needs", "E2-ES-EPC", 'E3-Full-EPC', 
                 'R1-Hist-emi', 'R2-Benefits', 'R3-C-budget', 'R4-RES-expansion', 
                 'R5-Budget/cap', 'Fit for 55')
  
  colorList <- c("#000000", '#000000', '#737373', '#bdbdbd', 
                 "#000000", "#000000", "#bdbdbd", 
                 '#000000', '#000000', '#737373', '#bdbdbd', '#C0C0C0', 
                 'darkgrey')
  
  shapeList <- c(22, 0, 0, 22, 
                 24, 2, 2, 
                 21, 1, 1, 1, 21, 
                 23)
  
  # panel A - countries mainly above the FitFor55 target:
  # BGR, ESP, MLT, PRT, SWE, ITA, CYP
  
  above_fit55 <- full_spec_df %>% filter(countries %in% c('ROU', 'BGR', 'ESP', 
                                                          'MLT', 'PRT', 'SWE',
                                                          'ITA', 'SVK', 'LTU')) 
  above_fit55plot <- ggplot(data = filter(above_fit55, type != 'esd'), aes(fac, value, 
                                                   fill = type, color = type,
                                                   shape = type)) + 
    geom_point(size = 3, stroke = 2, position=position_dodge(width = 0.8)) +
    xlab(NULL) + ylab(NULL) + 
    scale_fill_manual(values = colorList, name ="Interpretation", 
                      breaks = breakList, labels = labelList) + 
    scale_color_manual(values = colorList, name="Interpretation",
                       breaks = breakList, labels= labelList) +
    scale_shape_manual(values = shapeList, 
                       name="Interpretation", breaks = breakList,
                       labels = labelList) + 
    theme_bw() + geom_hline(data = filter(above_fit55, type == 'esd'), 
                            aes(yintercept = value), color = '#525252', 
                            linetype = 'dashed', size = 1) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme(legend.position = "none") + facet_grid(~countries)
  
  above_fit55plot
  # panel B - countries around equal / targets to either side: EST
  # DEU, GRC, HRV, LTV, LVA, POL, ROU, SVK, HUN, FRA
  equal_fit55 <- full_spec_df %>% filter(countries %in% c('DEU', 'GRC', 'HRV',
                                                          'FRA', 'LTV', 'LVA',
                                                          'POL', 'CYP', 'EST', 
                                                          'HUN'))
  equal_fit55plot <- ggplot(data = filter(equal_fit55, type != 'esd'), aes(fac, value, 
                                             fill = type, color = type,
                                             shape = type)) + 
    geom_point(size = 3, stroke = 2, position=position_dodge(width = 0.8)) +
    ylab(NULL) + xlab(NULL) + 
    scale_fill_manual(values = colorList, name ="Interpretation", 
                      breaks = breakList, labels = labelList) + 
    scale_color_manual(values = colorList, name="Interpretation",
                       breaks = breakList, labels= labelList) +
    scale_shape_manual(values = shapeList, 
                       name="Interpretation", breaks = breakList,
                       labels = labelList) + 
    theme_bw() + geom_hline(data = filter(equal_fit55, type == 'esd'), 
                            aes(yintercept = value), color = '#525252', 
                            linetype = 'dashed', size = 1) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme(legend.position = "none") + facet_grid(~countries)
  
  
  
  # panel C: countries below mainly:
  # AUT, BEL, DNK, FIN, IRL, LUX, NLD, CZE, FRA, SVN
  below_fit55 <- full_spec_df %>% filter(countries %in% c('AUT', 'BEL', 'DNK',
                                                             'FIN', 'IRL', 'LUX',
                                                             'NLD', 'CZE', 'SVN'))
  below_fit55plot <- ggplot(data = filter(below_fit55,  type != 'esd'), aes(fac, value, 
                                              fill = type, color = type,
                                              shape = type)) + 
    geom_point(size = 3, stroke = 2, position=position_dodge(width = 0.8)) +
    xlab('Country') + ylab(NULL) + 
    scale_fill_manual(values = colorList, name ="Interpretation", 
                      breaks = breakList, labels = labelList) + 
    scale_color_manual(values = colorList, name="Interpretation",
                       breaks = breakList, labels= labelList) +
    scale_shape_manual(values = shapeList, 
                       name="Interpretation", breaks = breakList,
                       labels = labelList) + 
    theme_bw() + geom_hline(data = filter(below_fit55, type == 'esd'), 
                            aes(yintercept = value), color = '#525252', 
                            linetype = 'dashed', size = 1) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme(legend.position = "bottom") + facet_grid(~countries)
  
  
  
  allPlots <- above_fit55plot / equal_fit55plot / below_fit55plot
  allPlots + plot_annotation(tag_levels = 'A')
  return(allPlots)
}

# Plot Figure 2 ----------------------------------------------------------------
# load the pre-calculated and saved data so we don't have to run python scripts
full_spec <- readRDS('./data/single_interpretations.rds')

# Plots the resulting country allocation given a 100% implementation of a 
# specified interpretation
plotInterpretations <- plot100(full_spec) + plot_annotation(tag_levels = 'A') 

plotInterpretations 

# FIGURE 3 --------------------------------------------------------------------


# See section "SI FIGURES" where method for generating any ternary chart
# figures are specified

# FIGURE 4 --------------------------------------------------------------------
# Figure 4 inset panel                                                     ----
# Ternary plot of where all optimals are --

# Add data points (data stored as N,C,R) but here need to be C, R, N)
points_MaxAlloc <- list( 
   c(maxResults[[1]][2] * 100, maxResults[[1]][1] * 100, 
     maxResults[[1]][3] * 100)
)

for( x in 2:length(maxResults)) {
  points_MaxAlloc <- c(points_MaxAlloc, list(c(maxResults[[x]][2] * 100, 
                                               maxResults[[x]][1] * 100, 
                                               maxResults[[x]][3] * 100)))
}

names(points_MaxAlloc) <- names(maxResults)

points_esdAlloc <- list(
  c(esdOptResults[[1]][2] * 100, esdOptResults[[1]][1] * 100, 
    esdOptResults[[1]][3] * 100)
)

for( x in 2:length(esdOptResults)) {
  points_esdAlloc <- c(points_esdAlloc, list(c(esdOptResults[[x]][2] * 100, 
                                               esdOptResults[[x]][1] * 100, 
                                               esdOptResults[[x]][3] * 100)))
}

names(points_esdAlloc) <- names(esdOptResults)


# try to change the sizes of points depending on how many there are
test <- do.call(rbind.data.frame, points_esdAlloc)
n_occur <- data.frame(table(test))
n_occur <- n_occur %>% subset(Freq != 0)

colnames(n_occur)[1] <- "C" 
colnames(n_occur)[2] <- "E" 
colnames(n_occur)[3] <- "R" 

test2 <- do.call(rbind.data.frame, points_MaxAlloc)
n_occur2 <- data.frame(table(test2))
n_occur2 <- n_occur2 %>% subset(Freq != 0)

colnames(n_occur2)[1] <- "C" 
colnames(n_occur2)[2] <- "E" 
colnames(n_occur2)[3] <- "R" 


TernaryPlot(alab = paste("Increasing level of capability"), 
            blab = paste("Increasing level of equality"), 
            clab = paste("Increasing level of responsibility"), 
            atip = 'C', btip = 'E', ctip = 'R')

AddToTernary(points, points_MaxAlloc, pch = 21, cex = sqrt(n_occur2$Freq), 
             bg = c(rep('#f46d43', length(maxResults))
))

AddToTernary(points, points_esdAlloc, pch = 21, cex = sqrt(n_occur$Freq), 
             bg = c(rep('#74add1', length(esdOptResults))
))




# Figure 4 main panel------------------------------------------------------------
# Box-whisker plot of country results for all optimal points

# First, get the relevant data needed for box whiskers (country reductions
# at optimal points)

getCorners <- function(inputName){
  # takes scenario name and splits into the four parameters it contains
  # (Cap, Equ, Res, and type - absolute or reduction)
  dfToText <- deparse(names(inputName))
  return(dfToText %>% str_replace_all("\"", "") %>% str_split_fixed("-", 4))
}

maxOptCountryAlloc <- getAllocation(getCorners(maxResults[1])[1],
                                    maxResults[[1]][2],
                                    getCorners(maxResults[1])[2],
                                    maxResults[[1]][1],
                                    getCorners(maxResults[1])[3],
                                    maxResults[[1]][3],
                                    all_shares[[1]])

for(x in 2:length(maxResults)) {
  tempMatrix <- getAllocation(getCorners(maxResults[x])[1],
                              maxResults[[x]][2],
                              getCorners(maxResults[x])[2],
                              maxResults[[x]][1],
                              getCorners(maxResults[x])[3],
                              maxResults[[x]][3],
                              all_shares[[x]])
  
  maxOptCountryAlloc <- bind_rows(maxOptCountryAlloc, tempMatrix)
}
maxOptCountryAlloc$type <- "Proximity to upper \nbound of equity-\ncompatible emissions"

esdOptCountryAlloc <- getAllocation(getCorners(esdOptResults[1])[1],
                                    esdOptResults[[1]][2],
                                    getCorners(esdOptResults[1])[2],
                                    esdOptResults[[1]][1],
                                    getCorners(esdOptResults[1])[3],
                                    esdOptResults[[1]][3],
                                    all_shares[[1]])

for(x in 2:length(esdOptResults)) {
  tempMatrix <- getAllocation(getCorners(esdOptResults[x])[1],
                              esdOptResults[[x]][2],
                              getCorners(esdOptResults[x])[2],
                              esdOptResults[[x]][1],
                              getCorners(esdOptResults[x])[3],
                              esdOptResults[[x]][3],
                              all_shares[[x]])
  
  esdOptCountryAlloc <- bind_rows(esdOptCountryAlloc, tempMatrix)
}
esdOptCountryAlloc$type <- "Minimal change\nfrom the 2018 ESR"

boxWhiskerFig2Data <- bind_rows(maxOptCountryAlloc, esdOptCountryAlloc)

boxWhiskerFig2Data <- left_join(boxWhiskerFig2Data, esdReduct, by = 'countries')

boxWhiskerFig2Data$countries <- factor(boxWhiskerFig2Data$countries)

legendText <- 'Possible emission reductions resulting \nfrom negotiation points when considering:'

boxWhiskFig2 <- ggplot(boxWhiskerFig2Data, 
                       aes(x = reorder(countries, -value.y), y = value.x, 
                           color = type.x, fill = type.x)) + 
  geom_boxplot(coef = 5, outlier.shape = NA) +
  geom_step(data = boxWhiskerFig2Data, 
             aes(y = value.y, x= countries, group = 1), color = 'black', 
            linetype = 'solid', size = 1.05) + 
  xlab('Country') + ylab('Change compared to 2005 emissions') +
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(legend.position = 'bottom') + 
  scale_color_manual(values = c('#4575b4', '#d73027')) +
  scale_fill_manual(values = c('#abd9e9', '#fdae61')) +
  labs(color = legendText, 
       fill = legendText) +
  geom_hline(aes(yintercept = -0.44), linetype = 'dotted', size = 1.05)

boxWhiskFig2

# SI FIGURES ------------------------------------------------------------------
# Functions to plot and save ternary charts for each country

# create a dictionary so I don't have to type in all the labels
library(hash)
labD <- hash()
labD[['EU_GDPPOP']] <- 'C1-EU-cap.'
labD[['GDPPOP_div']] <- 'C2-GDP/cap'
labD[['gee_n']] <- 'C3-Gov.'
labD[['RES_cap']] <- 'C4-RES-cap'
labD[['1995']] <- 'R1-Hist-emi'
labD[['BENEFITS']] <- 'R2-Benefits'
labD[['CBUDGET']] <- 'R3-C-budget'
labD[['RES']] <- 'R4-RES-exp.'
labD[['CBUD_stefan']] <- 'R5-Budget/cap'
labD[['MIN_WELF']] <- 'E1-Basic needs'
labD[['CPOP']] <- 'E2-ES-EPC'
labD[['CPOP_stefan']] <- 'E3-Full-EPC'

# function to plot a ternary chart, given 4 args:
# countryCode - 3 letter ISO code for country
# cVar, nVar, rVar - codes for the three equity interpretations, 
# as listed in the dictionary above (e.g. "EU_GDPPOP")
ternNew <- function(countryCode, cVar, nVar, rVar) {

  scen_index <- which(names(all_shares) == paste(cVar, nVar, rVar,
                                                 'reduction-_zTrue', sep = '-'))
  temp <- all_shares[[scen_index]]
  
  # Now select just the country of interest
  country <- countryCode
  
  temp <- temp %>% filter(countries == country)
  
  # Now just select the 100% values for each of c, e, and r
  # first get only the fully qualified versions
  temp <- findTotQual(temp, cVar, nVar,
                      rVar)
  
  temp <- temp %>% mutate(value = round(value, 4))
  cap <- (temp[temp[[cVar]] > 1.00, ])$value
  res <- (temp[temp[[rVar]] > 1.00, ])$value
  equ <- (temp[temp[[nVar]] > 1.00, ])$value
  
  # determine the optimals for max and esd versions
  
  scen_index_max <- which(names(maxResults) == paste(cVar, nVar,rVar,
                                                     'absolute-_zTrue', sep = '-'))
  scen_index_esd <- which(names(esdOptResults) == paste(cVar, nVar, rVar,
                                                        'absolute-_zTrue', sep = '-'))
  # get just the chosen combination's results
  max_reduct_scen <- maxResults[[scen_index_max]]
  esd_reduct_scen <- esdOptResults[[scen_index_esd]]
  
  # save as data_points
  data_points <- list(M = c(max_reduct_scen[2]*100, max_reduct_scen[1]*100,
                            max_reduct_scen[3]*100),
                      E = c(esd_reduct_scen[2]*100, esd_reduct_scen[1]*100,
                            esd_reduct_scen[3]*100))
  
  functionToContour <- function(c,n,r) {
    c * cap + n * equ + r * res
  }
  
  
  TernaryPlot(alab = paste('Increasing level of', labD[[cVar]]),
              blab = paste('Increasing level of', labD[[nVar]]),
              clab = paste('Increasing level of', labD[[rVar]]),
              atip = 'C', btip = 'E', ctip = 'R')
  
  values = TernaryPointValues(functionToContour, resolution = 24L)
  ColourTernary(values, spectrum = viridis::cividis(256L, alpha = 0.75))
  TernaryContour(functionToContour, resolution = 36L)
  
  AddToTernary(points, data_points, pch = 21, cex = 2.8,
               bg = c('#f46d43','#74add1'))

}

# Plot an example country ternary -----
ternNew('AUT', "GDPPOP_div", "MIN_WELF", '1995')



# function to plot only the cases where no variation occurs, due to 
# e.g. indifference between two corners
ternSpecial <- function(countryCode, cVar, nVar, rVar) {
  
  scen_index <- which(names(all_shares) == paste(cVar, nVar, rVar,
                                                 'reduction-_zTrue', sep = '-'))
  temp <- all_shares[[scen_index]]
  
  # Now select just the country of interest
  country <- countryCode
  
  temp <- temp %>% filter(countries == country)
  
  # Now just select the 100% values for each of c, e, and r
  # first get only the fully qualified versions
  temp <- findTotQual(temp, cVar, nVar,
                      rVar)
  
  temp <- temp %>% mutate(value = round(value, 4))
  cap <- (temp[temp[[cVar]] > 1.00, ])$value
  res <- (temp[temp[[rVar]] > 1.00, ])$value
  equ <- (temp[temp[[nVar]] > 1.00, ])$value
  
  # determine the optimals for max and esd versions
  
  scen_index_max <- which(names(maxResults) == paste(cVar, nVar,rVar,
                                                     'absolute-_zTrue', sep = '-'))
  scen_index_esd <- which(names(esdOptResults) == paste(cVar, nVar, rVar,
                                                        'absolute-_zTrue', sep = '-'))
  # get just the chosen combination's results
  max_reduct_scen <- maxResults[[scen_index_max]]
  esd_reduct_scen <- esdOptResults[[scen_index_esd]]
  
  # save as data_points
  data_points <- list(M = c(max_reduct_scen[2]*100, max_reduct_scen[1]*100,
                            max_reduct_scen[3]*100),
                      E = c(esd_reduct_scen[2]*100, esd_reduct_scen[1]*100,
                            esd_reduct_scen[3]*100))
  
  functionToContour <- function(c,n,r) {
    c * cap + n * equ + r * res
  }
  
  
  TernaryPlot(alab = paste('Increasing level of', labD[[cVar]]),
              blab = paste('Increasing level of', labD[[nVar]]),
              clab = paste('Increasing level of', labD[[rVar]]),
              atip = 'C', btip = 'E', ctip = 'R')
  
  values = TernaryPointValues(functionToContour, resolution = 24L)
  middle_triangle <- matrix(c(
    0, 0, 100,
    0, 100, 0,
    100, 0, 0
  ), ncol = 3, byrow = TRUE)
  TernaryPolygon(middle_triangle, col = '#fde725BF', border = 'grey')
  
  AddToTernary(points, data_points, pch = 21, cex = 2.8,
               bg = c('#f46d43','#74add1'))
  
}
# Save pdfs of special cases, e.g. where normal function would not plot
# the color gradients correctly, due to indifference between 2 corners
pdf(paste('./tern_figs/', 'BGR', '_c2.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('BGR', "GDPPOP_div", "MIN_WELF", '1995')
dev.off()

pdf(paste('./tern_figs/', 'BGR', '_c3.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('BGR', "gee_n", "MIN_WELF", '1995')
dev.off()

pdf(paste('./tern_figs/', 'MLT', '_c2.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('MLT', "GDPPOP_div", "MIN_WELF", '1995')
dev.off()

pdf(paste('./tern_figs/', 'MLT', '_c3.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('MLT', "gee_n", "MIN_WELF", '1995')
dev.off()

pdf(paste('./tern_figs/', 'MLT', '_c4.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('MLT', "RES_cap", "MIN_WELF", '1995')
dev.off()

pdf(paste('./tern_figs/', 'POL', '_c2.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('POL', "GDPPOP_div", "MIN_WELF", '1995')
dev.off()

pdf(paste('./tern_figs/', 'POL', '_c4.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('POL', "RES_cap", "MIN_WELF", '1995')
dev.off()

pdf(paste('./tern_figs/', 'ROU', '_c2.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('ROU', "GDPPOP_div", "MIN_WELF", '1995')
dev.off()

pdf(paste('./tern_figs/', 'ROU', '_c3.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('ROU', "gee_n", "MIN_WELF", '1995')
dev.off()

pdf(paste('./tern_figs/', 'ROU', '_c4.pdf', sep = ''), 
    width = 6, height = 6)
ternSpecial('ROU', "RES_cap", "MIN_WELF", '1995')
dev.off()

# function to save pdfs of all ternary charts of a country 
# for Supplementary Info
saveCharts <- function(countryCode) {
  pdf(paste('./tern_figs/', countryCode, '_main.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, cVar = 'EU_GDPPOP', rVar = '1995', nVar = 'MIN_WELF')
  dev.off()
  
  pdf(paste('./tern_figs/', countryCode, '_e2.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, 'EU_GDPPOP',  rVar ='1995', nVar ='CPOP')
  dev.off()
  
  pdf(paste('./tern_figs/', countryCode, '_r4.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, 'EU_GDPPOP',  rVar ='RES', nVar ='MIN_WELF')
  dev.off()
  
  pdf(paste('./tern_figs/', countryCode, '_r2.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, 'EU_GDPPOP',  rVar ='BENEFITS', nVar ='MIN_WELF')
  dev.off()
  
  pdf(paste('./tern_figs/', countryCode, '_r3.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, 'EU_GDPPOP',  rVar ='CBUDGET', nVar ='MIN_WELF')
  dev.off()
  
  pdf(paste('./tern_figs/', countryCode, '_c2.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, 'GDPPOP_div',  rVar ='1995', nVar ='MIN_WELF')
  dev.off()
  
  pdf(paste('./tern_figs/', countryCode, '_e3.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, 'EU_GDPPOP',  rVar ='1995', nVar ='CPOP_stefan')
  dev.off()
  
  pdf(paste('./tern_figs/', countryCode, '_r5.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, 'EU_GDPPOP',  rVar ='CBUD_stefan', nVar ='MIN_WELF')
  dev.off()
  
  pdf(paste('./tern_figs/', countryCode, '_c3.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, 'gee_n',  rVar ='1995', nVar ='MIN_WELF')
  dev.off()
  
  pdf(paste('./tern_figs/', countryCode, '_c4.pdf', sep = ''), 
      width = 6, height = 6)
  ternNew(countryCode, 'RES_cap',  rVar ='1995', nVar ='MIN_WELF')
  dev.off()
}


# Charts for Figure 3 -----
pdf(paste('./tern_figs/', 'DEU', 'fig3_a.pdf', sep = ''), 
    width = 6, height = 6)
test <- ternNew("DEU", cVar = 'EU_GDPPOP', nVar = 'MIN_WELF', rVar = '1995')
dev.off()

pdf(paste('./tern_figs/', 'DEU', 'fig3_b.pdf', sep = ''), 
    width = 6, height = 6)
test <- ternNew("DEU", cVar = 'EU_GDPPOP', nVar = 'CPOP_stefan', rVar = '1995')
dev.off()


# Ternary charts for all countries for SI ----
lapply(as.list(countryISOs), saveCharts)


# Other figures for SI ------------------------------------------------------------

# Comparing starting years for historical budgets
h_sens <- read_csv('./data/sensitivity/h_year_sens.csv')

interest <- c('CYP', 'CZE', 'DEU', 'ESP', 'FRA', 'GRC', 'ITA', 'SVN')

h_diff <- h_sens %>% subset(countries %in% interest)

diff_plot <- ggplot(h_diff, aes(h_year, value)) + geom_point() +
  facet_wrap(~countries) + labs(x = 'Starting year for historical budgets',
                                y = 'Percentage reduction compard to 2005 emissions') +
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(1990, 2000, 2))

diff_plot

# Comparing RES as a Cap or Resp interpretation
res_sens <- readRDS('./data/single_interpretations.rds')

res_sens <- res_sens %>% subset(type %in% c('cv1_all', 'rv3_all', 'cv3_all'))
res_plot <- ggplot(res_sens, aes(x = reorder(countries, value), y = value)) + 
  geom_col(position = position_dodge(width = 0.75), 
           aes(fill = type)) +
  scale_fill_manual(values = c('#4575b4', 'orange', 'black'), 
                    labels = c('GDP per capita', 'RES (capability)', 
                               'RES (responsibility)'),
                    name = 'Interpretation') + 
  theme_minimal() + labs(y = 'Emission reduction compared to 2005',
                         x = 'Country',
                         fill = 'Interpretation') + 
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  theme(legend.position = 'bottom')

res_plot

# comparing different governance indicators
gov_sens <- read_csv('./data/sensitivity/reductions_2030_new_gov_indic_with_zero_res.csv')

# reshape to long format
gov_long <- gather(gov_sens, index, value, !country)


gov_plot <- ggplot(gov_long, aes(x=index, y=value)) + geom_point() + 
  facet_wrap(~country) + theme_minimal() + 
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position = 'bottom') + 
  labs(x='Index used to generate capability interpretation',
       y='Percent reductions from 2005 emissions') + 
  scale_x_discrete(limits = c("GDP/cap", "log(GDP/cap)", 'State capacity',
                              'Gov. quality', 'Gov. effectiveness'))
  #scale_color_colorblind(name='Interpretation') 

gov_plot
