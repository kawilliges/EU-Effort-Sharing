library(shinythemes)
library(plyr)
library(shiny)
library(ggplot2)
library(ggthemes)
library(tidyverse)
#library(readxl)
library(plotly)
library(hash)
library(Ternary)
#library(shinyBS)
#library(esData)

#rsconnect::setAccountInfo(name='',
#                          token='',
#                          secret='')

# name of the tool on shinyapps: EU_effort_sharing_tool
# temporary -------------------------------------------------------------------
#wkDir <- 'C:/GitHub/python_effort_sharing_tool/r_analysis/'
#setwd(wkDir)
#getwd()
# temp -----------------------------------------------------------------------#

# 1. Load data -----
load('all_shares_with_RES.rda')
load('max_optimals_new.rda')
load('esd_optimals_new.rda')
load('fig4_boxData.rda')
load('popData.rda')

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

# Define the previous reduction targets
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

## Define a list of country ISOs
countryISOs <- c('AUT','BEL','BGR','CYP','CZE','DEU','DNK','ESP','EST','FIN',
                 'FRA','GRC','HRV','HUN','IRL','ITA','LTU','LUX','LVA','MLT',
                 'NLD','POL','PRT','ROU','SVK','SVN','SWE')

# define dictionaries for the dropdowns to correspond to the labels in
# datasets
c_dict <- hash()
c_dict[['C1-EU-capability']] <- 'EU_GDPPOP'
c_dict[['C2-GDP/cap']] <- 'GDPPOP_div'
c_dict[['C3-Governance']] <- 'gee_n'
c_dict[['C4-RES-cap.']] <- 'RES_cap'

e_dict <- hash()
e_dict[['E1-Basic needs']] <- 'MIN_WELF'
#TODO: Check that CPOP_stefan corresponds to ES_EPC and same for CPOP
e_dict[['E2-ES-EPC']] <- 'CPOP'
e_dict[['E3-Full-EPC']] <- 'CPOP_stefan'

r_dict <- hash()
r_dict[['R1-Hist-emi']] <- '1995'
r_dict[['R2-Benefits']] <- 'BENEFITS'
r_dict[['R3-C-budget']] <- 'CBUDGET'
r_dict[['R4-RES-expansion']] <- 'RES'
r_dict[['R5-Budget/cap']] <- 'CBUD_stefan'

iso_dict <- hash()
iso_dict[['Austria']] <- 'AUT'
iso_dict[['Belgium']] <- 'BEL'
iso_dict[['Bulgaria']] <- 'BGR'
iso_dict[['Cyprus']] <- 'CYP'
iso_dict[['Czechia']] <- 'CZE'
iso_dict[['Germany']] <- 'DEU'
iso_dict[['Denmark']] <- 'DNK'
iso_dict[['Spain']] <- 'ESP'
iso_dict[['Estonia']] <- 'EST'
iso_dict[['Finland']] <- 'FIN'
iso_dict[['France']] <- 'FRA'
iso_dict[['Greece']] <- 'GRC'
iso_dict[['Croatia']] <- 'HRV'
iso_dict[['Hungary']] <- 'HUN'
iso_dict[['Ireland']] <- 'IRL'
iso_dict[['Italy']] <- 'ITA'
iso_dict[['Lithuania']] <- 'LTU'
iso_dict[['Luxembourg']] <- 'LUX'
iso_dict[['Latvia']] <- 'LVA'
iso_dict[['Malta']] <- 'MLT'
iso_dict[['The Netherlands']] <- 'NLD'
iso_dict[['Poland']] <- 'POL'
iso_dict[['Portugal']] <- 'PRT'
iso_dict[['Romania']] <- 'ROU'
iso_dict[['Slovakia']] <- 'SVK'
iso_dict[['Slovenia']] <- 'SVN'
iso_dict[['Sweden']] <- 'SWE'

longdict <- hash()
longdict[['R1-Hist-emi']] <- "Reflects the use of fossil fuels since the year 1995, when countries were liable to know the impacts of GHG emissions on climate and had the ability to abate."
longdict[['C1-EU-capability']] <- "This interpretation of GDP per capita corresponds to the 2021 EU policy proposal that puts a cap on the relevance of countries per capita GDP differences for the required emission reductions and is derived empirically from the 2021 proposal"
longdict[['C2-GDP/cap']] <- "This interpretation weighs the relevance of all per capita GDP differences equally in specifying the countries' required emission reductions, increasing or reducing emissions allocations based on the percentage deviation from the EU average GDP per capita"
longdict[['C3-Governance']] <- "By expanding the interpretation of capability beyond macroeconomic measures, this interpretation takes into account additional factors which contribute to the ability of countries to reduce emissions. Similar to the GDP per capita interpretation, here a governance indicator (government effectiveness) derived by the World Bank is used in its place"
longdict[['C4-RES-cap.']] <- "As countries that have in the recent past (since 2005) more strongly developed their RES capacity are more likely to have the human capital, institutional framework, and first-mover advantages to enable ease of further growth in the future, we allocate greater emissions reduction burdens to those countries, and comparatively less reductions are required of countries with less capacity to quickly develop their RES."
longdict[['E1-Basic needs']] <- "In a first step, each member state is pre-allocated the emissions equivalent to energy use (at current emission intensities) required to meet basic needs. After this initial step, the remainder of the budget is distributed to states in an equal-per-capita manner, so that all Member States are assigned at least enough emissions to reach the basic needs threshold"
longdict[['E2-ES-EPC']] <- "Reflects a convergence to equal per capita emissions by 2030 (beginning at today's unequal levels of emissions), based on country emissions in effort-sharing sectors"
longdict[['E3-Full-EPC']] <- "Reflects convergence to equal per capita emissions by 2030 (beginning at today's unequal level of emissions), based on all sectors' emissions"
longdict[['R2-Benefits']] <- "Incorporates the benefits a country has obtained due to emissions prior to the year 1995, interpreted as the embodied emissions in national capital stock"
longdict[['R3-C-budget']] <- "The total emissions budget for the ES sector (calculated by a fictitious linear path from 2020 to 2030) is split among Member States according to population, without any convergence period, thus eliminating any aspect of grandfathering"
longdict[['R4-RES-expansion']] <- "Reflects the difference of the Member States in terms of their change in renewable share from 2005 to 2019 compared to the EU-27 total. Countries with a higher change over that time period compared to the EU average are allocated as a reward a larger emission budget"
longdict[['R5-Budget/cap']] <- "Proposes an alternative method to address historic emissions as compared to R1-Hist-emi, by scaling future emission allowances based on differences in historical cumulative emissions per capita"

# Section for defining explanatory texts!!!####################################
tab1_heading <- "Manual weighting"
tab1_text <- 'Below, you can select equity interpretations and their weights. 
All weights must sum to 100. The tool will make sure of it by automatically 
adjusting the sliders, starting with the bottom-most, if the sum of weights 
is beyond 100.'
# Functions -------------------------------------------------------------------
findTotQual <- function(inputData, capVar, needVar, respVar) {
  return(inputData %>% mutate(totQual = (as.numeric(inputData[[capVar]]) + 
                                           as.numeric(inputData[[needVar]]) + 
                                           as.numeric(inputData[[respVar]]))) %>% 
           subset(totQual < 1.04 & totQual > 0.96))
}

# function to find only the country preferred allocation, given inputs
getOnlyMaxes <- function(inputData, capVar, needVar, respVar) {
  return(findTotQual(inputData, capVar, needVar, respVar) %>% 
           group_by(countries) %>% slice(which.max(value)))
}



# Define drop down menu items
c_drop <- c('C1-EU-capability', 'C2-GDP/cap', 'C3-Governance', 'C4-RES-cap.')
e_drop <- c('E1-Basic needs', 'E2-ES-EPC', 'E3-Full-EPC')
r_drop <- c('R1-Hist-emi', 'R2-Benefits', 'R3-C-budget', 'R4-RES-expansion', 'R5-Budget/cap')

country_drop <- c('Austria', 'Belgium', 'Bulgaria', 'Cyprus', 'Czechia', 
                  'Germany', 'Denmark', 'Spain', 'Estonia', 'Finland', 'France',
                  'Greece', 'Croatia', 'Hungary', 'Ireland', 'Italy', 
                  'Lithuania', 'Luxembourg', 'Latvia', 'Malta', 
                  'The Netherlands', 'Poland', 'Portugal', 'Romania', 'Slovakia',
                  'Slovenia', 'Sweden')

ui <-  navbarPage(theme = shinytheme("simplex"),
                  "Effort sharing in EU member states", fluid = TRUE,
                  tabPanel("EU allocation (manual)",
                           fluidPage(
                             titlePanel('2030 emission reductions'),
                             
                             sidebarLayout(
                               sidebarPanel(
                                 tags$div(class = 'H2', checked = NA,
                                          tags$h2(tab1_heading),
                                          tags$p(tab1_text)),
                                 selectInput('r_var', 'Responsibility interpretation: ',
                                             r_drop),
                                 div(id='my_div', style = 'margin-top:-10px;'),
                                 textOutput("r_desc1"),
                                 div(id='my_div', style = 'margin-top:10px;'),
                                 sliderInput("r_choice", "Percentage",
                                             min = 0, max = 100, value = 30, step = 5),
                                 div(id='my_div',style='margin-top:-10px;'),
                                 selectInput('c_var', 'Capability interpretation: ',
                                             c_drop),
                                 div(id='my_div', style = 'margin-top:-10px;'),
                                 textOutput("c_desc1"),
                                 div(id='my_div', style = 'margin-top:10px;'),
                                 sliderInput("c_choice", "Percentage",
                                             min = 0, max = 100, value = 30, step = 5),
                                 div(id='my_div',style='margin-top:-10px;'),
                                 selectInput('e_var', 'Equality interpretation: ',
                                             e_drop),
                                 div(id='my_div', style = 'margin-top:-10px;'),
                                 textOutput("e_desc1"),
                                 div(id='my_div', style = 'margin-top:10px;'),
                                 sliderInput("e_choice", "Percentage",
                                             min = 0, max = 100, value = 40, step = 5),
                                 div(id='my_div',style='margin-top:-10px;')
                               ),
                               mainPanel(
                                 plotOutput('man_bar_out')
                               )
                             )
                           )
                  ),
                  tabPanel("Optimal-weighting allocation",
                           fluidPage(
                             titlePanel('EU allocations derived from negotiation points'),
                             sidebarLayout(
                               sidebarPanel(
                                 selectInput('r_var2', 'Responsibility interpretation: ',
                                             r_drop),
                                 div(id='my_div', style = 'margin-top:-10px;'),
                                 textOutput("r_desc2"),
                                 div(id='my_div', style = 'margin-top:10px;'),
                                 selectInput('c_var2', 'Capability interpretation: ',
                                             c_drop),
                                 div(id='my_div', style = 'margin-top:-10px;'),
                                 textOutput("c_desc2"),
                                 div(id='my_div', style = 'margin-top:10px;'),
                                 selectInput('e_var2', 'Equality interpretation: ',
                                             e_drop),
                                 div(id='my_div', style = 'margin-top:-10px;'),
                                 textOutput("e_desc2"),
                                 div(id='my_div', style = 'margin-top:10px;'),
                                 plotOutput('opt_tern_out', height = 300)
                               ),
                               mainPanel(
                                 plotOutput('opt_bar_out')
                               )
                             )
                           )
                  ),
                  tabPanel('Individual country ternary charts',
                           fluidPage(
                             titlePanel('Individual country ternary charts'),
                             sidebarLayout(
                               sidebarPanel(
                                 selectInput('r_var3', 'Responsibility interpretation: ',
                                             r_drop),
                                 div(id='my_div', style = 'margin-top:-10px;'),
                                 textOutput("r_desc3"),
                                 div(id='my_div', style = 'margin-top:10px;'),
                                 selectInput('c_var3', 'Capability interpretation: ',
                                             c_drop),
                                 div(id='my_div', style = 'margin-top:-10px;'),
                                 textOutput("c_desc3"),
                                 div(id='my_div', style = 'margin-top:10px;'),
                                 selectInput('e_var3', 'Equality interpretation: ',
                                             e_drop),
                                 div(id='my_div', style = 'margin-top:-10px;'),
                                 textOutput("e_desc3"),
                                 div(id='my_div', style = 'margin-top:10px;'),
                                 selectInput('country', 'Country: ', country_drop)
                               ),
                               mainPanel(
                                 imageOutput("ternary_out")
                               )
                             )
                           ))
)

server <- function(input, output, session) {
  
  observe({
    e_lvl <- 100 - (input$c_choice + input$r_choice) 
    updateSliderInput(session = getDefaultReactiveDomain(), 'e_choice', value = e_lvl)
    
    if(input$c_choice + input$r_choice > 100) {
      updateSliderInput(session = getDefaultReactiveDomain(), 'c_choice', value = (100 - input$r_choice))
    }
    
    if(input$e_choice + input$c_choice > 100) {
      updateSliderInput(session = getDefaultReactiveDomain(), 'r_choice', value = (100 - input$e_choice))
    }
  })
  output$r_desc1 <- renderText(longdict[[input$r_var]])
  
  output$c_desc1 <- renderText(longdict[[input$c_var]])
  
  output$e_desc1 <- renderText(longdict[[input$e_var]])
  
  output$r_desc2 <- renderText(longdict[[input$r_var2]])
  
  output$c_desc2 <- renderText(longdict[[input$c_var2]])
  
  output$e_desc2 <- renderText(longdict[[input$e_var2]])
  
  output$r_desc3 <- renderText(longdict[[input$r_var3]])
  
  output$c_desc3 <- renderText(longdict[[input$c_var3]])
  
  output$e_desc3 <- renderText(longdict[[input$e_var3]])
  
  output$man_bar_out <- renderPlot({
    
    # TEST VALUES ---------------------------------------------------------#
    #input <- NULL
    #input$c_var <- 'C1-EU-capability'
    #input$e_var <- 'E3-Full-EPC'
    #input$r_var <- 'R2-Benefits'
    #input$e_choice <- 0.05
    #input$c_choice <- 0.45
    #input$r_choice <- 0.55
    
    # scale the input weights to sum to 1, rounded to nearest 0.05 FOR NOW
    #if((input$e_choice + input$c_choice + 
    #    input$r_choice) > 100){
    e_scaled = round_any(input$e_choice / (input$e_choice + input$c_choice + 
                                             input$r_choice), 0.05)
    c_scaled = round_any(input$c_choice / (input$e_choice + input$c_choice + 
                                             input$r_choice), 0.05)
    r_scaled = round_any(input$r_choice / (input$e_choice + input$c_choice + 
                                             input$r_choice), 0.05)
    #} else{
    #  e_scaled = input$e_choice / 100
    #  c_scaled = input$c_choice / 100
    #  r_scaled = input$r_choice / 100
    #}
    
    # Identify the scenario (interpretations combination) chosen by the user
    scen_index <- which(names(all_shares_with_RES) == paste(c_dict[[input$c_var]], 
                                                       e_dict[[as.character(input$e_var)]], 
                                                       r_dict[[as.character(input$r_var)]],
                                                       'reduction-_zTrue', sep = '-'))
    
    # get just the chosen combination's results
    temp <- all_shares_with_RES[[scen_index]]
    
    # filter out only the results for the user-specified weighting
    temp <- temp[(temp[[e_dict[[input$e_var]]]] > e_scaled - 0.01 & 
                    temp[[e_dict[[input$e_var]]]] < e_scaled + 0.01), ]
    temp <- temp[(temp[[c_dict[[input$c_var]]]] > c_scaled - 0.01 & 
                    temp[[c_dict[[input$c_var]]]] < c_scaled + 0.01), ]
    vals_to_plot <- temp[(temp[[r_dict[[input$r_var]]]] > r_scaled - 0.01 & 
                            temp[[r_dict[[input$r_var]]]] < r_scaled + 0.01), ]
    
    # plot the results
    plot <- ggplot(vals_to_plot, aes(x = countries, y = value)) + 
      geom_col(fill = '#4575b4') + 
      theme_minimal() + xlab('Country') + 
      ylab('Change compared to 2005 emissions') +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) #+ 
    theme(legend.position = 'none')# + 
    #   scale_color_manual(values = c('#d73027', '#4575b4')) +
    #   scale_fill_manual(values = c('#fdae61', '#abd9e9')) 
    
    plot
  })
  
  output$opt_bar_out <- renderPlot({
    # Identify the scenario (interpretations combination) chosen by the user
    scen_index_max <- which(names(max_optimals_new) == paste(c_dict[[input$c_var2]], 
                                                       e_dict[[as.character(input$e_var2)]], 
                                                       r_dict[[as.character(input$r_var2)]],
                                                       'absolute-_zTrue', sep = '-'))
    scen_index_esd <- which(names(esd_optimals_new) == paste(c_dict[[input$c_var2]], 
                                                          e_dict[[as.character(input$e_var2)]], 
                                                          r_dict[[as.character(input$r_var2)]],
                                                          'absolute-_zTrue', sep = '-'))
    scen_index <- which(names(all_shares_with_RES) == paste(c_dict[[input$c_var2]],
                                                       e_dict[[as.character(input$e_var2)]], 
                                                       r_dict[[as.character(input$r_var2)]],
                                                       'reduction-_zTrue', sep = '-'))
    # get just the chosen combination's results
    max_reduct_scen <- max_optimals_new[[scen_index_max]]
    esd_reduct_scen <- esd_optimals_new[[scen_index_esd]]
    
    # get the country reduction results for maxReults optimals -----------
    # get just the chosen combination's results
    temp <- all_shares_with_RES[[scen_index]]
    
    # the results of max_ and esd are 3 values, for lvls of c, e, r (in order)
    temp <- temp[(temp[[e_dict[[input$e_var2]]]] > max_reduct_scen[2] - 0.01 & 
                    temp[[e_dict[[input$e_var2]]]] < max_reduct_scen[2] + 0.01), ]
    temp <- temp[(temp[[c_dict[[input$c_var2]]]] > max_reduct_scen[1] - 0.01 & 
                    temp[[c_dict[[input$c_var2]]]] < max_reduct_scen[1] + 0.01), ]
    vals_to_plot_max <- temp[(temp[[r_dict[[input$r_var2]]]] > max_reduct_scen[3] - 0.01 & 
                                temp[[r_dict[[input$r_var2]]]] < max_reduct_scen[3] + 0.01), ]
    
    vals_to_plot_max$type <- "Proximity to upper \nbound of equity-\ncompatible emissions"
    
    # get the country reduction results for maxReults optimals -----------
    # get just the chosen combination's results
    temp <- all_shares_with_RES[[scen_index]]
    
    # the results of max_ and esd are 3 values, for lvls of c, e, r (in order)
    temp <- temp[(temp[[e_dict[[input$e_var2]]]] > esd_reduct_scen[2] - 0.01 & 
                    temp[[e_dict[[input$e_var2]]]] < esd_reduct_scen[2] + 0.01), ]
    temp <- temp[(temp[[c_dict[[input$c_var2]]]] > esd_reduct_scen[1] - 0.01 & 
                    temp[[c_dict[[input$c_var2]]]] < esd_reduct_scen[1] + 0.01), ]
    vals_to_plot_esd <- temp[(temp[[r_dict[[input$r_var2]]]] > esd_reduct_scen[3] - 0.01 & 
                                temp[[r_dict[[input$r_var2]]]] < esd_reduct_scen[3] + 0.01), ]
    
    vals_to_plot_esd$type <- 'Minimal change\ncompared to 2018 ESR'
    
    # join the two dataframes
    df <- rbind(vals_to_plot_max, vals_to_plot_esd)
    
    legendText <- 'Possible emission reductions resulting \nfrom weighted combinations which \nminimize differences compared to:'
    
    plot <- ggplot(df, aes(x = countries, y = value, fill = type)) + 
      geom_col(position = 'dodge')  + 
      xlab('Country') + ylab('Change compared to 2005 emissions') +
      theme_minimal() + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      theme(legend.position = 'bottom') + 
      scale_color_manual(values = c('#d73027', '#4575b4')) +
      scale_fill_manual(values = c('#fdae61', '#abd9e9')) +
      labs(color = legendText, 
           fill = legendText)
    plot
  })
  
  output$opt_tern_out <- renderPlot({
    # Identify the scenario (interpretations combination) chosen by the user
    scen_index_max <- which(names(max_optimals_new) == paste(c_dict[[input$c_var2]], 
                                                       e_dict[[as.character(input$e_var2)]], 
                                                       r_dict[[as.character(input$r_var2)]],
                                                       'absolute-_zTrue', sep = '-'))
    scen_index_esd <- which(names(esd_optimals_new) == paste(c_dict[[input$c_var2]], 
                                                          e_dict[[as.character(input$e_var2)]], 
                                                          r_dict[[as.character(input$r_var2)]],
                                                          'absolute-_zTrue', sep = '-'))
    
    # get just the chosen combination's results
    max_reduct_scen <- max_optimals_new[[scen_index_max]]
    esd_reduct_scen <- esd_optimals_new[[scen_index_esd]]
    
    data_points <- list(M = c(max_reduct_scen[2]*100, max_reduct_scen[1]*100,
                              max_reduct_scen[3]*100),
                        E = c(esd_reduct_scen[2]*100, esd_reduct_scen[1]*100,
                              esd_reduct_scen[3]*100))
    
    TernaryPlot(alab = paste(input$c_var2),
                blab = paste(input$e_var2),
                clab = paste(input$r_var2),
                atip = 'C', btip = 'E', ctip = 'R')
    
    AddToTernary(points, data_points, pch = 21, cex = 2.8,
                 bg = c('#f46d43','#74add1'))
  })
  
  # Plotting individual country ternary charts
  output$ternary_out <- renderImage({
    # First isolate just the scenario selected by user
    scen_index <- which(names(all_shares_with_RES) == paste(c_dict[[input$c_var3]],
                                                       e_dict[[input$e_var3]], 
                                                       r_dict[[input$r_var3]],
                                                       'reduction-_zTrue', sep = '-'))
    temp <- all_shares_with_RES[[scen_index]]
    
    # Now select just the country of interest
    country <- iso_dict[[input$country]]
    
    temp <- temp %>% filter(countries == country)
    
    # Now just select the 100% values for each of c, e, and r
    # first get only the fully qualified versions
    temp <- findTotQual(temp, c_dict[[input$c_var3]], e_dict[[input$e_var3]],
                        r_dict[[input$r_var3]])
    cap <- (temp[temp[[c_dict[[input$c_var3]]]] > 1.00, ])$value
    res <- (temp[temp[[r_dict[[input$r_var3]]]] > 1.00, ])$value
    equ <- (temp[temp[[e_dict[[input$e_var3]]]] > 1.00, ])$value
    
    # determine the optimals for max and esd versions
    
    scen_index_max <- which(names(max_optimals_new) == paste(c_dict[[input$c_var3]], 
                                                       e_dict[[input$e_var3]], 
                                                       r_dict[[input$r_var3]],
                                                       'absolute-_zTrue', sep = '-'))
    scen_index_esd <- which(names(esd_optimals_new) == paste(c_dict[[input$c_var3]], 
                                                          e_dict[[input$e_var3]], 
                                                          r_dict[[input$r_var3]],
                                                          'absolute-_zTrue', sep = '-'))
    # get just the chosen combination's results
    max_reduct_scen <- max_optimals_new[[scen_index_max]]
    esd_reduct_scen <- esd_optimals_new[[scen_index_esd]]
    # save as data_points
    data_points <- list(M = c(max_reduct_scen[2]*100, max_reduct_scen[1]*100,
                              max_reduct_scen[3]*100),
                        E = c(esd_reduct_scen[2]*100, esd_reduct_scen[1]*100,
                              esd_reduct_scen[3]*100))
    
    functionToContour <- function(c,n,r) {
      c * cap + n * equ + r * res
    }
    
    # Read ternary_out's width and height. These are reactive values, so this
    # expression will re-run whenever they change.
    width  <- session$clientData$output_ternary_out_width*1.25
    height <- session$clientData$output_ternary_out_height*1.25
    
    # For high-res displays, this will be greater than 1
    pixelratio <- session$clientData$pixelratio
    
    
    outfile <- tempfile(fileext = '.png')
    
    png(outfile, width = width*pixelratio, height = height*pixelratio)
    TernaryPlot(alab = paste('Increasing level of', input$c_var3),
                blab = paste('Increasing level of', input$e_var3),
                clab = paste('Increasing level of', input$r_var3),
                atip = 'C', btip = 'E', ctip = 'R')
    
    values = TernaryPointValues(functionToContour, resolution = 24L)
    ColourTernary(values, spectrum = viridis::cividis(256L))
    TernaryContour(functionToContour, resolution = 36L)
    
    AddToTernary(points, data_points, pch = 21, cex = 2.8,
                 bg = c('#f46d43','#74add1'))
    dev.off()
    
    # Return a list containing the filename
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
  
}

# Run shiny app (locally) ----------------------------------------------------#
shinyApp(ui, server)

