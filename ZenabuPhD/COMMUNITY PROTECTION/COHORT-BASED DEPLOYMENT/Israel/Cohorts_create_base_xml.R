###############################
# Script for creating the base xml for the SSA commodity forecast analysis using
# the latest version of OpenMalariaUtilities available at:
# https://github.com/SwissTPH/r-openMalariaUtilities
#
# Varied components:
# - age structure
# - seasonality
# - vectors contribution to transmission
# - historical interventions (ITN, IRS, CM)
#
# Output: a base .xml file containing @placeholders@ for the varied parameters.
# These placeholders will be populated with country_name-specific values.
#
# 20.09.2022
# monica.golumbeanu@unibas.ch

# FC: This version was adapted for the calibration of Ghana as part of SNT 2023
###############################


#################################
###########################
####additional scripts needed for future bednet decay
#####################3
####################################

library(minpack.lm)
library(ggplot2)
library(tidyverse)

smoothcompact_function = function(time, k=2.14285714285714, L,init_cov=1){
  return(init_cov*(time<L)*exp( k - k/(1 - (time/L)^2) ))
}

fit_ITN_decay = function( distrib_year = 2011, distrib_month=1, distrib_day = 5
                          , coverages
                          , years
                          , months
                          , days=NULL
                          , k=2.14285714285714
                          , nls_control = list(maxiter = 100)) {
  #' fit half-life and initial coverage from ITN coverage values
  #' @param distrib_year year of initial distribution
  #' @param distrib_month month of initial distribution
  #' @param distrib_day day of initial distribution
  #' @param coverages vector of ITN coverage values
  #' @param years vector of years corresponding to coverage values
  #' @param months vector of months corresponding to coverage values
  #' @param days vector of days corresponding to coverage values
  #' @param nls_control maximal number of iterations for fitting
  #' @param k shape parameter, fixed in ITN decay parameterisations
  #' @export
  #' @importFrom minpack.lm nlsLM
  #' @examples fit_ITN_decay(coverages=c(.6,.3,.2),years=rep(2012,3),months=c(1,9,10))
  
  if(length(coverages)!=length(years)) stop("Coverages and years should be of the same length")
  if(length(coverages)!=length(months)) stop("Coverages and months should be of the same length")
  if(length(years)!=length(months)) stop("Years and months should be of the same length")
  
  tab=NULL
  tab$values=coverages
  tab$time=years-distrib_year+(months-distrib_month)/12
  if(!is.null(days)){tab$time=tab$time+(days-distrib_day)/365}
  tab = as.data.frame(tab)
  
  smoothcompact_fit = nlsLM(values ~ a * (c(1) - (time >= L)) * exp((c(1) - (time >= L)) *k* (c(1) - 1/(1 -(time/L)^2)))
                            , data = tab
                            , start = list(a = tab$values[1]
                                           #,k = 1
                                           ,L = 6)
                            , control = nls_control)
  return(list(init_cov=smoothcompact_fit$m$getPars()[1]
              ,L=smoothcompact_fit$m$getPars()[2]))
  
}

view_fitted_ITN_decay = function(distrib_year = 2011, distrib_month=1, distrib_day = 5
                                 , df_cov) {
  #' fit half-life and initial coverage from ITN coverage values
  #' @param distrib_year year of initial distribution
  #' @param distrib_month month of initial distribution
  #' @param distrib_day day of initial distribution
  #' @param df_cov data.frame with observed ITNcov, year and month
  #' @export
  #' @importFrom ggplot2 ggplot
  #' @examples view_fitted_ITN_decay(df_cov=data.frame(ITNcov=c(1,.9,.8),year=c(2011,2012,2013),month=c(1,5,12)))
  
  decay = NULL
  last_year = df_cov %>% filter(year==max(year)) %>% distinct(year) %>% pull()
  last_month = df_cov %>% filter(year==last_year) %>% filter(month==max(month)) %>% distinct(month) %>% pull()
  if(last_year-distrib_year>1){
    decay$month = c(distrib_month:12,rep(1:12,last_year-distrib_year-1),1:last_month)
    decay$year = c(rep(distrib_year,length(distrib_month:12))
                   ,rep((distrib_year+1):(last_year-1),each=12)
                   ,rep(last_year,length(1:last_month)))
  }
  if(last_year-distrib_year==1){
    decay$month = c(distrib_month:12,1:last_month)
    decay$year = c(rep(distrib_year,length(distrib_month:12))
                   ,rep(last_year,length(1:last_month)))
  }
  if(last_year==distrib_year){
    decay$month = distrib_month:last_month
    decay$year = rep(distrib_year,length(last_month-distrib_month))
  }
  
  fit = fit_ITN_decay(distrib_year = distrib_year
                      ,distrib_month = distrib_month
                      ,distrib_day = distrib_day
                      ,coverages = df_cov$ITNcov
                      ,years = df_cov$year
                      ,months = df_cov$month)
  
  decay = as.data.frame(decay)
  decay$time = decay$year-distrib_year+(decay$month-distrib_month)/12
  decay$L = fit$L
  decay$init_cov = fit$init_cov
  
  decay$ITNcov = smoothcompact_function(time=decay$time
                                        ,L=fit$L
                                        ,init_cov = fit$init_cov)
  
  g=ggplot(decay,aes(y = 100 * ITNcov, x = year+(month-1)/12)) +
    geom_line()+
    geom_point(data=df_cov)+
    geom_text(aes(x=(max(year)+min(year))/2,y=max(ITNcov)*100
                  ,label=paste("L=",round(L,2)
                               ,", Initial coverage=",round(init_cov*100,2),"%")))+
    theme_minimal() +
    labs(  x = "Year", y = "Coverage (%)")
  
  return(list(g,decay))
}

# view_fitted_ITN_decay(distrib_month = 4,df_cov=ITN_cov %>% filter(setting=="15.2"))
#
# view_fitted_ITN_decay(df_cov=data.frame(ITNcov=c(1,.9,.8)
#                                                           ,year=c(2011,2012,2013)
#                                                           ,month=c(1,5,12)))



#' Realistic net attrition for the future
#' Chemical barrier and survivorship: Comparative study of two brands of
#' polyester nets and one brand of polyethylene nets in different conditions of
#' used in Benin
#' Idelphonse B Ahogni, Rock Y Aïkpon, Jean-Fortuné Dagnon, Roseric Azondekon,
#'  Bruno Akinro, Germain G Padonou and Martin C Akogbeto
#'  https://www.dipterajournal.com/pdf/2020/vol7issue5/PartA/7-4-18-884.pdf

# source("/scicore/home/pothin/pothin/Gitrepo/examples_om/fit_ITN_initial_coverage_halflife.R")


distrib_year=2019
distrib_month=1
#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' Abtract, Findings

#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' Appendix page 4, explicit dates in main paper page 1230
#' "Proportion of participants reporting using a net the night before"
ITN_cov=data.frame(setting=rep(c("pyrethroid","pyriproxyfen","chlorfenapyr","PBO"),each=4),
                   year=rep(c(2019,2020,2020,2021),4),
                   month=rep(c(4,1,8,1),4),
                   ITNcov=c(868/1130,2848/4627,2606/4991,2488/5029,
                            764/1103,2646/4358,2124/4645,2087/5455,
                            752/1099,3157/4833,2676/5194,2585/5576,
                            771/1046,2506/4247,1885/4633,1534/5186)
)

ITNdecay_pyrethroid = fit_ITN_decay(distrib_year = distrib_year
                                    ,distrib_month = distrib_month
                                    ,coverages=ITN_cov %>% filter(setting=="pyrethroid") %>% pull(ITNcov)
                                    ,years=ITN_cov %>% filter(setting=="pyrethroid") %>% pull(year)
                                    ,months=ITN_cov %>% filter(setting=="pyrethroid") %>% pull(month))

ITNdecay_pyriproxyfen = fit_ITN_decay(distrib_year = distrib_year
                                      ,distrib_month = distrib_month
                                      ,coverages=ITN_cov %>% filter(setting=="pyriproxyfen") %>% pull(ITNcov)
                                      ,years=ITN_cov %>% filter(setting=="pyriproxyfen") %>% pull(year)
                                      ,months=ITN_cov %>% filter(setting=="pyriproxyfen") %>% pull(month))

ITNdecay_chlorfenapyr = fit_ITN_decay(distrib_year = distrib_year
                                      ,distrib_month = distrib_month
                                      ,coverages=ITN_cov %>% filter(setting=="chlorfenapyr") %>% pull(ITNcov)
                                      ,years=ITN_cov %>% filter(setting=="chlorfenapyr") %>% pull(year)
                                      ,months=ITN_cov %>% filter(setting=="chlorfenapyr") %>% pull(month))

ITNdecay_PBO= fit_ITN_decay(distrib_year = distrib_year
                            ,distrib_month = distrib_month
                            ,coverages=ITN_cov %>% filter(setting=="PBO") %>% pull(ITNcov)
                            ,years=ITN_cov %>% filter(setting=="PBO") %>% pull(year)
                            ,months=ITN_cov %>% filter(setting=="PBO") %>% pull(month))

view_fitted_ITN_decay(distrib_year=distrib_year
                      ,distrib_month = distrib_month
                      ,df_cov=ITN_cov %>% filter(setting=="pyrethroid"))
view_fitted_ITN_decay(distrib_year=distrib_year
                      ,distrib_month = distrib_month
                      ,df_cov=ITN_cov %>% filter(setting=="pyriproxyfen"))
view_fitted_ITN_decay(distrib_year=distrib_year
                      ,distrib_month = distrib_month
                      ,df_cov=ITN_cov %>% filter(setting=="chlorfenapyr"))
view_fitted_ITN_decay(distrib_year=distrib_year
                      ,distrib_month = distrib_month
                      ,df_cov=ITN_cov %>% filter(setting=="PBO"))



# To install OpenMalariaUtilities:
# devtools::install_github("SwissTPH/r-openMalariaUtilities", force = TRUE)

# Load the necessary packages
library(devtools)
library(openMalariaUtilities)
library(tidyverse)
library(omucompat)

# Function that creates a list with all the elements which are specific for the
# base xml (placeholders, interventions, etc.)
create_baseList = function(country_name="GHA", 
                           sim_start="1918-01-01", 
                           versionnum=44L) {
  
  # Load the script that computes attrition for different LLIN types
  #source("ITN_netAttrition.R")
  
  ## Basic xml skeleton
  baseList = list(
    # Mandatory
    expName = country_name,
    # Mandatory
    OMVersion = versionnum,
    # Mandatory
    demography = list(),
    monitoring = list(),
    interventions = list(),
    healthSystem = list(),
    entomology = list(),
    # These are optional for OM
    # parasiteGenetics = list(),
    # pharmacology = list(),
    # diagnostics = list(),
    model = list()
  )
  
  baseList = defineDemography(
    baseList,
    name = country_name,
    popSize = "@pop@",
    maximumAgeYrs = 94,
    lowerbound = 0,
    poppercent = GHA$poppercent,
    upperbound = GHA$upperbound
  )
  
  ## Create monitoring snippet
  baseList[["monitoring"]] = list(
    name = "Surveys",
    ## Mandatory, different from OM schema
    startDate = sim_start,
    continuous = monitoringContinuousGen(period = 1,
                                         options = list(
                                           name = c("input EIR", "simulated EIR", "human infectiousness", "N_v0",
                                                    "immunity h", "immunity Y", "new infections",
                                                    "num transmitting humans", "ITN coverage", "GVI coverage", "alpha",
                                                    "P_B", "P_C*P_D"),
                                           value = c("true", "true", "true", "true", "true", "true", "true", "true",
                                                     "true", "true", "true", "false", "false")
                                         )
    ),
    SurveyOptions = monitoringSurveyOptionsGen(
      options = list(
        name = c("nHost", "nPatent", "nUncomp", "nSevere", "nDirDeaths",
                 "inputEIR", "simulatedEIR","nTreatments1","nTreatments2","nTreatments3","expectedDirectDeaths"),
        value = c("true", "true", "true", "true", "true", 
                  "true", "true", "true", "true", "true", "true")
      )
    ),
    # We are setting a survey on the 5th day of each month
    surveys = monitoringSurveyTimesGen(detectionLimit = 100, startDate = "2015-01-01", #sim_start changed from 1999 to 2015
                                       endDate = "2035-01-01",
                                       interval = list(days = c(5), months = c(1:12), years = c(2015:2035)),
                                       simStart = sim_start),
    ## surveyAgeGroupsGen will write thirdDimension table to cache, important for postprocessing
    ageGroup = surveyAgeGroupsGen(lowerbound = 0, upperbounds = c(1, 2, 5, 10, 100)),
    # add cohorts in monitoring
    cohorts = list(subPop = list(id="LLINusers", number="1"))
  )
  
  # Mosquito species and contribution to transmission
  # We generically define the contribution of the three dominant species
  # and assign 0 for the settings where they do not occur
  contrib = c("@conteir_arabiensis_indoor@", "@conteir_arabiensis_outdoor@",
              "@conteir_gambiae_indoor@", "@conteir_gambiae_outdoor@",
              "@conteir_funestus_indoor@", "@conteir_funestus_outdoor@")
  mosqs = c("arabiensis_indoor", "arabiensis_outdoor",
            "gambiae_indoor", "gambiae_outdoor",
            "funestus_indoor", "funestus_outdoor")
  
  
  ### ALL HUMAN INTERVENTIONS NEED TO BE REVIEWED
  # Begin interventions for humans
  ########################################################################################################
  ## Definition section
  ########################################################################################################
  
  
  # # # Define IRS snippets
  baseList = define_IRS_compat(baseList = baseList, mosqs = mosqs,
                               component = "Actellic300CS")
  
  # # Define historical ITN and PBO
  baseList = define_ITN_compat(baseList = baseList, component = "histITNresist2015",
                               mosqs = mosqs, hist = TRUE, resist = TRUE,
                               versionnum = versionnum)
  baseList = define_ITN_compat(baseList = baseList, component = "histITNpre2015",
                               mosqs = mosqs, hist = TRUE, resist = FALSE,
                               versionnum = versionnum)
  baseList = define_ITN_compat(baseList = baseList, component = "histITNresistweak2015",
                               mosqs = mosqs, hist = TRUE, resist = TRUE,
                               versionnum = versionnum)
  
  #comment halflife as a variable out
  # baseList = define_ITN_compat(baseList = baseList, component = "PBO",
  #                              mosqs = mosqs, hist = TRUE,
  #                              resist = FALSE, halflife = "@ITNhalflife@", versionnum = versionnum)
  # # Define ITN and PBO starting in 2020
  # baseList = define_ITN_compat(baseList = baseList, component = "futITN",
  #                              mosqs = mosqs, hist = FALSE, resist = TRUE, halflife = "@ITNhalflife@",
  #                              versionnum = versionnum)
  # baseList = define_ITN_compat(baseList = baseList, component = "futITNweak",
  #                              mosqs = mosqs, hist = FALSE, resist = TRUE, halflife = "@ITNhalflife@",
  #                              versionnum = versionnum)
  # baseList = define_ITN_compat(baseList = baseList, component = "futPBO",
  #                              mosqs = mosqs, hist = FALSE,
  #                              resist = FALSE, halflife = "@ITNhalflife@", versionnum = versionnum)
  # baseList = define_ITN_compat(baseList = baseList, component = "futPBO3",
  #                              mosqs = mosqs, hist = FALSE,
  #                              resist = FALSE, halflife = "@ITNhalflife@", versionnum = versionnum)
  # baseList = define_ITN_compat(baseList = baseList, component = "futIG2",
  #                              mosqs = mosqs, hist = FALSE,
  #                              resist = FALSE, halflife = "@ITNhalflife@", versionnum = versionnum)
  
  baseList = define_ITN_compat(baseList = baseList, component = "PBO",
                               mosqs = mosqs, hist = TRUE,
                               resist = FALSE, halflife = 2, versionnum = versionnum)
  # Define ITN and PBO starting in 2020
  baseList = define_ITN_compat(baseList = baseList, component = "futITN",
                               mosqs = mosqs, hist = FALSE, resist = TRUE, halflife = ITNdecay_pyrethroid$L/2,
                               versionnum = versionnum)
  baseList = define_ITN_compat(baseList = baseList, component = "futITNweak",
                               mosqs = mosqs, hist = FALSE, resist = TRUE, halflife = ITNdecay_pyrethroid$L/2,
                               versionnum = versionnum)
  baseList = define_ITN_compat(baseList = baseList, component = "futPBO",
                               mosqs = mosqs, hist = FALSE,
                               resist = FALSE, halflife = ITNdecay_PBO$L/2, versionnum = versionnum)
  baseList = define_ITN_compat(baseList = baseList, component = "futPBO3",
                               mosqs = mosqs, hist = FALSE,
                               resist = FALSE, halflife = 3, versionnum = versionnum)
  baseList = define_ITN_compat(baseList = baseList, component = "futIG2",
                               mosqs = mosqs, hist = FALSE,
                               resist = FALSE, halflife = ITNdecay_chlorfenapyr$L/2, versionnum = versionnum)
  
  
  
  ### To make it a weak resistance, we set postprandialKillingEffect and preprandialKillingEffect to zero ###
  #Adapted from Jeanne's code
  for (net_name in c("histITNresistweak2015", "futITNweak")){
    indexNet = as.numeric(which(unlist(lapply(baseList$interventions$human, function(x) x[["id"]] == net_name))))
    
    lenFields = length(baseList$interventions$human[indexNet]$component$ITN$anophelesParams$postprandialKillingEffect)
    #assign 0 values for all tags inside pre and post
    #to act on the 6 mosquitoes, index 8:13 (CAREFUL WITH THIS PART IN THE FUTURE!!!) 
    indexMosq=as.numeric(extractList(baseList$interventions$human[indexNet]$component$ITN,name="anophelesParams",onlyIndex=T))
    for (j in indexMosq){
      
      for (i in 1:lenFields){
        baseList$interventions$human[indexNet]$component$ITN[j]$anophelesParams$postprandialKillingEffect[[i]] <- 0
        baseList$interventions$human[indexNet]$component$ITN[j]$anophelesParams$preprandialKillingEffect[[i]] <- 0
      }
      
    }
    
  }
  
  
  ## Define empty ITN/IRS   #### This is to have something when coverage is zero??
  baseList <- define_nothing_compat(
    baseList = baseList, component = "nothing", mosqs = mosqs)
  
  ## Define SMC
  baseList <- define_treatSimple(baseList, component = "SMC",
                                 durationBlood = "30d")
  
  ## Define PMC
  baseList <- define_treatSimple(baseList, component = "PMC",
                                 durationBlood = "10d")
  
  
  ## Define RTS,S
  #coming from work from Melissa Penny
  vaccineParameterization=list(RTSSVaccine=list(PEV=list(decay=list(
    L="223d",`function`="weibull", k="0.84"),
    efficacyB=list(value="10.0"),
    initialEfficacy=list(value="0.91"))
  ))
  
  
  #this defaults to a component id= "RTSSVaccine", same as in vaccineParametrization
  baseList <- defineVaccine(baseList, vaccineParameterization, verbatim=T)
  
  
  
  ########################################################################################################
  ## Deployment section
  ########################################################################################################
  
  
  #######################
  ## Historical
  #######################
  
  # Historical IRS 
  baseList = deploy_it_compat(baseList = baseList, component = "Actellic300CS",
                              coverage = "@histIRScov@", byyear = T, y1 = 2012, 
                              y2 = 2022, every = 1, interval = "year",
                              m1 = 4, m2 = 4, d1 = 5, d2 = 5,
                              SIMSTART = sim_start)
  
  
  # Historical ITN campaign, we assume no resistance before 2015
  baseList = deploy_it_compat(baseList = baseList, component = "histITNpre2015", 
                              coverage = "@histITNcov@",byyear = TRUE,
                              y1 = 2000, y2 = 2014, every = 1, interval = "year",
                              m1 = 4, m2 = 4, d1 = 5, d2 = 5, SIMSTART = sim_start
  )
  baseList = deploy_it_compat(baseList = baseList, component = "@histITNtype2015@", 
                              coverage = "@histITNcov@",byyear = TRUE,
                              y1 = 2015, y2 = 2021, every = 1, interval = "year",
                              m1 = 4, m2 = 4, d1 = 5, d2 = 5, SIMSTART = sim_start
  )
  
  
  ## Past SMC
  baseList <- deploy_it_compat(
    baseList = baseList, component = "SMC",
    coverage = "@histSMCcov@",
    byyear = TRUE,
    ## minAge, maxAge
    minAge = 0.25, maxAge = "@maxSMCage@",
    y1 = 2012, y2 = 2022, every = 1, interval = "month",
    m1 = 6, m2 = 9, d1 = 5, d2 = 5, SIMSTART = SIMSTART
  )
  
  
  ###Past LSM
  DatesLarviciding_unimodal_past<- expand.grid(2021:2022,paste(str_pad(c(5:6,9:10),2,pad="0"),"01",sep="-"))
  DatesLarviciding_unimodal_past<-sort(paste(DatesLarviciding_unimodal_past$Var1,DatesLarviciding_unimodal_past$Var2,sep="-"))
  
  
  baseList <- defineLarv(baseList,    mosquitos=mosqs,
                         coverage = "@histLSMcov_unimodal@",
                         dates=DatesLarviciding_unimodal_past
  )
  
  
  
  ############################
  ## Future  (TO BE REVIEWED)
  ############################
  
  ############## Future IRS ############
  baseList <- deployIT(
    baseList = baseList, component = "Actellic300CS",
    coverage = "@futIRScov@",
    interval="1 year",
    startDate="2023-01-01",
    endDate="2030-12-31"
  )
  
  ############## Future nets campaign ############
  
  ##starting 2021
  baseList <- deployIT(
    baseList = baseList, component = "@futITNtypecampaign1@", coverage = "@futNetcovstart2021@",
    dates="2021-05-01"
  )
  
  baseList <- deployIT(
    baseList = baseList, component = "@futITNtypecampaign1@", coverage = "@futNetcovstart2021@",
    dates="2024-05-01"
  )
  
  baseList <- deployIT(
    baseList = baseList, component = "@futITNtypecampaign1@", coverage = "@futNetcovstart2021@",
    dates="2027-05-01"
  )
  
  
  # ###starting 2022
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypecampaign1@", coverage = "@futNetcovstart2022@",
  #   dates="2022-05-01"
  # )
  # 
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypecampaign2@", coverage = "@futNetcovstart2022@",
  #   dates="2025-05-01"
  # )
  # 
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypecampaign3@", coverage = "@futNetcovstart2022@",
  #   dates="2028-05-01"
  # )
  # 
  # 
  # ###starting 2023
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypecampaign1@", coverage = "@futNetcovstart2023@",
  #   dates="2023-05-01"
  # )
  # 
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypecampaign2@", coverage = "@futNetcovstart2023@",
  #   dates="2026-05-01"
  # )
  # 
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypecampaign3@", coverage = "@futNetcovstart2023@",
  #   dates="2029-05-01"
  # )
  # 
  # 
  # ###starting 2024
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypecampaign1@", coverage = "@futNetcovstart2024@",
  #   dates="2024-05-01"
  # )
  # 
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypecampaign2@", coverage = "@futNetcovstart2024@",
  #   dates="2027-05-01"
  # )
  # 
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypecampaign3@", coverage = "@futNetcovstart2024@",
  #   dates="2030-05-01"
  # )
  
  
  
  ############## Future nets continuous ############
  # baseList = deploy_it_compat(baseList = baseList, component = "@futITNtypeContinuous@",
  #                             cumulative = TRUE, coverage = "@futITNContinuouscov@", byyear = TRUE,
  #                             y1 = 2023, y2 = 2030, every = 1, interval = "year",
  #                             SIMSTART = sim_start)
  # 
  # baseList <- deployIT(
  #   baseList = baseList, component = "@futITNtypeContinuous@", coverage = "@futITNContinuouscov@",
  #   interval="1 year",
  #   cumulative = TRUE,
  #   startDate="2021-01-01",
  #   endDate="2030-12-31"
  # )
  
  
  
  ############## Future SMC ############
  # ##2021
  # ## July-October 2021
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 5
  #                      , startDate = "2021-07-15", endDate = "2021-10-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to5JulOct@")
  # 
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 10
  #                      , startDate = "2021-07-15", endDate = "2021-10-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to10JulOct@")
  # 
  # 
  # ## June-September 2021
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 5
  #                      , startDate = "2021-06-15", endDate = "2021-09-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to5JunSept@")
  # 
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 10
  #                      , startDate = "2021-06-15", endDate = "2021-09-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to10JunSept@")
  # 
  # ## June-October 2021
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 5
  #                      , startDate = "2021-06-15", endDate = "2021-10-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to5JunOct@")
  # 
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 10
  #                      , startDate = "2021-06-15", endDate = "2021-10-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to10JunOct@")
  # 
  # 
  # ##2022
  # ## July-October 2022
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 5
  #                      , startDate = "2022-07-15", endDate = "2022-10-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to5JulOct@")
  # 
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 10
  #                      , startDate = "2022-07-15", endDate = "2022-10-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to10JulOct@")
  # 
  # 
  # ## June-September 2022
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 5
  #                      , startDate = "2022-06-15", endDate = "2022-09-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to5JunSept@")
  # 
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 10
  #                      , startDate = "2022-06-15", endDate = "2022-09-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to10JunSept@")
  # 
  # ## June-October 2022
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 5
  #                      , startDate = "2022-06-15", endDate = "2022-10-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to5JunOct@")
  # 
  # baseList <- deployIT(baseList, component = "SMC"
  #                      , minAge = 0.25
  #                      , maxAge = 10
  #                      , startDate = "2022-06-15", endDate = "2022-10-30"
  #                      , interval = "1 month"
  #                      , coverage = "@futcovSMC0to10JunOct@")
  
  ###2023
  ## July-October 2023
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2023-07-15", endDate = "2023-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JulOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2023-07-15", endDate = "2023-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JulOct@")
  
  
  ## June-September 2023
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2023-06-15", endDate = "2023-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunSept@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2023-06-15", endDate = "2023-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunSept@")
  
  ## June-October 2023
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2023-06-15", endDate = "2023-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2023-06-15", endDate = "2023-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunOct@")
  
  ## 2024
  
  ## July-October 2024
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2024-07-15", endDate = "2024-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JulOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2024-07-15", endDate = "2024-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JulOct@")
  
  ## June-September 2024
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2024-06-15", endDate = "2024-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunSept@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2024-06-15", endDate = "2024-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunSept@")
  
  ## June-October 2024
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2024-06-15", endDate = "2024-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2024-06-15", endDate = "2024-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunOct@")
  
  
  ### 2025
  
  ## July-October 2025
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2025-07-15", endDate = "2025-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JulOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2025-07-15", endDate = "2025-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JulOct@")
  
  
  ## June-September 2025
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2025-06-15", endDate = "2025-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunSept@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2025-06-15", endDate = "2025-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunSept@")
  
  ## June-October 2025
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2025-06-15", endDate = "2025-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2025-06-15", endDate = "2025-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunOct@")
  
  
  ### 2026
  
  ## July-October 2026
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2026-07-15", endDate = "2026-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JulOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2026-07-15", endDate = "2026-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JulOct@")
  
  ## June-September 2026
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2026-06-15", endDate = "2026-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunSept@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2026-06-15", endDate = "2026-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunSept@")
  
  ## June-October 2026
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2026-06-15", endDate = "2026-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2026-06-15", endDate = "2026-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunOct@")
  
  
  ### 2027
  
  ## July-October 2027
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2027-08-15", endDate = "2027-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JulOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2027-07-15", endDate = "2027-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JulOct@")
  
  ## June-September 2027
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2027-06-15", endDate = "2027-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunSept@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2027-06-15", endDate = "2027-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunSept@")
  
  ## June-October 2027
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2027-06-15", endDate = "2027-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2027-06-15", endDate = "2027-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunOct@")
  
  
  ### 2028
  
  ## July-October 2028
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2028-07-15", endDate = "2028-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JulOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2028-07-15", endDate = "2028-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JulOct@")
  
  ## June-September 2028
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2028-06-15", endDate = "2028-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunSept@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2028-06-15", endDate = "2028-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunSept@")
  
  ## June-October 2028
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2028-06-15", endDate = "2028-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2028-06-15", endDate = "2028-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunOct@")
  
  
  ### 2029
  
  ## July-October 2029
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2029-07-15", endDate = "2029-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JulOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2029-07-15", endDate = "2029-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JulOct@")
  
  ## June-September 2029
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2029-06-15", endDate = "2029-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunSept@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2029-06-15", endDate = "2029-09-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunSept@")
  
  ## June-October 2029
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 5
                       , startDate = "2029-06-15", endDate = "2029-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to5JunOct@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 10
                       , startDate = "2029-06-15", endDate = "2029-10-30"
                       , interval = "1 month"
                       , coverage = "@futcovSMC0to10JunOct@")
  
  
  
  ############## Future PMC ############
  ## Future PMC: 10 and 14 weeks, 6, 9, 12, 15, 18 and 24 months
  
  ### Communes starting in January 2024
  baseList <- deploy_cont_compat(baseList, component = "PMC"
                                 , begin = "2024-03-01"
                                 , end = "2030-12-31"
                                 , targetAgeYrs = c(10/52, 14/52, 0.5, 0.75,
                                                    1 ,1.25, 1.5, 2)
                                 , coverage = rep("@futPMCMar2024age1cov@",8))
  ### Communes starting in January 2023
  baseList <- deploy_cont_compat(baseList, component = "PMC"
                                 , begin = "2023-03-01"
                                 , end = "2030-12-31"
                                 , targetAgeYrs = c(10/52, 14/52, 0.5, 0.75,
                                                    1, 1.25, 1.5, 2)
                                 , coverage = rep("@futPMCMar2023age1cov@",8))
  
  
  
  
  ############## Future IPTsc ############
  ## Future IPTsc: age group 5-16 
  
  #2024-2030
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 5
                       , maxAge = 16
                       , dates = c("2024-03-01","2024-07-01","2024-11-01",
                                   "2025-03-01","2025-07-01","2025-11-01",
                                   "2026-03-01","2026-07-01","2026-11-01",
                                   "2027-03-01","2027-07-01","2027-11-01",
                                   "2028-03-01","2028-07-01","2028-11-01",
                                   "2029-03-01","2029-07-01","2029-11-01",
                                   "2030-03-01","2030-07-01","2030-11-01")
                       , coverage = "@futcovIPTsc_uni@")
  
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 5
                       , maxAge = 16
                       , dates = c("2024-03-01","2024-07-01","2024-11-01",
                                   "2025-03-01","2025-07-01","2025-11-01",
                                   "2026-03-01","2026-07-01","2026-11-01",
                                   "2027-03-01","2027-07-01","2027-11-01",
                                   "2028-03-01","2028-07-01","2028-11-01",
                                   "2029-03-01","2029-07-01","2029-11-01",
                                   "2030-03-01","2030-07-01","2030-11-01")
                       , coverage = "@futcovIPTsc_bim@")
  
  ### IPTi, validated by Branwen
  baseList <- deploy_cont_compat(baseList, component = "PMC"
                                 , begin = "2023-01-01"
                                 , end = "2030-12-31"
                                 , targetAgeYrs = c(6/52, 10/52,14/52,0.75,
                                                    1.5)
                                 , coverage = rep("@futIPTicov@",5))
  
  
  
  
  # ############## Future RTSS ############
  # 
  #3rd dose at 9months, 9/12 = 0.75
  baseList <- deploy_cont_compat(baseList, component = "RTSSVaccine"
                                 , begin = "2021-01-01"
                                 , end = "2022-12-31"
                                 #, targetAgeYrs = c(0.75, 2.25)
                                 , targetAgeYrs = c(0.75, 1.5)
                                 , coverage = c("@futRTSS3covStart2021@", "@futRTSS4covStart2021@")
                                 , vaccMaxCumDoses = c(1,2)
                                 , vaccMinPrevDoses = c(0,1))
  
  baseList <- deploy_cont_compat(baseList, component = "RTSSVaccine"
                                 , begin = "2023-01-01"
                                 , end = "2030-12-31"
                                 #, targetAgeYrs = c(0.75, 2.25)
                                 , targetAgeYrs = c(0.75, 1.5)
                                 , coverage = c("@futRTSS3covStart2023@", "@futRTSS4covStart2023@")
                                 , vaccMaxCumDoses = c(1,2)
                                 , vaccMinPrevDoses = c(0,1))
  
  
  
  
  #####Larviciding
  
  ### Hard code the dates for LSM
  DatesLarviciding_unimodal<- expand.grid(2023:2030,paste(str_pad(c(5:6,9:10),2,pad="0"),"01",sep="-"))
  DatesLarviciding_unimodal<-sort(paste(DatesLarviciding_unimodal$Var1,DatesLarviciding_unimodal$Var2,sep="-"))
  
  # DatesLarviciding_bimodal<- expand.grid(2023:2030,paste(str_pad(c(2:3,6:7,9:10),2,pad="0"),"01",sep="-"))
  # DatesLarviciding_bimodal<-sort(paste(DatesLarviciding_bimodal$Var1,DatesLarviciding_bimodal$Var2,sep="-"))
  #   
  baseList <- defineLarv(baseList,    mosquitos=mosqs,
                         coverage = "@futLSMcov_unimodal@",
                         dates=DatesLarviciding_unimodal
  )
  
  # baseList <- defineLarv(baseList,    mosquitos=mosqs,
  #                        coverage = "@futLSMcov_bimodal@",
  #                        dates=DatesLarviciding_bimodal
  # )
  # 
  
  ###future MDA
  baseList <- deployIT(baseList, component = "SMC"
                       , minAge = 0.25
                       , maxAge = 100
                       , startDate = "2023-01-01", endDate = "2030-01-01"
                       , interval = "1 month"
                       , coverage = "@futMDAcov@")
  
  
  ############################
  ## End of Future deployment
  ############################
  
  
  
  #### HUMAN INTERVENTIONS STOP HERE
  
  
  # Importation
  baseList = define_importedInfections_compat(baseList = baseList, 10, time = 0)
  
  # Health system changes (we set the coverage of access to treatment, the effective coverage, each year) for past
  baseList = define_changeHS_compat(baseList = baseList, access = "histAccess",
                                    y1 = 2003, y2 = 2020,
                                    use_at_symbol = TRUE,
                                    ## Default values
                                    pSelfTreatUncomplicated = 0.01821375,
                                    pSeekOfficialCareSevere = .48,
                                    SIMSTART = sim_start)
  # Health system changes (we set the coverage of access to treatment, the effective coverage, each year) for futre
  baseList = define_changeHS_compat(baseList = baseList, access = "futCMcov",
                                    y1=2023,y2=2023,
                                    use_at_symbol = TRUE,
                                    ## Default values
                                    pSelfTreatUncomplicated = 0.01821375,
                                    pSeekOfficialCareSevere = .48,
                                    SIMSTART = sim_start)
  
  
  # Write a "dummy" health system (the simulation will use the previously defined health system changes)
  baseList = write_healthsys_compat(baseList = baseList, access = 0)
  
  ## Entomology section
  baseList = make_ento_compat(baseList = baseList, mosqs, contrib, EIR = "@EIR@",
                              seasonality = paste0("@m", 1:12, "@"), ## NEED TO MATCH with
                              propInfected = .078, propInfectious = .021)
  
  ## Specify seed and finish XML file
  baseList = write_end_compat(baseList = baseList,
                              seed = "@seed@", modelname = "base")
  
  return(baseList)
}

