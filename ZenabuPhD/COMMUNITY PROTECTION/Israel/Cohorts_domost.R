#################################
# Main initialization script of the workflow, creates for each country
# an analysis directory with all the files and scripts
# consisting of the following steps:
# 1. Defining the base xml
# 1. Defining the simulation settings (scenarios)
# 2. Creating the scripts for generating the simulation scenarios
# 3. Creating the scripts for running OpenMalaria simulations
# 4. Creating the scripts for postprocessing
#
#################################

# History cleanup
rm(list=ls())

# Load the necessary packages
library(devtools)
library(openMalariaUtilities)
library(OMAddons)
library(omuslurm)
library(omucompat)
library(dplyr)

#####################################
# Initialization
#####################################
git_repo_path="/scicore/home/pothin/suboi0000/GitRepositories/om-ghana"
#git_repo_path="/scicore/home/pothin/suboi0000/GitRepositories/om-ghana"

setwd(file.path(git_repo_path,"emulator_training/OMsimulations/nsp2023/COMMUNITY_PROTECTION/COHORTS/"))

# Load the base xml list setup function and auxiliary functions
# If the script is executed outside of RStudio, the full path to this file needs to be provided:
source(file.path(git_repo_path,"emulator_training/OMsimulations/nsp2023/COMMUNITY_PROTECTION/COHORTS/00_create_base_xml_CP_cohort.R"))

# Define root directory with all the experiments according to the user
#root_dir_path = "/scicore/home/pothin/GROUP/ZenabuPhD/OMOutputs/COMMUNITY_PROTECTION/COHORTS"
root_dir_path = "/scicore/home/pothin/GROUP/ZenabuPhD/OMOutputs/COMMUNITY_PROTECTION/COHORTS/LLIN_Only/"
countryDat_file = "/scicore/home/pothin/GROUP/ZenabuPhD/OMinput/final_countrydat_ghana.csv"
calibration_file ="/scicore/home/pothin/GROUP/ZenabuPhD/OMinput/calibrated_ghana_2015.csv"
scenario_file="/scicore/home/pothin/GROUP/ZenabuPhD/OMinput/future_scenarios_vnew.csv"

# Load the scenario specifications (countryDat file) and the country data
countryDat = read.csv(countryDat_file,sep=";")
# Adjust column names
names(countryDat) <- gsub(pattern = "^ITN", replacement = "histITNcov", x = names(countryDat))
names(countryDat) <- gsub(pattern = "^EffCov14d", replacement = "histAccess", x = names(countryDat))
names(countryDat) <- gsub(pattern = "^X([0-9]{4})_IRS", replacement = "histIRScov\\1", x = names(countryDat))
names(countryDat) <- gsub(pattern = "^X([0-9]{4})_SMC", replacement = "histSMCcov\\1", x = names(countryDat))

calibrationEIR = read.csv(calibration_file)

iso_code <- "GHA"

#-------------------------------
#####################################
# Access conversion
#####################################
# Convert values
## -- need to convert Access to 5-day time steps ( current access to care, 85% eff. treat rate)
# katya = FALSE, scale 1
countryDat%>%
  mutate(across(starts_with("histAccess"),function(x)x/100))%>%
  convert_access(countryDat, pattern = "histAccess", katya = FALSE, scale = 1)->countryDat

## according to NMCP MAP values are too low
countryDat$histITNcov2020<-0.45

#####################################
# Defining the simulation scenarios
#####################################
# Define the list with all the scenario variations per country
# (in this example: population size, seed, and EIR vary)
full = list()
full$pop = 10000L
full$seed = c(1:5)
full$setting = unique(calibrationEIR$setting)
full$maxSMCage <- 5
full$futcovSMC0to5JulOct <- 0 #seq(0,1,0.2)
full$futNetcovstart2021 <-  c(seq(0,0.6,0.2), 0.7, 0.8, 0.82, 0.85, 0.87, 0.9) # c(0, 0.5)
full$histITNtype2015 <- "histITNresistweak2015"
# full$futIRScov <- seq(0,1,0.2)
# full$futRTSS3covStart2023 <- seq(0,1,0.2)
# FC: EIRs are needed for calibration, not for forward simulation where they'll be an input part of countryDat
full$EIR = c(8, 16, 32, 48, 60, 80, 100)

#### 'scens1' will contain all possible combinations of these scenario variations
scens = expand.grid(full)
scens = left_join(scens, countryDat, by="setting")

scens = scens %>%
  mutate(across(c("setting"),tolower))%>%
  left_join(calibrationEIR%>%
              select(-c("EIR_lci","EIR_uci","EIR"))%>%
              # pivot_longer(cols = c("EIR_lci","EIR_uci","EIR"))%>%
              select(setting,sub)%>% #,name,value
              # rename(EIR="value",EIR_CI="name")%>%
              mutate(sub=gsub(" ","_",tolower(sub)))
            ,
            by="setting")#%>%
# mutate(EIR=ifelse(EIR_CI=="EIR_lci"&is.na(EIR),1,EIR))%>%
# mutate(EIR=ifelse(EIR_CI=="EIR_uci"&is.na(EIR),max(EIR,na.rm=T),EIR))

future_scenarios<-read_csv2(scenario_file)

##before you join, check if the sub names are the same, there is a correspondence table to fix this...
anti_join(scens,future_scenarios,by=c("setting","sub"))%>%distinct(sub)

read_csv(file.path("/scicore/home/pothin/GROUP/ZenabuPhD/OMinput/countryDat2futureScenarios_ST.csv"))%>%
  rename(sub=sub_futureScenarios)->countryDat2futureScenarios
future_scenarios%>%
  left_join(countryDat2futureScenarios,by="sub")%>%
  mutate(sub=ifelse(is.na(sub_countryDat),sub,sub_countryDat))%>%
  mutate(sub_countryDat=NULL) ->future_scenarios #%>%
# mutate(futITNtypecampaign1=recode(futITNtypecampaign1,`futPBO`="futITN"),
#        futITNtypecampaign2=recode(futITNtypecampaign2,`futPBO`="futPBO"),
#        futITNtypecampaign3=recode(futITNtypecampaign3,`futPBO`="futIG2"))->future_scenarios

future_scenarios <- future_scenarios%>%
  select(-futNetcovstart2021, -futcovSMC0to5JulOct)

##################################### keep only CM
scens = right_join(scens,future_scenarios,by=c("setting","sub"))%>%
  filter(sub=="biakoye"&scenario_name=="bau")


# filter for 5 scenarios or simulations
scens <- scens%>%
  filter(sub=="biakoye"&scenario_name=="bau")%>% #&futcovSMC0to5JulOct==0&futNetcovstart2021==0.8
  mutate(histSMCcov2021=0,histSMCcov2022=0,
         futRTSS3covStart2021=0,futRTSS3covStart2022=0,futRTSS3covStart2023=0,
         futRTSS4covStart2021=0,futRTSS4covStart2022=0,futRTSS4covStart2023=0,
         futITNcov2023=0, futcovSMC0to5JulOct=0)#%>% futNetcovstart2021=0.5,
#head(5)



cols_NA<-colnames(scens)[which(scens%>%summarize(across(everything(),function(x) any(is.na(x))))%>%as.vector==T)]
if(length(cols_NA)!=0){
  print("The following scens columns have NAs, please revise:")
  print(cols_NA)
}
print(paste(nrow(scens), "scenarios created for", iso_code))

#--------------------end of scens-----------------------------------


sub="biakoye"
expName<-paste0("comm_protect_llin_CMOnly_", sub)

#####################################
# Experiment setup
#####################################
# Definition of the main country folder where all results will be stored
setupDirs(experimentName = expName, rootDir = root_dir_path, replace = TRUE)

# Initialize cache and create the base xml file
baseList_country = create_baseList(country_name = iso_code,
                                   sim_start = "1918-01-01",
                                   versionnum = 44L)

createBaseXml(baseList_country, replace = TRUE)
baseXML<-file.path(root_dir_path,"comm_protect_llin_CMOnly_biakoye_base.xml")
file.copy(baseXML, getCache(x = "experimentDir"), overwrite = TRUE )

# list Cache
listCache()
thirdDimension<-getCache("thirdDimension")
data.frame(number=1001:1005,
           id=paste0("LLINusers_",thirdDimension$id[1:5]))->LLINusers
thirdDimension<-rbind(thirdDimension,
                      LLINusers)

putCache("thirdDimension",thirdDimension)

## createBaseXml, otherwise the cache is not set up.
setupOM()
#--------------------end of OM setup

scens = finalizeScenarios(scens)

# Store the scenarios in the cache folder
storeScenarios(scens) # 54810 scenarios/simulations

#####################################
# Prepare scripts for creating, running and postprocessing all the scenarios and simulations
#####################################
#####################################
# Validate the xml
#####################################
if (validateXML(xmlfile = getCache(x = "baseXml"),
                schema = file.path(getCache(x = "experimentDir"), "/scenario_44.xsd"),
                scenarios = scens)) {
  print ("XML definition is valid.")
}

# Prepare SLURM scripts
slurmPrepareScenarios(
  expName = expName,
  scenarios = scens,
  nCPU = 10, 
  memCPU = "1GB", 
  time = "00:50:00", 
  qos = "6hours",
  bSize = 200, 
  rModule = "R/4.2.1-foss-2022a"
)

slurmPrepareSimulations(
  expName = paste0(expName, "_", sub),
  scenarios = scens, 
  memCPU = "5GB", 
  nCPU = 50,
  time = "5:55:00", 
  qos = "6hours",
  bSize = 200,
  omModule = "OpenMalaria/44.0-intel-compilers-2023.1.0",
  rModule = "R/4.2.1-foss-2022a"
)


## 3. Prepare the scripts for post-processing the OpenMalaria outputs
# Define the age groups of interest for the outputs (including aggregations)
age_groups_list = list(all=c("0-5","0-10","0-100"),
                       LLINusers=c("0-5","0-10","0-100"))

# Define the OpenMalaria outputs of interest; these will correspond to the
# columns of the results table to be stored in the database
results_columns = c("scenario_id",
                    "date", "age_group", "date_aggregation",
                    "nTreatments1", "nTreatments2", "nTreatments3",
                    "nHost", "nUncomp", "nSevere","nPatent",
                    "tUncomp", "tSevere","nDirDeaths",
                    "incidenceRate", "prevalenceRate")

# Remove the results database if it already exists
# Overwriting an existing database with the same scenario IDs will not work
db_file = file.path(paste0(getCache("rootDir"), expName, ".sqlite"))
if (file.exists(db_file)) {
  print(paste0("A database for ", expName, " exists already and will be removed."))
  file.remove(db_file)
}

# Generate the postprocessing script
# run the caliouptut function below before here
# Make sure to adjust the nCPU, memCPU, time and qos if you run larger experiments
# FC: Increase memory and time, according to https://git.scicore.unibas.ch/idm/countrymodelling/country-modelling-workflow-doc/-/wikis/sciCORE-Resources

slurmPrepareResults(
  expDir = getCache("experimentDir"), dbName = expName,
  resultsName = "om_results", resultsCols = results_columns,
  aggrFun = CalcEpiOutputsCohort,
  aggrFunArgs = list(
    indicators = results_columns,
    aggregateByAgeGroup = age_groups_list,
    aggregateByDate = "year"
  ),
  ntasks = 1, 
  mem = "80G", 
  nCPU = 32, 
  time = "05:00:00", 
  strategy = "batch", 
  qos = "6hours", 
  indexOn = NULL, 
  rModule = "R/4.2.1-foss-2022a"
)







# CalcEpiOutputsCohort function

CalcEpiOutputsCohort <- function(df, indicators = c("incidenceRate", "prevalenceRate"),
                                 aggregateByAgeGroup = NULL, aggregateByDate = NULL,
                                 timeHorizon = NULL, use.gc = TRUE) {
  ## aggregateByAgeGroup needs to be a list of cohorts
  ## aggregateByAgeGroup=list(all=c("0-5","0-10"),LLINusers=c("0-5","0-10"))
  
  ## Make sure input is a data.table
  df <- data.table::as.data.table(df)
  
  ## Check if the required measures are available
  reqMeasures <- list(
    prevalenceRate = c("nPatent", "nHost"),
    incidenceRate = c("nUncomp", "nSevere", "nHost"),
    incidenceRatePerThousand = c("nUncomp", "nSevere", "nHost"),
    tUncomp = c("nTreatments1", "nTreatments2"),
    tSevere = c("nTreatments3"),
    nHosp = c("nHospitalDeaths", "nHospitalRecovs", "nHospitalSeqs"),
    edeath = c("expectedDirectDeaths", "expectedIndirectDeaths", "nHost"),
    edeathRatePerHundredThousand = c(
      "expectedDirectDeaths",
      "expectedIndirectDeaths", "nHost"
    ),
    edirdeath = c("expectedDirectDeaths", "nHost"),
    edirdeathRatePerHundredThousand = c("expectedDirectDeaths", "nHost"),
    ddeath = c("nIndDeaths", "nDirDeaths", "nHost"),
    ddeathRatePerHundredThousand = c("nIndDeaths", "nDirDeaths", "nHost")
  )
  
  ## These columns are not measures and should be ignored
  ## TODO Add cohorts, genotypypes, etc.
  fixedCols <- c(
    "scenario_id", "date_aggregation", "date", "age_group",
    "cohort"
  )
  indicators <- indicators[!indicators %in% fixedCols]
  
  ## Select needed measures for calculations and check if they are in the input.
  ## If not, stop.
  missingMeasures <- unlist(
    reqMeasures[indicators]
  )[!unlist(reqMeasures[indicators]) %in% unique(df[["measure"]])]
  if (length(missingMeasures > 0)) {
    stop(paste("The following measures are missing in your OM output:",
               paste0(unique(missingMeasures), collapse = "\n"),
               "Aborting.",
               sep = "\n"
    ))
  }
  
  ## Collect measures for calculation
  neededMeasures <- unique(unlist(reqMeasures[indicators]))
  ## Add other selected measures
  neededMeasures <- union(neededMeasures, indicators)
  ## We only keep the output columns as defined in indicators.
  dropCols <- neededMeasures[!neededMeasures %in% indicators]
  
  ## Narrow table to selected dates and measures
  df <- df[measure %in% neededMeasures]
  df <- df[, survey_date := as.Date(survey_date)]
  if (!is.null(timeHorizon)) {
    df <- df[survey_date >= timeHorizon[1] & survey_date <= timeHorizon[2]]
  }
  
  ## Split and rename third dimension column and aggregate
  ##
  ## TODO Implement cohorts, which requires the creation of an extra column.
  ##      | third_dimension |  ->  | cohort | age_roup |
  ##      |     AB:0-1      |  ->  |   AB   |    0-1   |
  ##
  ## Don't forget the "none" cohort, which should contain individuals not in any
  ## cohort.
  ## Also take special care with the "none"/0 age group!
  data.table::setnames(df, old = "third_dimension", new = "age_group")
  
  ## NOTE This needs to be extended if we also use genotypes, etc.
  
  
  ## NOTE All of the following aggregations are bottlenecks. So make sure that
  ##      these are as fast as possible. Our best bet is to make sure that
  ##      data.table's GForce optimizes the calls. These are fucking fast.
  ##      Verify and benchmark with options(datatable.verbose = TRUE).
  
  ## Age aggregation
  if (!is.null(aggregateByAgeGroup)) {
    dg<-df[FALSE,]
    for (k in names(aggregateByAgeGroup)) {
      pattern <- if(k=="all"){""}else{paste0(k,"_")}
      groups <- aggregateByAgeGroup[[k]]
      ## Translate "All" to a huge age range. Nobody should be older than 200,
      ## right?
      groups <- data.table::data.table(x = replace(groups, groups == "All", "0-200"))
      ages <- unique(df[["age_group"]])
      
      ## Remove "none" group, if present. This suppresses an annoying warning in
      ## the following age group selection.
      #ages <- ages[!ages %in% "none"]
      ages <- if(k=="all"){ages[!grepl("_",ages)]}else{ages[grepl(k,ages)]}
      ages <- gsub(pattern,"",ages)
      
      tmp <- data.table::data.table(x = ages)
      
      ## Split the age group string (e.g. "0-1") into lo and hi columns (e.g. 1 0)
      groups <- groups[, c("lo", "hi") := data.table::tstrsplit(
        x, "-",
        fixed = TRUE
      )][, c("lo", "hi") := list(as.numeric(lo), as.numeric(hi))]
      
      tmp <- tmp[, c("lo", "hi") := data.table::tstrsplit(
        x, "-",
        fixed = TRUE
      )][
        ,
        c("lo", "hi") := list(as.numeric(lo), as.numeric(hi))
      ]
      
      ## Now check which monitoring age groups are in the range of the requested
      ## age groups. We will use this information later to sum up the values.
      selAges <- list()
      
      for (i in seq_len(nrow(groups))) {
        nn <- unlist(strsplit(as.character(ages[which(tmp[["lo"]] >= groups[[i, "lo"]] & tmp[["hi"]] <= groups[[i, "hi"]])]), "-", fixed = TRUE))
        #nn <- ifelse(length(nn) > 0, paste0(nn[1], "-", nn[length(nn)]), NA)
        nn <- NA
        selAges[[groups[[i, "x"]]]][["groups"]] <- as.character(ages[which(tmp[["lo"]] >= groups[[i, "lo"]] & tmp[["hi"]] <= groups[[i, "hi"]])])
        selAges[[groups[[i, "x"]]]][["adj_name"]] <- nn
      }
      ## Change 0-200 back to All to avoid confusion
      names(selAges)[names(selAges) == "0-200"] <- "All"
      ## Also add back the "none" group
      selAges <- c(selAges, list(none = list(
        groups = "none",
        adj_name = NA
      )))
      
      ## Now we sum up the values for the desired age groups
      dff <- data.table::rbindlist(
        lapply(seq_along(selAges), function(i, list) {
          if (length(list[[i]][["groups"]]) == 0) {
            stop(paste0(
              "Data for age group ", names(list[i]),
              " could not retrieved. Check boundaries.\n"
            ))
          }
          
          # if (!is.na(list[[i]][["adj_name"]]) & list[[i]][["adj_name"]] != names(list[i])) {
          #   warning(
          #     paste0(
          #       "Requested age group ", names(list[i]),
          #       " does not fit monitored boundaries and was adjusted to ",
          #       list[[i]][["adj_name"]], "\n"
          #     )
          #   )
          # }
          
          df[age_group %in% paste0(pattern,list[[i]][["groups"]]),
             .(value = sum(value)),
             by = .(scenario_id, survey_date, measure)
          ][
            , age_group := if (!is.na(list[[i]][["adj_name"]])) {
              list[[i]][["adj_name"]]
            } else {
              if (k=="all"){
                names(list[i])
              }else{
                paste0(k,"_",names(list[i]))
              }
            }
          ]
        }, list = selAges)
      )
      ## Run garbage collector
      dg <- rbind(dg,dff)
      if (use.gc) {
        rm(groups, ages, tmp, selAges, dff)
        gc()
      }
    }
  }
  df<-dg
  if (use.gc) {
    rm(dg)
    gc()
  }
  
  ## Date aggregation
  if (!is.null(aggregateByDate)) {
    ## Join the aggregated column from the dictionary so can apply the correct
    ## modifications.
    df <- df[omOutputDict(),
             aggregated := i.aggregated,
             on = c(measure = "measure_name")
    ]
    ## We create the resulting data table by creating the individual sub-data
    ## tables and then using rbindlist.
    df <- data.table::rbindlist(l = list(
      if ("year" %in% aggregateByDate) {
        data.table::rbindlist(list(
          df[aggregated == TRUE, .(value = sum(value)),
             by = .(scenario_id,
                    date = data.table::year(as.Date(survey_date)),
                    measure, age_group
             )
          ],
          df[aggregated == FALSE, .(value = mean(value)),
             by = .(scenario_id,
                    date = data.table::year(as.Date(survey_date)),
                    measure, age_group
             )
          ]
        ))[
          , c("date", "date_aggregation") := .(
            paste0(date, "-12-31"), rep("year", times = nrow(.SD))
          )
        ]
      },
      if ("month" %in% aggregateByDate) {
        data.table::rbindlist(list(
          df[aggregated == TRUE, .(value = sum(value)),
             by = .(scenario_id,
                    year = data.table::year(as.Date(survey_date)),
                    month = data.table::month(as.Date(survey_date)),
                    measure, age_group
             )
          ],
          df[aggregated == FALSE, .(value = mean(value)),
             by = .(scenario_id,
                    year = data.table::year(as.Date(survey_date)),
                    month = data.table::month(as.Date(survey_date)),
                    measure, age_group
             )
          ]
        ))[
          , c(
            "date", "date_aggregation", "year", "month"
          ) := .(
            paste0(year, "-", sprintf("%02d", month), "-15"),
            rep("month", times = nrow(.SD)), NULL, NULL
          )
        ]
      },
      if (is.null(aggregateByDate) || "survey" %in% aggregateByDate) {
        data.table::rbindlist(list(
          df[aggregated == TRUE, .(value = sum(value)),
             by = .(scenario_id, survey_date, measure, age_group)
          ],
          df[aggregated == FALSE, .(value = mean(value)),
             by = .(scenario_id, survey_date, measure, age_group)
          ]
        ))[
          , c(
            "date", "date_aggregation", "survey_date"
          ) := .(
            as.character(survey_date), rep("survey", times = nrow(df)), NULL
          )
        ]
      }
    ), use.names = TRUE)
    ## Run garbage collector
    if (use.gc) {
      gc()
    }
  }
  
  ## Transform to wide format
  ## TODO Cohorts need to be added
  df <- data.table::dcast(
    df,
    scenario_id + date_aggregation + date + age_group ~ measure,
    value.var = "value"
  )
  
  ## Calculate indicators
  if ("incidenceRate" %in% indicators) {
    df[, incidenceRate := (nUncomp + nSevere) / nHost]
  }
  
  if ("incidenceRatePerThousand" %in% indicators) {
    df[, incidenceRatePerThousand := ((nUncomp + nSevere) / nHost) * 1000]
  }
  
  if ("prevalenceRate" %in% indicators) {
    df[, prevalenceRate := nPatent / nHost]
  }
  
  if ("tUncomp" %in% indicators) {
    df[, tUncomp := nTreatments1 + nTreatments2]
  }
  
  if ("tSevere" %in% indicators) {
    df[, tSevere := nTreatments3]
  }
  
  if ("nHosp" %in% indicators) {
    df[, nHosp := nHospitalDeaths + nHospitalRecovs + nHospitalSeqs]
  }
  
  if ("edeath" %in% indicators) {
    df[, edeath := expectedDirectDeaths + expectedIndirectDeaths]
  }
  
  if ("edeathRatePerHundredThousand" %in% indicators) {
    df[, edeathRatePerHundredThousand := (expectedDirectDeaths + expectedIndirectDeaths) / nHost * 1e5]
  }
  
  if ("edirdeath" %in% indicators) {
    df[, edirdeath := expectedDirectDeaths]
  }
  
  if ("edirdeathRatePerHundredThousand" %in% indicators) {
    df[, edirdeathRatePerHundredThousand := expectedDirectDeaths / nHost * 1e5]
  }
  
  if ("ddeath" %in% indicators) {
    df[, ddeath := nIndDeaths + nDirDeaths]
  }
  
  if ("ddeathRatePerHundredThousand" %in% indicators) {
    df[, ddeathRatePerHundredThousand := (nIndDeaths + nDirDeaths) / nHost * 1e5]
  }
  
  ## Drop non-requested columns
  df <- df[, (dropCols) := NULL]
  ## Make sure the date is a string (Needed for SQLite)
  df <- df[, date := as.character(date)]
  
  return(df)
}

