data(immigrationconjoint)
data(immigrationdesign)

##-----------------------------------------------------
## Define Functions
##-----------------------------------------------------
## Additive, no arguments explicit
amce_mod <- function(varlist, data){
  myformula <- as.formula(paste("Chosen_Immigrant ~ ", paste(varlist, collapse = "+")))
  amce(myformula, data = data)
}

## only additive, no cluster
amce_additive_mod <- function(varlist, data){
  myformula <- as.formula(paste("Chosen_Immigrant ~ ", paste(varlist, collapse = "+")))
  amce(myformula, data = data,
       cluster = FALSE,
       respondent.id = NULL,
       design = "uniform",
       respondent.varying = NULL,
       subset = NULL,
       baselines = NULL,
       weights = NULL,
       na.ignore = FALSE)
}

## Additive and clustering
amce_additive_mod_cluster <- function(varlist, data, clustervar){
  myformula <- as.formula(paste("Chosen_Immigrant ~ ", paste(varlist, collapse = "+")))
  amce(myformula, data = data,
       cluster = TRUE,
       respondent.id = clustervar,
       design = "uniform",
       respondent.varying = NULL,
       subset = NULL,
       baselines = NULL,
       weights = NULL,
       na.ignore = FALSE)
}

## Additive, with cluster and weights
amce_additive_mod_cluster_weights <- function(varlist, data, clustervar, weightvar){
  myformula <- as.formula(paste("Chosen_Immigrant ~ ", paste(varlist, collapse = "+")))
  amce(myformula, data = data,
       cluster = TRUE,
       respondent.id = clustervar,
       design = "uniform",
       respondent.varying = NULL,
       subset = NULL,
       baselines = NULL,
       weights = weightvar,
       na.ignore = FALSE)
}

##Additive with interactions all variables, no cluster
amce_add_interact_mod <- function(varlist, data){
  myformula <- as.formula(paste("Chosen_Immigrant ~ ", paste(varlist, collapse = "*")))
  amce(myformula, data = data,
       cluster = FALSE,
       respondent.id = NULL,
       design = "uniform",
       respondent.varying = NULL,
       subset = NULL,
       baselines = NULL,
       weights = NULL,
       na.ignore = FALSE)
}

## Respondent varying
amce_add_mod_respvar <- function(varlist, data, respvar){
  myformula <- as.formula(paste("Chosen_Immigrant ~ ", paste(varlist, collapse = "+")))
  amce(myformula, data = data,
       cluster = FALSE,
       respondent.id = NULL,
       design = "uniform",
       respondent.varying = respvar,
       subset = NULL,
       baselines = NULL,
       weights = NULL,
       na.ignore = FALSE)
}

##-----------------------------------------------------
## Define Inputs
##-----------------------------------------------------
varlist1 <- c('Gender', 'Education', '`Language Skills`', '`Country of Origin`', 'Job',
              '`Job Experience`', '`Job Plans`', '`Reason for Application`', '`Prior Entry`')
varlist2 <- c('Gender', 'Education')

# create weights in data
weights <- runif(nrow(immigrationconjoint))
immigrationconjoint$weights <- weights
##-----------------------------------------------------
## Run Models
##-----------------------------------------------------
mod1 <- amce_additive_mod(varlist1, immigrationconjoint)
mod2 <- amce_additive_mod_cluster(varlist1, data = immigrationconjoint, clustervar = 'CaseID')
mod3 <- amce_additive_mod_cluster_weights(varlist1, data = immigrationconjoint, clustervar = 'CaseID',
                                          weightvar = 'weights')
##-----------------------------------------------------
## Run Tests
##-----------------------------------------------------
test_that("Output class is correct",{
  expect_is(mod1,"amce")
  expect_is(mod2, "amce")
  expect_is(mod3, "amce")
})

test_that("Warnings are issued",{
  # setting cluster to FALSE
  expect_warning(amce_mod(varlist1, immigrationconjoint),
                 "respondent.id is NULL - setting cluster to FALSE. Please specify a respondent.id variable if you want to estimate clustered standard errors")
})

# test_that("Interactions are calculated"){
#   
# })

#make sure no duplicates after spaces and special characters removed
varlist3 <- c(varlist1, '`Educ ation`')
tempdata <- immigrationconjoint
tempdata$`Educ ation` <- tempdata$Education

test_that("error for duplicate variables",{
  expect_error(amce_additive_mod(varlist3, tempdata),
               'Error: Variable names must be unique when whitespace and meta-characters are removed. Please rename.')
})

# add in any missing base terms for interactions
test_that("warning for missing base terms added",{
  expect_warning(amce(Chosen_Immigrant ~ Gender:Education, data = immigrationconjoint),
                 'Missing base terms for interactions added to formula')
})

## test for missing variables
test_that("error missing variables",{
  var <- 'Education'
  expect_error(amce_additive_mod(varlist1,(immigrationconjoint[, !names(immigrationconjoint) %in% var])),
               paste("Error:", var, "not in 'data'"))
})

# Make sure non-respondent varying are factors
test_that("warning change to factor", {
  temp_num <- immigrationconjoint
  temp_num$Gender <- as.numeric(temp_num$Gender)
  expect_warning(amce(Chosen_Immigrant ~ Education + Gender,
                      respondent.varying = NULL,
                      data = temp_num,
                      cluster = FALSE),
                 'Warning: Gender changed to factor')
})

## missing data
cat("\nMISSING VALUES TEST")
test_that("error missing values in data",{
  tempdf <- immigrationconjoint
  missing <- sample(x = 1:nrow(tempdf), size = .1*nrow(tempdf), replace = FALSE)
  tempdf$Education[missing] <- NA
  expect_error(amce(Chosen_Immigrant ~ Education + Gender,
                      respondent.varying = NULL,
                      data = tempdf,
                      cluster = FALSE),
                 "Error: Education has missing values in 'data'")
})

## missing respondent var in formula
cat("\nMISSING RESPONDENT VAR TEST")
test_that("error missing respondent var in formula",{
  expect_error(amce(Chosen_Immigrant ~ Education + Gender,
                    respondent.varying = "Job",
                    data = immigrationconjoint,
                    design = immigrationdesign,
                    cluster = FALSE),
               "Error: Job is specified in respondent.varying, but is not in the formula")
})



## outcome variable is numeric
cat("\nTEST NUMERIC OUTCOME")
test_that("error missing respondent var in formula",{
  tempdf <- immigrationconjoint
  tempdf$Chosen_Immigrant <- as.factor(tempdf$Chosen_Immigrant)
  expect_error(amce(Chosen_Immigrant ~ Education + Gender,
                      respondent.varying = NULL,
                      data = tempdf,
                      cluster = FALSE),
               "Error: ChosenImmigrant is not numeric or integer")
})


## user supplied baselines
cat("\nTEST BASELINE ERRORS")
test_that("error supplied baseline is not a level",{
  baselines <- list()
  baselines$Gender <- "other"
  expect_error(amce(Chosen_Immigrant ~ Gender + Education +  `Country of Origin`+ `Country of Origin` + Job,
                    data=immigrationconjoint, design=immigrationdesign, respondent.varying=NULL,
                    subset=!is.na(immigrationconjoint$ethnocentrism), respondent.id="CaseID",
                    cluster=TRUE, na.ignore=TRUE, weights=NULL, baselines = baselines),
               "Error: user supplied baseline other is not a level of Gender")
})

## sanity checks clustering
cat("\nTEST CLUSTERING")
test_that("Warning respondent id NULL setting cluster to false",{
  expect_warning(amce(Chosen_Immigrant ~ Gender + Education +  `Country of Origin`+ `Country of Origin` + Job,
                    data=immigrationconjoint, design=immigrationdesign, respondent.varying=NULL,
                    respondent.id=NULL,
                    cluster=TRUE, na.ignore=TRUE, weights=NULL),
               "respondent.id is NULL - setting cluster to FALSE. Please specify a respondent.id variable if you want to estimate clustered standard errors")
})

cat("\nTEST RESPONDENT ID")
test_that("Error respondent.id not in data",{
  expect_error(amce(Chosen_Immigrant ~ Gender + Education +  `Country of Origin`+ `Country of Origin` + Job,
                      data=immigrationconjoint, design=immigrationdesign, respondent.varying=NULL,
                      respondent.id="idvar",
                      cluster=TRUE, na.ignore=TRUE, weights=NULL),
                 "Error: respondent.id not in data")
})

## weights not in data
cat("\nTEST WEIGHTS")
test_that("Error: weights not in data",{
  expect_error(amce(Chosen_Immigrant ~ Gender + Education +  `Country of Origin` + Job,
                      data=immigrationconjoint, design=immigrationdesign, respondent.varying=NULL,
                      respondent.id=NULL,
                      cluster=FALSE, na.ignore=TRUE, weights="wts"),
                 "Error: weights not in data")
})

cat("\nTEST RESPONDENT VARYING")
test_that("Error:respondent.varying var not in the formula",{
  expect_error(amce(Chosen_Immigrant ~ Gender+  `Country of Origin` + Job,
                    data=immigrationconjoint, design=immigrationdesign, respondent.varying="Education",
                    respondent.id="CaseID",
                    cluster=TRUE, na.ignore=TRUE, weights=NULL, baselines = NULL),
               "Error: Education is specified in respondent.varying, but is not in the formula")
})

cat("\nTEST DESIGN OBJECT")
test_that("Error:not in 'design' object",{
  attribute_list <- list()
  attribute_list[["Gender"]] <- c("female","male")
  attribute_list[["Country of Origin"]] <- c("Germany","France","Mexico",
                                             "Philippines","Poland","India",
                                             "China","Sudan","Somalia","Iraq")
  temp_design <- makeDesign(type='constraints', attribute.levels=attribute_list)

  expect_error(amce(Chosen_Immigrant ~ Gender + Education + Job,
                    data=immigrationconjoint, design=temp_design, respondent.varying="Education",
                    respondent.id="CaseID",
                    cluster=TRUE, na.ignore=FALSE, weights=NULL, baselines = NULL),
               "Error: Job not in 'design' object")
})

cat("\nTEST DESIGN OBJECT IN DATA")
test_that("Error:'design' object is not in 'data'",{
  temp_conjoint <- immigrationconjoint[,-4]
  attribute_list <- list()
  attribute_list[["Gender"]] <- c("female","male")
  attribute_list[["Country of Origin"]] <- c("Germany","France","Mexico",
                                             "Philippines","Poland","India",
                                             "China","Sudan","Somalia","Iraq")
  temp_design <- makeDesign(type='constraints', attribute.levels=attribute_list)
  expect_error(amce(Chosen_Immigrant ~ Education,
                    data=temp_conjoint, design=temp_design, respondent.varying="Education",
                    respondent.id="CaseID",
                    cluster=TRUE, na.ignore=FALSE, weights=NULL, baselines = NULL),
               "Error: attribute Gender in 'design' object is not in 'data'")
})

cat("\nTEST FACTOR LEVEL OF ATTRIBUTE")
test_that("Error:factor level of attribute not in 'design' object",{
  attribute_list <- list()
  attribute_list[["Gender"]] <- c("female","male")
  attribute_list[["Country of Origin"]] <- c("Germany","France","Mexico",
                                             "Philippines","Poland","India")
  temp_design <- makeDesign(type='constraints', attribute.levels=attribute_list)

  expect_error(amce(Chosen_Immigrant ~ Gender + Education + `Country of Origin`,
                    data=immigrationconjoint, design=temp_design, respondent.varying="Education",
                    respondent.id="CaseID",
                    cluster=TRUE, na.ignore=FALSE, weights=NULL, baselines = NULL),
               "Error: factor level China of attribute CountryofOrigin not in 'design' object")
})

cat("\nTEST VALID STRING")
test_that("Error:'design' must be a valid character string ('uniform') or a conjointDesign object",{
  expect_error(amce(Chosen_Immigrant ~ Gender + Education + `Country of Origin`,
                    data=immigrationconjoint, design="temp", respondent.varying="Education",
                    respondent.id="CaseID",
                    cluster=TRUE, na.ignore=FALSE, weights=NULL, baselines = NULL),
               "Design object must be a valid character string 'uniform' or a conjointDesign object")
})

## subsetting checks
cat("\nTEST SUBSETTING")
test_that("Warning: invalid argument to 'subset' - must be a logical",{
  sub_data <- immigrationconjoint[!is.na(immigrationconjoint$ethnocentrism),]
  expect_warning(amce(Chosen_Immigrant ~ Gender + Education +  `Country of Origin` + Job,
                      data=immigrationconjoint,
                      design=immigrationdesign,
                      cluster = FALSE,
                      subset=sub_data),
               "Warning: invalid argument to 'subset' - must be a logical")
})


## subsetting checks
cat("\nTEST SUBSET LENGTH")
test_that("Warning: subset must be same length as data",{
  expect_warning(amce(Chosen_Immigrant ~ Gender + Education +  `Country of Origin` + Job,
                      data=immigrationconjoint,
                      design=immigrationdesign,
                      cluster = FALSE,
                      subset=(1 == 1)),
                 "Warning: invalid argument to 'subset' - must be the same length as the number of rows in data")
})


## read qualtrics
cat("\nTEST READ QUALTRICS")
test_that("Either responses or ranks must be non-NULL",{
  expect_error(conjoint_data <- read.qualtrics("CandidateConjointQualtrics.csv",
                                               responses=NULL,
                                               covariates=c("Q6.6"), respondentID="V1"),
                 "Either responses or ranks must be non-NULL")
})

cat("\nTEST NUMBER OF RESPONSES")
test_that("Error: Number of response columns doesn't equal number of tasks in data",{
  expect_error(conjoint_data <- read.qualtrics("CandidateConjointQualtrics.csv",
                                               responses=c("Q2.3", "Q2.7", "Q2.10"),
                                               covariates=c("Q6.6"), respondentID="V1"),
               "Error: Number of response columns doesn't equal number of tasks in data")
})

cat("\nTEST NUMBER OF COLUMNS")
test_that("Error: Number of rank columns doesn't equal number of tasks times number of profiles in data",{
  expect_error(conjoint_data <- read.qualtrics("CandidateConjointQualtrics.csv",
                                               ranks=c("Q2.3", "Q2.7", "Q2.10"),
                                               covariates=c("Q6.6"), respondentID="V1"),
               "Error: Number of rank columns doesn't equal number of tasks times number of profiles in data")
})

cat("\nTEST LEVELS OF VARIABLES")
test_that("Error: levels of variables",{
  df <- immigrationconjoint
  levels(df$Education) <- c(levels(df$Education), "college degree+")
  df$Education[df$Education == "graduate degree"] <- "college degree+"
  #table(df$Education,exclude = NULL)
  expect_error(amce(Chosen_Immigrant ~ Gender + Education + `Language Skills` +
                      `Country of Origin` + Job + `Job Experience` + `Job Plans` +
                      `Reason for Application` + `Prior Entry`, data=df,
                    cluster=TRUE, respondent.id="CaseID", design=immigrationdesign),
               "Error: levels of variable Education are not unique when whitespace and meta-characters are removed. Please rename.")
  expect_error(amce(Chosen_Immigrant ~ Gender + `Language Skills` +
                      `Country of Origin` + Job + `Job Experience` + `Job Plans` +
                      `Reason for Application` + `Prior Entry`, data=df,
                    cluster=TRUE, respondent.id="CaseID", design=immigrationdesign),
               "Error: levels of variable Education used as a dependency, are not unique when whitespace and meta-characters are removed. Please rename.")
})
