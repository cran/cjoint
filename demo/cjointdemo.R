################ Demo for cjoint package - Strezhnev, Hainmueller, Hopkins, Yamamoto ####################

# Load data from immigration conjoint data (from the replication data to "Causal Inference in Conjoint Analysis")
data("immigrationconjoint")
data("immigrationdesign")

##  You can construct the conjoint design manually in in R
attribute_list <- list()
attribute_list[["Education"]] <-c("no formal","4th grade",
                                  "8th grade","high school",
                                  "two-year college","college degree",
                                  "graduate degree")
attribute_list[["Gender"]] <- c("female","male")
attribute_list[["Country of Origin"]] <-  c("Germany","France","Mexico",
                                            "Philippines","Poland","India",
                                            "China","Sudan","Somalia","Iraq")
attribute_list[["Reason for Application"]] <- c("reunite with family",
                                                "seek better job",
                                                "escape persecution")
attribute_list[["Job"]] <- c("janitor","waiter","child care provider",
                             "gardener","financial analyst",
                             "construction worker","teacher",
                             "computer programmer","nurse",
                             "research scientist","doctor")
attribute_list[["Job Experience"]] <- c("none","1-2 years",
                                        "3-5 years","5+ years")
attribute_list[["Job Plans"]] <- c("contract with employer",
                                   "interviews with employer", "will look for work",
                                   "no plans to look for work")
attribute_list[["Prior Entry"]] <- c("never","once as tourist",
                                     "many times as tourist","six months with family",
                                     "once w/o authorization")
attribute_list[["Language Skills"]] <- c("fluent English",
                                         "broken English",
                                         "tried English but unable",
                                         "used interpreter")

# Randomization constraints in the conjoint design
constraint_list <- list()
constraint_list[[1]] <- list()
constraint_list[[1]][["Education"]] <- c("no formal","4th grade","8th grade","high school")
constraint_list[[1]][["Job"]] <- c("financial analyst","computer programmer","research scientist","doctor")
constraint_list[[2]] <- list()
constraint_list[[2]][["Reason for Application"]] <- c("escape persecution")
constraint_list[[2]][["Country of Origin"]] <- c("Germany","France","Mexico","Philippines","Poland","India")

immigrationdesign <- makeDesign(type='constraints', attribute.levels=attribute_list, constraints=constraint_list)

# Run AMCE estimator - this calls the main amce() function in the cjoint package
results <- amce(Chosen_Immigrant ~  Gender + Education + `Language Skills` + `Country of Origin` + Job + `Job Experience` + `Job Plans` + `Reason for Application` + `Prior Entry`, data=immigrationconjoint, cluster=TRUE, respondent.id="CaseID", design=immigrationdesign)
# Print summary
summary(results)
# Plot results
plot(results, xlab="Change in Pr(Immigrant Preferred for Admission to U.S.)", xlim=c(-.3,.3), breaks=c(-.2, 0, .2), labels=c("-.2","0",".2"), text.size=13)

# You can specify interactions in the formula using : or *
interaction_results <- amce(Chosen_Immigrant ~ Gender + Education + Job + Education*Job, data=immigrationconjoint, cluster=TRUE, respondent.id="CaseID", design=immigrationdesign) 
# Print summary
summary(interaction_results)
# Plot results
plot(interaction_results, xlab="Change in Pr(Immigrant Preferred for Admission to U.S.)", xlim=c(-.3,.3), breaks=c(-.2, 0, .2), labels=c("-.2","0",".2"), text.size=13)
