## cjoint: An R Package for estimating Average Marginal Component-specific Effects from conjoint survey experiments
## July, 2015
## Anton Strezhnev, Elissa Berwick, Jens Hainmueller, Daniel Hopkins, Teppei Yamamoto

#############################
# Clustered standard errors
#############################

cluster_se_glm <- function(model, cluster){
    
    #  Drop unused cluster indicators, if cluster var is a factor
    if (class(cluster) == "factor") {
        cluster <- droplevels(cluster)
    }
    
    if (nrow(model.matrix(model)) != length(cluster)) {
        stop("check your data: cluster variable has different N than model - you may have observations with missing data") 
    }
    
    M <- length(unique(cluster))
    N <- length(cluster)           
    K <- model$rank
    
    ## if(M<50) {
    ##     warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
    ## }
   
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
    rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
    return(rcse.cov)
}

##########################################
## function for removing ALL whitespaces
## from elements of vector "vec"
## trailing and leading stripped
## middle replaced with "repl"
## also removes backticks
##########################################

#version for normal vectors
space.begone <- function(vec,repl="_") {   
    sapply(vec,USE.NAMES = F,function(x) {
        #remove leading and trailing whitespaces
        x <- gsub("^\\s+|\\s+$", "", x)
        #also remove backticks
        x <- gsub("`", "", x)
        #replace middle with repl
        gsub("\\s+",repl,x)
    })
}

#########################
## amce function
#########################

amce <- function(formula, data, design="uniform", respondent.varying = NULL, subset=NULL, respondent.id=NULL, cluster=TRUE, na.ignore=FALSE, weights=NULL, baselines = NULL) {

###### Formula and Variables

    #we will split the formula into separate lists
    #unique_vars = all VARIABLES in formula
    #respondent_vars = VARIABLES varying by respondent
    #profile_vars = VARIABLES varying by profile
    #orig_effects = all EFFECTS in formula 
    #inter_effects = higher order EFFECTS in formula
    
    #extract UNIQUE variables and remove ALL whitespaces 
    formula_char <- all.vars(formula)
    formula_char <- space.begone(formula_char)
    respondent.varying <- space.begone(respondent.varying)
    #if this makes for non-unique names, stop
    if(length(unique(formula_char)) != length(formula_char)) {
        stop("Variable names must be unique without whitespace")
    }
    #separate dependent and independent variables, respondent and profile varying
    y_var <-formula_char[1]
    unique_vars <- space.begone(rownames(attr(terms(formula),"factor"))[-1])
    #base respondent variables
    respondent_vars <- unlist(sapply(respondent.varying,function(x) unique_vars[grepl(x,unique_vars)]))
    #base profile variables
    profile_vars <- unique_vars[!is.element(unique_vars,respondent_vars)]

    #fix whitespaces in formula proper as well
    vars <- attr(terms(formula),"term.labels")
    space.vars <- vars[grep("`+",vars)]
    if(length(space.vars) > 0) {
        nospace.vars <- space.begone(space.vars)
        vars <- c(vars[-grep("`+",vars)],nospace.vars)
    }
    #go through all terms and sort within interactions
    vars <- sapply(vars,function(x) paste(sort(strsplit(x,":")[[1]]),collapse = ":"))
    #sort all terms
    vars <- paste(sort(vars),collapse = " + ")
    form <- formula(paste(c(y_var,vars),collapse = "~"))
    
    #and identify ALL original effects; will add in missing base terms automatically
    orig_effects <- attr(terms(form),"term.labels")
    ## add in any missing base terms for interactions
    base_terms <- unique_vars[!is.element(unique_vars,orig_effects)]
    if (length(base_terms > 0)) {
        orig_effects <- c(orig_effects,base_terms)
        form <- formula(paste(c(form,base_terms),collapse="+"))
        warning("Missing base terms for interactions added to formula")
    }

    #identify the requested profile effects
    if (length(respondent.varying) > 0) {
        #identify profile only effects
        profile_effects <- unlist(sapply(respondent.varying,function(x) {
            orig_effects[!grepl(x,orig_effects)]
        }))
    } else {
        profile_effects <- orig_effects
    }

    #remove spaces from column names in data    
    colnames(data) <- space.begone(colnames(data))

#######  Sanity Checks Re: Data

    # Are variables in data?
    for(var in formula_char) {
        if(!(var %in% colnames(data))) {
            stop(paste("Error:", var, "not in 'data'"))
        }
    }

    # Is there missing data?
    if(na.ignore == FALSE){
      for(variab in formula_char){
        if (sum(is.na(data[[variab]])) != 0 ){
          stop(paste("Error:", variab, "has missing values in 'data'"))
        }
      }
    }

    # Is the respondent varying characteristic even in the formula obj?
    if (!is.null(respondent.varying)){
      for (var in respondent.varying){
        found <- 0
        for (formulavars in formula_char){
          if (var == formulavars){
            found <- 1
          }
        }
        if (found == 0){
          stop(paste("Error:", var, "is specified in respondent.varying, but is not in the formula"))
        }
      }
      
    }

    # Check whether outcome variable is a binary 0-1 or numeric
    if (!is.numeric(data[[y_var]]) & !is.integer(data[[y_var]])) {
        stop(paste("Error:", y_var, "is not numeric or integer"))
    }

    # Are the user-supplied desired baselines in the data?
    if (!is.null(baselines)) {
        for(var in names(baselines)) {
            if(!(baselines[[var]] %in% data[[var]])) {
                stop(paste("Error: user supplied baseline",baselines[[var]],"is not a level of",var))
            }
        }      
    }

#######  Sanity Checks Re: clustering

    # If respondent.id is NULL make sure cluster is FALSE
    if(is.null(respondent.id) & cluster == TRUE) {
        warning("respondent.id is NULL - setting cluster to FALSE. Please specify a respondent.id variable if you want to estimate clustered standard errors")
        cluster <- FALSE
    }

    # If (cluster = TRUE), make sure respondent.id is specified and in data
    if (cluster == TRUE) {
        if (is.null(respondent.id)) {
            stop('Error: Must specify a respondent.id if cluster = TRUE')
        } else if (!(respondent.id %in% colnames(data))) {
            stop('Error: respondent.id not in data')
        }
    }

######## Sanity checks re: CMD
    
    #Make R CMD check happy
    J_baseline <- NULL
    J_effect <- NULL

##### Sanity Checks re: design matrix
    
    # If design is already conjointDesign object, proceed to relevant sanity checks
    if (class(design) == "conjointDesign") {
        # Remove whitespaces etc from dimension names of design array 
        names(dimnames(design$J)) <- space.begone(names(dimnames(design$J)))
        #and design dependencies
        names(design$depend) <- space.begone(names(design$depend))
        design$depend <- lapply(design$depend,function(x) space.begone(x))  
        #Now check to make sure profile varying attributes are in conjointDesign
        for (eff in profile_vars) {   
            if (!(eff %in% names(dimnames(design$J)))) {
                stop(paste("Error:", var, "not in 'design' object"))
            }
        }      
        #Check to make sure conjointDesign attributes are in data and level names match
        for (eff in names(dimnames(design$J))) {
            if (!(eff %in% colnames(data))){
                stop(paste("Error: attribute", eff, "in 'design' object is not in 'data'"))
            } else {
        # Check all level names for the attribute in dataset appear in design
                for (lev in levels(as.factor(data[[eff]]))) {
                    if (!(lev %in% dimnames(design$J)[[eff]])) {
                        stop(paste("Error: factor level", lev, "of attribute", eff, "not in 'design' object"))
                    }
                }
            }
        }    
    } else if (design == "uniform") {    
        # else if design == "uniform", create J-dimensional array 
        design <- list()        
        # Determine dimensions
        # And create J matrix with uniform probabilities across all vars
        design.dim <- vector(length=length(profile_vars))
        dim_list <- list()
        for (i in 1:length(profile_vars)) {
            design.dim[i] <- length(unique(data[[profile_vars[i]]]))
            dim_list[[i]] <- levels(factor(data[[profile_vars[i]]]))
        }
        names(dim_list) <- profile_vars
        design$J <- array(1/prod(design.dim), dim=design.dim, dimnames=dim_list)
        design$depend <- compute_dependencies(design$J)       
    } else {
         #if neither uniform nor conjointDesign, error
        stop('Error: argument \'design\' must be a valid character string ("uniform") or a conjointDesign object')   
    }
   
####### Subsetting data    

    if (is.null(subset)) {
        data <- data 
    } else {
        if (class(subset) == "logical") {
            if (length(subset) == nrow(data)) {
                data <- subset(data, subset) 
            } else {
                warning("Warning: invalid argument to 'subset' - must be the same length as the number of rows in data")
            }
        } else {
            warning("Warning: invalid argument to 'subset' - must be a logical")
        }
    }

######### Set baselines if manually assigning
######### And re-factor all of the variables that appear in the design matrix

    attrib <- names(dimnames(design$J))
    for (i in 1:length(names(dimnames(design$J)))) {
        #if a manual baseline was set use that
        if (attrib[i] %in% names(baselines)) {
            base_level <- baselines[[attrib[i]]]
        } else {
            #otherwise use first level in data
            base_level <- levels(data[[attrib[i]]])[1] 
        }
        data[[attrib[i]]] <- factor(data[[attrib[i]]])
        data[[attrib[i]]] <- relevel(data[[attrib[i]]], base_level)
    }
    
####### Adding relevant interaction terms to model

    #If there are any dependencies-- only for profile-varying!
    if(any(profile_vars %in% names(design$depend))) {
        #initialize full interaction set
        depend_vars <- c()     
        #loop over effects with dependencies
        for(eff in profile_vars[profile_vars %in% names(design$depend)]) {
            #initialize interaction set for given variable
            inter <- c()
            #identify higher order occurences of variable
            #make sure it's just that variable, not followed or begun by "_"
            eff_all <- grep(paste(c(":",eff,"(?!_)","|",eff,":(?<!_)"),collapse=""),
                            orig_effects,value=T,perl=T)
            #if you find some, break up, sort and replace ":" with "*"
            if (length(eff_all) > 0) {
                eff_all <- sapply(strsplit(eff_all,":"),function(x) paste(sort(x),collapse="*"))
            }
            #combine with lower order (main effect)
            eff_all <- c(eff,eff_all)
            #for each occurrence, create interaction
            inter <- sapply(eff_all,USE.NAMES = F,function(x) {
                #get conditioning set
                T_r <- design$depend[[eff]]
                #make factors
                for (t in T_r){
                    data[[t]] <- as.factor(data[[t]])
                }
                #combine name and dependency
                T_r_d <- c(x,T_r)
                #sort alphabetically
                T_r_d <- T_r_d[sort.list(T_r_d)]
                #make interaction term
                paste(T_r_d,collapse="*")
            })
            #add to list
            depend_vars <- c(depend_vars,inter)
        }      
        #drop repeats
        depend_vars <- unique(depend_vars)
        #add to formula
        form.full <- formula(paste(c(form,sort(depend_vars)),collapse = " + "))
    }else{
      form.full <- form
    }
    
    #keep everything in alphabetical order!!!!!!!!!!!!
    full.terms <- attr(terms(form.full),"term.labels")
    full.terms <- sapply(full.terms,function(x) paste(sort(strsplit(x,":")[[1]]),collapse=":"))
    form.full <- formula(paste(c(y_var,paste(sort(full.terms),collapse= " + ")),collapse = " ~ "))
    
####### If there are respondent varying terms, split into two formulas
######## One contains only profile effects
######## Second contains all effects
    
    if (length(respondent.varying) > 0) {
        #all variables to be run
        all_run_vars <- attr(terms(form.full),"term.labels")
        #remove those involving respondent things
        prof_only <- sapply(respondent.varying,function(x) {
            all_run_vars[!grepl(x,all_run_vars)] 
        })
        prof_only <- paste(prof_only,collapse = " + ")
        #formula with profile only
        form.prof <- paste(all.vars(form.full)[1],prof_only,sep=" ~ ")
    } else {
        #otherwise use full formula
        form.prof <- form.full
    }
    form.prof <- formula(form.prof)
    
####### Running OLS

    #run model(s)-- if using weights, use tools from "survey" package
    #otherwise run usual lm function
    if (is.null(weights)) {
        lin.mod.prof <- lm(form.prof, data=data)
        if (length(respondent.varying) > 0) {
            lin.mod.full <- lm(form.full, data=data)
        } else {
            lin.mod.full <- NULL
        }
    } else {
        if (cluster) {
            out.design <- svydesign(ids = data[[respondent.id]], weights = data[[weights]], data=data)
        } else {
            out.design <- svydesign(ids = ~0, weights = data[[weights]], data=data)
        }
        lin.mod.prof <- svyglm(form.prof, data=data, design = out.design)
        if (length(respondent.varying) > 0) {
            lin.mod.full <- svyglm(form.full, data=data, design = out.design)
        } else {
            lin.mod.full <- NULL
        }
    }

    #If there's missing data - flag it
    #if (na.ignore == TRUE & !is.null(lin.mod.full$na.action)) {
     #   stop(paste("Error: Observations with missing data in 'data'"))
    #}
  
    #Get sample size
    sample_size_prof <- length(lin.mod.prof$residuals)
    if (length(respondent.varying) > 0) {
        sample_size_full <- length(lin.mod.full$residuals)
    } else {
        sample_size_full <- NULL
    }
    
    #Compute vcov of OLS
    if (is.null(weights) && cluster == TRUE) {
        vcov_mat_prof <- cluster_se_glm(lin.mod.prof, data[[respondent.id]])
        if (length(respondent.varying) > 0) {
            vcov_mat_full <- cluster_se_glm(lin.mod.full, data[[respondent.id]])
        } else {
            vcov_mat_full <- NULL
        }
    } else {
    #Not clustered
        vcov_mat_prof <- vcovHC(lin.mod.prof,type="HC2")
        if (length(respondent.varying) > 0) {
            vcov_mat_full <- vcovHC(lin.mod.full,type="HC2")
        } else {
            vcov_mat_full <- NULL
        }
    }

######### Extract Effects from the profile-vars only linear model

# proposed nomenclature here:
# effect = attribute in question, which has "effect levels"
# depends = attributes it depends on, each of which has "depend levels"

    #before we start, make a blank call to the design array J
    J_call <- Quote(design$J[])
    J_call <- J_call[c(1, 2, rep(3, length(dim(design$J))))]

##############loop over unique profile vars only (AMCE); interactions etc. below

    estimates <- list()
    for(i in 1:length(profile_effects)) {

        #split into sections if it's an interaction
        substrings <- strsplit(profile_effects[i], "[:*]", perl=TRUE)[[1]]

        #administrative loop to find levels
        all_levels <- list()
        all_levels_coefs <- list()
        for(effect in substrings) {
            #get all level names and coefficient names-- sans baseline!!!
            all_levels[[effect]] <- levels(data[[effect]])[-1]
            all_levels_coefs[[effect]] <- sapply(all_levels[[effect]], function(x) {
                paste(effect, x, sep="")
            })
        }
  
        #find all combinations of level names-- add as FIRST column
        levels <- expand.grid(all_levels, stringsAsFactors = FALSE)
        #make level combos in first column
        levels <- cbind(apply(levels,1,function(x) paste(x,collapse=":")),levels)
        colnames(levels) <- c("name",substrings)

        #and all combinations of actual coefficient names
        coefs <- expand.grid(all_levels_coefs, stringsAsFactors = FALSE)
        coefs <- apply(coefs,1,function(x) paste(x,collapse=":"))
       
        # Initialize the results 
        results <- matrix(nrow=2, ncol = nrow(levels))
        if (length(substrings) > 1) {
            rownames(results) <- c("ACIE", "Std. Error")
        } else {
            rownames(results) <- c("AMCE", "Std. Error")
        }
        colnames(results) <- as.character(levels[,1])

         #### find extra times when this effect is mentioned
        #exactly this, not preceded or followed by an "_"
        regexp_name <- paste("(?<!_)","(?<![0-9A-Za-z])",profile_effects[i], "(?![0-9A-Za-z])","(?!_)",sep="")
        all_depends <- grep(regexp_name, attr(terms(form.prof),"term.labels"),value=T, perl=T)
        # remove the actual term
        all_depends <- all_depends[-is.element(all_depends,profile_effects[i])]

 #### loop over every combination of levels of component effects
        for(j in 1:nrow(levels)) {
                
            #figure out which level of inter we're doing
            effect_level <- as.character(levels[j,1])
            effect_level_coef <- coefs[j]

            #get its beta and var-cov matrix
            initial_beta <- coefficients(lin.mod.prof)[effect_level_coef]
            if (effect_level_coef %in% colnames(vcov_mat_prof)) {
                initial_var <- vcov_mat_prof[effect_level_coef, effect_level_coef]
            } else {
                initial_var <- NA
            }

            #if interaction,make sure there is baseline support for this level combination
            if (!is.na(initial_beta) & !is.na(initial_var) & length(substrings) > 1) {
                for (effect1 in substrings) {
                    #get effect base
                    effect_base1 <- levels(data[[effect1]])[1]
                    #subset data to that level
                    base.subset <- data[which(data[[effect1]] == effect_base1),]                       
                    #loop over other profile-varying vars in interaction to subset further
                    for(effect in substrings[!(substrings %in% effect1)]) {
                        base.subset <- base.subset[which(base.subset[[effect]] == as.character(levels[j,effect])),]
                    }                        
                    #if there's no support left, change beta and var to NA
                    if (nrow(base.subset) == 0) {
                        initial_beta <- initial_var <- NA
                        #and give a warning that you had to do it
                        warning(paste("Warning: level",effect_level,"lacks support at",effect1,"baseline",effect_base1, ", effect undefined unless alternative baseline is provided."))
                    }
                }
            }
               
            # If initial_beta and initial_variance are not NA (are valid level combination)
            # and there are dependent variables to incorporate
            if (!is.na(initial_beta) & !is.na(initial_var) & length(all_depends) > 0) {

                #get the slice of design array J associated with baseline and inter level
                J_effect_call <- J_base_call <-  J_call
                for(effect in substrings) {
                    #identify its baseline and modify baseline call accordingly
                    base <- levels(data[[effect]])[1]
                    effect_index <- which(names(dimnames(design$J)) == effect)
                    J_base_call[effect_index + 2] <- base                            
                    #identify level of each effect and modify inter call accordingly
                    level <- levels[j,effect]
                    J_effect_call[effect_index + 2] <- level
                }
                eval(call("<-", Quote(J_baseline), J_base_call))
                eval(call("<-", Quote(J_effect), J_effect_call))

                # Initialize some vectors to store interactions and probabilities
                interaction_probabilities <- c()
                interaction_names <- c()

                #### loop over dependencies for all components of interaction
                for(k in 1:length(all_depends)) {

                    #attribute effect is dependent on
                    depend <- all_depends[[k]]
                    #figure out what levels of what variables are involved
                    substrings_d <- strsplit(depend,":")[[1]]
                    substrings_d <- substrings_d[!is.element(substrings_d,substrings)]
                    all_depend_coefs <- list()
                    for (sub in substrings_d) {
                        all_depend_coefs[[sub]] <-  sapply(levels(data[[sub]]), function(x) paste(sub, x,sep=""))
                    }
                    all_depend_levels <- expand.grid(all_depend_coefs)
                    substrings_l <- strsplit(effect_level_coef,":")[[1]]
                    for (l in length(substrings_l):1) {
                        all_depend_levels <- cbind(rep(substrings_l[l], nrow(all_depend_levels)), all_depend_levels)
                    }
                    colnames(all_depend_levels)[1:length(substrings_l)] <- substrings

                    ####put terms together in proper order
                    regexp_terms <- c(substrings,substrings_d)
                    #split up depends; terms should already be in "proper" order
                    subterms <- strsplit(depend,":")[[1]]
                    #figure out order of just our terms
                    coef.order <- order(sapply(regexp_terms,function(x) grep(x,subterms)))
                    #now put together interaction terms
                    all_depend_level_coefs <- apply(all_depend_levels,1,function(x) {
                        paste(x[coef.order],collapse=":")
                    })

                    #baseline support for depend attribute level 
                    if (!(is.null(dim(J_baseline)))){
                        baseline_support <- apply(J_baseline,substrings_d,sum)
                    } else {
                        baseline_support <- J_baseline
                    }
                    baseline_support[baseline_support != 0] <- 1

                    #probs for depend attribute levels WITH baseline support
                    if (!is.null(dim(J_effect))) {
                        joint_prob <- apply(J_effect, substrings_d, sum)*baseline_support
                    } else {
                        joint_prob <- J_effect*baseline_support
                    }
                    #make it a vector
                    joint_prob <- as.vector(joint_prob)
                    names(joint_prob) <- all_depend_level_coefs
                            
                    ##### loop over levels of depends attribute
                    #if present baselines will omit automatically because coefs are NA
                    for (z in 1:length(all_depend_level_coefs)) {

                        #coefficient name that goes with this effect level & depend level
                        depend_level_coef <- all_depend_level_coefs[z] 
                        #calculate probabilities for this effect and depend level 
                        var_prob <- joint_prob[depend_level_coef]
                        var_prob <- as.numeric(var_prob)/as.numeric(sum(joint_prob))
                                
                        #now add interaction beta and variance to initial 
                        if (!is.na(lin.mod.prof$coefficients[depend_level_coef])) {
                            if (!is.na(vcov_mat_prof[depend_level_coef,depend_level_coef]) & !is.na(vcov_mat_prof[effect_level_coef, depend_level_coef])) {
                                # add weighted beta to initial_beta
                                initial_beta <- initial_beta + var_prob*lin.mod.prof$coefficients[depend_level_coef]
                                # add weighted variance + covariance terms too
                                initial_var <- initial_var + (var_prob^2)*vcov_mat_prof[depend_level_coef, depend_level_coef] +  2*(var_prob)*vcov_mat_prof[effect_level_coef, depend_level_coef]

                                # add probabilities and names to compute covariances
                                interaction_probabilities <- c(interaction_probabilities, var_prob)
                                interaction_names <- c(interaction_names, depend_level_coef)
                            }
                        } #end if that added beta & var, cov
                        
                    } #end loop over levels of dependent attribute
                } #end for loop over different dependent attributes

                # after going through all levels of the depends and all dependent attributes
                # add remaining covariance terms to the parameter variance 
                if (length(interaction_probabilities) > 1) {
                    #loop over all depend attributes 1 to N-1
                    for (x in 1:(length(interaction_probabilities) - 1)) {
                        #loop over depend attributes one ahead of previous to the end
                        for (y in (x+1):(length(interaction_probabilities))) {
                            initial_var <- initial_var + 2*interaction_probabilities[x]*interaction_probabilities[y]*vcov_mat_prof[interaction_names[x], interaction_names[y]]
                        }
                    }
                } #end if has more than 1 depend levels and/or depend attributes

            } #end if has valid beta, var, dependencies
                
            # Store effect and standard error estimates
            results[1,j] <- initial_beta
            if (!is.na(initial_var)) {
                results[2,j] <- sqrt(initial_var)
            } else {
                results[2,j] <- NA
            }
                
        } #end for loop over all level combinations

        # drop elements of the estimates matrix that are NAs
        #keep_vector <- as.vector(which(!is.na(results[1,])))
        #results <- as.matrix(results[c(1,2),keep_vector])

        # combine estimates + SEs into single matrix - store in list
        estimates[[profile_effects[i]]] <- results
            
    } #end for loop over profile effects      

######### Extract Effects from the full model (if have respondent interactions)

# proposed nomenclature here:
# effect = attribute in question, which has "effect levels"
# depends = attributes it depends on, each of which has "depend levels"
# inters = attributes in interaction terms each of which has "inter levels"

 #if there are any respondent effects
    if (length(respondent.varying) > 0) {

        conditional.estimates <- list()
        #get anything from full formula that has respondent vars
        all_run_vars <- attr(terms(form.full),"term.labels")
        resp_effects <- unlist(sapply(respondent.varying,function(x) {
            all_run_vars[grepl(x,all_run_vars)]
        }))
        #find missing base terms
        missing_base <- unlist(sapply(resp_effects,function(x) {
            x1 <- strsplit(x,":")[[1]]
            x2 <- x1[!is.element(x1,respondent_vars)]
            if (length(x2) > 0) paste(x2,collapse=":")
        }))
        resp_effects <- c(missing_base,resp_effects)
        resp_effects_plus <- paste(sort(resp_effects),collapse=" + ")
        form.resp <- formula(paste(c(y_var,resp_effects_plus),collapse=" ~ "))
        all_resp <- attr(terms(form.resp),"term.labels")
        if (any(!is.element(all_resp,all_run_vars))) {
            warning("Warning: mismatch of term names between formulas")
        }
        
        #run admin mod
        admin.mod <- lm(form.resp,data=data)
        all_coef_names <- colnames(model.matrix(admin.mod))
        #ones that don't appear in original model matrix
        all_coef_names <- all_coef_names[is.element(all_coef_names,colnames(vcov_mat_full))]
        #make blank vcov matrix to modify using dependencies
        vcov_new <- matrix(nrow=length(all_coef_names),ncol=length(all_coef_names))
        colnames(vcov_new) <- rownames(vcov_new) <- all_coef_names
        #fill it with the initial values
        vcov_new[all_coef_names,all_coef_names] <- vcov_mat_full[all_coef_names,all_coef_names]

        #initialize list for weighted cross-terms
        covariance_list <- list()

        #loop over respondent-related effects
        for (i in 1:length(all_resp)) {
            
            #split into component effects, if interaction
            substrings <- strsplit(all_resp[i], "[:*]", perl=TRUE)[[1]]

            ## start by finding levels
            #administrative loop over components of interaction
            all_levels <- list()
            all_levels_coefs <- list()
            for(effect in substrings) {
                #if it's not a factor, only has the 1 "level" and coefficient name stays
                if (class(data[[effect]]) != "factor") {
                    all_levels[[effect]] <- effect
                    all_levels_coefs[[effect]] <- effect
                } else {
                    #if it is a factor, get all level names and coefficient names-- sans baseline!!!
                    all_levels[[effect]] <- levels(data[[effect]])[-1]
                    all_levels_coefs[[effect]] <- sapply(all_levels[[effect]],
                                                         function(x) {
                                                             paste(effect, x, sep="")
                                                         })
                }
            }
                
            #find all combinations of level names-- add as FIRST column
            levels <- expand.grid(all_levels, stringsAsFactors = FALSE)
            levels <- cbind(apply(levels,1,function(x) paste(x,collapse=":")),levels)
            colnames(levels) <- c("name",substrings)

            #and all combinations of coefficient names
            coefs <- expand.grid(all_levels_coefs, stringsAsFactors = FALSE)
            coefs <- apply(coefs,1,function(x) paste(x,collapse=":"))
       
            # Initialize the results 
            results <- matrix(nrow=2, ncol = nrow(levels))
            rownames(results) <- c("Conditional Estimate", "Std. Error")
            colnames(results) <- levels[,1]

            #### find extra times when this effect is mentioned in full formula
            # only if anything related to profile var is involved
            if (any(substrings %in% profile_vars)) {
                regexp_name <- paste(sapply(substrings,function(x) paste(c("(?<!_)(?=.*",x,")(?!_)"), collapse="")),collapse="")
                all_depends <- grep(regexp_name,attr(terms(form.full),"term.labels"),value=T, perl=T)
                # remove the actual term
                all_depends <- all_depends[!is.element(all_depends,all_resp[i])]
                # remove any that involve any other respondent varying terms
                resp.other <- respondent_vars[!is.element(respondent_vars,substrings)]
                all_depends <- unlist(sapply(all_depends,function(x) {   
                    sub_depends <- strsplit(x,":")[[1]]
                    if (all(!is.element(sub_depends,resp.other))) x
                }))
            } else {
                #no profile vars, no depends
                all_depends <- c()
            }
                              
            #### loop over every combination of levels of component effects
            for(j in 1:nrow(levels)) {
                
                #figure out which level of inter we're doing
                effect_level <- as.character(levels[j,1])
                effect_level_coef <- coefs[j]

                #get its beta and var-cov matrix
                initial_beta <- coefficients(lin.mod.full)[effect_level_coef]
                if (effect_level_coef %in% colnames(vcov_mat_full)) {
                    initial_var <- vcov_mat_full[effect_level_coef, effect_level_coef]
                } else {
                    initial_var <- NA
                }

                #make sure there is baseline support for this level combination
                if (!is.na(initial_beta) & !is.na(initial_var)) {                   
                    for (effect1 in substrings[substrings %in% profile_vars]) {
                        #get effect base
                        effect_base1 <- levels(data[[effect1]])[1]
                        #subset data to that level
                        base.subset <- data[which(data[[effect1]] == effect_base1),]                       
                        #loop over other profile-varying vars in interaction to subset further
                        for(effect in substrings[substrings %in% profile_vars][!(substrings[substrings %in% profile_vars] %in% effect1)]) {
                            base.subset <- base.subset[which(base.subset[[effect]] == as.character(levels[j,effect])),]
                        }                        
                        #if there's no support left, change beta and var to NA
                        if (nrow(base.subset) == 0) {
                            initial_beta <- initial_var <- NA
                            #and give a warning that you had to do it
                            warning(paste("Warning: level",effect_level,"lacks support at",effect1,"baseline",effect_base1, ", effect undefined unless alternative baseline is provided."))
                        }
                    }
                }
               
                # If initial_beta and initial_variance are not NA and there are depends
                # proceed to add to beta and var
                if (!is.na(initial_beta) & !is.na(initial_var) & length(all_depends) > 0) {

                    #get the slice of design array J associated with baseline and inter level
                    #profile variables only!
                    J_effect_call <- J_base_call <-  J_call
                    for(effect in substrings[substrings %in% profile_vars]) {
                        #identify its baseline and modify baseline call accordingly
                        base <- levels(data[[effect]])[1]
                        effect_index <- which(names(dimnames(design$J)) == effect)
                        J_base_call[effect_index + 2] <- base                           
                        #identify level of each effect and modify inter call accordingly
                        level <- levels[j,effect]
                        J_effect_call[effect_index + 2] <- level
                    }
                    eval(call("<-", Quote(J_baseline), J_base_call))
                    eval(call("<-", Quote(J_effect), J_effect_call))

                    # Initialize some vectors to store interactions and probabilities
                    interaction_probabilities <- c()
                    interaction_names <- c()
                    covariance_names <- c()
                    covariance_probs <- c()
  
                    #### loop over dependencies for all components of effect
                    for(k in 1:length(all_depends)) {

                        #attribute effect is dependent on
                        depend <- all_depends[[k]]
                       #figure out what levels of what variables are involved
                        substrings_d <- strsplit(depend,":")[[1]]
                        substrings_d <- substrings_d[!is.element(substrings_d,substrings)]
                        all_depend_coefs <- list()
                        for (sub in substrings_d) {
                            all_depend_coefs[[sub]] <-  sapply(levels(data[[sub]]), function(x) paste(sub, x,sep=""))
                        }
                        all_depend_levels <- expand.grid(all_depend_coefs)
                        substrings_l <- strsplit(effect_level_coef,":")[[1]]
                        for (l in length(substrings_l):1) {
                            all_depend_levels <- cbind(rep(substrings_l[l], nrow(all_depend_levels)), all_depend_levels)
                        }
                        colnames(all_depend_levels)[1:length(substrings_l)] <- substrings

                        ####put terms together in proper order
                                        #fix escape characters
                        regexp_terms <- unlist(sapply(c(substrings, substrings_d), function(x) {
                            if (grepl("\\(|\\^",x)) {
                                x <- gsub("\\(","\\\\(",x)
                                x <- gsub("\\)","\\\\)",x)
                                gsub("\\^","\\\\^",x)
                            } else {
                                x
                            }
                        }))
                        #split up depends; terms should already be in "proper" order
                        subterms <- strsplit(depend,":")[[1]]
                        #figure out order of just our terms
                        coef.order <- order(sapply(regexp_terms,function(x) grep(x,subterms)))
                        #now put together interaction terms
                        all_depend_level_coefs <- apply(all_depend_levels,1,function(x) {
                            paste(x[coef.order],collapse=":")
                        })

                       #baseline support for depend attribute level in inter
                        if (!(is.null(dim(J_baseline)))){
                            baseline_support <- apply(J_baseline,substrings_d,sum)
                        } else {
                            baseline_support <- J_baseline
                        }
                        baseline_support[baseline_support != 0] <- 1

                        #support for depend attribute levels WITH baseline support
                        if (!is.null(dim(J_effect))) {
                            joint_prob <- apply(J_effect, substrings_d, sum)*baseline_support
                        } else {
                            joint_prob <- J_effect*baseline_support
                        }
                        #make it a vector
                        joint_prob <- as.vector(joint_prob)
                        names(joint_prob) <- all_depend_level_coefs
                            
                        ##### loop over levels of depends attribute
                        #baselines will omit automatically because coefs are NA
                        for (z in 1:length(all_depend_level_coefs)) {

                            #coefficient name that goes with this effect level & depend level
                            depend_level_coef <- all_depend_level_coefs[z] 
                            #calculate probabilities for this effect and depend level 
                            var_prob <- joint_prob[depend_level_coef]
                            var_prob <- as.numeric(var_prob)/as.numeric(sum(joint_prob))
                                
                            #now add interaction beta and variance to initial 
                            if (!is.na(lin.mod.full$coefficients[depend_level_coef])) {
                                if (!is.na(vcov_mat_full[depend_level_coef,depend_level_coef]) & !is.na(vcov_mat_full[effect_level_coef, depend_level_coef])) {
                                    # add probabilities to initial_beta
                                    initial_beta <- initial_beta +
                                        var_prob*lin.mod.full$coefficients[depend_level_coef]
                                    # add variance + covariance terms too
                                    initial_var <- initial_var + (var_prob^2)*vcov_mat_full[depend_level_coef,depend_level_coef] +  2*(var_prob)*vcov_mat_full[effect_level_coef, depend_level_coef]

                                    #modify vcov matrix with terms that exist between
                                    #all variable and variable with depends 
                                    for (coef_name in all_coef_names) {
                                        if (!is.na(vcov_mat_full[coef_name, depend_level_coef])) {
                                            vcov_new[coef_name,effect_level_coef] <- vcov_new[coef_name,effect_level_coef] + (var_prob)*vcov_mat_full[coef_name, depend_level_coef]
                                                #same for symmetrical term!
                                            vcov_new[effect_level_coef,coef_name] <- vcov_new[effect_level_coef,coef_name] + (var_prob)*vcov_mat_full[coef_name, depend_level_coef]
                                        }
                                    }
                                        
                                    # add probabilities and names to compute covariances
                                    interaction_probabilities <- c(interaction_probabilities, var_prob)
                                    interaction_names <- c(interaction_names, depend_level_coef)
                                    #and across different variables
                                    covariance_probs <- c(covariance_probs,var_prob)
                                    covariance_names <- c(covariance_names,depend_level_coef)
                                        
                                }
                            } #end if that added beta & var, cov       
                        } #end loop over levels of dependent attribute
                    } #end for loop over different dependent attributes

                    # after going through all levels of the depends and all dependent attributes
                    # add remaining covariance terms to the parameter variance 
                    if (length(interaction_probabilities) > 1) {
                            #loop over all depend attributes 1 to N-1
                        for (x in 1:(length(interaction_probabilities) - 1)) {
                                #loop over depend attributes one ahead of previous to the end
                            for (y in (x+1):(length(interaction_probabilities))) {
                                initial_var <- initial_var + 2*interaction_probabilities[x]*interaction_probabilities[y]*vcov_mat_full[interaction_names[x], interaction_names[y]]
                            }
                        }
                    } #end if has more than 1 depend levels and/or depend attributes

                    #add names of depend levels and their var probs to list
                    covariance_list[[effect_level_coef]] <- data.frame(covariance_names,covariance_probs)
                        
                } #end if initial beta and var are NA, has depends
                
                # Store effect and standard error estimates
                results[1,j] <- initial_beta
                if (!is.na(initial_var)) {
                    results[2,j] <- sqrt(initial_var)
                } else {
                    results[2,j] <- NA
                }
                
            } #end for loop over all level combinations

            # drop elements of the estimates matrix that are NAs
            #keep_vector <- as.vector(which(!is.na(results[1,])))
            #results <- as.matrix(results[c(1,2),keep_vector])

            # combine estimates + SEs into single matrix - store in list
            conditional.estimates[[all_resp[i]]] <- results
            
        } #end for loop over respondent related effects      
    
       #final modifications for var-cov matrix
       #these only exist when both variables have depends terms
       #so only modify previously modified variables
        if (length(covariance_list) > 1) {
            #loop over each modified coefficient
            for (x in 1:length(covariance_list)) {
                var1 <- names(covariance_list)[x]
                names1 <- as.character(covariance_list[[var1]][,1])
                probs1 <- covariance_list[[var1]][,2]
                #loop over all other modified coefficients, so i != j
                for (y in 1:length(covariance_list)) {
                    var2 <- names(covariance_list)[y]
                    names2 <- as.character(covariance_list[[var2]][,1])
                    probs2 <- covariance_list[[var2]][,2]
                    #loop over each one of interaction names for var1
                    for (z in 1:length(names1)) {
                        vcov_new[var1,var2] <- vcov_new[var1,var2] + sum(probs1[z]*probs2*vcov_mat_full[names1[z],names2])  
                    }
                } 
            }
        }

    } #end if there are any respondent related effects
    
    
############  create conjoint object for output
    
    output <- list()
    class(output) <- c("amce")
    
    #saving things for unconditional estimates
    output$estimates <- estimates
    #saving profile attributes
    output$attributes <- dimnames(design$J)
    #save original profile-only vcov matrix
    output$vcov.coefs <- vcov_mat_prof
    #save sample size used for unconditional estimates
    output$samplesize_prof <- sample_size_prof
    #save style edited formula (no depends)
    output$formula <- form
    
    #saving things for conditional estimates
    if (length(respondent.varying) > 0) {     
        output$cond.estimates <- conditional.estimates
        output$vcov.effects <- vcov_new
        output$samplesize_full <- sample_size_full
        #save style edited formula (no depends)
        output$cond.formula <- form.resp
    }
    
    # Save baselines of unique (main) effects (if factor) to "baselines"
    # If continuous save summary information to "continuous"
    output$baselines <- list()
    output$continuous <- list()
    for (k in unique_vars) {
        if (class(data[[k]]) == "factor") {
            output$baselines[[k]] <- levels(data[[k]])[1]
        } else {
             output$continuous[[k]] <- quantile(model.matrix(form,data)[,k], probs=c(0.25,0.5,0.75), na.rm=T)
            #output$continuous[[k]] <- quantile(data[[k]],  probs=c(0.25,0.5,0.75),na.rm=T)
        }
    }

    #save number of respondents if ID given
    if (!is.null(respondent.id)) {
        output$numrespondents <- length(unique(data[[respondent.id]]))
    } else {
        output$numrespondents <- NULL
    }

    #save respondent variables if given
    if (!is.null(respondent.varying)) {
        output$respondent.varying <- respondent.varying
    } else {
        output$respondent.varying <- NULL
    }

    #save weights as output (if any)
    if (!is.null(weights)) {
        output$weights <- subset(data, select = weights)
    } else {
        output$weights <- NULL
    }

    #save the original data
    output$data <- data
    return(output)
}

############################################################
## summary function for results of main AMCE function
############################################################

# Function for summarizing output from main amce function
# LIST given to "interaction.values" contains VECTORS at which ...
# ... interaction effects will be calculated; default is quantiles...
# ... in the case of continuous, levels otherwise
# ... can be given manual names by naming entries
# Note that both must be NAMED LISTS, name is the respondent.varying effect
# LOGICAL given to show.all indicates whether or not to show ...
# ...  non-interaction terms when reporting interaction effects

summary.amce <- function(object, interaction.values=NULL, show.all=FALSE, ...) {
    amce_obj <- object

######################### administrative section
    
    # Initialize a list to store summary object
    summary_results <- list()
    baselines.prof <- c()
    baselines.resp <- c()

    # Create header of data.frame
    header <- c("Attribute", "Level", "Estimate", "Std. Err", "z value", "Pr(>|z|)", " ")
    names.prof <- c()
    names.resp <- c()
    
    # Extract elements of the formula
    formula_char <- all.vars(amce_obj$formula)
    formula_char <- space.begone(formula_char)

    y_var <-formula_char[1]

##################### reporting AMCE profile-varying results

    #all attribute estimates
    all.prof <- names(amce_obj$estimates)
    nprof <- sum(sapply(all.prof,function(x) length(unlist(amce_obj$estimates[[x]]))/2))

    # Create results matrix (for profile varying only)
    summary_results[["estimates"]] <- matrix(nrow= nprof, ncol=length(header))
    colnames(summary_results[["estimates"]]) <- header
    summary_results[["estimates"]] <- as.data.frame(summary_results[["estimates"]])
    index <- 1

    # Loop over non-respondent varying attributes, which are all factors by assumption
    for (effect in all.prof) {
        
    # Figure out the baseline levels    
        variates <- strsplit(effect, ":")[[1]]
        lev_list <- c()
        for (var in variates) {
            lev_list <- c(lev_list, amce_obj$baselines[[var]])
        }
        baselines.prof <- c(baselines.prof, paste(lev_list,sep="",collapse=":"))
        names.prof <- c(names.prof, effect)
                
        # Append results to the estimates dataframe
        for (p in 1:ncol(amce_obj$estimates[[effect]])) {
                
            summary_results[["estimates"]][index,1] <- effect
            summary_results[["estimates"]][index,2] <- colnames(amce_obj$estimates[[effect]])[p]
            summary_results[["estimates"]][index,3] <- amce_obj$estimates[[effect]][1,p]
            summary_results[["estimates"]][index,4] <- amce_obj$estimates[[effect]][2,p]
            zscr <- amce_obj$estimates[[effect]][1,p]/amce_obj$estimates[[effect]][2,p]
            summary_results[["estimates"]][index,5] <- zscr
            pval <- 2*pnorm(-abs(zscr))
            summary_results[["estimates"]][index,6] <- pval
            
            # Stars!
            if (!is.na(pval)) {
                if (pval < .001) {
                    summary_results[["estimates"]][index,7] <- "***"
                } else if (pval < .01) {
                    summary_results[["estimates"]][index,7] <- "**"
                } else if (pval < .05) {
                    summary_results[["estimates"]][index,7] <- "*"
                } else {
                    summary_results[["estimates"]][index,7] <- ""
                }
            } else {
                summary_results[["estimates"]][index,7] <- ""
            }
            index <- index + 1               
        }
    }

        #save as data frame and save baselines
    summary_results[["estimates"]] <- as.data.frame(summary_results[["estimates"]])
    summary_results[["baselines.prof"]] <- data.frame("Attribute" = names.prof, "Level" = baselines.prof)

################ reporting respondent-varying results

########### any other appearances by the interaction term, like "Education*ethnocentrism"
########### have already had the beta and var, cov for any other appearances added in
########### (ex: "Education:ethnocentrism:Job", "Education:ethnocentrism:Prior_Entry","Education:ethnocentrism:Job:Prior_Entry" )
########### So the only things that need to be grabbed are the profile term (Education)
########### and the collection of things that need "ethnocentrism", in "Education:ethnocentrism"    
    
    # If there are respondent varying attributes, add estimates at levels/quantiles
    if (length(amce_obj$cond.estimates) > 0) {

        # Extract vectors of betas for all coefficients
        # add in a 1 for the unreported intercept
        all.coef <- unlist(lapply(amce_obj$cond.estimates,function(x) colnames(x)))
        beta.vector <- c(1,do.call(cbind,amce_obj$cond.estimates)[1,])
        #also get levels for all factors in the cond.formula
        xlevs <- sapply(all.vars(amce_obj$cond.formula)[-1] [all.vars(amce_obj$cond.formula)[-1] %in% names(amce_obj$baselines)],function(x) levels(amce_obj$data[[x]]), simplify = F)
        
        # Loop over respondent-varying characteristics
        for (effect in amce_obj$respondent.varying) {

            #if it's a factor add its baseline to output
            if (effect %in% names(amce_obj$baselines)) {
                baselines.resp <- c(baselines.resp, amce_obj$baselines[[effect]])
                names.resp <- c(names.resp, effect)
            }

######### now deal with modifying interaction terms          

            #### identify all REQUESTED terms involving effect
            all_req_vars <- attr(terms(amce_obj$formula),"term.labels")
            resp_effects <- gsub(":","*",all_req_vars[grepl(effect,all_req_vars)])
            #use formula to get any missing base terms, just in case
            form.resp <- formula(paste(c(all.vars(amce_obj$formula)[1],paste(resp_effects,collapse=" + ")),collapse=" ~ "))
            all_resp <- attr(terms(form.resp),"term.labels")
            #go through all terms and sort within interactions
            all_resp <- sapply(all_resp,function(x) paste(sort(strsplit(x,":")[[1]]),collapse = ":"))
            #figure out profile attributes these refer to
            all.mod <- unlist(sapply(all_resp,function(x) {
                subs <- strsplit(x,":")[[1]]
                subs <- subs[which(is.element(subs,all.prof))]
                if (length(subs) > 0) paste(subs,collapse=":")
            }))
            #just unique ones
            all.mod <- unique(all.mod)
            #make sure there are some
            if (length(all.mod) == 0) {
                stop(paste(c("respondent characteristic",effect,"not interacted with profile attributes, no interpretation"),collapse=" "))
            }
            #how many terms are affected by this respondent characteristic?
            neffect <- sum(sapply(all.mod,function(x) length(unlist(amce_obj$estimates[[x]]))/2))
            
            # if there are any interactions, proceed to loop over them
            if (length(all.mod) > 0) {
                
                #empty vectors for table key
                tab_name <- c()
                tab_var <- c()
                tab_val <- c()

                #get level names default if no values given or install given names
                if (is.null(interaction.values)) {
                    # if it's a factor get levels
                    if (effect %in% names(amce_obj$baselines)) {
                        lev_list <- colnames(amce_obj$estimates[[effect]])
                    } else {
                    # otherwise get summary information from "continuous"
                        lev_list <- names(amce_obj$continuous[[effect]])
                    }
                } else if (is.null(names(interaction.values[[effect]]))) {
                    # if there are levels given but no names, just number
                    lev_list <- seq(1,length(interaction.values[[effect]]))
                } else {
                    # if the entries are named use those names
                    lev_list <- names(interaction.values[[effect]])
                }
                
               # loop over levels/quantiles of "effect"
                for (i in 1:length(lev_list)) {
                    
                    # make new list entries
                    entry.name <- paste(c(effect,i),collapse="")
                    #if show all is true copy full results matrix
                    if (show.all == TRUE) {
                        summary_results[[entry.name]] <- summary_results[["estimates"]]
                    } else {
                        #if show all is false, make empty results matrix
                        summary_results[[entry.name]] <- matrix(NA,nrow = neffect, ncol=length(header))
                        colnames(summary_results[[entry.name]]) <- header
                        summary_results[[entry.name]] <- as.data.frame(summary_results[[entry.name]])
                        index <- 1
                    }
                    
                    # get the appropriate model matrix for prediction
                    # first, make a data.dummy data matrix
                    data.dummy <- amce_obj$data
                    # if continuous, setting value of "effect" using fake data matrix
                    # also name of level; if continuous just effect name
                    if (effect %in% names(amce_obj$continuous)) {
                        #set value of continuous var
                        if (!is.null(interaction.values)) {
                            data.dummy[[effect]] <- interaction.values[[effect]][i]
                        } else {
                            data.dummy[[effect]] <- amce_obj$continuous[[effect]][i]
                        }
                        #name of "level" is just effect name
                        resp.lev <- effect
                        #as is variable name
                        orig.resp <- effect
                    } else {
                        #if resp var is a factor set level
                        resp.lev <- lev_list[i]
                        # make original var name
                        orig.resp <- paste(c(effect,resp.lev),collapse="")
                         #set level value in data
                        data.dummy[[effect]] <- lev_list[i]
                    }

                    #edit table key
                    tab_name <- c(tab_name,entry.name)
                    tab_var <- c(tab_var, effect)
                    tab_val <- c(tab_val, lev_list[i])
                  
                    # loop over all interactions with "effect"
                    for (prof.var in all.mod) {

                        #loop over the associated betas 
                        for (p in 1:ncol(amce_obj$cond.estimates[[prof.var]])) {

                            #name of level we're modifying
                            prof.lev <- colnames(amce_obj$cond.estimates[[prof.var]])[p]
                            #and where it is in original results, if showing all 
                            if (show.all == TRUE) {
                                index <- which(summary_results[["estimates"]][,1] == prof.var & summary_results[["estimates"]][,2] == prof.lev)
                            }

                            #### set levels within data dummy for prediction
                            # and paste profile name,level together 
                            prof.vars <- strsplit(prof.var, ":")[[1]]
                            prof.levs <- strsplit(prof.lev, ":")[[1]]
                            orig.profs <- c()
                            for (v in 1:length(prof.vars)) {
                                #put together var and level for original name
                                orig.profs[v] <- paste(c(prof.vars[v], prof.levs[v]), collapse="")
                                data.dummy[[prof.vars[v]]] <- prof.levs[v]
                            }
                            #if interaction, put all back together
                            if (length(orig.profs > 1)) {
                                orig.prof <- paste(orig.profs,collapse=":")
                            } else {
                                orig.prof <- orig.profs
                            }
                            
                            #use modified data.dummy to make model matrix, preserving old levels
                            pred.mat <- model.matrix(amce_obj$cond.formula,data.dummy,xlev = xlevs)
                            
                            #all column names containing original profile var name
                            pred.cols <- grep(paste(c("(?=.*", orig.resp, ")(?=.*", orig.prof, ")"), collapse=""),colnames(pred.mat) ,value=T,perl=T)                            
                            #make sure it contains no unrelated vars
                            pred.cols <- unlist(sapply(pred.cols,function(x) {
                                subs <- strsplit(x,":")[[1]]
                                subs <- subs[!is.element(subs,orig.profs)]
                                subs <- subs[!grepl(effect,subs)]
                                if(length(subs) < 1) x
                            }))
                            #add in base
                            pred.cols <- c(orig.prof,pred.cols)
                          
                            #version appearing in results only has level
                            pred.cols2 <- grep(paste(c("(?=.*", resp.lev, ")(?=.*", prof.lev, ")"), collapse=""), names(beta.vector),value=T,perl=T)
                           #make sure it contains ONLY the profile var and respondent var
                           pred.cols2 <- unlist(sapply(pred.cols2,function(x) {
                                subs <- strsplit(x,":")[[1]]
                                subs <- subs[!is.element(subs,prof.levs)]
                                subs <- subs[!grepl(effect,subs)]
                                if(length(subs) < 1) x
                            }))
                            pred.cols2 <- c(prof.lev,pred.cols2)

                            #calculation for coefficient
                            beta.inter <- pred.mat[1,pred.cols] %*% beta.vector[pred.cols2]
                            if (!is.na(beta.inter)) {
                            #gather covariance terms; other occurences dealt with in main fxn
                            all.cov <- c()
                            for (a in 1:length(pred.cols)) {
                                for (b in 1:length(pred.cols)) {
                                    all.cov <- c(all.cov,pred.mat[1,pred.cols[a]]*pred.mat[1,pred.cols[b]]*amce_obj$vcov.effects[pred.cols[a], pred.cols[b]])
                                }
                            }
                            var.inter <- sum(all.cov)
                            #modified zscr and pval
                                se.inter <- sqrt(var.inter)
                                zscr <- beta.inter/se.inter
                                pval <- 2*pnorm(-abs(zscr))
                            } else {
                                se.inter <- zscr <- pval <- NA
                            }

                            #and re-write entries
                            summary_results[[entry.name]][index,1] <- prof.var
                            summary_results[[entry.name]][index,2] <- prof.lev
                            summary_results[[entry.name]][index,3] <- beta.inter
                            summary_results[[entry.name]][index,4] <- se.inter
                            summary_results[[entry.name]][index,5] <- zscr
                            summary_results[[entry.name]][index,6] <- pval

                            # Stars!
                            if (!is.na(beta.inter)) {
                                if (pval < .001) {
                                    summary_results[[entry.name]][index,7] <- "***"
                                } else if (pval < .01) {
                                    summary_results[[entry.name]][index,7] <- "**"
                                } else if (pval < .05) {
                                    summary_results[[entry.name]][index,7] <- "*"
                                } else {
                                    summary_results[[entry.name]][index,7] <- ""
                                }
                            } else {
                                summary_results[[entry.name]][index,7] <- ""
                            }
                            if (show.all == FALSE) {
                                index <- index + 1
                            }
                            
                        } #end loop over levels of profile var
                    } #end loop over modified profile vars
                    
                    # Convert results to data frame
                    summary_results[[entry.name]] <- as.data.frame(summary_results[[entry.name]])

                } #end loop over levels of effect  
            } #end if there are modified profile vars
        } #end loop over respondent varying characteristics

        #save results as data frame and save baselines
         summary_results[["table_values"]] <- data.frame("Table Name" = tab_name, "Level Name" = tab_var, "Level Value" = tab_val)
        summary_results[["table_values"]] <- apply(summary_results[["table_values"]],2,function(x) as.character(x))
        summary_results[["baselines.resp"]] <- data.frame("Attribute" = names.resp, "Level" = baselines.resp)
        
    } else {
        summary_results[["table_values"]] <- NULL
        summary_results[["baselines.resp"]] <- NULL
    }
  
    # Save sample size(s)
    summary_results[["samplesize.estimates"]] <- amce_obj$samplesize_prof
    if (!is.null(amce_obj$samplesize_full)) {
        summary_results[["samplesize.resp"]] <- amce_obj$samplesize_full
    } else {
        summary_results[["samplesize.resp"]] <- NULL
    }

    # If there's a respondent number, add that as well
    if (!is.null(amce_obj$numrespondents)) {
        summary_results[["respondents"]] <- amce_obj$numrespondents
    } else {
        summary_results[["respondents"]] <- NULL
    }

    # Set class
    class(summary_results) <- c("summary.amce")
  
    # Return
    return(summary_results)
}

############################################################
## print summary function for results of main AMCE function
############################################################

print.summary.amce <- function(x, digits=5, ...) {
    summary_result <- x

    #basic print
    cat("-----------------------------------\n")
    cat("Unconditional Component Effects:\n")
    cat("-----------------------------------\n")
    print(summary_result$estimates, digits=digits, row.names=F)
    cat("---\n")
    cat(paste("Number of Obs. = ", summary_result$samplesize.estimates, sep=""))
    cat("\n")
    cat("---\n")   
    if (!is.null(summary_result$respondents)) {
        cat(paste("Number of Respondents = ", summary_result$respondents, sep=""))
        cat("\n")
        cat("---\n")
    }
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05")
    cat("\n")
    cat("\n")
    cat("------------------\n")
    cat("Baseline Levels:\n")
    cat("------------------\n")
    print(summary_result$baselines.prof, row.names=F)
    cat("\n")
    cat("\n")

    #add extra tables for interactions with respondent varying
    if (!is.null(summary_result$table_values)) {
        for (i in 1:nrow(summary_result$table_values)) {
            cat("------------------------------------------------------------\n")
            cat(paste(c("Component Effects Conditional on Respondent-Varying Characteristics","(",summary_result$table_values[i,2],"=", summary_result$table_values[i,3],"):\n"),collapse=" "))
            cat("------------------------------------------------------------\n")
            print(summary_result[[summary_result$table_values[i,1]]], digits=digits, row.names=F)
            cat("---\n")
            cat(paste("Number of Obs. = ", summary_result$samplesize.resp, sep=""))
            cat("\n")
            cat("---\n")   
            if (!is.null(summary_result$respondents)) {
                cat(paste("Number of Respondents = ", summary_result$respondents, sep=""))
                cat("\n")
                cat("---\n")
            }
            cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05")
            cat("\n")
            cat("\n")
            if(!is.null(summary_result$baselines.resp[[summary_result$table_values[i,2]]])) {
                cat("----------------\n")
                cat("Baseline Levels:\n")
                cat("----------------\n")
                print(summary_result$baselines.resp[[summary_result$table_values[i,2]]], row.names=F)
            }
        }   
    }
}

########################################################
# plot amce function
#######################################################

# default will return single plot with point estimates and CI's
# to facet plots, give "facet.name" the name of the variable to facet by...
# ... by default variable facets will be ALL level combinations...
# ... or if a continuous variable is given, will use quantiles...
# ... to customize, directly give "facet.levels" a named LIST ...
# ... with desired values of each variable entered as DATA FRAME with desired names
# show.all tells you what to do with the non-interacted variables; default is false
# show.interaction.base allows you to show base term of interaction or not; default false

plot.amce <- function(x, main="", xlab="Change in E[Y]", ci=.95, colors=NULL, xlim=NULL, breaks=NULL, labels=NULL, attribute_names = NULL, level_names = NULL, label.baseline = TRUE, text.size=11, text.color = "black", point.size = .6, plot.theme = NULL, facet.name = NULL, facet.levels = NULL, show.all = FALSE, show.interaction.base = FALSE, ...) {
    
    # You need ggplot2
    amce_obj <- x
    ylim <- xlim
    
    # Make R CMD check happy
    pe <- NULL
    group <- NULL
    lower <- NULL
    upper <- NULL

################## basic set-up: unconditional effects

    # Extract content from formula object
    formula_char <- all.vars(amce_obj$formula)
    formula_char <- space.begone(formula_char)
  
    y_var <-formula_char[1]

    # Extract raw attribute names from the amce_obj$estimates object
    raw_attributes <- names(amce_obj$estimates)
    # Extract raw levels 
    raw_levels <- list()
    for(m in names(amce_obj$estimates)) {
        raw_levels[[m]] <- colnames(amce_obj$estimates[[m]])
    }
    all.prof <- names(amce_obj$estimates)
    # Determine baseline level for each effect estimate in raw_levels and append to beginning of each vector in raw_levels
    for (effect in names(raw_levels)) {
    # If it's an AMCE
        if (!grepl(":",effect)) {
            #if effect is a factor variable, get actual baseline
            if (effect %in% names(amce_obj$baselines)) {
                raw_levels[[effect]] <- c(amce_obj$baselines[[effect]], raw_levels[[effect]])
            } #otherwise no change     
    # If it's an ACIE  
        } else {
            effect_elements <- strsplit(effect, ":")[[1]]
            baseline_interactions <- c()
            for (elem in effect_elements) {
                #if effect element is a factor variable get baseline
                if (elem %in% names(amce_obj$baselines)) {
                    baseline_interactions <- c(baseline_interactions, amce_obj$baselines[[elem]])
                } else {
                    #otherwise just add name of continuous variable
                    baseline_interactions <- c(baseline_interactions,elem)
                }
            }
            interaction_str <- paste(baseline_interactions,sep="",collapse=":")
            raw_levels[[effect]] <- c(interaction_str, raw_levels[[effect]]) 
        }
    } #matches for names of raw_levels
    
################ facet set-up; unconditional or conditional

    if (length(facet.name) > 0) {
        #is the facet name in unconditional estimates?
        if (is.element(facet.name,names(amce_obj$estimates))) {
            # identify correct ACIE term(s)
            #where it's mentioned
            all.mod <- grep(paste(c("(?<!_)",facet.name,"(?!_)"),collapse=""), names(amce_obj$estimates), perl=T,value=T)
            #remove from "other"
            # identify other terms
            all.other <- names(amce_obj$estimates)[!is.element(names(amce_obj$estimates),all.mod)]
            #figure out what's being modified
            all.mod <- unlist(sapply(all.mod,function(x) {
                subs <- strsplit(x,":")[[1]]
                subs <- subs[!is.element(subs,facet.name)]
                if (length(subs) > 0) paste(subs,collapse=":")
            }))
            #just unique ones
            all.facet <- unique(all.mod)
            #make sure there are some
            if (length(all.facet) == 0) {
                stop(paste(c("Stop: Facet variable",facet.name,"not interacted with other variables"),collapse=" "))
            }
            # remove modified from other too
            all.other <- all.other[!is.element(all.other,all.facet)]
        } else {
            #### identify all REQUESTED terms involving effect
            all_req_vars <- attr(terms(amce_obj$formula),"term.labels")
            resp_effects <- gsub(":","*",all_req_vars[grepl(facet.name,all_req_vars)])
            #use formula to get any missing base terms, just in case
            form.resp <- formula(paste(c(y_var,paste(resp_effects,collapse=" + ")),collapse=" ~ "))
            all_resp <- attr(terms(form.resp),"term.labels")
            #go through all terms and sort within interactions
            all_resp <- sapply(all_resp,function(x) paste(sort(strsplit(x,":")[[1]]),collapse = ":"))
            #figure out profile attributes these refer to
            all.mod <- unlist(sapply(all_resp,function(x) {
                subs <- strsplit(x,":")[[1]]
                subs <- subs[which(is.element(subs,all.prof))]
                if (length(subs) > 0) paste(subs,collapse=":")
            }))
             #make sure there are some
            if (length(all.mod) == 0) {
                stop(paste(c("Stop: Facet variable",facet.name,"not interacted with other variables"),collapse=" "))
            }
            #just unique ones
            all.facet <- unique(all.mod)
             #if we're showing the interaction base, add that on
            if(show.interaction.base == T) {
                all.mod <- c(all.mod,facet.name)
            }
            all.other <- names(amce_obj$estimates)[!is.element(names(amce_obj$estimates),all.facet)]
        }
    } else {
        all.facet <- NULL
        all.other <- names(amce_obj$estimates)
    }

    #if we're faceting but user didn't input levels, establish defaults
    if (length(facet.name) > 0 & length(facet.levels) == 0) {
        facet.levels <- list()      
        #if it's a factor, default facet levels are all levels (sans baseline)
        if (facet.name %in%  names(amce_obj$baselines)) {
            facet.levels[[facet.name]] <- colnames(amce_obj$estimates[[facet.name]])
            #names are levels too
            names(facet.levels[[facet.name]]) <- colnames(amce_obj$estimates[[facet.name]])
        } else if (facet.name %in% names(amce_obj$continuous)) {
             #if it's continuous, default is quantiles
            facet.levels[[facet.name]] <-  amce_obj$continuous[[facet.name]]
            #names stay as values, eg 25%
            names(facet.levels[[facet.name]]) <- names(amce_obj$continuous[[facet.name]])
        } else {
          stop("Error - facet.name in neither baselines nor continuous")
        }
    }
    
############# Adjust attribute and level names
    
    # Sanity check attribute_names and level_names against amce_obj$attributes
    if (is.null(attribute_names)) {
        attribute_names <- raw_attributes
        attribute_names <- sub(":", " X ", attribute_names)
    } else {
        if (length(attribute_names) != length(raw_attributes)) {
            cat(paste("Error: The number of elements in attribute_names ", length(attribute_names), " does not match the number of attributes in amce object for which estimates were obtained", length(amce_obj$attributes), "\n", sep=""))
            cat("Defaulting attribute_names to attribute names in AMCE object\n")
            attribute_names <- raw_attributes
        }
    }
    
  # level_names
  if (is.null(level_names)) {
      level_names <- raw_levels
  } else { 
      for (name in names(raw_levels)) {
          if (length(level_names[[name]]) != length(raw_levels[[name]])) {
              cat(paste("Error: level_names lengths do not match levels for attribute ", name, "\n",sep=""))
              cat(paste("Defaulting level_names for attribute ", name, " to level names in AMCE object", "\n",sep=""))
        level_names[[name]] <- raw_levels[[name]]
          }
      }
  }
  
  ## Fix non-unique attribute and level labels
    attr_uniques <- unique(attribute_names)
    for (i in 1:length(attr_uniques)) {
        attr <- attr_uniques[i]
        num <- 0
        for (j in 1:length(attribute_names)) {
            if (attribute_names[j] == attr){
                leading_spaces <- paste(rep(" ", num), sep="",collapse="")
                attribute_names[j] <- paste(attr, leading_spaces,  sep="")
                num <- num + 1
            }
        }
    }
    
  # Adjust attributes
    uniques <- unique(unlist(level_names))
  # Level labels
    for (i in 1:length(uniques)) {
        level <- uniques[i]
        num <- 0
        for (attr in names(level_names)) {
            for (j in 1:length(level_names[[attr]])) {
                if (level_names[[attr]][j] == level) {
                    leading_spaces <- paste(rep(" ", num), sep="",collapse="")
                    level_names[[attr]][j] <- paste(level, leading_spaces, sep="")
                    num <- num + 1
                }
            }
        }
    }
  
   # If label.baseline = T, label the first level of the attribute as the baseline
   # Unless it's a continious variable; in that case leave it alone
    if (label.baseline) {
        for (attr_ind in names(level_names)) {
            if (!attr_ind %in% names(amce_obj$continuous)) {
                level_names[[attr_ind]][1] <- paste("(Baseline = ", level_names[[attr_ind]][1], ")", sep="")
            }
        }
    }
  
  # Convert ci to z-score
    if (ci < 1 & ci > 0) {
        zscr <- qnorm(1- ((1-ci)/2))
    } else {
        cat("Invalid confidence interval -- Defaulting to 95%")
        zscr <- qnorm(1- ((1-.95)/2))
    }

########################### Compile estimates into plottable objects
    d <- data.frame(pe=c(), se=c(), upper=c(), lower=c(), var=c(), group=c(), facet=c())
  
########## Iterate over estimates that are NOT part of facet calculations

    #if there are any
    if (length(all.other) > 0) {
        
        for (i in 1:length(all.other)) {

            #set up basic group header
            attr_name <- all.other[i]
            print_attr_name <- attribute_names[which(is.element(names(amce_obj$estimates), all.other))][i]
            d_head <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste(print_attr_name, ":", sep=""), group="<NA>",facet="<NA>")

            #if there are facets and show.all=T, repeat same for each
            if (!is.null(facet.name) & show.all == T) {
                for (k in 1:length(facet.levels[[facet.name]])) {
                    #set facet name
                    d_head$facet <- paste(c(facet.name, names(facet.levels[[facet.name]])[k]), collapse = " = ")
                    #add line to plot data
                    d <- rbind(d, d_head)  
                }           
            } else if (is.null(facet.name)) {
                #if no facets, just add once
                d <- rbind(d,d_head)
            } #and nothing if show.all = F and there are facets
       
            #iterate over levels
            for (j in 1:length(level_names[[attr_name]])) {
                
                # If it's the baseline level (all factors here)
                if (j == 1) {  
                    #get the baseline and print a blank line
                    level_name <- level_names[[attr_name]][j]
                    d_lev <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste("   ", level_name,sep=""), group=print_attr_name, facet=NA)    
                } else {
                    #actual level name
                    level <- raw_levels[[attr_name]][j]
                    #level name to print
                    level_name <- level_names[[attr_name]][j]
                    #retrieve estimate and SE
                    val_pe <- amce_obj$estimates[[attr_name]][1,level]
                    val_se <- amce_obj$estimates[[attr_name]][2,level]       
                    #calculate bounds
                    upper_bnd <- val_pe + zscr*val_se
                    lower_bnd <- val_pe - zscr*val_se 
                    #make line to add to plot data
                    d_lev <- data.frame(pe=val_pe, se=val_se, upper=upper_bnd, lower=lower_bnd, var=paste("   ", level_name,sep=""), group=print_attr_name, facet=NA)       
                } #end if a baseline

                #if there are facets and show.all=T, repeat for each
                if (!is.null(facet.name) & show.all == T) {                        
                    for (k in 1:length(facet.levels[[facet.name]])) {
                    #set facet name
                        d_lev$facet <- paste(c(facet.name, names(facet.levels[[facet.name]])[k]),collapse = " = ")
                    #add line to plot data
                        d <- rbind(d, d_lev)  
                    }
                } else if (is.null(facet.name)) {
                #if no facets, just add once
                    d <- rbind(d,d_lev)
                } #and if there are facets and show.all = F, nothing
                
            } #end loop over levels
        
        } #end loop over non-facet related attribute names

    } #end if there are non-facet vars

########## Iterate over estimates that ARE part of facet calculations
    #duplicate d_lev for however many levels/quantiles we have,
    #modifying values and "facet" accordingly 
    if (!is.null(facet.name)) {

        #for each facet level make new set of plot data
        for (k in 1:length(facet.levels[[facet.name]])) {

            # get the appropriate model matrix for prediction
            # first, make a data.dummy data matrix
            data.dummy <- amce_obj$data
            # if facet is continuous, setting value of "effect" using fake data matrix
            # also set name of level; if continuous level is just effect name
            if (facet.name %in% names(amce_obj$continuous)) {
                #set value of continuous var
                data.dummy[[facet.name]] <- facet.levels[[facet.name]][k]
                #name of "level" is name above the facet value
                resp.lev <- facet.name
                #as is variable name
                orig.facet <- facet.name
            } else {
                # if facet is a factor, get the level name
                resp.lev <- as.character(facet.levels[[facet.name]][k])
                # make original var name
                orig.facet <- paste(c(facet.name,resp.lev),collapse="")
                data.dummy[[effect]] <- resp.lev  
            }

            #loop over variables to be modified
            for (mod.var in all.facet) {

               #set up header to reflect base (non-facet) category
                print_attr_name <- mod.var
                d_head <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste(print_attr_name, ":", sep=""), group="<NA>",facet="<NA>")
                d_head$facet <- paste(c(facet.name, names(facet.levels[[facet.name]])[k]),collapse = " = ")
                #add new header
                d <- rbind(d, d_head)

                if (is.element(facet.name,names(amce_obj$estimates))) {
                    #original interaction term
                    attr_name <- grep(paste(c("(?=.*",facet.name,")(?=.*",mod.var,")"), collapse=""),names(amce_obj$estimates),perl=T,value=T)
                    #iterate over levels of mod.var
                    for (p in 1:length(level_names[[mod.var]])) {
                        # If it's the baseline level (all factors here)
                        if (p == 1) {  
                            #get the baseline and print a blank line
                            level_base <- level_names[[mod.var]][p]
                            d_lev <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste("   ",  level_base,sep=""), group=print_attr_name, facet=NA)          
                        } else {
                           #figure out what level we're at
                            mod.lev <- raw_levels[[mod.var]][p]
                            #its name
                            level_base <- level_names[[mod.var]][p]

                            #put together interaction level
                            coef.order <- sapply(c(facet.name,mod.var),function(x) grep(x,strsplit(attr_name,":")[[1]]))
                            level <- paste(c(resp.lev,mod.lev)[coef.order],collapse = ":")
                            
                            #retrieve estimate and SE
                            val_pe <- amce_obj$estimates[[attr_name]][1,level]
                            val_se <- amce_obj$estimates[[attr_name]][2,level]       
                            #calculate bounds
                            upper_bnd <- val_pe + zscr*val_se
                            lower_bnd <- val_pe - zscr*val_se 
                            #make line to add to plot data
                            d_lev <- data.frame(pe=val_pe, se=val_se, upper=upper_bnd, lower=lower_bnd, var=paste("   ", level_base,sep=""), group=print_attr_name, facet=NA)       
                        } #end if a baseline
                        #set facet name
                        d_lev$facet <- paste(c(facet.name, names(facet.levels[[facet.name]])[k]), collapse = " = ")
                        #add level data to plot data
                        d  <- rbind(d, d_lev)
                    } #end loop over levels
                } else { #if not for ACIE: 
  
                    # Extract vectors of betas for all coefficients
                    # add in a 1 for the unreported intercept
                    all.coef <- unlist(lapply(amce_obj$cond.estimates,function(x) colnames(x)))
                    beta.vector <- c(1,do.call(cbind,amce_obj$cond.estimates)[1,])
                    #also get levels for all factors in the cond.formula
                    xlevs <- sapply(all.vars(amce_obj$cond.formula)[-1] [all.vars(amce_obj$cond.formula)[-1] %in% names(amce_obj$baselines)],function(x) levels(amce_obj$data[[x]]), simplify=F)
  
                #iterate over levels of mod.var
                    for (p in 1:length(level_names[[mod.var]])) {

                        # If it's the baseline level (mod.var must be factor)
                        if (p == 1) {
                            #get the baseline and make a blank line
                            level_base <- level_names[[mod.var]][p]
                            d_lev <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste("   ",  level_base,sep=""), group=print_attr_name, facet=NA)        
                        } else {
                            #figure out what level we're at
                            mod.lev <- raw_levels[[mod.var]][p]
                            #its name
                            level_base <- level_names[[mod.var]][p]
                        
                            #paste profile name,level together to get original profile var name
                            #set levels within data dummy for prediction
                            mod.vars <- strsplit(mod.var, ":")[[1]]
                            mod.levs <- strsplit(mod.lev, ":")[[1]]
                            orig.mods <- c()
                            for (v in 1:length(mod.vars)) {
                                #put together var and level
                                orig.mods[v] <- paste(c(mod.vars[v], mod.levs[v]), collapse="")
                                data.dummy[[mod.vars[v]]] <- mod.levs[v]
                            }
                            #if interaction, put all back together
                            if (length(orig.mods > 1)) {
                                orig.mod <- paste(orig.mods,collapse=":")
                            } else {
                                orig.mod <- orig.mods
                            }
                            #use modified data.dummy to make model matrix, preserving old levels
                            pred.mat <- model.matrix(amce_obj$cond.formula,data.dummy,xlev = xlevs)

                            #all column names containing original profile var name
                            pred.cols <- grep(paste(c("(?=.*", orig.mod, ")(?=.*", orig.facet, ")"), collapse=""),colnames(pred.mat) ,value=T,perl=T)                            
                            #make sure it contains no unrelated vars
                            pred.cols <- unlist(sapply(pred.cols,function(x) {
                                subs <- strsplit(x,":")[[1]]
                                subs <- subs[!is.element(subs,orig.mods)]
                                subs <- subs[!grepl(facet.name,subs)]
                                if(length(subs) < 1) x
                            }))
                            #add in base
                            pred.cols <- c(orig.mod,pred.cols)
                          
                            #version appearing in results only has level
                            pred.cols2 <- grep(paste(c("(?=.*", resp.lev, ")(?=.*", mod.lev, ")"), collapse=""), names(beta.vector),value=T,perl=T)
                           #make sure it contains ONLY the profile var and respondent var
                           pred.cols2 <- unlist(sapply(pred.cols2,function(x) {
                                subs <- strsplit(x,":")[[1]]
                                subs <- subs[!is.element(subs,mod.levs)]
                                subs <- subs[!grepl(facet.name,subs)]
                                if(length(subs) < 1) x
                            }))
                            pred.cols2 <- c(mod.lev,pred.cols2)

                            ##### do actual calculation for beta and variance
                            beta.inter <- pred.mat[1,pred.cols] %*% beta.vector[pred.cols2]
                            if (!is.na(beta.inter)) {
                                # for SE, first obtain relevant covariance
                                all.cov <- c()
                                for (a in 1:length(pred.cols)) {
                                    for (b in 1:length(pred.cols)) {
                                        all.cov <- c(all.cov, pred.mat[1,pred.cols[a]]*pred.mat[1,pred.cols[b]]*amce_obj$vcov.effects[pred.cols[a], pred.cols[b]])
                                    }
                                }
                                var.inter <- sum(all.cov)
                                #to report
                                val_pe <- beta.inter
                                val_se <- sqrt(var.inter)
                                #calculate bounds
                                upper_bnd <- val_pe + zscr*val_se
                                lower_bnd <- val_pe - zscr*val_se
                                        
                            #make line to add to plot data
                                d_lev <- data.frame(pe=val_pe, se=val_se, upper=upper_bnd, lower=lower_bnd, var=paste("   ", level_base,sep=""), group=print_attr_name, facet=NA)
                            } else {
                             #if coefficient doesn't exist print NA row
                                d_lev <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste("   ", level_base,sep=""), group=print_attr_name, facet=NA)
                            }
                        } #end if p==1 or continuous

                        #set facet name
                        d_lev$facet <- paste(c(facet.name, names(facet.levels[[facet.name]])[k]), collapse = " = ")
                        #add level data to plot data
                        d  <- rbind(d, d_lev)
                    } #end loop over levels
                } #end if ACIE
   
            } #end loop over all modified vars
        } #end loop over level of facetted variable        
    } else {
        #if there are no facets, remove that column
        d <- d[,1:6]
    }
    
#################    format "d" dataframe
    
  # Set Y bounds
    if (is.null(ylim)) {
        max_upper <- max(d$upper, na.rm=T) + .05
        min_lower <- min(d$lower, na.rm=T) - .05
        ylim <- c(min_lower, max_upper)
        d[is.na(d)] <- max_upper + 100
    } else {
        d[is.na(d)] <- max(ylim) + 100
    }

  # Make group factors <NA> actually NA
    d$group[d$group == "<NA>"] <- NA
    #same with facet
    if(!is.null(facet.name)) d$facet[d$facet == "<NA>"] <- NA
  
  # Reverse factor ordering
    d$var <- factor(d$var,levels=unique(d$var)[length(d$var):1])
    #make facet into factor, if it exists
    if (!is.null(facet.name)) {
        d$facet <- factor(d$facet,levels=unique(d$facet))
    }
    
########## plot output
                                        
    p = ggplot(d, aes(y=pe, x=var, colour=group))
    p = p + coord_flip(ylim = ylim)  
    p = p + geom_hline(yintercept = 0,size=.5,colour="black",linetype="dotted") 
    p = p + geom_pointrange(aes(ymin=lower,ymax=upper),position="dodge",size=point.size)

    #add facetting
    if (!is.null(facet.name)) {
        p = p + facet_wrap(~ facet)
    } 
    
  # If breaks and labels Null, use default
    if (is.null(breaks) & is.null(labels)) {
        p = p + scale_y_continuous(name=xlab)
    } else if (is.null(breaks) & !is.null(labels)) {
        p = p + scale_y_continuous(name=xlab, labels=labels)
    } else if (!is.null(breaks) & is.null(labels)) {
        p = p + scale_y_continuous(name=xlab, breaks=breaks)
    } else if (!is.null(breaks) & !is.null(labels)) {
        p = p + scale_y_continuous(name=xlab, breaks=breaks, labels=labels)
    }
  
    p = p + scale_x_discrete(name="")
  # If there's a title,add it
    if (!is.null(main)) {
        if (main != "") {
            p = p + ggtitle(main)
        }
    }
  # If no colors specified, use default
    if (is.null(colors)) {
        p = p + scale_colour_discrete(" ") 
    } else if (is.vector(colors)) {
    # Make manual palette
        cPal <- rep(colors, ceiling(length(unique(d$group))/length(colors)))
    # Use manual palette
        p = p + scale_colour_manual(values=cPal)
    } else {
        cat("Error: 'colors' must be a vector. Using default colors\n")
        p = p + scale_colour_discrete(" ") 
    }

  # colour scheme
    # if no theme specified, use default
    if (is.null(plot.theme)){
        
        theme_bw1 <- function(base_size = text.size, base_family = "") {
            theme_grey(base_size = base_size, base_family = base_family) %+replace%
            theme(axis.text.x = element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),axis.text.y = element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), axis.ticks = element_line(colour = "grey50"),axis.title.y =  element_text(size = base_size,angle=90,vjust=.01,hjust=.1),plot.title = element_text(face = "bold"),legend.position = "none")
        }
        
        p = p + theme_bw1()
        print(p)
        
    } else if (is.null(class(plot.theme)))  {
      
      cat("Error: 'plot.theme' is not a valid ggplot theme object. Using default theme\n")
      theme_bw1 <- function(base_size = text.size, base_family = "") {
          theme_grey(base_size = base_size, base_family = base_family) %+replace%
          theme(axis.text.x = element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),axis.text.y = element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), axis.ticks = element_line(colour = "grey50"),axis.title.y =  element_text(size = base_size,angle=90,vjust=.01,hjust=.1),plot.title = element_text(face = "bold"),legend.position = "none")
      }
      
      p = p + theme_bw1()
      print(p)
      
  } else if (class(plot.theme)[1] != "theme") {
      
      cat("Error: 'plot.theme' is not a valid ggplot theme object. Using default theme\n")
      theme_bw1 <- function(base_size = text.size, base_family = "") {
          theme_grey(base_size = base_size, base_family = base_family) %+replace%
          theme(axis.text.x = element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),axis.text.y = element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), axis.ticks = element_line(colour = "grey50"),axis.title.y =  element_text(size = base_size,angle=90,vjust=.01,hjust=.1),plot.title = element_text(face = "bold"),legend.position = "none")
      }
      
      p = p + theme_bw1()
      print(p)
        
    # otherwise use the user-passed theme   
    } else {
        p = p + plot.theme
        print(p)
    }
      
}

#####################
# dependencies
######################

compute_dependencies <- function(J, tol=1e-14){
  # Get attribute names
  attribute_names <- names(dimnames(J))
  # If only one attribute, no dependence
  if (length(attribute_names) == 1){
    dependency_list <- list()
    dependency_list[[attribute_names[1]]] <- c()
    return(dependency_list)
  }else{
  # Create list for each attribute_name
  dependency_list <- list()
  for (k in attribute_names){
    dependency_list[[k]] <- c()
  }
  # Loop over each pair of attributes - figure out if they're independent
  for (i in 1:(length(attribute_names)-1)){
    for(j in (i+1):length(attribute_names)){
      attr1 <- attribute_names[i]
      attr2 <- attribute_names[j]
      cross_tab <- apply(J, c(attr1,attr2), sum)
      # Standardize
      sums <- apply(cross_tab, 1, sum)
      cross_tab_std <- cross_tab/sums
      # Compute similarities
      is_equal = TRUE
      r_1 <- cross_tab_std[1,]
      if (nrow(cross_tab_std) > 1){
        for (m in 2:nrow(cross_tab_std)){
          if (any(as.vector(r_1) - as.vector(cross_tab_std[m,]) > tol)){
            is_equal <- FALSE 
          }
        }
      }
      
      # If not the same, append to dependency dictionary
      if (!is_equal){
        dependency_list[[attr1]] <- c(dependency_list[[attr1]], attr2)
        dependency_list[[attr2]] <- c(dependency_list[[attr2]], attr1)
      } 
    }
  }
  return(dependency_list)
  }
}

##############################
# Make design matrix 
#############################

makeDesign <- function(type="file", J=NULL, filename=NULL, attribute.levels = NULL, constraints=NULL, level.probs=NULL, tol=1e-14){
  ## Initialize conjointDesign object
  design.obj <- NULL
  ## If type="file", then try to load from file generated by Conjoint SDT
  if (type == 'file'){
    if (is.null(filename)){
      cat("Error: Must provide a valid filename argument for type = 'file'\n")
      stop()
    }
    # Make a connection and read the lines
    connection <- file(filename, open="r")
    file_lines <- readLines(connection)
    close(connection)
    attr_index <- which(file_lines == "Attributes")
    weight_index <- which(file_lines == "Weights")
    restriction_index <- which(file_lines == "Restrictions")
    
    attributes <- file_lines[(attr_index+1):(weight_index-1)] 
    weight <- file_lines[(weight_index+1):(restriction_index-1)]
    if (restriction_index+1 != length(file_lines)){
      constr <- file_lines[(restriction_index+1):length(file_lines)]
    }else{
      constr <- NULL
    }
    attribute.levels <- list()
    for (attrstr in attributes){
      attributename <- strsplit(attrstr, ":")[[1]][1]
      levels <- strsplit(strsplit(attrstr, ":")[[1]][2], ",")[[1]]
      attribute.levels[[attributename]] <- levels
    }
    level.probs <- list()
    for (probstr in weight){
      attributename <- strsplit(probstr, ":")[[1]][1]
      weights <- strsplit(strsplit(probstr, ":")[[1]][2], ",")[[1]]
      level.probs[[attributename]] <- as.numeric(weights)
    }
    if (is.null(constr) != TRUE){
      constraints <- list()
      for (i in 1:length(constr)){
        allconstraints <- strsplit(constr[i], ";")[[1]]
        constraints[[i]] <- list()
        for (m in allconstraints){
          attributename <- strsplit(m, ":")[[1]][1]
          constrained_levels <- strsplit(strsplit(m, ":")[[1]][2], ",")[[1]]
          constraints[[i]][[attributename]] <- constrained_levels
        }
      }
    }else{
      constraints <- NULL
    }
  }
  ## if type = "array", check whether J is a valid array, then create conjointDesign object  
  if (type == 'array'){
    if (sum(J) != 1){
      cat("Error: Profile assignment probability array invalid: Does not sum to 1\n")
    }else{
      design.obj$J <- J
      design.obj$dependence <- compute_dependencies(J)
    }
  }else if (type == 'constraints' | type == 'file'){
    ## if type = "constraints"
    if (is.null(attribute.levels) | is.list(attribute.levels) != TRUE){
     cat("Error: Must provide a valid list() object in attribute.levels argument for type = 'constraints'\n")
    }
    # Calculate number of dimensions
    dimensions <- c()
    for (attr in names(attribute.levels)){
      dimensions <- c(dimensions, length(attribute.levels[[attr]]))
    }
    # Initialize 
    J_mat <- array(NA, dim=dimensions, dimnames = attribute.levels)

    # Fill in constrained cells with 0 probability
    for (cstr in constraints){
      # Save the names of the constrained attributes
      constraint_names <- names(cstr)
      # Construct a call to select relevant rows
      select_call <- Quote(J_mat[])
      select_call <- select_call[c(1, 2, rep(3, length(dim(J_mat))))]
      for (f in 1:length(constraint_names)){
        name <- constraint_names[f]
        index <- which(names(dimnames(J_mat)) == name)
        select_call[index+2] <- cstr[f]
      }
      # Make a Call
      eval(call("<-", select_call, 0))
    }
    
    # If no randomization weights, then uniform marginal randomization
    if (is.null(level.probs)){
      cell_prob <- 1/sum(is.na(J_mat))
      J_mat[is.na(J_mat)] <- cell_prob
    }else{
      # Normalize level.probs to sum to 1
      for (attr in names(level.probs)){
        level.probs[[attr]] <- level.probs[[attr]]/sum(level.probs[[attr]])
      }
      
      # If no names in level.probs, assume they're in order
      for (attr in names(level.probs)){
        if (is.null(names(level.probs[[attr]]))){
          names(level.probs[[attr]]) <- attribute.levels[[attr]]
        }
      }
      # If randomization weights specified: more calculations
      unconstrained_probs <- J_mat
      unconstrained_probs[TRUE] <- 1
      # For each Attribute
      for (attr in names(dimnames(J_mat))){
        # For each level in each Attribute
        for (level in attribute.levels[[attr]]){
          # Get a marginal probability
          marg_prob <- level.probs[[attr]][[level]]
        
          # Extract the rows pertaining to that marginal probability
          select_call <- Quote(unconstrained_probs[])
          select_call <- select_call[c(1, 2, rep(3, length(dim(J_mat))))]
          
          index <- which(names(dimnames(J_mat)) == attr)
          select_call[index+2] <- level
          # Make a Call
          
          eval(call("<-", select_call, eval(call("*",select_call, marg_prob))))
        }
      }
      missing_prob <- sum(unconstrained_probs[is.na(J_mat)==FALSE])
      increase_prob <- unconstrained_probs*1/(1-missing_prob)

      J_mat[is.na(J_mat)] <- increase_prob[is.na(J_mat)]
    }

    
    
    design.obj$J <- J_mat
    design.obj$dependence <- compute_dependencies(J_mat, tol)
    
  }else{
    cat("Invalid type argument: Must be either 'file', 'array', or 'constraints")
  }

  # Return design.obj if valid
  ## Make it a conjointDesign object 
  class(design.obj) <- "conjointDesign"
  return(design.obj)
}

###########################
# read in qualtrics data
###########################

read.qualtrics <- function(filename, responses, covariates=NULL, respondentID=NULL){
  # If nothing in "responses"

  if (is.null(responses)){
    cat("responses cannot be NULL")
    return(NULL)
  }
  
  # Load CSV Results
  qualtrics_results <- read.csv(filename, stringsAsFactors=F)

  # Extract variable names/question names
  var_names <- as.character(qualtrics_results[1,])
  q_names <- colnames(qualtrics_results)
  qualtrics_results <- qualtrics_results[2:nrow(qualtrics_results),]
  colnames(qualtrics_results) <- var_names
  
  # Find the attribute names
  attr_name_cols <- var_names[grep("F-[0-9]+-[0-9]+(?!-)", var_names, perl=TRUE)]

  # If no attributes fit the description
  if (length(attr_name_cols) == 0){
    cat("Error: Cannot find any columns designating attributes and levels. Please make sure the input file originated from a Qualtrics survey designed using the Conjoint SDT")
    return(NULL)
  }
  
  # Parse to matrix
  attr_name_matrix <- matrix(unlist(strsplit(attr_name_cols,"-")),nrow=3,ncol=length(attr_name_cols))
  colnames(attr_name_matrix) <- attr_name_cols
  attr_name_matrix <- attr_name_matrix[2:nrow(attr_name_matrix),]
  attr_name_matrix <- as.data.frame(t(attr_name_matrix))

  num_tasks <-unique(as.integer(attr_name_matrix[,1]))
  
  # If number of responses doesn't match num_tasks
  if (length(num_tasks) != length(responses)){
    cat("Error: Number of response columns doesn't equal number of tasks in data")
    return(NULL)
  }
  
  # Find the level names
  level_name_cols <- var_names[grep("F-[0-9]+-[0-9]+-[0-9]", var_names, perl=TRUE)]

  # Convert to matrix
  level_name_matrix <- matrix(unlist(strsplit(level_name_cols,"-")),nrow=4,ncol=length(level_name_cols))
  colnames(level_name_matrix) <- level_name_cols
  level_name_matrix <- level_name_matrix[2:nrow(level_name_matrix),]
  level_name_matrix <- as.data.frame(t(level_name_matrix))

  # If respondentID is not null
  if (!is.null(respondentID)){
    respondent_index <- qualtrics_results[,which(q_names %in% respondentID)]

  }else{
    respondent_index <- 1:nrow(qualtrics_results)
  }
  
  # Get the response rows
  if (is.character(responses[1])){
    response_vars <- which(q_names %in% responses)
  }else{
    response_vars <- responses
  }

  # Initialize output dataframe
  out_data_dataset <- NULL
  
  # Parse each row (this is probably really slow...)
  for (r in 1:nrow(qualtrics_results)){
    # Get the attributes
    skipRow = FALSE
    # If no attributes, skip the row
    attribute_refs_1 <- rownames(attr_name_matrix[attr_name_matrix[,1] == 1,])
    attribute_vector_1 <- qualtrics_results[r,attribute_refs_1]

    if (is.na(attribute_vector_1[1])){
      skipRow = TRUE
    }else if (attribute_vector_1[1] == ""){
      skipRow = TRUE
    }
    if (skipRow != TRUE){
    # Extract a covariate vector
    if (!is.null(covariates)){
      covariate_index <- which(q_names %in% covariates)
      covnames <- q_names[covariate_index]
      unit_cov <- qualtrics_results[r,covariate_index]

    }else{
      unit_cov <- c()
    }
    # For each question
    for (k in num_tasks){
       attribute_refs <- rownames(attr_name_matrix[attr_name_matrix[,1] == k,])
       
       attribute_vector <- qualtrics_results[r,attribute_refs]
       
       num_profiles <- as.integer(unique(level_name_matrix[,2]))


       selec_num <- qualtrics_results[r,response_vars[k]]

       # For each profile
       for (j in num_profiles){
         profile_ref <- rownames(level_name_matrix[level_name_matrix[,1] == k&level_name_matrix[,2] == j,])
         profile_levels <- qualtrics_results[r,profile_ref]


         names(profile_levels) <- attribute_vector
         
         if (is.na(as.integer(selec_num))){
           selec <- NA 
         }else if (as.integer(selec_num) == as.integer(j)){
           selec <- 1
         }else{
           selec <- 0
         }

         row_vec <- data.frame(r,respondent_index[r], k, j, profile_levels, selec, unit_cov)

         header <- as.vector(unlist(c("respondentIndex", "respondent","task","profile",attribute_vector, "selected", covnames)))

         colnames(row_vec) <- header

         if (is.null(out_data_dataset)){
           out_data_dataset <- row_vec
         }else{
           out_data_dataset <- rbind(out_data_dataset, row_vec)
         }

       }
    }
    }
  }
  # Do some post-processing
  for (m in attribute_vector){
    out_data_dataset[[m]] <- as.factor(out_data_dataset[[m]])
  }
  out_data_dataset$respondentIndex <- as.integer(out_data_dataset$respondentIndex)
  out_data_dataset$selected <- as.integer(out_data_dataset$selected)
  out_data_dataset$task <- as.integer(out_data_dataset$task)
  out_data_dataset$profile <- as.integer(out_data_dataset$profile)
  
  # Return dataset
  return(out_data_dataset)
}




