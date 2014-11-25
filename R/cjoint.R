## cjoint: An R Package for estimating Average Marginal Component-specific Effects from conjoint survey experiments
## January, 2014
## Anton Strezhnev, Jens Hainmueller, Daniel Hopkins, Teppei Yamamoto

# Clustered standard errors
cluster_se_glm <- function(model, cluster){

  if(nrow(model.matrix(model))!=length(cluster)){
    stop("check your data: cluster variable has different N than model - you may have observations with missing data")
  }
  M <- length(unique(cluster))
  N <- length(cluster)           
  K <- model$rank   
  if(M<50){
    warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
  }
  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
  rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
  return(rcse.cov)
}


## amce function
amce <- function(formula, data, design="uniform", subset=NULL, respondent.id=NULL, cluster=TRUE, na.ignore=FALSE){
  ## Extract formula objects
  formula_char <- as.character(terms(formula))
  effects <- strsplit(formula_char[3], "\\+")[[1]]
  # Strip whitespace from names
  for (i in 1:length(effects)){
    effects[i] <- gsub("^\\s+|\\s+$", "", effects[i])
  }
  y_var <- formula_char[2]
  
  # If respondent.id is NULL make sure cluster is FALSE
  if (is.null(respondent.id)&cluster==TRUE){
    warning("respondent.id is NULL - setting cluster to FALSE. Please specify a respondent.id variable if you want to estimate clustered standard errors")
    cluster <- FALSE
  }
  
  #### Sanity Checks
  ## Ensure all variables in formula are also in data
  if (!(y_var %in% colnames(data))){
    ## If the name contains backticks
    if (gsub("`","",y_var) %in% colnames(data)){
      colnames(data)[colnames(data)==gsub("`","",y_var)] <- y_var
    }else{
      stop(paste("Error:", y_var, "not in 'data'"))
    }
  }
  for (obj in effects){
    vars <- strsplit(obj,"[:*]", perl=TRUE)[[1]]
    # Strip whitespace
    for (i in 1:length(vars)){
      vars[i] <- gsub("^\\s+|\\s+$", "", vars[i])
    }
    for (var in vars){
      
      if (!(var %in% colnames(data))){
        if (gsub("`", "", var) %in% colnames(data)){
          colnames(data)[colnames(data)==gsub("`","",var)] <- var
        }else{
          stop(paste("Error:", var, "not in 'data'"))
        }
      }
    }
  }
  
  ### For interactions, append to effects interactions of * if not in effects
  for (obj in effects){
    vars <- strsplit(obj,"[*]", perl=TRUE)[[1]]
    # Strip whitespace
    for (i in 1:length(vars)){
      vars[i] <- gsub("^\\s+|\\s+$", "", vars[i])
    }
    if (length(vars) > 1){
      for (var in vars){
        if (!(var %in% effects)){
          effects <- c(effects, var)
        } 
      }
    }
  }

  ## Make R CMD check happy
  J_baseline <- NULL
  J_mat <- NULL
  ## Check whether outcome variable is a binary 0-1 
  if (!is.numeric(data[[y_var]])&!is.integer(data[[y_var]])){
    stop(paste("Error:", y_var, "is not numeric or integer"))
  }
  
  
  # Check whether design is valid character string or conjointDesign object
  # Create J-dimensional array if design == "uniform"
  if (class(design) == "conjointDesign"){
    # Check to make sure formula attributes are in conjointDesign
    for (eff in effects){
      vars <- strsplit(eff,"[:*]", perl=TRUE)[[1]]
      # Strip whitespace
      for (i in 1:length(vars)){
        vars[i] <- gsub("^\\s+|\\s+$", "", vars[i])
      }
      for (var in vars){

        if (!(var %in% names(dimnames(design$J)))){
          if (gsub("`", "", var) %in% names(dimnames(design$J))){
            names(dimnames(design$J))[names(dimnames(design$J))==gsub("`", "", var)] <- var
            # Fix dependencies
            for (jm in 1:length(names(design$depend))){
              if (names(design$depend)[jm] == gsub("`", "", var)){
                names(design$depend)[jm] <- var
              }
              for (qf in 1:length(design$depend[[jm]])){
                if (design$depend[[jm]][qf] == gsub("`", "", var)){
                  design$depend[[jm]][qf] <- var
                }
              } 
            }
          }else{
            stop(paste("Error:", var, "not in 'design' object"))
          }
        }
      }
    }

    # Check to make sure conjointDesign attributes are in data and level names match
    for (eff in names(dimnames(design$J))){
      if (!(eff %in% colnames(data))){
        stop(paste("Error: attribute", eff, "in 'design' object is not in 'data'"))
      }else{
        # Check all level names for the attribute in dataset appear in design
        for (lev in levels(as.factor(data[[eff]]))){
          if (!(lev %in% dimnames(design$J)[[eff]])){
            stop(paste("Error: factor level", lev, "of attribute", eff, "not in 'design' object"))
          }
        }
      }
    }    
  }else if (design == "uniform"){
    design <- list()
    # Extract attributes from formula
    attrib <- c()
    for (eff in effects){
      all_eff <- strsplit(eff, "[:*]", perl=TRUE)[[1]]
      # Strip whitespace
      for (i in 1:length(all_eff)){
        all_eff[i] <- gsub("^\\s+|\\s+$", "", all_eff[i])
      }
      for (at in all_eff){
        if (!(at %in% attrib)){
          attrib <- c(attrib, at)
        }
      }
    }

    # Determine dimensions
    design.dim <- vector(length=length(attrib))
    for (i in 1:length(attrib)){
      design.dim[i] <- length(unique(data[[attrib[i]]]))
    }
    
    # Create J matrix with uniform probabilities across all vars
    dim_list <- list()
    for (i in 1:length(attrib)){
      dim_list[[i]] <- levels(factor(data[[attrib[i]]]))
    }

    names(dim_list) <- attrib
    design$J <- array(1/prod(design.dim), dim=design.dim, dimnames=dim_list)
    design$depend <- compute_dependencies(design$J)
  }  
  else{
    stop('Error: argument \'design\' must be a valid character string ("uniform") or a conjointDesign object')
  }

  # If (cluster = TRUE), make sure respondent.id is specified and in data
  if (cluster == TRUE){
    if (is.null(respondent.id)){
      stop('Error: Must specify a respondent.id if cluster = TRUE')
    }else if(!(respondent.id %in% colnames(data))){
      stop('Error: respondent.id not in data')
    }
  }
  
  ## Subset the data
  if (is.null(subset)){
    data <- data 
  }else{
    if (class(subset) == "logical"){
      if (length(subset) == nrow(data)){
        data <- subset(data, subset)
      }else{
        warning("Warning: invalid argument to 'subset' - must be the same length as the number of rows in data")
      }
    }else{
      warning("Warning: invalid argument to 'subset' - must be a logical")
    }
  }
  
  ## Re-factor all of the variables
  for (i in 1:length(names(dimnames(design$J)))){
    attrib <- names(dimnames(design$J))
    base_level <- levels(data[[attrib[i]]])[1]
    data[[attrib[i]]] <- factor(data[[attrib[i]]])
    data[[attrib[i]]] <- relevel(data[[attrib[i]]], base_level)
  }
  
  ## Initialize the interaction set
  regression_vars <- c()
  ### Store effect estimates in a list
  estimates <- list()
  ## Store the unique effect variables in a list
  main_effect_vars <- c()
  
  ### Loop over every effect and add relevant interaction terms to the model
  for (k in 1:length(effects)){
    effect <- effects[k]
    # Check whether it's an amce or acie
    substrings <- strsplit(effect, "[:*]",perl=TRUE)[[1]]
    # Strip whitespace
    for (qr in 1:length(substrings)){
      substrings[qr] <- gsub("^\\s+|\\s+$", "", substrings[qr])
    }
   
    if (length(substrings)==1){
      # It's an AMCE
      # Add to main_effect_vars
      if (!(effect %in% main_effect_vars)){
        main_effect_vars <- c(main_effect_vars, effect)
      }
      # Make it a factor
      
      data[[effect]] <- as.factor(data[[effect]])
      # Figure out the conditioning set
      T_r <- design$depend[[effect]] 
      # Make them factors
      for (t in T_r){
        data[[t]] <- as.factor(data[[t]])
      }
      # If there's anything, append the interaction to the interactions list
      if (length(T_r) > 0){
        T_r_d <- c()
        # If there's a space in T_r, add backticks
        for (T_r_elem in T_r){
          if (grepl("[[:space:]]", T_r_elem)&!grepl("`", T_r_elem)){
            T_r_d <- c(T_r_d, paste("`",T_r_elem, "`",sep=""))
          }else{
            T_r_d <- c(T_r_d, T_r)
          } 
        }
        regression_vars <- c(regression_vars, paste(substrings, T_r_d, sep="*", collapse="*"))
      }else{
        regression_vars <- c(regression_vars, effect)
      }
    }
    if (length(substrings) > 1){
      # If it's an ACIE
      T_r <- c()
      for (effectvar in substrings){
        # Add to main_effect_vars
        if (!(effectvar %in% main_effect_vars)){
          main_effect_vars <- c(main_effect_vars, effectvar)
        }
        # Make the variable a factor 
        data[[effectvar]] <- as.factor(data[[effectvar]])
        # Append to conditioning set
        T_r <- c(T_r, design$depend[[effectvar]])
      }
      # Make the conditioning set factor variables
      for (r in T_r){
        data[[r]] <- as.factor(data[[r]])
      }
      # Append interactions to regression
      if (length(T_r) > 0){
        T_r_d <- c()
        # If there's a space in T_r, add backticks
        for (T_r_elem in T_r){
          if (grepl("[[:space:]]", T_r_elem)&!grepl("`", T_r_elem)){
            T_r_d <- c(T_r_d, paste("`",T_r_elem, "`",sep=""))
          }else{
            T_r_d <- c(T_r_d, T_r)
          } 
        }
        regression_vars <- c(regression_vars, paste(substrings, T_r_d, sep="*", collapse="*"))
      }else{
        regression_vars <- c(regression_vars, paste(substrings, collapse="*"))
      }
    }
  }
  

  form <- paste(y_var, " ~ ", sep = "")
  reg_sum <- paste(regression_vars, sep="", collapse=" + ")
  form <- paste(form, reg_sum, sep=" ")
  form <- as.formula(form)

  old_colnames <- colnames(data)
  # Replace all column names briefly
  for (k in 1:length(colnames(data))){
    if (grepl("`", colnames(data)[k])){
      colnames(data)[k] <- gsub("`","",colnames(data)[k])
    }
  }

  # Run OLS
  lin.mod <- lm(form, data=data)

  # If there's missing data - flag it
  if (na.ignore == TRUE&!is.null(lin.mod$na.action)){
    stop(paste("Error: Observations with missing data in 'data'"))
  }
  
  # Get sample size
  sample_size <- length(lin.mod$residuals)
  # Compute vcov of OLS
  if (cluster == TRUE){
    vcov_mat <- cluster_se_glm(lin.mod, data[[respondent.id]])
  }else{
    # Not clustered
    vcov_mat <- vcov(lin.mod)
  }
  #print(vcov_mat)
  colnames(data) <- old_colnames
  ## Extract Effects from the linear model
  for (k in 1:length(effects)){
    effect <- effects[k]
    # Check whether it's an amce or aice
    substrings <- strsplit(effect, "[:*]", perl=TRUE)[[1]]
    # Strip Whitespace
    for (qr in 1:length(substrings)){
      substrings[qr] <- gsub("^\\s+|\\s+$", "", substrings[qr])
    }

    # If AMCE
    if (length(substrings) == 1){
        results <- matrix(nrow=2, ncol=length(levels(data[[effect]]))-1)
        rownames(results) <- c("AMCE", "Std. Error")
        colnames(results) <- levels(data[[effect]])[2:length(levels(data[[effect]]))]
        
        # Get the support of the baseline
        baseline_lev <- as.character(levels(data[[effect]])[1])
        baseline_level_name <- paste(effect, as.character(levels(data[[effect]])[1]),sep="")
        # Select the subset of J related to the level of beta
        J_call <- Quote(design$J[])
        J_call <- J_call[c(1, 2, rep(3, length(dim(design$J))))]
        index <- which(names(dimnames(design$J)) == effect)
        J_call[index+2] <- sub(effect, "", baseline_lev)
        # Make a Call - store resulting matrix in J_mat
        eval(call("<-", Quote(J_baseline), J_call))       
        
        # For each level (excluding baseline)
        for (l in 2:length(levels(data[[effect]]))){
       
          # Define the name of the level
          lev <- as.character(levels(data[[effect]])[l])
          # Define the combined name of the factor and level
          level_name <- paste(effect, as.character(levels(data[[effect]])[l]), sep="")

          # Get the beta and variance estimates
          initial_beta <- coefficients(lin.mod)[level_name]
          initial_variance <- vcov_mat[level_name, level_name]
          # Make a regular expression to find all other terms containing level_name
          regexp_level_name <- paste(level_name, "(?![ 0-9A-Za-z])",sep="")
          # Find all interaction term coefficients
          extra_effects <- names(coefficients(lin.mod))[grep(regexp_level_name, names(coefficients(lin.mod)), perl=T)]
          # Select only effects with an interaction :
          extra_effects <- extra_effects[grep(":", extra_effects)]

          # If there are interaction terms
          if (length(extra_effects) > 0){
            # Select the subset of J related to the level of beta
            J_call <- Quote(design$J[])
            J_call <- J_call[c(1, 2, rep(3, length(dim(design$J))))]
            index <- which(names(dimnames(design$J)) == effect)
            J_call[index+2] <- sub(effect, "", lev)
            # Make a Call - store resulting matrix in J_mat
            eval(call("<-", Quote(J_mat), J_call))
            
            # Initialize some vectors to store interactions and probabilities
            interaction_probabilities <- c()
            interaction_names <- c()
            
            # Loop over every extra_effect and append to beta/variance estimates
            for (term in extra_effects){
             
              # If the coefficient of that term is not NA
              if (!is.na(lin.mod$coefficients[term])){
                # Split up the terms of the interaction
                with_interaction_terms <- strsplit(term, ":", perl=TRUE)[[1]]
                with_interaction_terms <- with_interaction_terms[with_interaction_terms != level_name]
                # If term is an interaction then figure out what the component levels are and add to beta
                if (length(with_interaction_terms) > 0){
                  # Match attributes to the interaction terms
                  find_prob <- with_interaction_terms
                  matched_attr <- sapply(names(dimnames(design$J)), grepl, find_prob)

                  ### If more matches than length of with_interaction_terms, we need to fix the problem
                  if (sum(as.integer(matched_attr)) != length(with_interaction_terms)){
                    ## If we have more than 1 
                    if (length(with_interaction_terms) > 1){
                      attr_names_match <- colnames(matched_attr)
                      for (m in 1:nrow(matched_attr)){
                        match_fix <- names(matched_attr[m,matched_attr[m,]==TRUE])
                        for (testname in match_fix){
                          matched_attr[m,] <- matched_attr[m,]*grepl(testname, names(matched_attr[m,]))
                        }
                      }
                      matched_attr[,] <- !(matched_attr %in% c(0, FALSE))
                      #matched_attr <- sapply(matched_attr, as.logical)

                      colnames(matched_attr) <- attr_names_match
                    }else{
                      attr_names_match <- names(matched_attr)
                      match_fix <- names(matched_attr[matched_attr == TRUE])
                      for(testname in match_fix){
                        matched_attr <- matched_attr*grepl(testname, names(matched_attr))

                      }
                      matched_attr <- as.logical(matched_attr)
                      names(matched_attr) <- attr_names_match
                    }
                  }

                  interact_attrib <- c()
                  for (m in 1:nrow(as.matrix(matched_attr))){
                    interact_attrib <- c(interact_attrib, as.matrix(matched_attr)[m,][as.matrix(matched_attr)[m,] == TRUE])
                  }
                  interact_attrib <- names(interact_attrib)
                  # Figure out the level for each attribute
                  interact_levels <- c()
                  for (m in 1:length(interact_attrib)){
                      if (grepl("`", find_prob[m])){
                        templev <- sub(interact_attrib[m],"",with_interaction_terms[m])
                        interact_levels <- c(interact_levels, sub("``","",templev))
                      }else{
                        interact_levels <- c(interact_levels, sub(interact_attrib[m],"",with_interaction_terms[m]))
                      }
                  }
                  names(interact_levels) <- interact_attrib
                  # Compute randomization probabilities
                  # Figure out the support of the baseline vars

                  if (!(is.null(dim(J_baseline)))){
                   baseline_support <- apply(J_baseline, interact_attrib, sum)
                  }else{
                    baseline_support <- J_baseline
                  }

                  baseline_support[baseline_support != 0] <- 1

                  
                  # Subset the part of the frame 
                  if (!is.null(dim(J_mat))){
                    joint_prob <- apply(J_mat, interact_attrib, sum)*baseline_support
                  }else{
                    joint_prob <- J_mat*baseline_support
                  }
                  
                  # Make a call to calculate probabilities
                  select_call <- Quote(joint_prob[])
                  select_call <- select_call[c(1, 2, rep(3, length(interact_levels)))]
                  for (leve in 1:length(interact_levels)){
                    select_call[leve+2] <- interact_levels[leve]
                  } 
                  var_prob <- 0

                  eval(call("<-",Quote(var_prob), select_call))

                  var_prob <- as.numeric(var_prob)/as.numeric(sum(joint_prob))
                  if (!is.na(lin.mod$coefficients[term])&!is.na(vcov_mat[term,term])&!is.na(vcov_mat[level_name, term])){
                    # Add probabilities to initial_beta
                    initial_beta <- initial_beta + var_prob*lin.mod$coefficients[term]
                    # Add Variance + Covariance with Beta term to
                    initial_variance <- initial_variance + (var_prob^2)*vcov_mat[term,term] + 2*(var_prob)*vcov_mat[level_name, term]
                  
                    # Add probabilities and names to lists to compute remainder of covariances
                    interaction_probabilities <- c(interaction_probabilities, var_prob)
                    interaction_names <- c(interaction_names, term)
                  }
                }
              }
            }
            # Append remaining covariance terms to the parameter variance
            if (length(interaction_probabilities) > 1){
              for (i in 1:(length(interaction_probabilities)-1)){
                for (j in (i+1):(length(interaction_probabilities))){
                  initial_variance <- initial_variance + 2*interaction_probabilities[i]*interaction_probabilities[j]*vcov_mat[interaction_names[i], interaction_names[j]]                 
                }
              }
            }
          }
          
          # Store effect and standard error estimates
          results[1,lev] <- initial_beta
          results[2,lev] <- sqrt(initial_variance)
        }
        # Combine estimates + SEs into single matrix - store in list
        estimates[[effect]] <- results 
      ### If ACIE
      }else{
        # Denote effect variable
        effectvar <- substrings[1]
        # Get interaction variables
        interaction_vars <- substrings[2:length(substrings)]
        # Confirm baseline has non-zero probability
        baseline_levs <- c(as.character(levels(data[[effectvar]])[1]))
        baseline_level_names <- c(paste(effectvar, as.character(levels(data[[effectvar]])[1]),sep=""))
        for (interact_var in interaction_vars){
          baseline_levs <- c(baseline_levs, as.character(levels(data[[interact_var]])[1]))
          baseline_level_names <- c(baseline_level_names, paste(interact_var, as.character(levels(data[[interact_var]])[1]),sep=""))
        }
        baseline_value <- paste(baseline_levs, collapse=":", sep="")
        
        # Select the subset of J related to the baseline levels of the interaction variables
        J_call <- Quote(design$J[])
        J_call <- J_call[c(1, 2, rep(3, length(dim(design$J))))]
        for (p in 1:length(substrings)){
          index <- which(names(dimnames(design$J)) == substrings[p])
          J_call[index+2] <- sub(substrings[p], "", baseline_levs[p])
        }
        # Make a Call - store resulting matrix in J_baseline
        eval(call("<-", Quote(J_baseline), J_call))         

        # Compute column names for ACIEs
        acie_list <- list()

        for (k in 1:length(substrings)){
          string <- substrings[k]
          len_lev <- length(levels(data[[string]]))
          level_string <- levels(data[[string]])[2:len_lev]
          acie_list[[string]] <- level_string 
        }
        
        acie_matrix <- expand.grid(acie_list, stringsAsFactors = FALSE)
        colnames(acie_matrix) <- substrings
        # Compute number of ACIEs
        num_acie <- nrow(acie_matrix)

        # Determine interaction level names
        interaction_level_names <- c()
        for (p in 1:nrow(acie_matrix)){
          interaction_level_names <- c(interaction_level_names, paste(as.vector(acie_matrix[p,]), collapse=":",sep=""))
        }
        
        # Initialize the results 
        results <- matrix(nrow=2, ncol=num_acie)
        rownames(results) <- c("ACIE", "Std. Error")
        colnames(results) <- c(interaction_level_names)    
        
        # Loop over every row in the interaction list
        for (m in 1:nrow(acie_matrix)){
          
          j <- m
          # Attributes
          attr <- as.vector(colnames(acie_matrix))
          # Get the names of all of the levels
          lev <- as.vector(acie_matrix[m,])
          # Paste levels to attributes
          attribute_levels <- paste(attr, lev, sep="")
          # Find all items that contain the interaction elements
          all_coefficients <- names(coefficients(lin.mod))

          for (attrlev in attribute_levels){
            
            all_coefficients <- all_coefficients[grep(attrlev, all_coefficients, fixed=TRUE)]
            
          }
         
          # Find the "main effect" level (it has length(attribute_levels) - 1 colons)
          main_eff <- NULL
          for (coef in all_coefficients){
            sub_str <- gsub(":","",coef)
            if (nchar(coef) - nchar(sub_str) == length(attribute_levels) - 1){
              main_eff <- coef
            }
          }
          # If we for some reason can't find the 'main effect'
          if (is.null(main_eff)){
            stop("Something went wrong: Couldn't find the main interaction term: This shouldn't happen")
          }
          
          # Get the beta and variance estimates
          initial_beta <- coefficients(lin.mod)[main_eff]
          if (main_eff %in% colnames(vcov_mat)){
            initial_variance <- vcov_mat[main_eff, main_eff]
          }else{
            initial_variance <- NA
          }
          
          
          # If initial_beta and initial_variance are not NA (invalid level combination)
          if (!is.na(initial_beta)&!is.na(initial_variance)){
          
            # Pull out the remaining effect estimates
            extra_effects <- all_coefficients[all_coefficients != main_eff]
            
            # If there are some interaction terms
            if (length(extra_effects) > 0){
              # Select the subset of J related to the level of interaction terms
              J_call <- Quote(design$J[])
              J_call <- J_call[c(1, 2, rep(3, length(dim(design$J))))]
              for (p in 1:length(substrings)){
                index <- which(names(dimnames(design$J)) == substrings[p])
                J_call[index+2] <- sub(substrings[p], "", acie_matrix[m,substrings[p]])
              }
              # Make a Call - store resulting matrix in J_mat
              eval(call("<-", Quote(J_mat), J_call))
              
              # Initialize some vectors to store interactions and probabilities
              interaction_probabilities <- c()
              interaction_names <- c()
              
              # Loop over every extra_effect and append to beta/variance estimates
              for (term in extra_effects){
                # If the coefficient of that term is not NA
                if (!is.na(lin.mod$coefficients[term])){
                  # Split up the terms of the interaction
                  with_interaction_terms <- strsplit(term, ":")[[1]]
                  with_interaction_terms <- with_interaction_terms[!(with_interaction_terms %in% attribute_levels)]
                  # If term is an interaction then figure out what the component levels are and add to beta
                  if (length(with_interaction_terms) > 0){
                    # Match attributes to the interaction terms
                    find_prob <- with_interaction_terms
                    matched_attr <- sapply(names(dimnames(design$J)), grepl, find_prob)
                    ### If more matches than length of with_interaction_terms, we need to fix the problem
                    if (sum(as.integer(matched_attr)) != length(with_interaction_terms)){
                      ## If we have more than 1 
                      if (length(with_interaction_terms) > 1){
                        attr_names_match <- colnames(matched_attr)
                        for (m in 1:nrow(matched_attr)){
                          match_fix <- names(matched_attr[m,matched_attr[m,]==TRUE])
                          for (testname in match_fix){
                            matched_attr[m,] <- matched_attr[m,]*grepl(testname, names(matched_attr[m,]))
                          }
                        }
                        matched_attr[,] <- !(matched_attr %in% c(0, FALSE))
                        #matched_attr <- sapply(matched_attr, as.logical)
                        
                        colnames(matched_attr) <- attr_names_match
                      }else{
                        attr_names_match <- names(matched_attr)
                        match_fix <- names(matched_attr[matched_attr == TRUE])
                        for(testname in match_fix){
                          matched_attr <- matched_attr*grepl(testname, names(matched_attr))
                          
                        }
                        matched_attr <- as.logical(matched_attr)
                        names(matched_attr) <- attr_names_match
                      }
                    }
                    interact_attrib <- c()
                    for (w in 1:nrow(as.matrix(matched_attr))){
                      interact_attrib <- c(interact_attrib, as.matrix(matched_attr)[w,][as.matrix(matched_attr)[w,] == TRUE])
                    }
                    interact_attrib <- names(interact_attrib)
                    # Figure out the level for each attribute
                    interact_levels <- c()
                    for (w in 1:length(interact_attrib)){
                      if (grepl("`", find_prob[w])){
                        templev <- sub(interact_attrib[w],"",with_interaction_terms[w])
                        interact_levels <- c(interact_levels, sub("``","",templev))
                      }else{
                        interact_levels <- c(interact_levels, sub(interact_attrib[w],"",with_interaction_terms[w]))
                      }
                    }
                    names(interact_levels) <- interact_attrib
                    # Compute randomization probabilities

                    # Figure out the support of the baseline vars
                    if (!(is.null(dim(J_baseline)))){
                      baseline_support <- apply(J_baseline, interact_attrib, sum)
                    }else{
                      baseline_support <- J_baseline
                    }
                    baseline_support[baseline_support != 0] <- 1
                    
                    # Subset the part of the frame 
                    if (!is.null(dim(J_mat))){
                      joint_prob <- apply(J_mat, interact_attrib, sum)*baseline_support
                    }else{
                      joint_prob <- J_mat*baseline_support
                    }
                    
                    # Make a call to calculate probabilities
                    select_call <- Quote(joint_prob[])
                    select_call <- select_call[c(1, 2, rep(3, length(interact_levels)))]
                    for (leve in 1:length(interact_levels)){
                      select_call[leve+2] <- interact_levels[leve]
                    } 
                    var_prob <- 0
                    eval(call("<-",Quote(var_prob), select_call))
                    var_prob <- as.numeric(var_prob)/as.numeric(sum(joint_prob))
                    if (!is.na(lin.mod$coefficients[term])&!is.na(vcov_mat[term,term])){
                      # Add probabilities to initial_beta
                      initial_beta <- initial_beta + var_prob*lin.mod$coefficients[term]
                      # Add Variance + Covariance with Beta term to
                      initial_variance <- initial_variance + (var_prob^2)*vcov_mat[term,term] + 2*(var_prob)*vcov_mat[main_eff, term]
                      
                      # Add probabilities and names to lists to compute remainder of covariances
                      interaction_probabilities <- c(interaction_probabilities, var_prob)
                      interaction_names <- c(interaction_names, term)
                    }
                  }
                }
              }
              # Append remaining covariance terms to the parameter variance
              if (length(interaction_probabilities) > 1){
                for (q in 1:(length(interaction_probabilities)-1)){
                  for (z in (q+1):(length(interaction_probabilities))){
                    initial_variance <- initial_variance + 2*interaction_probabilities[q]*interaction_probabilities[z]*vcov_mat[interaction_names[q], interaction_names[z]]                 
                  }
                }
              }
            }
          }
          
          # Store effect and standard error estimates
          results[1,m] <- initial_beta
          if (!is.na(initial_variance)){
            results[2,m] <- sqrt(initial_variance)
          }else{
            results[2,m] <- NA
          }
        }
        
        # Drop elements of the estimates matrix that are NAs
        #print(results)
        keep_vector <- as.vector(which(!is.na(results[1,])))
        results <- as.matrix(results[c(1,2),keep_vector])
        # Reintroduce colnames if results used to be a vector
        if (length(interaction_level_names) == 1){
          colnames(results) <- c(interaction_level_names)  
        }
          
        # Combine estimates + SEs into single matrix - store in list
        estimates[[paste(colnames(acie_matrix),sep="",collapse=":")]] <- results 
      }
    }
  # Create conjoint object
  output <- list()
  class(output) <- c("amce")
  output$formula <- formula
  output$estimates <- estimates
  # Get rid of backquotes
  for (k in 1:length(names(output$estimates))){
    if (grepl("`", names(output$estimates)[k])){
      names(output$estimates)[k] <- gsub("`","",names(output$estimates)[k])
    } 
  }
  # Save attribute levels
  output$attributes <- dimnames(design$J)
  for (k in 1:length(names(output$attributes))){
    if (grepl("`", names(output$attributes)[k])){
      names(output$attributes)[k] <- gsub("`","",names(output$attributes)[k])
    } 
  }
  # Save relevant baselines
  output$baselines <- list()
  for (k in main_effect_vars){
    output$baselines[[k]] <- levels(data[[k]])[1]
  }
  for(j in 1:length(names(output$baselines))){
    if (grepl("`", names(output$baselines)[j])){
      names(output$baselines)[j] <- gsub("`","",names(output$baselines)[j])
    } 
  }
  output$samplesize <- sample_size
  if (!is.null(respondent.id)){
    output$numrespondents <- length(unique(data[[respondent.id]]))
  }else{
    output$numrespondents <- NULL
  }
  return(output)
}


## Print summary of results
summary.amce <- function(object, ...){
  amce_obj <- object
  # Initialize a list to store summary object
  summary_results <- list()
  baselines <- c()
  # Create header of data.frame
  header <- c("Attribute", "Level", "Estimate", "Std. Err", "z value", "Pr(>|z|)", " ")
  names_vec <- c()
  # Create results matrix 
  summary_results[["estimates"]] <- matrix(nrow=length(unlist(amce_obj$estimates))/2, ncol=length(header))
  colnames(summary_results[["estimates"]]) <- header
  summary_results[["estimates"]] <- as.data.frame(summary_results[["estimates"]])
  index <- 1
  
  # Loop over every attribute
  for (effect in names(amce_obj$estimates)){
    # Figure out the baseline level
    
    # If it's an AMCE, easy
    if (!grepl(":", effect)){
      baselines <- c(baselines , amce_obj$baselines[[effect]])
      names_vec <- c(names_vec, effect)
    # If it's not an AMCE, slightly less easy  
    }else{
      variates <- strsplit(effect, ":")[[1]]
      lev_list <- c()
      for (var in variates){
        lev_list <- c(lev_list, amce_obj$baselines[[var]])
      }
      baselines <- c(baselines, paste(lev_list,sep="",collapse=":"))
      names_vec <- c(names_vec, effect)
    }

    # Append results to the estimates dataframe
    for (p in 1:ncol(amce_obj$estimates[[effect]])){
      summary_results[["estimates"]][index,1] <- effect
      summary_results[["estimates"]][index,2] <- colnames(amce_obj$estimates[[effect]])[p]
      summary_results[["estimates"]][index,3] <- amce_obj$estimates[[effect]][1,p]
      summary_results[["estimates"]][index,4] <- amce_obj$estimates[[effect]][2,p]
      zscr <- amce_obj$estimates[[effect]][1,p]/amce_obj$estimates[[effect]][2,p]
      summary_results[["estimates"]][index,5] <- zscr
      pval <- 2*pnorm(-abs(zscr))
      summary_results[["estimates"]][index,6] <- pval
      # Stars!
      if (pval < .001){
        summary_results[["estimates"]][index,7] <- "***"
      }else if (pval < .01){
        summary_results[["estimates"]][index,7] <- "**"
      }else if (pval < .05){
        summary_results[["estimates"]][index,7] <- "*"
      }else{
        summary_results[["estimates"]][index,7] <- ""
      }
      index <- index + 1
    }
  }
    
  # Convert results to data frame
  summary_results[["estimates"]] <- as.data.frame(summary_results[["estimates"]])
  summary_results[["baseline"]] <- data.frame("Attribute" = names_vec, "Level" = baselines)
  # Save sample size
  summary_results[["samplesize"]] <- amce_obj$samplesize
  # If there's a respondent number, add that as well
  if (!is.null(amce_obj$numrespondents)){
    summary_results[["respondents"]] <- amce_obj$numrespondents
  }else{
    summary_results[["respondents"]] <- NULL
  }
  # Set class
  class(summary_results) <- c("summary.amce")
  
  # Return
  return(summary_results)
}

print.summary.amce <- function(x, digits=5, ...){
  summary_result <- x
  cat("----------------\n")
  cat("Estimates:\n")
  cat("----------------\n")
  print(summary_result$estimates, digits=digits, row.names=F)
  cat("---\n")
  cat(paste("Number of Obs. = ", summary_result$samplesize, sep=""))
  cat("\n")
  cat("---\n")
  if (!is.null(summary_result$respondents)){
    cat(paste("Number of Respondents = ", summary_result$respondents, sep=""))
    cat("\n")
    cat("---\n")
  }
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05")
  cat("\n")
  cat("\n")
  cat("----------------\n")
  cat("Baseline Levels:\n")
  cat("----------------\n")
  print(summary_result$baseline, row.names=F)

}

plot.amce <- function(x, main="", xlab="Change in E[Y]", ci=.95, colors=NULL, xlim=NULL, breaks=NULL, labels=NULL, attribute_names = NULL, level_names = NULL, label.baseline = TRUE, text.size=11, text.color = "black", point.size = .6, ...){
  # You need ggplot2
  amce_obj <- x
  ylim <- xlim

  # Make R CMD check happy
  pe <- NULL
  group <- NULL
  lower <- NULL
  upper <- NULL
  
  # Extract raw attribute names from the amce_obj$estimates object
  raw_attributes <- names(amce_obj$estimates)
  # Extract raw levels 
  raw_levels <- list()
  for(m in names(amce_obj$estimates)){
    raw_levels[[m]] <- colnames(amce_obj$estimates[[m]])
  }
  
  # Determine baseline level for each effect estimate in raw_levels and append to beginning of each vector in raw_levels
  for (effect in names(raw_levels)){
    # If it's an AMCE
    if (!grepl(":",effect)){
      raw_levels[[effect]] <- c(amce_obj$baselines[[effect]], raw_levels[[effect]])
    # If it's an ACIE  
    }else{
      effect_elements <- strsplit(effect, ":")[[1]]
      baseline_interactions <- c()
      for (elem in effect_elements){
        baseline_interactions <- c(baseline_interactions, amce_obj$baselines[[elem]])
      }
      interaction_str <- paste(baseline_interactions,sep="",collapse=":")
      raw_levels[[effect]] <- c(interaction_str, raw_levels[[effect]])
    }
  }
  
  # Sanity check attribute_names and level_names against amce_obj$attributes
  if (is.null(attribute_names)){
    attribute_names <- raw_attributes
    attribute_names <- sub(":", " X ", attribute_names)
  }else{
    if (length(attribute_names) != length(raw_attributes)){
      cat(paste("Error: The number of elements in attribute_names ", length(attribute_names), " does not match the number of attributes in amce object for which estimates were obtained", length(amce_obj$attributes), "\n", sep=""))
      cat("Defaulting attribute_names to attribute names in AMCE object\n")
      attribute_names <- raw_attributes
    }
  }
  
  # level_names
  if (is.null(level_names)){
    level_names <- raw_levels
  }else{
    for (name in names(raw_levels)){
      if (length(level_names[[name]]) != length(raw_levels[[name]])){
        cat(paste("Error: level_names lengths do not match levels for attribute ", name, "\n",sep=""))
        cat(paste("Defaulting level_names for attribute ", name, " to level names in AMCE object", "\n",sep=""))
        level_names[[name]] <- raw_levels[[name]]
      }
    }
  }
  
  ## Fix non-unique attribute and level labels
  attr_uniques <- unique(attribute_names)
  for (i in 1:length(attr_uniques)){
    attr <- attr_uniques[i]
    num <- 0
    for (j in 1:length(attribute_names)){
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
  for (i in 1:length(uniques)){
    level <- uniques[i]
    num <- 0
    for (attr in names(level_names)){
        for (j in 1:length(level_names[[attr]])){
          if (level_names[[attr]][j] == level){
            leading_spaces <- paste(rep(" ", num), sep="",collapse="")
            level_names[[attr]][j] <- paste(level, leading_spaces, sep="")
            num <- num + 1
          }
        }
    }
  }
  
  # If label.baseline = T, label the first level of the attribute as the baseline
  if (label.baseline){
    for (attr_ind in names(level_names)){
      level_names[[attr_ind]][1] <- paste("(Baseline = ", level_names[[attr_ind]][1], ")", sep="")
    }
  }
  
  # Convert ci to z-score
  if (ci < 1 & ci > 0){
    zscr <- qnorm(1- ((1-ci)/2))
  }else{
    cat("Invalid confidence interval -- Defaulting to 95%")
    zscr <- qnorm(1- ((1-.95)/2))
  }

  # Compile estimates into plottable objects
  d <- data.frame(pe=c(), se=c(), upper=c(), lower=c(), var=c(), group=c())
  
  # Iterate over the amce object
  for (i in 1:length(amce_obj$estimates)){
    attr_name = names(amce_obj$estimates)[i]
    print_attr_name <- attribute_names[i]
    d_head <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste(print_attr_name, ":", sep=""), group="<NA>")
    d <- rbind(d, d_head)
    for (j in 1:length(level_names[[attr_name]])){
      # If it's the baseline level
      if (j == 1){
        level_name <- level_names[[attr_name]][j]
        d_lev <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste("   ", level_name,sep=""), group=print_attr_name)
        d <- rbind(d, d_lev)
      # if it's an effect level
      }else{
        level <- colnames(amce_obj$estimates[[attr_name]])[j-1]
        level_name <- level_names[[attr_name]][j]
        val_pe <- amce_obj$estimates[[attr_name]][1,level]
        val_se <- amce_obj$estimates[[attr_name]][2,level]
        upper_bnd <- val_pe + zscr*val_se
        lower_bnd <- val_pe - zscr*val_se
        d_lev <- data.frame(pe=val_pe, se=val_se, upper=upper_bnd, lower=lower_bnd, var=paste("   ", level_name,sep=""), group=print_attr_name)
        d <- rbind(d, d_lev)
      }
    }
    
  }
  # Set Y bounds
  if (is.null(ylim)){
    max_upper <- max(d$upper, na.rm=T) + .05
    min_lower <- min(d$lower, na.rm=T) - .05
    ylim <- c(min_lower, max_upper)
    d[is.na(d)] <- max_upper + 100
  }else{
    d[is.na(d)] <- max(ylim) + 100
  }

  # Make group factors <NA> actually NA
  d$group[d$group == "<NA>"] <- NA
  
  # Reverse factor ordering
  d$var <- factor(d$var,levels=unique(d$var)[length(d$var):1])
  # Output to ggplot
  p = ggplot(d, aes(y=pe, x=var, colour=group))
  p = p + coord_flip(ylim = ylim)  
  p = p + geom_hline(yintercept = 0,size=.5,colour="black",linetype="dotted") 
  p = p + geom_pointrange(aes(ymin=lower,ymax=upper),position="dodge",size=point.size)
  # If breaks and labels Null, use default
  if (is.null(breaks) & is.null(labels)){
    p = p + scale_y_continuous(name=xlab)
  }else if (is.null(breaks) & !is.null(labels)){
    p = p + scale_y_continuous(name=xlab, labels=labels)
  }else if (!is.null(breaks) & is.null(labels)){
    p = p + scale_y_continuous(name=xlab, breaks=breaks)
  }else if (!is.null(breaks) & !is.null(labels)){
    p = p + scale_y_continuous(name=xlab, breaks=breaks, labels=labels)
  }
  
  p = p + scale_x_discrete(name="")
  # If there's a title,add it
  if (!is.null(main)){
    if (main != ""){
      p = p + ggtitle(main)
    }
  }
  # If no colors specified, use default
  if (is.null(colors)){
    p = p + scale_colour_discrete(" ") 
  }else if (is.vector(colors)){
    # Make manual palette
    cPal <- rep(colors, ceiling(length(unique(d$group))/length(colors)))
    # Use manual palette
    p = p + scale_colour_manual(values=cPal)
  }else{
    cat("Error: 'colors' must be a vector. Using default colors\n")
    p = p + scale_colour_discrete(" ") 
  }
  
  
  # colour scheme
  theme_bw1 <- function(base_size = text.size, base_family = "") {
    
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
      theme(
        axis.text.x =       element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),
        axis.text.y =       element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), # changes position of X axis text
        
        axis.ticks =        element_line(colour = "grey50"),
        axis.title.y =      element_text(size = base_size,angle=90,vjust=.01,hjust=.1),
        plot.title =             element_text(face = "bold"),
        legend.position = "none"
      )
  }
  p = p + theme_bw1()
  print(p)
  
}

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
## make design

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

         row_vec <- as.vector(unlist(c(r,respondent_index[r], k, j, profile_levels, selec, unit_cov)))

         header <- as.vector(unlist(c("respondentIndex", "respondent","task","profile",attribute_vector, "selected", covariates)))

         names(row_vec) <- header

         if (is.null(out_data_dataset)){
           out_data_dataset <- as.data.frame(t(row_vec), stringsAsFactors = F)
         }else{
           out_data_dataset <- rbind(out_data_dataset, as.data.frame(t(row_vec), stringsAsFactors = F))
         }

       }
    }
    }
  }
  # Do some post-processing
  out_data_dataset$respondentIndex <- as.integer(out_data_dataset$respondentIndex)
  out_data_dataset$task <- as.integer(out_data_dataset$task)
  out_data_dataset$profile <- as.integer(out_data_dataset$profile)
  
  # Return dataset
  return(out_data_dataset)
}




