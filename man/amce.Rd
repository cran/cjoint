\name{amce}
\alias{amce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating Causal Effects in Conjoint Experiments
}
\description{
This function takes a dataset and a conjoint design and returns Average Marginal Component Effects (AMCEs) and Average Component Interaction Effects (ACIE) for the attributes specified in the formula. By default, this function assumes uniform randomization of attribute levels and no profile restrictions. If your design incorporates weighted randomization or restrictions on displayable profiles, first generate a design object using \code{makeDesign}. Interactions with respondent-level characteristics are handled by identifying relevant variables as respondent-varying.}
\usage{
amce(formula, data, design = "uniform",
	      respondent.varying = NULL, subset = NULL,
	      respondent.id = NULL, cluster = TRUE, na.ignore=FALSE,
	      weights = NULL, baselines = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  A \code{formula} object specifying the name of the outcome variable on the left-hand side and the attributes for which effects are to be estimated on the right-hand side. RHS attributes should be separated by \code{+} signs. Interaction effects can be specified using standard interaction syntax - joining attribute names using either \code{:} or \code{*}.  However using the \code{:} syntax will produce the same results as \code{*} since missing base terms are automatically added to the formula. For example \code{Y ~ X1 + X2} will return AMCEs for X1 and X2. \code{Y ~ X1 + X2 + X1:X2} will return AMCEs for X1 and X2 along with an ACIE for X1/X2. \code{Y ~ X1*X2} and \code{Y ~ X1:X2} will produce identical results to \code{Y ~ X1 + X2 + X1:X2}. Note that you can place backticks around a variable name containing spaces in order to have \code{\link{formula}} interpret it as a single variable name. Any respondent characteristics must be designated as such in \code{redpondent.varying}. 
}
  \item{data}{
  A dataframe containing the outcome variable, attributes, respondent identifiers, respondent covariate data and sampling weights from a conjoint experiment.
}
  \item{design}{
  Either the character string \code{"uniform"} or a \code{conjointDesign} object created by the \code{\link{makeDesign}} function. If a \code{conjointDesign} is not passed, the function will assume all attribute levels have an equal probability of being presented to a respondent and that no profiles are restricted. Defaults to \code{"uniform"}.
}
\item{respondent.varying}{
A vector of character strings giving the names of any respondent-varying characteristics being interacted with AMCEs or ACIEs in the \code{formula}.
}
  \item{subset}{
  A logical vector with length \code{nrow(data)} denoting which rows in \code{data} should be included in estimation. This can for example be used to subset the data along respondent-level covariates. Defaults to \code{NULL}.
}
  \item{respondent.id}{
  A character string indicating the column of \code{data} containing a unique identifier for each respondent. Defaults to \code{NULL}.
}
  \item{cluster}{
  A logical indicating whether estimated standard errors should be clustered on \code{respondent.id}. Defaults to \code{TRUE}.
}
  \item{na.ignore}{
  A logical indicating whether the function should ignore missing rows in \code{data}. If \code{FALSE}, amce() will raise an error if there are rows with missing values. Defaults to \code{FALSE}.
  }
\item{weights}{
 A character string giving the name of the column in the data containing any survey weights. See documentation for \code{survey} package for more information.
}
\item{baselines}{
Manually adjust the baselines of select factor variables (either attributes or respondent varying) by supplying a list. Names of list entries should correspond with variable names. The content of each entry should be a character string giving the new baseline.
}
}

\value{
An object of class "amce" containing:
\item{attributes}{A list containing the names of attributes.}
\item{baselines}{Baseline levels for each attribute in \code{estimates}. Baselines determined using the first element of \code{levels()}. If a different baseline level is desired for an attribute, use the \code{relevel()} function on the variable prior to calling the \code{amce()} routine or supply an alternative baseline in \code{baselines} argument.}
\item{continuous}{List of quantiles for any non-factor variables, whether attributes or respondent varying.}
\item{data}{The original data.}
\item{estimates}{A list containing AMCE and ACIE estimates for each attribute in \code{formula}. Each element of \code{estimates} corresponds to a single attribute or interaction.}
\item{formula}{The \code{formula} passed to the \code{amce()} routine.}
\item{samplesize_prof}{The number of valid profiles (rows) in the dataset}
\item{user.names}{A vector with the original user supplied names for any attributes. These may differ from the attribute names in \code{estimates} if the original names contain spaces.}
\item{vcov.prof}{The modified variance-covariance matrix for AMCE and ACIE estimates. Incorporates cluster corrections as well as attribute dependencies. Profile varying attributes only.}
\item{numrespondents}{The number of respondents in the dataset (if \code{respondent.id} is not \code{NULL}). }
\item{respondent.varying}{Names of respondent-varying variables, if any.}
\item{cond.formula}{The formula used for calculating estimates conditional on respondent varying characteristics. Only returned when respondent-varying characteristics are present.}
\item{cond.estimates}{A list containing estimated effects of respondent-varying characteristics conditional on attribute values. Each element of \code{cond.estimates} corresponds to a single attribute or interaction. Only returned when respondent-varying characteristics are present. To obtain AMCE and ACIE estimates conditional on the values of the respondent-varying characteristics see \code{\link{summary.amce}}}
\item{samplesize_full}{The number of valid profiles (rows) in the dataset when respondent varying characteristics are included. Only returned when respondent-varying characteristics are present. }
\item{vcov.resp}{The modified variance-covariance matrix for effect estimates conditional on respondent-varying characteristics, where dependent relationships have been incorporated into variances and covariances. Only returned when respondent-varying characteristics are present. }
}
\references{
Hainmueller, J., Hopkins, D., and Yamamoto T. (2014) Causal Inference in Conjoint Analysis: Understanding Multi-Dimensional Choices via Stated Preference Experiments. Political Analysis 22(1):1-30
}
\seealso{
\code{\link{summary.amce}} for summaries and \code{\link{plot.amce}} for generating a coefficient plot using \code{ggplot2}.

\code{\link{makeDesign}} to create \code{conjointDesign} objects.
}
\examples{
# Immigration Choice Conjoint Experiment Data from Hainmueller et. al. (2014).
data("immigrationconjoint")
data("immigrationdesign")

# Run AMCE estimator using all attributes in the design
results <- amce(Chosen_Immigrant ~  Gender + Education + `Language Skills` +
                `Country of Origin` + Job + `Job Experience` + `Job Plans` +
                `Reason for Application` + `Prior Entry`, data=immigrationconjoint,
                cluster=TRUE, respondent.id="CaseID", design=immigrationdesign)
# Print summary
summary(results)

\dontrun{
# Run AMCE estimator using all attributes in the design with interactions
interaction_results <- amce(Chosen_Immigrant ~  Gender + Education + `Language Skills` +
                `Country of Origin` + Job + `Job Experience` + `Job Plans` +
                `Reason for Application` + `Prior Entry` + Education:`Language Skills` +
		Job: `Job Experience` + `Job Plans`:`Reason for Application`,
		data=immigrationconjoint, cluster=TRUE, respondent.id="CaseID",
		design=immigrationdesign)
# Print summary
summary(interaction_results)

# create weights in data
weights <- runif(nrow(immigrationconjoint))
immigrationconjoint$weights <- weights
# Run AMCE estimator using weights
results <- amce(Chosen_Immigrant ~  Gender + Education + `Language Skills` +
                `Country of Origin` + Job + `Job Experience` + `Job Plans` +
                `Reason for Application` + `Prior Entry`, data=immigrationconjoint,
                cluster=TRUE, respondent.id="CaseID", design=immigrationdesign,
		weights = "weights")
# Print summary
summary(results)

# Include a respondent-varying interaction
results <- amce(Chosen_Immigrant ~ Gender + Education + Job +
	   	ethnocentrism:Job + Education:Job,
		data=immigrationconjoint, na.ignore = TRUE,
		cluster=FALSE,design=immigrationdesign,
		respondent.varying = "ethnocentrism")
# Print summary
summary(results)

# Change the baseline for "Education"
baselines <- list()
baselines$Education <- "graduate degree"

results <- amce(Chosen_Immigrant ~ Gender + Education + Job +
		 Education:Job, data=immigrationconjoint,
		 cluster=FALSE,design=immigrationdesign,
		 baselines=baselines)
# Print summary
summary(results)
}

}

