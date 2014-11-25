\name{amce}
\alias{amce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating Causal Effects in Conjoint Experiments
}
\description{
This function takes a dataset and a conjoint design and returns Average Component-specific Causal Effects (AMCEs) and Average Component Interaction Effects (ACIE) for the attributes specified in the formula. By default, this function assumes uniform randomization of attribute levels and no profile restrictions. If your design incorporates weighted randomization or restrictions on displayable profiles, first generate a design object using \code{makeDesign}. Note that interactions with respondent-level characteristics are currently unsupported but are forthcoming in future versions.
}
\usage{
amce(formula, data, design = "uniform", subset = NULL, respondent.id = NULL, 
      cluster = TRUE, na.ignore=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  A \code{formula} object specifying the name of the outcome variable on the left-hand side and the attributes for which effects are to be estimated on the right-hand side. RHS attributes should be separated by \code{+} signs. Interaction effects can be specified using standard interaction syntax - joining attribute names using either \code{:} or \code{*}. For example \code{Y ~ X1 + X2} will return AMCEs for X1 and X2. \code{Y ~ X1 + X2 + X1:X2} will return AMCEs for X1 an X2 along with an ACIE for X1/X2. Note that you can place backticks around a variable name containing spaces in order to have \code{\link{formula}} interpret it as a single variable name.
}
  \item{data}{
  A dataframe containing the outcome variable, attributes and respondent covariate data from a conjoint experiment.
}
  \item{design}{
  Either the character string \code{"uniform"} or a \code{conjointDesign} object created by the \code{\link{makeDesign}} function. If a \code{conjointDesign} is not passed, the function will assume all attribute levels have an equal probability of being presented to a respondent and that no profiles are restricted. Defaults to \code{"uniform"}.
}
  \item{subset}{
  A logical vector with length \code{nrow(data)} denoting which rows in \code{data} should be included in estimation. This can be used to subset the data along respondent-level covariates and estimate separate effects for diferent sub-groups. Defaults to \code{NULL}.
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
}

\value{
An object of class "amce" containing:
\item{formula }{The \code{formula} passed to the \code{amce()} routine}
\item{estimates }{A list containing AMCE and ACIE estimates for each attribute in \code{formula}. Each element of \code{estimates} corresponds to a single attribute or interaction.}
\item{attributes }{A list containing the names }
\item{baselines }{Baseline levels for each attribute in \code{estimates}. Baselines determined using the first element of \code{levels()}. Use the \code{relevel()} function to change the baseline levels for attributes.}
\item{samplesize }{The number of profiles (rows) in the dataset}
\item{numrespondents }{The number of respondents in the dataset (if \code{respondent.id} is not \code{NULL}) }

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

}

