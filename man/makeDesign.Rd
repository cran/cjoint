\name{makeDesign}
\alias{makeDesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create conjoint design
}
\description{
Generates a "conjointDesign" object to be passed to \code{\link{amce}}. Average Marginal Component Effects (AMCEs) are defined relative to the distribution of potential choice profiles. When the probability of each profile being presented to a respondent is not constant (e.g. some profiles are restricted from appearing), some attributes will not be independent and simple difference-in-means estimates will not be unbiased for the AMCE. This function allows users to specify non-uniform profile assignment schemes to be used by \code{\link{amce}}.
}
\usage{
makeDesign(type="file", J=NULL, filename=NULL, attribute.levels = NULL,
            constraints=NULL, level.probs=NULL, tol=1e-14)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
  A character string - either "file", "constraints", or "array". Each option requires a different set of additional arguments. If "file", the user must specify a filename of a ".dat" file exported via the Conjoint Survey Design Tool. If "constraints," the user must provide a list of attributes and levels in "attribute.levels" along with any constraints in "constraints" and attribute randomization weights in "level.probs". If "array," the user must pass an array to J.
}
  \item{J}{
  If type="array", \code{makeDesign} requires a d-dimensional array to be passed to J where d is the number of attributes in a single profile. . Each dimension should have a number of indices corresponding to the number of levels for that attribute. Attribute and level names are taken from the \code{dimnames} of the array (level names are taken from the list elements, attribute names are taken from the names of the elements). Each cell of the array is the joint probability of assigning that particular profile (combination of attribute-levels).
}
  \item{filename}{
  If type="file", this is a character string giving the name of a design file exported from the Conjoint Survey Design Tool using the "Export design to R" option.
}
  \item{attribute.levels}{
  If type="constraints", attribute.levels is a required argument. This takes a named list with each element containing a character vector with the level names of a single attribute. The names of the elements should be the attribute names.
}
  \item{constraints}{
  If type = "constraints", this is an optional argument. A list of lists. Each element list is a "restriction" and describes a set of restricted profiles (attribute-level groupings that cannot be displayed to a respondent) with character vectors containing the attribute-levels that cannot be displayed together. If your design prevents certain profiles from being seen by respondents (for example, if the particular profile would be highly impluasible), this argument allows you to specify those restrictions. Each element of a restriction is a character vector of levels for a particular attribute. That element must have the name of one of the attributes. For each specified restriction, \code{makeDesign} will treat any profile containing any combination of the restricted levels as having zero probability of being displayed. For example, if a restriction contains two elements: A = c("1","2") and B = c("4", "5"), then any profiles containing levels 1 or 2 of attribute A and levels 4 or 5 of attribute B will restricted from the design.
}
  \item{level.probs}{
  If type = "constraints", this is an optional argument. A named list containing numeric vector elements. Each list element should have the name of one of the attributes passed in attribute.levels and the numeric vector it contains should have the names of the corresponding levels. Each element of the numeric vector specifies the marginal probability of that level appearing in a profile before any restricted profiles are eliminated. Each of the vectors should sum to 1.  If NULL, then uniform randomization is assumed (excluding constrained profiles).
}
  \item{tol}{
  A very small non-zero number used in the routine for determining which attributes are dependent (the distribution of one attribute depends on the value of another attribute). Two attributes are independent if the distribution of one attribute does not change conditional on the other - this is calculated automatically through the J matrix Do not change unless you are having issues with \code{makeDesign} incorrectly labeling independent attributes as dependent (in which case, lower the value) or incorrectly labeling dependent attributes as independent (in which case, raise the value).
}
}
\details{
See Hainmueller et. al. (2014) for details on the AMCE estimator in the presence of conditionally independent attributes. If you have no restricted profiles in your design and uniform randomization across profiles is used, then creating a design object is not necessary (the regression coefficients will give you AMCEs). However, if some profile combinations are restricted or some profiles are made to appear more often than others, it is necessary to specify them using \code{makeDesign} in order to obtain the correct AMCE estimates using \code{\link{amce}}.

You can create designs with assignment probabilities different from those of your original experiment. This allows you to estimate AMCEs for another tailored user-defined distribution of attribute-levels, such as the empirical distribution of the attributes in the population of interest. Note that the distribution of the original experiment must have a support that covers the entire support of the custom-made distribution. See Hainmueller et. al. (2014) for more details about the AMCE estimator and its relationship to the attribute-level distribution.
}
\value{
Returns an object of class "conjointDesign" containing:
\item{J }{
  A d-dimensional array where d is the number of attributes in a single profile. Each dimension has a number of indices corresponding to the number of levels for that attribute. Each entry is the probability of that profile being displayed to a survey respondent. The \code{dimnames} of the array contain the attribute and level names.
}
\item{dependence }{
  A list of character vectors denoting which attributes are marginally dependent - used in the \code{\link{amce}} estimation routine.
}
}
\references{
Strezhnev, A., Hainmueller, J., Hopkins, D., and Yamamoto, T. (2014) Conjoint Survey Design Tool. \url{https://github.com/astrezhnev/conjointsdt/}

Hainmueller, J., Hopkins, D., and Yamamoto T. (2014) Causal Inference in Conjoint Analysis: Understanding Multi-Dimensional Choices via Stated Preference Experiments. Political Analysis 22(1):1-30
}

\seealso{
\code{\link{amce}}
}
\examples{
data("immigrationconjoint")
\dontrun{
## You can load a design from a .dat file from the Conjoint SDT
immigrationdesign <- makeDesign(type="file", filename="immigrant.dat")
}

## Or you can construct the conjoint design manually in in R
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
# Constraints on Education and Job attributes
# Cannot have "doctor" with "no formal" - "high school" education
constraint_list[[1]] <- list()
constraint_list[[1]][["Education"]] <- c("no formal","4th grade",
                                          "8th grade","high school")
constraint_list[[1]][["Job"]] <- c("financial analyst",
                                  "computer programmer","research scientist",
                                  "doctor")
# Constraints on Reason for Application/Country of Origin
constraint_list[[2]] <- list()
constraint_list[[2]][["Reason for Application"]] <- c("escape persecution")
constraint_list[[2]][["Country of Origin"]] <- c("Germany","France",
                                                "Mexico","Philippines",
                                                "Poland","India")

immigrationdesign <- makeDesign(type='constraints', attribute.levels=attribute_list,
                                constraints=constraint_list)

## You can set your own user-defined distribution over the attribute-levels
marginal_weights <- list() # Uniform for all except education which is bimodal
marginal_weights[["Education"]] <- c(1/20, 1/20, 2/10, 2/10, 1/10, 2/10, 2/10)
marginal_weights[["Gender"]] <- rep(1/length(attribute_list[["Gender"]]),
                                length(attribute_list[["Gender"]]))
marginal_weights[["Country of Origin"]] <- rep(1/length(
                                      attribute_list[["Country of Origin"]]),
                                      length(attribute_list[["Country of Origin"]]))
marginal_weights[["Reason for Application"]] <- rep(1/length(
                              attribute_list[["Reason for Application"]]),
                              length(attribute_list
                              [["Reason for Application"]]))
marginal_weights[["Job"]] <- rep(1/length(attribute_list[["Job"]]),
                              length(attribute_list[["Job"]]))
marginal_weights[["Job Experience"]] <- rep(1/length(attribute_list[["Job Experience"]]),
                                      length(attribute_list[["Job Experience"]]))
marginal_weights[["Job Plans"]] <- rep(1/length(attribute_list[["Job Plans"]]),
                                      length(attribute_list[["Job Plans"]]))
marginal_weights[["Prior Entry"]] <- rep(1/length(attribute_list[["Prior Entry"]]),
                                      length(attribute_list[["Prior Entry"]]))
marginal_weights[["Language Skills"]] <- rep(1/length(attribute_list[["Language Skills"]]),
                                      length(attribute_list[["Language Skills"]]))

\dontrun{
immigrationdesign_reweight <- makeDesign(type='constraints', attribute.levels=attribute_list,
                                constraints=constraint_list, level.probs=marginal_weights)

## Note that estimated AMCEs can depend on the randomization distribution

results <- amce(Chosen_Immigrant ~  Gender + Education + Job, data=immigrationconjoint,
                cluster=TRUE, respondent.id="CaseID", design=immigrationdesign)
summary(results)

results_wt <- amce(Chosen_Immigrant ~  Gender + Education + Job, data=immigrationconjoint,
                cluster=TRUE, respondent.id="CaseID", design=immigrationdesign_reweight)
summary(results_wt)
}
}

