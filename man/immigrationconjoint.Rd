\name{immigrationconjoint}
\alias{immigrationconjoint}
\docType{data}
\title{
Immigration Conjoint Experiment Dataset from Hainmueller et. al. (2014)
}
\description{
A dataset containing the results of a conjoint survey of a representative sample of American adults who were asked to choose which hypothetical immigrants they think should be admitted into the United States. Each row corresponds to a single profile presented to the respondent.
}
\usage{data("immigrationconjoint")}
\format{
  A data frame with 13,960 observations on the following 16 variables.
  \describe{
    \item{\code{CaseID}}{a numeric vector indicating the respondent to which the particular profile corresponds}
    \item{\code{contest_no}}{a numeric vector indicating the number of the task to which the profile corresponds}
    \item{\code{Education}}{a factor with levels \code{no formal}, \code{4th grade}, \code{8th grade}, \code{high school}, \code{two-year college}, \code{college degree}, \code{graduate degree}}
    \item{\code{Gender}}{a factor with levels \code{female}, \code{male}}
    \item{\samp{Country of Origin}}{a factor with levels \code{India}, \code{Germany}, \code{France}, \code{Mexico}, \code{Philippines}, \code{Poland}, \code{China}, \code{Sudan}, \code{Somalia}, \code{Iraq}}
    \item{\samp{Reason for Application}}{a factor with levels \code{reunite with family}, \code{seek better job}, \code{escape persecution}}
    \item{\code{Job}}{a factor with levels \code{janitor}, \code{waiter}, \code{child care provider}, \code{gardener}, \code{financial analyst}, \code{construction worker}, \code{teacher}, \code{computer programmer}, \code{nurse}, \code{research scientist}, \code{doctor}}
    \item{\samp{Job Experience}}{a factor with levels \code{none}, \code{1-2 years}, \code{3-5 years}, \code{5+ years}}
    \item{\samp{Job Plans}}{a factor with levels \code{will look for work}, \code{contract with employer}, \code{interviews with employer}, \code{no plans to look for work}}
    \item{\samp{Prior Entry}}{a factor with levels \code{never}, \code{once as tourist}, \code{many times as tourist}, \code{six months with family}, \code{once w/o authorization}}
    \item{\samp{Language Skills}}{a factor with levels \code{fluent English}, \code{broken English}, \code{tried English but unable}, \code{used interpreter}}
    \item{\code{Chosen_Immigrant}}{a numeric vector denoting whether the immigrant profile was selected}
    \item{\code{ethnocentrism}}{a numeric vector}
    \item{\code{profile}}{a numeric vector giving the profile number}
    \item{\code{LangPos}}{a numeric vector}
    \item{\code{PriorPos}}{a numeric vector}
  }
}

\source{
Hainmueller, J., Hopkins, D., and Yamamoto T. (2014) Causal Inference in Conjoint Analysis: Understanding Multi-Dimensional Choices via Stated Preference Experiments. Political Analysis 22(1):1-30
}

\examples{
\dontrun{
data("immigrationconjoint")
data("immigrationdesign")

# Run AMCE estimator using all attributes in the design
results <- amce(Chosen_Immigrant ~  Gender + Education + `Language Skills` + 
                `Country of Origin` + Job + `Job Experience` + `Job Plans` + 
                `Reason for Application` + `Prior Entry`, data=immigrationconjoint, 
                cluster=TRUE, respondent.id="CaseID", design=immigrationdesign)
}
}
\keyword{datasets}
