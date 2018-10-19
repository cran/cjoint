\encoding{UTF-8}
\name{japan2014conjoint}
\alias{japan2014conjoint}
\docType{data}
\title{
Japan 2014 Conjoint Experiment Dataset from Horiuchi et. al. (2014)
}
\description{
A dataset containing the results of a conjoint survey of a representative sample of Japanese adults who were asked to choose which of two parties (with specified policy manifestos) they would support in the 2014 House of Representatives general election. Each row corresponds to a single profile presented to the respondent. See the original replication archive (Horiuchi et. al. 2014) for more detailed descriptions.
}
\usage{data("japan2014conjoint")}
\format{
  A data frame with 20,360 observations on the following 11 variables.
  \describe{
    \item{respondent}{a character vector uniquely identifying the respondent}
    \item{respondentIndex}{a numeric vector uniquely indexing the respondent}
    \item{task}{a numeric vector indexing the task presented to the respondent}
    \item{profile}{a numeric vector indexing the profile presented to the respondent}
    \item{selected}{a numeric vector indicating whether the profile was selected}
    \item{Consumption tax}{a factor indicating the profile position on `Consumption Tax`}
    \item{Consumption tax.rowpos}{a numeric vector indicating which row the attribute `Consumption Tax` appeared in the given profile}
    \item{Employment}{a factor indicating the profile position on `Employment`}
    \item{Employment.rowpos}{a numeric vector indicating which row the attribute `Employment` appeared in the given profile}
    \item{Monetary and fiscal policy}{a factor indicating the profile position on `Monetary and Fiscal Policy`}
    \item{Monetary and fiscal policy.rowpos}{a numeric vector indicating which row the attribute `Monetary and Fiscal Policy` appeared in the given profile}
    \item{Economic growth strategy}{a factor indicating the profile position on `Economic Growth Strategy`}
    \item{Economic growth strategy.rowpos}{a factor indicating which row the attribute `Economic Growth Strategy` appeared in the given profile}
    \item{Nuclear power}{a factor indicating the profile position on `Nuclear Power`}
    \item{Nuclear power.rowpos}{a factor indicating which row the attribute `Nuclear Power` appeared in the given profile}
    \item{TPP}{a factor indicating the profile position on `Trans-Pacific Partnership (TPP)`}
    \item{TPP.rowpos}{a factor indicating which row the attribute `Trans-Pacific Partnership (TPP)` appeared in the given profile}
    \item{Collective self-defense}{a factor indicating which row the attribute `Collective Self-Defense` appeared in the given profile}
    \item{Collective self-defense.rowpos}{a factor indicating which row the attribute `Collective Self-Defense` appeared in the given profile}
    \item{Constitutional revision}{a factor indicating the profile position on `Constitutional Revision`}
    \item{Constitutional revision.rowpos}{a factor indicating which row the attribute `Constitutional Revision` appeared in the given profile}
    \item{National assembly seat reduction}{a factor indicating the profile position on `National Assembly Seat Reduction`}
    \item{National assembly seat reduction.rowpos}{a factor indicating which row the attribute `National Assembly Seat Reduction` appeared in the given profile}
    \item{wgt}{post-stratification weights to map the survey sample to the census population}
  }
}
\source{
Horiuchi, Yusaku; Smith, Daniel M.; Yamamoto, Teppei, 2017, "Replication Data for: Measuring Voters' Multidimensional Policy Preferences with Conjoint Analysis: Application to Japan's 2014 Election", https://doi.org/10.7910/DVN/KUMMUJ, Harvard Dataverse, V1
}

\section{Note}{
The origional survey was conducted in Japanese. To comply with CRAN's policy, we translate the variable names into English.
}

\examples{
data("japan2014conjoint")

# Run AMCE estimator using all attributes and uniform design
results <- amce(selected ~ `Consumption tax` + `Employment` + `Monetary and fiscal policy` + 
                  `Economic growth strategy` + `Nuclear power` + `TPP` + 
                  `Collective self-defense` + `Constitutional revision` + 
                  `National assembly seat reduction`,  data=japan2014conjoint, cluster=TRUE, 
                  respondent.id="respondentIndex", weights="wgt", design="uniform")
}
\keyword{datasets}
