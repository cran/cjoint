\name{immigrationdesign}
\alias{immigrationdesign}
\docType{data}
\title{
Conjoint Design for the Immigration Experiment in Hainmueller et. al. (2014)
}
\description{
A "conjointDesign" object for the randomization scheme used in the immigration conjoint experiment in "Causal Inference in Conjoint Analysis." See \code{\link{immigrationconjoint}} for the accompanying dataset. 
}
\usage{data("immigrationdesign")}
\format{
A "conjointDesign" object. See \code{\link{makeDesign}} for more information about the structure.

}

\source{
Hainmueller, J., Hopkins, D., and Yamamoto T. (2014) Causal Inference in Conjoint Analysis: Understanding Multi-Dimensional Choices via Stated Preference Experiments. Political Analysis 22(1):1-30
}

\examples{
# Loads the immigrationconjoint dataset and the immigrationdesign object
data("immigrationconjoint")
data("immigrationdesign")

# immigrationdesign is passed to the amce() function through the "design" argument.
results <- amce(Chosen_Immigrant ~  Gender + Education + `Language Skills` + 
                `Country of Origin` + Job + `Job Experience` + `Job Plans` + 
                `Reason for Application` + `Prior Entry`, data=immigrationconjoint, 
                cluster=TRUE, respondent.id="CaseID", design=immigrationdesign)

}
\keyword{datasets}
