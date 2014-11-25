\name{read.qualtrics}
\alias{read.qualtrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Data from a Conjoint Qualtrics Experiment
}
\description{
Converts the raw .CSV data file downloaded from an online conjoint experiment run using the Qualtrics survey software into a data frame usable by the \code{\link{amce}} routine. Each row of the Qualtrics .CSV file is a single survey respondent. The rows of the resulting dataframe correspond to individual profile choices. Currently, this function only works for designs with a binary outcome variable (selected/not selected) that requires a single choice among profiles.
}
\usage{
read.qualtrics(filename, responses, covariates = NULL, respondentID = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
  A .CSV file containing responses exported from a Qualtrics survey experiment. The first row contains question and variable identifiers (e.g. V1, Q1). The second row contains variable names and question texts. Subsequent rows contain the answers of each respondent. This .CSV file can be exported through the "Download Data" panel in the Qualtrics web interface. Note that answers should be shown as "coded values" and not as choice text.
}
  \item{responses}{
  A character or integer vector with the identifiers of the CSV columns that contain the conjoint responses. The first element corresponds to the identifier of the first question, the second element corresponds to the second question and so on. The length of this vector must be equal to the number of conjoint tasks answered by each respondent. Identifiers are typically in the form "Q#" - e.g. \code{c("Q1", "Q3", "Q5")} would represent a three question survey where the conjoint questions are Q1, Q3 and Q5. If specified as an integer vector, the elements are the column numbers corresponding to each question response.
}
  \item{covariates}{
  A character vector denoting the column names of any other respondent-level characteristics measured in the survey that should be included in the resulting dataframe.
}
  \item{respondentID}{
  A character string denoting a column containing a unique identifier for each respondent (e.g. an IP address). This identifier will be carried over into the output. If NULL, each respondent will be given an arbitrary identifier in the output dataframe. Leave as NULL if you do not want responses to be linked back to a known respondent identifier.
}
}
\details{
This function currently only works with experiments that generate profiles using .PHP scripts created by the Conjoint Survey Design Tool. It also is only able to handle standard conjoint designs (binary outcome variable/forced choice).\cr
(\url{http://scholar.harvard.edu/astrezhnev/conjoint-survey-design-tool}). \cr

For each respondent in the .CSV file, attribute and level names are stored using the following naming convention:\cr
\bold{Level Name}: F-[task number]-[profile number]-[attribute number]\cr
\bold{Attribute Name}: F-[task number]-[attribute number]\cr
\bold{Example}: F-1-3-2 denotes the level corresponding to Task 1, Profile 3, Attribute 2\cr
F-3-3 denotes the attribute name corresponding to Task 3, Attribute 3
}


\value{
A dataframe in which each row corresponds to a single profile. The column "selected" denotes whether that profile was selected by the respondent. The columns "respondent" and "task" denote the respondent and task numbers to which the profile was assigned. Respondent-level covariates are appended to each row.
}
\references{
Strezhnev, A., Hainmueller, J., Hopkins, D., and Yamamoto, T. (2014) Conjoint Survey Design Tool. \url{http://scholar.harvard.edu/astrezhnev/conjoint-survey-design-tool}
}
\examples{
\dontrun{
# An example file with 5 conjoint tasks per respondent and a single covariate
# You can demonstrate this function's output using the CandidateConjointQualtrics.csv 
# demonstration file in the 'inst' subdirectory
conjoint_data <- read.qualtrics("CandidateConjointQualtrics.csv", 
responses=c("Q2.3", "Q2.7", "Q2.10", "Q2.13", "Q2.16"), 
covariates=c("Q6.6"), respondentID="V1")
}
}


