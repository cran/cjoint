\name{read.with.qualtRics}
\alias{read.with.qualtRics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Read Data from a Conjoint Qualtrics Experiment
}
\description{
	Converts the raw dataframe downloaded with R Package \code{"qualtRics"} from an online conjoint experiment into a data frame usable by the \code{\link{amce}} routine. Each row of the raw dataframe is a single survey respondent. The rows of the resulting dataframe correspond to individual profile choices. Currently, the function supports various types of outcomes in a conjoint experiment, with details explained below.
}
\usage{
read.with.qualtRics(filename, responses = NULL, covariates = NULL,
	respondentID = NULL, letter = "F",  new.format=FALSE, ranks = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
  A dataframe containing responses download from a Qualtrics survey experiment with R package \code{"qualtRics"}. Each row contains the answers of each respondent. Note that answers should be shown as "coded values" and not as choice text.
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
\item{new.format}{
  Indicator for whether the .csv file is from the new Qualtrics export format with three title rows (TRUE) or from the old format (FALSE) with two title rows. Defaults to FALSE.
}
\item{ranks}{
	An integer vector with the identifiers of the CSV columns that contain the conjoint rankings or ratings. 
}
\item{letter}{
	The beginning letter used in the naming convention of levels and attributes.
}
}
\details{
This function currently only works with experiments that generate profiles using .PHP scripts created by the Conjoint Survey Design Tool. It also is only able to handle standard conjoint designs (binary outcome variable/forced choice).\cr
(\url{http://scholar.harvard.edu/astrezhnev/conjoint-survey-design-tool}). \cr

For each respondent in the dataframe, attribute and level names are stored using the following naming convention:\cr
\bold{Level Name}: F-[task number]-[profile number]-[attribute number]\cr
\bold{Attribute Name}: F-[task number]-[attribute number]\cr
\bold{Example}: F-1-3-2 denotes the level corresponding to Task 1, Profile 3, Attribute 2\cr
F-3-3 denotes the attribute name corresponding to Task 3, Attribute 3

\bold{East Asian Language Support:} The read.qualtrics function relies on the read.csv function in R-core. The read.csv function only works well for some of character encoding, but not others, for East Asian languages in some Operation 
systems. In Windows, .csv files containing East Asian languages such as Chinese or Japanese should be stored in the ANSI encoding rather than UTF-8 for the read.csv function to work.

\bold{Different types of responses}: This function supports various types of responses commonly used in conjoint analyses. Here are some illustrations on some typical types. 

1. The respondent is asked to fill in the profile she prefers the most, and her choice is restored in one response variable. In this case, set the argument responses=the response variable.

2. The respondent is asked to give each profile a rank within each task, and her ranks for each profile within each task are restored in J response variables, suppose there are J profiles within each task. In this case, set the argument ranks as a vector restoring the variables names of these responses variables, in the order of rank of profile 1 in task 1, rank of profile 2 in task 1, ... rank of profile J in task 1, rank of profile 1 in task 2, rank of profile 2 in task 2, ... rank of profile J in task 2, ..., rank of profile J in task K.

3. The respondent is asked to rate each profile within each task, and her ratings for each profile within each task are restored in J response variables, suppose there are J profiles within each task. In this case, set the argument ratings as a vector restoring the variables names of these responses variables, in the order of rating of profile 1 in task 1, rating of profile 2 in task 1, ... rating of profile J in task 1, rating of profile 1 in task 2, rating of profile 2 in task 2, ... rating of profile J in task 2, ..., rating of profile J in task K.

4. The respondent is asked to select the top L profiles she prefers within each task. L colud be smaller than the number of profiles J available within each task. Her choices are recorded with L variables for each task, indicating the first choice among J profiles in task 1, the second choice among J profiles in task 1, ..., the L-th choice among J profiles in task 1, the first choice among J profiles in task 2, the second choice among J profiles in task 2, ..., the L-th choice among J profiles in task 2, ..., the L-th choice among J profiles in task K. In this case, the read.qualtrics function should be applied in the following step:

(a) Set the argument response as a vector restoring the variables names of the first choice among J profiles in all tasks, in the order of first choice among J profiles in task 1, first choice among J profiles in task 2, ..., first choice among J profiles in task K. Save the respective data output.

(b) Set the argument response as a vector restoring the variables names of the second choice among J profiles in all tasks, in the order of second choice among J profiles in task 1, second choice among J profiles in task 2, ..., second choice among J profiles in task K. Save the respective data output.

(c) Repeat until you have completed the above steps for all top L choices. And then merge all data according to their respondentID. However, you may need to change the variable name of the responses a little bit for a successful merge.

5. The respondent is asked to select the top L profiles she prefers within each task. L colud be smaller than the number of profiles J available within each task. The questionnaire is designed in such format that the respondent fills in the ranking of top L profiles that she prefers, while leave the ranking of the other J-L profiles blank. In this case, set the argument ranks as a vector restoring the variables names of these J responses variables, in the order of rank of profile 1 in task 1, rank of profile 2 in task 1, ... rank of profile J in task 1, rank of profile 1 in task 2, rank of profile 2 in task 2, ... rank of profile J in task 2, ..., rank of profile J in task K. 

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
# Suppose "CandidateConjointQualtrics" is the dataframe exported by R Package qualtRics
conjoint_data <- read.with.qualtRics(CandidateConjointQualtrics, 
responses=c("Q2.3", "Q2.7", "Q2.10", "Q2.13", "Q2.16"), 
covariates=c("Q6.6"), respondentID="V1")
}
}


