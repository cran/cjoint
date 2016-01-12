\name{summary.amce}
\alias{summary.amce}
\alias{print.summary.amce}

\title{
Summarizing AMCE estimates
}
\description{
\code{summary} method for class "amce"
}
\usage{
\method{summary}{amce}(object, covariate.values=NULL, ...)

\method{print}{summary.amce}(x, digits=5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class "amce", a result of a call to \code{\link{amce}}.
  }
  \item{covariate.values}{
   An optional list containing a vector at which conditional effects will be calculated in the case of AMCE and ACIE's conditional on respondent-varying characteristics. The class of the values in the vector must match the class of the respondent-varying characteristic in question. If the "amce" object contains respondent varying characteristics, when set to NULL (default) interaction effects will be reported at quantiles in the case of a continuous variable and levels in the case of a factor. Names of list entries must correspond to variable names. Multiple respondent-varying characteristics cannot currently be varied at the same time within the same calcuation. Instead conditional effects will be calculated for one respondent-characteristic at a time while others are held at 0 (for continuous variables) or baseline (for factors).
  }
  \item{x}{
  An object of class "summary.amce", a result of a call to \code{summary.amce}.
  }
  \item{digits}{
  The number of significant digits to use when printing.
  }
  \item{\dots}{
  Further arguments from other methods.
}
}
\value{
The function \code{summary.amce} computes and returns formatted data frames of effect estimates returned by \code{\link{amce}}
\item{amce }{A dataframe containing AMCE estimates and standard errors. Each row corresponds to a single attribute-level effect. }
\item{baselines_amce }{Baseline levels for each attribute relative to which the AMCEs are calculated. }
\item{acie }{A dataframe containing ACIE estimates and standard errors, if any. Each row corresponds to a single attribute-level effect. }
\item{baselines_acie }{Baseline levels for each attribute relative to which the ACIEs are calculated. }
\item{baselines_amce_resp }{Baseline levels for conditional AMCE estimates, if any, relative to which interactions are calculated. }
\item{baselines_acie_resp }{Baseline levels for conditional ACIE estimates, if any, relative to which interactions are calculated. }
\item{samplesize_estimates }{The number of valid profiles (rows) in the dataset when only effects of profile varying attributes are calculated. }
\item{samplesize_resp }{The number of valid profiles (rows) in the dataset when respondent-varying characteristics are incorporated. }
\item{numrespondents }{The number of respondents in the dataset (if a \code{respondent.id} argument was passed to \code{amce})}
\item{table_values_amce}{A dataframe giving the names of additional tables of AMCE estimates conditional on respondent-varying characteristics. A separate table is produced for each level of each respondent varying characteristic. }
\item{table_values_acie}{A dataframe giving the names of tables of ACIE estimates conditional on respondent-varying characteristics. A separate table is produced for each level of each respondent varying characteristic. }
}

\seealso{
The estimation function \code{\link{amce}}.
}

\examples{

\dontrun{

#Results with respondent-varying characteristics
results <-amce(Chosen_Immigrant ~ Gender + Education + Education:ethnocentrism +
	`Country of Origin`+ `Country of Origin`:ethnocentrism + Job +
	Job:ethnocentrism + `Job Experience` + `Job Experience`:ethnocentrism,
	data=immigrationconjoint, design=immigrationdesign,cluster=FALSE,
	respondent.varying="ethnocentrism", na.ignore=TRUE)

#Calculate conditional estimates at user-supplied levels
int.vals<-list()
int.vals[["ethnocentrism"]]<-c(60,77,88,99,45)
summary(results, covariate.values = int.vals)

}
}




















