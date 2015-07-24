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
\method{summary}{amce}(object, interaction.values=NULL, show.all=FALSE, ...)

\method{print}{summary.amce}(x, digits=5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class "amce", a result of a call to \code{\link{amce}}.
  }
  \item{interaction.values}{
   An optional list containing a vector at which interaction effects will be calculated in the case of interactions with respondent-varying characteristics. If the "amce" object contains respondent varying characteristics, by default interaction effects will be reported at quantiles in the case of a continuous variable and levels in the case of a factor. Names of list entries must correspond to variable names. Multiple respondent-varying characteristics cannot currently be varied at the same time within the same calcuation.
  }
  \item{show.all}{
   A logical asking if attributes that are not part of respondent-varying interactions should also have their estimates reported at different levels of the respondent-varying variable.
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
\item{estimates }{A dataframe containing AMCE estimates and standard errors. Each row corresponds to a single attribute-level effect. }
\item{acie }{A dataframe containing ACIE estimates and standard errors. Each row corresponds to a single attribute-level effect. Separate data frames are produced for each level of interactions with respondent-varying characteristics.}
\item{baselines.amce }{Baseline levels for each attribute relative to which the AMCEs are calculated.}
\item{baselines.acie }{Baseline levels for each attribute relative to which the ACIEs are calculated.}
\item{baselines.resp }{Baseline levels for each respondent-varying factor variable, if any, relative to which interactions are calculated.}
\item{samplesize }{The number of profiles (rows) in the dataset}
\item{numrespondents }{The number of respondents in the dataset (if a \code{respondent.id} argument was passed to \code{amce})}
}

\seealso{
The estimation function \code{\link{amce}}.
}
