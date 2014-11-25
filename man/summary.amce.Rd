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
\method{summary}{amce}(object, ...)

\method{print}{summary.amce}(x, digits=5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object of class "amce", a result of a call to \code{\link{amce}}
}
  \item{x}{
  an object of class "summary.amce", a result of a call to \code{summary.amce}
  }
  \item{digits}{
  number of significant digits to use when printing
  }
  \item{\dots}{
  further arguments from other methods
}
}
\value{
The function \code{summary.amce} computes and returns a formatted data frame of effect estimates returned by \code{\link{amce}}
\item{estimates }{A dataframe containing AMCE and ACIE estimates and standard errors. Each row corresponds to a single attribute-level effect.}
\item{baselines }{Baseline levels for each attribute relative to which the AMCEs are calculated.}
\item{samplesize }{The number of profiles (rows) in the dataset}
\item{numrespondents }{The number of respondents in the dataset (if a \code{respondent.id} argument was passed to \code{amce})}
}

\seealso{
The estimation function \code{\link{amce}}.
}
