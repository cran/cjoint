\name{plot.amce}
\alias{plot.amce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot AMCE Estimates
}
\description{
\code{plot} method for "amce" objects
}
\usage{
\method{plot}{amce}(x, main = "", xlab = "Change in E[Y]", ci = 0.95, 
                    colors = NULL, xlim = NULL, breaks = NULL, 
                    labels = NULL, attribute_names = NULL, level_names = NULL, 
                    label.baseline = TRUE, text.size = 11, text.color = "black", 
                    point.size = 0.6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  an object of class "amce", a result of a call to \code{\link{amce}}
}
  \item{main}{
  title of the plot.
}
  \item{xlab}{
  Label of the x-axis of the plot (AMCE or ACIE). Default is "Change in Pr(Y=1)"
}
  \item{ci}{
  Size of confidence intervals to plot around point estimates. Must be between 0 and 1. Default is .95
}
  \item{colors}{
  Vector of color names to be used in the plot. The \code{plot} function will alternate between the colors in the vector for each attribute being plotted. If NULL, \code{plot} will use a default \code{ggplot2} color scheme.
}
  \item{xlim}{
 Numeric vector denoting the upper and lower bounds of the x-axis in the plot. If \code{NULL} the plot function will automatically set a range that includes all effect estimates.
}
  \item{breaks}{
  Numeric vector denoting where x-axis tick marks should be placed. If \code{NULL} plot will use ggplot2 defaults
}
  \item{labels}{
  Vector denoting how x-axis tick marks should be labeled. If \code{NULL} plot will use ggplot2 defaults
}
  \item{attribute_names}{
  Character vector of attribute names to be plotted as labels. By default \code{plot.amce} will use the attribute names in the "amce" object passed to it. 
}
  \item{level_names}{
  A list containing character vector elements with names in \code{attribute_names}. Each character vector in the list contains the level names to be plotted as labels beneath the corresponding attribute. By default \code{plot.amce} will use the level names in the "amce" object passed to it.
}
  \item{label.baseline}{
  If \code{TRUE}, the baseline levels for each attribute will be labeled as such. Defaults to \code{TRUE}.
}
  \item{text.size}{
  Size of text. Defaults to 11.
}
  \item{text.color}{
  Color of text in plot. Defaults to "black"
}
  \item{point.size}{
  Size of points in the plot. Defaults to 0.6
}
  \item{\dots}{
  Other graphical parameters passed to \code{plot}
}
}

\value{
A ggplot object containing a dotplot of estimated AMCEs
}
\references{
Hainmueller, J., Hopkins, D., and Yamamoto T. (2014) Causal Inference in Conjoint Analysis: Understanding Multi-Dimensional Choices via Stated Preference Experiments. Political Analysis 22(1):1-30
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{amce}} for the main estimation routine.
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

# Plot results
plot(results, xlab="Change in Pr(Immigrant Preferred for Admission to U.S.)", 
    ylim=c(-.3,.3), breaks=c(-.2, 0, .2), labels=c("-.2","0",".2"), text.size=13)

}

