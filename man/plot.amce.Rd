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
                    labels = NULL, attribute.names = NULL, level.names = NULL,
                    label.baseline = TRUE, text.size = 11, text.color = "black",
                    point.size = 0.6, plot.theme = NULL, facet.name = NULL, 
                    facet.levels = NULL, show.all =	FALSE, 
                    show.interaction.base = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{
  An object of class "amce", a result of a call to \code{\link{amce}}
}
  \item{main }{
  Title of the plot.
}
  \item{xlab }{
  Label of the x-axis of the plot (AMCE or ACIE). Default is "Change in E[Y]"
}
  \item{ci }{
  Levels for confidence intervals to plot around point estimates. Must be between 0 and 1. Default is .95
}
  \item{colors }{
  Vector of color names to be used in the plot. The \code{plot} function will alternate between the colors in the vector for each attribute being plotted. If NULL, \code{plot} will use a default \code{ggplot2} color scheme.
}
  \item{xlim }{
 Numeric vector denoting the upper and lower bounds of the x-axis in the plot. If \code{NULL} the plot function will automatically set a range that includes all effect estimates.
}
  \item{breaks }{
  Numeric vector denoting where x-axis tick marks should be placed. If \code{NULL} plot will use \code{ggplot2} defaults.
}
  \item{labels }{
  Vector denoting how x-axis tick marks should be labeled. If \code{NULL} plot will use ggplot2 defaults.
}
  \item{attribute.names }{
  Character vector of attribute names to be plotted as labels. By default \code{plot.amce} will use the attribute names in the "amce" object passed to it.
}
  \item{level.names }{
  A list containing character vector elements with names in \code{attribute_names}. Each character vector in the list contains the level names to be plotted as labels beneath the corresponding attribute. By default \code{plot.amce} will use the level names in the "amce" object passed to it.
}
  \item{label.baseline }{
  If \code{TRUE}, the baseline levels for each attribute will be labeled as such. Defaults to \code{TRUE}.
}
  \item{text.size }{
  Size of text. Defaults to 11.
}
  \item{text.color }{
  Color of text in plot. Defaults to "black".
}
  \item{point.size }{
  Size of points in the plot. Defaults to 0.6.
}
  \item{plot.theme }{
  A ggplot2 'theme' object to be added to the plot. If NULL, defaults to black-and-white theme. Note that passing a theme object will override text and point color/size options. 
  }
  \item{facet.name }{
  To facet plots (i.e., make separate plots for each value of a variable) give "facet.name" the name of the variable (either profile attribute or respondent-varying) to facet by. Unless given specific levels in "facet.levels", the plotted levels will consist of all levels of a factor variable or the quantiles of a continuous variable. Currently there is only support for a single facet.
}
  \item{facet.levels }{
  To manually set facet levels, provide a list to "facet.levels". Names of list entries should correspond with variable names. The content of each entry should be a vector giving the desired levels, whether factors or continuous, with the rows named with the desired level names.
}
  \item{show.all }{
   A logical asking if variables that do not vary by facet level should also have their estimates reported in all of the facetted plots.
}
  \item{show.interaction.base }{
  A logical asking if the base term of an interaction sbould be included in facetted plots.
}
  \item{\dots }{
  Other graphical parameters passed to \code{ggplot}.
}
}

\value{
A \code{ggplot} object containing a dotplot of estimated AMCEs.
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

\dontrun{
# Run AMCE estimator with an interaction with a respondent-varying characteristic
interaction_results <- amce(Chosen_Immigrant ~ Gender + Education
		    + Job + ethnocentrism*Job,
		    data = subset(immigrationconjoint, !is.na(immigrationconjoint$ethnocentrism)),
        cluster=TRUE, respondent.id="CaseID",
		    design = immigrationdesign,
		    respondent.varying = c("ethnocentrism"))

# Plot results with different plots for quantiles of the respondent-varying characteristic
# Do not show output for variables that do not vary with facetted levels
plot(interaction_results,facet.name = c("ethnocentrism"))

# RUN AMCE estimator with an interaction between two factor variables
interaction_results <- amce(Chosen_Immigrant ~ Gender + Education + Job
		    + Education*Job, data = immigrationconjoint,
		    cluster = FALSE, design = immigrationdesign)

# Plot results with different plots for all levels of ACIE
plot(interaction_results,facet.name = c("Education"))

# Plot results with different plots for only two levels of one interacted variable
facet.levels1 <- list()
facet.levels1[["Education"]] <- data.frame(c("college degree","graduate degree"))
rownames(facet.levels1[["Education"]]) <- c("college degree","graduate degree")

plot(interaction_results,facet.name = c("Education"),facet.levels = facet.levels1)
}
}

