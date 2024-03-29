% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_density.R
\name{plot_density}
\alias{plot_density}
\title{Generate density plot}
\usage{
plot_density(model, data, save_to, filename, save_pdf = FALSE)
}
\arguments{
\item{model}{a \code{runjags} object containing model results}

\item{data}{a \code{list} containing the data input for the runjags model,
generated with \code{jags_data()}}

\item{save_to}{a \code{string} specifying the location in which the output
(an \code{rds} file containing the density plot) should be saved}

\item{filename}{a \code{string} specifying the name of the output (an
\code{rds} file containing the density plot)}

\item{save_pdf}{a \code{boolean} specifying whether or not a pdf should be
saved}
}
\value{
Returns a density plot and saves an rds file
}
\description{
Generates a violin plot where violins represent the posterior probability of
resistance and points represent the calculated probability of resistance to
each antibiotic class. If a file already exists, it is read into memory.
}
