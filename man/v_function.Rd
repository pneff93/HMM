% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General_functions.R
\name{v_function}
\alias{v_function}
\title{Normalized v-function EM}
\usage{
v_function(m, N, beta, p, weight, Gamma, set, u)
}
\arguments{
\item{m}{number of likelihoods}

\item{N}{length of the supplied dataset}

\item{beta}{beta matrix}

\item{p}{vector of likelihood probabilities of dataset}

\item{weight}{weights provided by the alpha calculation}

\item{Gamma}{Gamma matrix}

\item{set}{index vector to align the p vector}

\item{u}{output matrix of the u_function}
}
\value{
returns the coresponding values, but due to further calculation steps
the m x m matrix for each timepoint of v-values is returned as one row in 
the output matrix.
}
\description{
This function calculates the v-function which is part of the 
estimation step. Each element in v aligns with the conditional expectation 
to reach state k from the previous state j given your datapoint x(t).
Due to the normalizing of the factors we also need to include the
corresponding weights into our function.
}
