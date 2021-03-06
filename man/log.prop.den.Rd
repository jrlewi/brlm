% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_prop_den.R
\name{log.prop.den}
\alias{log.prop.den}
\alias{log.prop.den2}
\title{Caluclate the proposal dens on the log scale
Computes proposal density by computing the necessary Jacobian of the transformation from y* to y**}
\usage{
\method{log}{prop.den}(y.prop, X, proj, l1obs, s1obs, fn.psi, fn.chi, n, p)

\method{log}{prop.den2}(y.prop, X, proj, Qt, l1obs, s1obs, fn.psi, fn.chi, n, p)
}
\arguments{
\item{y.prop}{proposed data vector satisfying T(y.prop)=T(y_obs)}

\item{X}{design matrix}

\item{proj}{the projection matrix onto the deviation space (aka the least squares residual space)}

\item{l1obs, }{s1obs observed statsitics}

\item{fn.psi}{custom psi function returning psi(x) or psi'(x) defining M-estimator for location ; default is fn.psi.huber}

\item{fn.chi}{custom chi function returning chi(x) or chi'(x) defining M-estimator for scale; default is fn.chi.prop2}

\item{n}{\code{length(y_obs)=length(y.prop)}}

\item{p}{number of regression coefficients}

\item{Qt}{Q transpose where Q is the orthonormalized X}
}
\value{
The proposal density for y.prop on the log scale (assuming original vector sampled uniformly on unit sphere)
}
\description{
Caluclate the proposal dens on the log scale
Computes proposal density by computing the necessary Jacobian of the transformation from y* to y**
}
\details{
Designed for use within \code{fn.one.rep.y}. Two equivalent versions, the second uses \code{fn.attenuation2} and is faster.
}
