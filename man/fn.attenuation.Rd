% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_attenuation.R
\name{fn.attenuation}
\alias{fn.attenuation}
\alias{fn.attenuation2}
\title{Computes attenuation factor}
\usage{
fn.attenuation(
  ystst,
  X,
  proj,
  l1obs,
  s1obs,
  fn.psi = fn.psi.huber,
  fn.chi = fn.chi.prop2
)

fn.attenuation2(
  ystst,
  X,
  proj,
  Qt,
  l1obs,
  s1obs,
  fn.psi = fn.psi.huber,
  fn.chi = fn.chi.prop2
)
}
\arguments{
\item{ystst}{a data vector such that T(ystst)=T(y_obs)}

\item{X}{design matrix}

\item{proj}{the projection matrix onto the deviation space (aka the least squares residual space)}

\item{l1obs}{s1obs: observed statsitics}

\item{Qt}{Q transpose where Q is the orthonormalized X}

\item{fn.psi:}{custom psi function returning psi(x) or psi'(x) defining M-estimator for location; default is fn.psi.huber}

\item{fn.chi:}{custom chi function returning chi(x) or chi'(x) defining M-estimator for scale; default is fn.chi.prop2}
}
\value{
Attenuations of the transformation to ystst. These attenuations make up part of the Jacobian.
}
\description{
Compute the attenuation factor of the Jacobian in the restricted likelihood MCMC step
}
\details{
Two equivalent functions computing the attenuation piece. The second one uses the transpose of the orthonormalized \code{X} as an input. This makes it somewhat quicker.
}
