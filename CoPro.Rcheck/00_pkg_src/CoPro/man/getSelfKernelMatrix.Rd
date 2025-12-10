\name{getSelfKernelMatrix}
\alias{getSelfKernelMatrix}
\title{Get Self-Kernel Matrix}
\usage{
getSelfKernelMatrix(object, sigma, cellType, slide = NULL, verbose = TRUE)
}
\arguments{
\item{object}{A CoPro object}

\item{sigma}{Sigma value}

\item{cellType}{Cell type name}

\item{slide}{Slide ID (for CoProMulti objects)}

\item{verbose}{Whether to print error messages}
}
\value{
Self-kernel matrix for the specified parameters
}
\description{
Convenience function to retrieve self-kernel matrices computed by computeSelfKernel.
}
