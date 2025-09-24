\name{getSelfDistMat}
\alias{getSelfDistMat}
\title{Get Self-Distance Matrix}
\usage{
getSelfDistMat(object, cellType, slide = NULL, verbose = TRUE)
}
\arguments{
\item{object}{A CoPro object}

\item{cellType}{Cell type name}

\item{slide}{Slide ID (for CoProMulti objects)}

\item{verbose}{Whether to print error messages}
}
\value{
Self-distance matrix for the specified cell type
}
\description{
Convenience function to retrieve self-distance matrices computed by computeSelfDistance.
}
