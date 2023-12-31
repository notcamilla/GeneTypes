% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classes.R
\name{name}
\alias{name}
\alias{name,Gene-method}
\title{Accessor function for the name of a Gene object.}
\usage{
name(x)

\S4method{name}{Gene}(x)
}
\arguments{
\item{x}{A Gene object.}
}
\value{
The name of the Gene.

The updated Gene object.
}
\description{
Accessor function for the name of a Gene object.

Accessor function for the name of a Gene object.
}
\examples{
gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3
metabolism regulating signaling molecule C',description = 'encodes a secreted
protein with a GG domain', structure = GRanges("chr7",IRanges(121348878,
121396364), strand = '-'), proteinID = 'Q92520')
name(gene1)
}
