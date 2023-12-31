% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classes.R
\name{anticodon}
\alias{anticodon}
\alias{anticodon,tRNAGene-method}
\title{Accessor function for the anticodon of a transfer RNAGene object.}
\usage{
anticodon(x)

\S4method{anticodon}{tRNAGene}(x)
}
\arguments{
\item{x}{A Gene object.}
}
\value{
The anticodon of the transfer RNAGene

The updated Gene object.
}
\description{
Accessor function for the anticodon of a transfer RNAGene object.

Accessor function for the anticodon of a transfer RNA Gene object.
}
\examples{
gene5 <- tRNAGene(id=100189346, symbol= 'TRA-AGC1-1',name='tRNA-Ala
(anticodon AGC) 1-1', description = 'tRNA-Ala (anticodon AGC) 1-1',
structure = GRanges("chr6", IRanges(28795964, 28796035), strand = "-"),
tRNAID='URS000063E4FD_9606', anticodon='AGC')
anticodon(gene5)
}
