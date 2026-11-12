\name{microfun_palette}
\alias{microfun_palette}
\title{MicroFun categorical palettes}
\usage{microfun_palette(side = c("taxa", "function"))}
\arguments{\item{side}{Taxa uses the light muted nine-colour palette; function uses the dark muted nine-colour palette.}}
\value{A character vector of nine hexadecimal colours.}
\description{Returns the categorical palettes used by chord and Sankey plots. Effect direction colours are controlled separately by each plotting function.}
\examples{microfun_palette("taxa")}
