\name{save_chord}
\alias{save_chord}
\title{Export a chord diagram}
\usage{save_chord(x, file = "microfun_chord.html", width = 18, height = 16,
  formats = NULL, ...)}
\arguments{
\item{x}{A microfun_chord object.}
\item{file}{Output path ending in pdf, png, svg, html or rds; defaults to microfun_chord.html.}
\item{formats}{Optional additional formats to save alongside the primary output: html, pdf, png, svg, rds.}
\item{width,height}{Figure dimensions in inches.}
\item{...}{Graphical arguments passed to the chord renderer; see \code{tax2fun_chord}.}
}
\description{Save a chord diagram and its label mapping. HTML includes zoom controls and a searchable, clickable sector legend. PDF legends continue on separate pages.}
\details{HTML uses a separate legend panel. Selecting an entry highlights its sector, including sectors with hidden labels. Labels that do not fit are omitted from the chart, and their aliases are omitted from static legends. PNG and SVG write a companion legend PDF when necessary. Rendered exports write a label mapping CSV. RDS stores the chord object without rendering.}
\value{Invisibly, the normalized output path or named vector of paths for multi-format export.}
\seealso{\code{\link{tax2fun_chord}}}
