\name{redraw_sankey}
\alias{redraw_sankey}
\title{Redraw saved tree-integrated Sankey results}
\usage{redraw_sankey(x, tree_options = list(), width = NULL, height = NULL)}
\arguments{
\item{x}{An existing tree-enabled MicroFun widget or its RDS path.}
\item{tree_options}{Named appearance options. Defaults enable aggregate_links, collapse_clades and bundle_mapping, and disable title, subtitle and show_footer.}
\item{width,height}{Minimum drawing dimensions in pixels.}
}
\description{Reuses saved associations, the retained original tree and taxon ordering. No differential model is fitted.}
\value{A tree-enabled htmlwidget with original edges, aggregated ribbons, clade membership and SVG geometry.}
