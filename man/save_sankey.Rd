\name{save_sankey}
\alias{save_sankey}
\title{Export a Sankey (standard or tree-integrated)}
\usage{
save_sankey(x, file, width = NULL, height = NULL,
  pdf_size = c(11.7, 16.5), pdf_mode = c("fit", "tiles"), selfcontained = TRUE)
}
\arguments{
\item{x}{A Sankey htmlwidget (ordinary or tree-integrated).}
\item{file}{Output path ending in html, svg, pdf, png or rds.}
\item{width,height}{Optional dimensions: pixels for PNG, points for PDF. PDF values override pdf_size.}
\item{pdf_size}{Maximum PDF page width and height in inches.}
\item{selfcontained}{For standard Sankey HTML, create a single self-contained HTML file (requires Pandoc). FALSE uses a dependency folder.}
\item{pdf_mode}{Fit the complete vector figure on one bounded page, or split it into consecutive vertical tiles. Tiles require qpdf.}
}
\description{Exports fixed SVG geometry. Standalone HTML retains hover highlighting without embedding the raw R data or dependency bundles. PDF and PNG require rsvg.}
\value{Invisibly, the normalized output path.}
\details{RDS preserves original statistics and the tree. PDF tiles are continuation viewports and may split labels or ribbons. No statistical analysis is recomputed. Ordinary Sankey (tree_mode="none") supports HTML and RDS directly through save_sankey; SVG, PDF and PNG remain tree-only.}
