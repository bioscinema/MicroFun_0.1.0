\name{tax2fun_chord}
\alias{tax2fun_chord}
\alias{print.microfun_chord}
\alias{plot.microfun_chord}
\title{Chord diagram of taxon-resolved functional contributions}
\usage{
tax2fun_chord(pathway_sdaa_result, stratified_path, physeq, group,
  taxon_level = NULL, direction = NULL,
  function_level = c("auto", "KEGG", "KO", "EC", "MetaCyc"),
  contrast = NULL, taxa_list = NULL, functions = NULL,
  padj_cutoff = 0.05, min_abs_effect = 0,
  max_taxa_per_function = Inf, max_functions_per_taxon = Inf,
  contribution = c("balanced_mean", "mean"), function_table = NULL,
  output_pdf = NULL, output_png = NULL, draw = TRUE, output_html = NULL, ...)
\method{print}{microfun_chord}(x, ...)
\method{plot}{microfun_chord}(x, ...)
}
\arguments{
\item{pathway_sdaa_result}{A uniquely named list from \code{pathway_sdaa()}.}
\item{stratified_path}{Path to the PICRUSt2 contribution TSV or TSV.GZ.}
\item{physeq}{The phyloseq object used for analysis.}
\item{group}{Name of the sample metadata grouping variable.}
\item{taxon_level}{Taxonomic rank. Inferred from recorded settings, otherwise Genus.}
\item{direction}{Result orientation, taxa or function. Inferred when recorded.}
\item{function_level}{Functional namespace. KEGG uses PICRUSt2 pathway contributions directly or sums raw contributions of member KOs using the recorded or package mapping. Auto infers the namespace from recorded settings or feature IDs.}
\item{contrast}{Two labels, \code{c(reference, target)}, matching result group1/group2.}
\item{taxa_list}{Optional exact, case-insensitive taxon names.}
\item{functions}{Optional functional IDs.}
\item{padj_cutoff}{Maximum adjusted P value; defaults to 0.05.}
\item{min_abs_effect}{Minimum absolute effect.}
\item{max_taxa_per_function,max_functions_per_taxon}{Positive integer or Inf. Deterministic top-N by contribution, applied in this order.}
\item{contribution}{balanced_mean gives equal weight to group means; mean weights groups by sample size. Values retain the original input units.}
\item{function_table}{Optional unique function_id and function_group columns, with optional function_name. Missing annotations otherwise use Unclassified function.}
\item{output_pdf,output_png,output_html}{Optional export paths.}
\item{draw}{Draw on the current device if no export path is supplied.}
\item{x}{A microfun_chord object.}
\item{...}{Renderer settings: start_degree=180, side_gap=10, sector_gap=0.9, link_alpha=0.32, bar_alpha=0.95, effect_scale=1, effect_cap=Inf, outer_bar_track_height=0.135, label_track_height=0.086, taxon_label_cex=0.9, function_label_cex=0.85, show_legend=TRUE, show_title=FALSE, title_text=NULL, direction_cols (named Up, Down, Neutral), phylum_base_cols and function_group_palette (colour vectors, optionally named by category), taxonomy_tint_step=0 (additional whitening per ancestor track), show_category_legend=show_legend, legend_cex=0.72, show_outer_border=TRUE, outer_border_col="#333333", outer_border_lwd=0.8.}
}
\description{Connects differential results with actual taxon-resolved functional contributions. Taxonomy tracks and function groups use independent contribution and effect encodings.}
\details{
Signed colours default to Up = "#F5BB7C", Down = "#B4CBE0"; ribbons use link_alpha = 0.32. Functional labels use name[ID], resolved from an explicit function_name, result annotations, then bundled references. Unknown identifiers remain visible. Terminal prefixes are omitted.

Input columns are sample, taxon, function and taxon_function_abun.
Missing edge/sample records within observed samples count as zero; a sample entirely absent from the contribution file causes an error.
Widths use the contribution column recorded by SDAA or taxon_function_abun when available.
Non-KEGG inputs can also use norm_taxon_function_contrib when no abundance column is present.
KEGG requires taxon_function_abun. KO input uses the sum of raw contributions of the mapped pathway members.
KEGG pathways overlap, so widths summed across pathways do not form a partition of total microbial function.
Maaslin2 AST coefficients are labelled as such and are not log fold changes.
Bars use absolute effect times effect_scale, capped at effect_cap. Their scale is stated in the legend.
Only ranks above the selected taxon_level are drawn as taxonomy rings. Duplicated taxon names with conflicting ancestral assignments cause an error.
A function maps to one displayed group. Hierarchical annotations use the text before the first semicolon. Available non-missing annotations are reconciled across taxa before Top-N filtering; known KEGG category case and whitespace are normalised. Genuine conflicts produce a warning and the display group Ambiguous annotation. The returned function_table records observed_groups and annotation_status for automatic mappings. All-missing annotations use Unclassified function, with KO reference fallback where available. Supply function_table to override this display mapping.
The default taxon palette is light muted nine-colour and the function palette is dark muted nine-colour.
Ancestor tracks inherit the root category colour without extra whitening by default.
Named palettes fix category assignments across plots; excess categories cycle with a warning.
Category legends use separate pages by default, or a right panel with legend_position="right": Taxa first, then Function. Taxa colours refer to the highest available ancestor rank, or terminal taxa when ancestry is absent. Swatches match the actual displayed category colours. Direction and encoding legends follow below. Legend text fits the selected layout. Set show_legend=FALSE to hide all legends; show_category_legend can override category visibility independently. A continuous border traces the outer edge of the effect-bar track; its colour and width are configurable. These settings also apply when redrawing saved chord objects.
No statistics are rerun by the plotting function. A backend without a usable effect estimate cannot be plotted.
This is a base/circlize plot, not a ggplot; use save_chord(), output_pdf/output_png/output_html or a graphics device.
}
\value{Invisibly returns a microfun_chord list: links, all_links (significant positive-width edges before top-N), tax_table, function_table, settings and diagnostics. Print or plot redraws it.}
\examples{
\dontrun{
p <- tax2fun_chord(result, "pred_metagenome_contrib.tsv.gz", ps1,
  group = "group", taxon_level = "Class", contrast = c("STS", "LTS"),
  max_functions_per_taxon = 8)
p
plot(p, show_legend = FALSE)
p$links
}
}

\section{Label options}{Graphical options include narrow_labels ("auto", "letters", "full"), label_codes (a uniquely named vector), label_fit (0.85), label_mapping_file and legend_min_cex (0.65). Defaults include legend_position="separate", taxonomy_label_cex=0.7, narrow_sector_degrees=0.5 and prefer_full_taxa=TRUE. Auto mode uses a letter code when the full name does not fit. With prefer_full_taxa enabled, terminal taxa narrower than narrow_sector_degrees also use a code. A label that exceeds the available arc is hidden. Letter codes are unique across taxon ranks and functions. The mapping CSV retains hidden labels. HTML lists every sector for click-to-highlight. PDF legends include visible codes. Reuse returned label_codes for consistent mappings across different plots.}
