
#' @importFrom magrittr %>%
NULL

#' Differential Abundance Analysis for Predicted Functional Pathways
#'
#' @description Performs differential abundance analysis on predicted functional pathway data using various statistical methods.
#' This function supports multiple methods for analyzing differences in pathway abundance between groups,
#' including popular approaches like ALDEx2, DESeq2, edgeR, and others.
#'
#' @param group Character string specifying the column name in metadata that contains group information
#'        for differential abundance analysis.
#'
#' @param daa_method Character string specifying the method for differential abundance analysis.
#'        Available choices are:
#'        \itemize{
#'          \item \code{"ALDEx2"}: ANOVA-Like Differential Expression tool
#'          \item \code{"DESeq2"}: Differential expression analysis based on negative binomial distribution
#'          \item \code{"edgeR"}: Exact test for differences between groups using negative binomial model
#'          \item \code{"limma voom"}: Limma-voom framework for RNA-seq analysis
#'          \item \code{"metagenomeSeq"}: Zero-inflated Gaussian mixture model
#'          \item \code{"LinDA"}: Linear models for differential abundance analysis
#'          \item \code{"Maaslin2"}: Multivariate Association with Linear Models
#'          \item \code{"Lefser"}: Linear discriminant analysis effect size
#'        }
#'        Default is "ALDEx2".
#'
#' @param select Vector of sample names to include in the analysis.
#'        If NULL (default), all samples are included.
#'
#' @param p.adjust Character string specifying the method for p-value adjustment.
#'        Choices are:
#'        \itemize{
#'          \item \code{"BH"}: Benjamini-Hochberg procedure (default)
#'          \item \code{"holm"}: Holm's step-down method
#'          \item \code{"bonferroni"}: Bonferroni correction
#'          \item \code{"hochberg"}: Hochberg's step-up method
#'          \item \code{"fdr"}: False Discovery Rate
#'          \item \code{"none"}: No adjustment
#'        }
#'
#' @references \itemize{
#'   \item ALDEx2: Fernandes et al. (2014) Unifying the analysis of high-throughput sequencing datasets:
#'         characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by
#'         compositional data analysis. Microbiome.
#'   \item DESeq2: Love et al. (2014) Moderated estimation of fold change and dispersion for
#'         RNA-seq data with DESeq2. Genome Biology.
#'   \item edgeR: Robinson et al. (2010) edgeR: a Bioconductor package for differential expression
#'         analysis of digital gene expression data. Bioinformatics.
#'   \item limma-voom: Law et al. (2014) voom: precision weights unlock linear model analysis tools
#'         for RNA-seq read counts. Genome Biology.
#'   \item metagenomeSeq: Paulson et al. (2013) Differential abundance analysis for microbial
#'         marker-gene surveys. Nature Methods.
#'   \item Maaslin2: Mallick et al. (2021) Multivariable Association Discovery in Population-scale
#'         Meta-omics Studies.
#' }
#'
#' @export
#'
#' @param abundance Numeric matrix or data frame with feature IDs as row names and unique sample IDs as column names. At least four selected samples and two observed groups are required. Choose input units appropriate for the selected backend; several backends round inputs to counts.
#'
#' @param metadata Data frame or tibble with a unique sample-identifier column matching the full set of abundance column names before sample selection. A matching `sample` column is preferred; otherwise other columns are searched. Row names alone are not used. Metadata is reordered to match the selected abundance columns.
#'
#' @param confounder Optional character vector of metadata covariate names. Only the MaAsLin2 backend uses this argument in the current implementation.
#'
#' @param reference Optional reference group. Used by MaAsLin2 and limma-voom; other backends use their own factor-level handling. For LinDA, supply an explicit reference matching the first fitted factor level; this wrapper does not relevel its model.
#'
#' @param ... Reserved for compatibility; currently not forwarded to backend functions.
#'
#' @param maaslin2_transform Transformation used only with `daa_method = "Maaslin2"`: `"AST"` (default), `"LOG"`, `"LOGIT"` or `"NONE"`. MaAsLin2 uses TSS normalization before this transformation. Supplying this argument with another method raises an error.
#'
#' @return A data frame whose columns depend on the backend. P-value backends return `feature`, `method`, group labels, `p_values`, `p_adjust` and `adj_method`. MaAsLin2 additionally returns `logFC` (the fitted coefficient) and `effect_type`; LinDA uses `log2FoldChange`. Lefser returns `effect_scores` without p-values or adjusted p-values. Some backends can return `NULL` when no result is produced.
#'
#' @details Backend names are case-sensitive: `ALDEx2`, `DESeq2`, `edgeR`, `limma voom`, `metagenomeSeq`, `LinDA`, `Maaslin2` and `Lefser`. Adjustment via `p.adjust()` is applied to all returned p-values together, including multiple contrasts or test types. MaAsLin2 coefficients are transformation-specific model effects, not automatically log fold changes. Other backends do not all return an effect estimate suitable for stratified plots. Multi-group contrast behavior depends on the backend; DESeq2 currently reports the second versus first factor level.
#'
#' @seealso [pathway_sdaa()], [pathway_annotation()]
#'
#' @md
pathway_daa <- function(abundance, metadata, group, confounder = NULL, daa_method = "ALDEx2",
                        select = NULL, p.adjust = "BH", reference = NULL, ..., maaslin2_transform = "AST") {
  if (daa_method == "Maaslin2") {
    maaslin2_transform <- .microfun_maaslin2_transform(maaslin2_transform)
  } else if (!missing(maaslin2_transform)) {
    stop("maaslin2_transform applies only when daa_method is 'Maaslin2'.")
  }
  # Check if the requested DAA method package is available
  method_packages <- list(
    "ALDEx2" = "ALDEx2",
    "DESeq2" = "DESeq2",
    "edgeR" = "edgeR",
    "limma voom" = "limma",
    "metagenomeSeq" = "metagenomeSeq",
    "Maaslin2" = "Maaslin2"
  )

  if (daa_method %in% names(method_packages)) {
    pkg_name <- method_packages[[daa_method]]
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for method '%s'. Please install it using BiocManager::install('%s')",
                   pkg_name, daa_method, pkg_name))
    }
  }
  # Input validation
  if (!is.data.frame(abundance) && !is.matrix(abundance)) {
    stop("abundance must be a data frame or matrix")
  }

  if (ncol(abundance) < 4) {
    stop("At least 4 samples are required for differential abundance analysis")
  }

  # Convert metadata to tibble
  if (!tibble::is_tibble(metadata)) {
    metadata <- tibble::as_tibble(metadata)
  }

  # Extract sample names from abundance data
  abundance_samples <- colnames(abundance)
  if (is.null(abundance_samples) || anyNA(abundance_samples) ||
      any(!nzchar(abundance_samples)) || anyDuplicated(abundance_samples)) {
    stop("Abundance data must have unique, non-empty sample IDs")
  }

  # Identify the column in metadata that matches sample names
  sample_col <- NULL
  for (col in unique(c(intersect("sample", names(metadata)), names(metadata)))) {
    ids <- as.character(metadata[[col]])
    if (!anyNA(ids) && !anyDuplicated(ids) && setequal(ids, abundance_samples)) {
      sample_col <- col
      break
    }
  }

  # Check if a matching column was found
  if (is.null(sample_col)) {
    stop("No column in metadata matches the sample names in abundance data")
  }

  message(sprintf("Using column '%s' as sample identifier", sample_col))

  # Verify if group column exists
  if (!group %in% colnames(metadata)) {
    stop(sprintf("group column '%s' not found in metadata", group))
  }

  # Handle sample selection
  if (!is.null(select)) {
    select <- as.character(select)
    if (!all(select %in% abundance_samples)) {
      stop("Some selected samples are not present in the abundance data")
    }
    if (anyDuplicated(select)) {
      stop("Selected sample IDs must be unique")
    }
    abundance <- abundance[, select, drop = FALSE]
  }

  # Align metadata to the final abundance columns, including the order of select.
  sample_ids <- colnames(abundance)
  if (length(sample_ids) < 4) {
    stop("At least 4 samples are required for differential abundance analysis")
  }
  metadata <- as.data.frame(metadata[match(sample_ids, as.character(metadata[[sample_col]])), ])
  rownames(metadata) <- sample_ids

  # Prepare data
  abundance_mat <- as.matrix(abundance)
  Group <- factor(metadata[[group]])
  Level <- levels(Group)
  length_Level <- length(Level)
  if (anyNA(Group)) {
    stop("Selected samples must have non-missing group assignments")
  }
  if (length_Level < 2) {
    stop("At least two groups are required for differential abundance analysis")
  }
  if (!is.null(reference) && !reference %in% Level) {
    stop("reference is not a group level in the selected samples")
  }

  # Perform differential analysis
  result <- switch(
    daa_method,
    "ALDEx2" = perform_aldex2_analysis(abundance_mat, Group, Level, length_Level),
    "DESeq2" = perform_deseq2_analysis(abundance_mat, metadata, group, Level),
    "LinDA" = perform_linda_analysis(abundance, metadata, group, reference, Level, length_Level),
    "limma voom" = perform_limma_voom_analysis(abundance_mat, Group, reference, Level, length_Level),
    "edgeR" = perform_edger_analysis(abundance_mat, Group, Level, length_Level),
    "metagenomeSeq" = perform_metagenomeseq_analysis(abundance_mat, metadata, group, Level),
    "Maaslin2" = perform_maaslin2_analysis(abundance_mat, metadata, group, confounder, reference, Level,
      length_Level, maaslin2_transform),
    "Lefser" = perform_lefser_analysis(abundance_mat, metadata, group, Level)
  )

  # Add multiple testing correction
  if (!is.null(result) && "p_values" %in% colnames(result)) {
    result$p_adjust <- p.adjust(result$p_values, method = p.adjust)
    result$adj_method <- p.adjust
  }

  return(result)
}

# Helper function: Perform ALDEx2 analysis
perform_aldex2_analysis <- function(abundance_mat, Group, Level, length_Level) {
  message("Running ALDEx2 analysis...")

  # Check if ALDEx2 is available
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    stop("Package 'ALDEx2' is required for ALDEx2 analysis. Please install it using BiocManager::install('ALDEx2')")
  }

  # Round the abundance data
  abundance_mat <- round(abundance_mat)

  # Save the original Group factor and convert to numeric for ALDEx2
  original_Group <- Group
  Group <- as.numeric(Group)

  # Perform different analyses based on the number of groups
  if (length_Level == 2) {
    message("Running ALDEx2 with two groups. Performing t-test...")

    # Create ALDEx2 object with numeric Group
    ALDEx2_object <- ALDEx2::aldex.clr(
      abundance_mat,
      Group,
      mc.samples = 256,
      denom = "all",
      verbose = FALSE
    )

    # Get t-test results
    results <- ALDEx2::aldex.ttest(
      ALDEx2_object,
      paired.test = FALSE,
      verbose = FALSE
    )

    # Build result dataframe using original Level names
    return(data.frame(
      feature = rep(rownames(results), 2),
      method = c(
        rep("ALDEx2_Welch's t test", nrow(results)),
        rep("ALDEx2_Wilcoxon rank test", nrow(results))
      ),
      group1 = rep(Level[1], 2 * nrow(results)),
      group2 = rep(Level[2], 2 * nrow(results)),
      p_values = c(results$we.ep, results$wi.ep),
      stringsAsFactors = FALSE
    ))

  } else {
    message("Running ALDEx2 with multiple groups. This might take some time...")

    # Create ALDEx2 object for multiple groups with numeric Group
    ALDEx2_object <- ALDEx2::aldex.clr(
      abundance_mat,
      Group,
      mc.samples = 256,
      denom = "all",
      verbose = FALSE
    )

    # Get Kruskal-Wallis and GLM test results
    results <- ALDEx2::aldex.kw(ALDEx2_object)

    # Build initial result dataframe
    result_df <- data.frame(
      feature = rep(rownames(results), 2),
      method = c(
        rep("ALDEx2_Kruskal-Wallace test", nrow(results)),
        rep("ALDEx2_glm test", nrow(results))
      ),
      p_values = c(results$kw.ep, results$glm.ep),
      stringsAsFactors = FALSE
    )

    # Add all group information using original Level names
    group_cols <- data.frame(matrix(
      NA,
      nrow = nrow(result_df),
      ncol = length_Level,
      dimnames = list(NULL, paste0("group", 1:length_Level))
    ))

    for (i in 1:length_Level) {
      group_cols[, i] <- Level[i]
    }

    # Merge results
    result_df <- cbind(
      result_df[, c("feature", "method")],
      group_cols,
      result_df[, "p_values", drop = FALSE]
    )

    message("ALDEx2 analysis with multiple groups complete.")
    return(result_df)
  }
}

# Helper function: Perform DESeq2 analysis
perform_deseq2_analysis <- function(abundance_mat, metadata, group, Level) {
  message("Running DESeq2 analysis...")

  # Convert to integer matrix
  message("converting counts to integer mode")
  counts <- round(as.matrix(abundance_mat))

  # Ensure group is a factor.
  metadata[[group]] <- factor(metadata[[group]], levels = Level)

  # Create DESeqDataSet object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = metadata
  )

  # Use try-catch to handle possible errors
  result <- tryCatch({
    suppressWarnings({
      # Create DESeqDataSet
      dds <- DESeq2::DESeqDataSet(se, design = as.formula(paste0("~", group)))

      # Choose the fitting method based on sample size.
      fitType <- if(ncol(abundance_mat) < 6) "mean" else "local"

      # Run DESeq2 pipeline
      dds <- DESeq2::estimateSizeFactors(dds)
      dds <- DESeq2::estimateDispersions(dds, fitType = fitType)
      dds <- DESeq2::nbinomWaldTest(dds)

      # Extract results
      res <- DESeq2::results(dds, contrast = c(group, Level[2], Level[1]))

      data.frame(
        feature = rownames(abundance_mat),
        method = "DESeq2",
        group1 = Level[1],
        group2 = Level[2],
        p_values = res$pvalue,
        stringsAsFactors = FALSE
      )
    })
  }, error = function(e) {
    stop("DESeq2 analysis failed: ", e$message)
  })

  if (is.null(result)) {
    stop("DESeq2 analysis failed to produce results")
  }

  return(result)
}

# Helper function: Perform limma voom analysis
perform_limma_voom_analysis <- function(abundance_mat, Group, reference, Level, length_Level) {
  message("Running limma voom analysis...")

  # Check if limma is available
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for limma voom analysis. Please install it using BiocManager::install('limma')")
  }

  # Ensure Group is a factor type
  Group <- factor(Group)
  if (!is.null(reference)) {
    Group <- relevel(Group, ref = reference)
  }

  # Create design matrix
  design <- model.matrix(~Group, contrasts.arg = list(Group = "contr.treatment"))

  # Perform voom transformation and analysis
  dge <- edgeR::DGEList(counts = abundance_mat, group = Group)
  dge <- edgeR::calcNormFactors(dge)
  v <- limma::voom(dge, design)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)

  # Extract results
  group_levels <- levels(Group)
  results <- do.call(rbind, lapply(seq.int(2L, ncol(design)), function(i) {
    data.frame(
      feature = rownames(abundance_mat),
      method = "limma voom",
      group1 = group_levels[1L],
      group2 = group_levels[i],
      p_values = fit$p.value[, i]
    )
  }))

  return(results)
}

# Helper function: Perform edgeR analysis
perform_edger_analysis <- function(abundance_mat, Group, Level, length_Level) {
  message("Running edgeR analysis...")

  # Create DGEList object
  dge <- edgeR::DGEList(counts = round(abundance_mat), group = Group)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateCommonDisp(dge, verbose = TRUE)

  if (length_Level == 2) {
    # Two-group comparison
    et <- edgeR::exactTest(dge, pair = c(1, 2))
    results <- data.frame(
      feature = rownames(abundance_mat),
      method = "edgeR",
      group1 = Level[1],
      group2 = Level[2],
      p_values = et$table$PValue
    )
  } else {
    # Multi-group comparison
    results_list <- list()
    combinations <- utils::combn(seq_along(Level), 2)

    for (i in 1:ncol(combinations)) {
      et <- edgeR::exactTest(dge, pair = combinations[,i])
      results_list[[i]] <- data.frame(
        feature = rownames(abundance_mat),
        method = "edgeR",
        group1 = Level[combinations[1,i]],
        group2 = Level[combinations[2,i]],
        p_values = et$table$PValue
      )
    }
    results <- do.call(rbind, results_list)
  }

  return(results)
}

# Helper function: Perform metagenomeSeq analysis
perform_metagenomeseq_analysis <- function(abundance_mat, metadata, group, Level) {
  message("Running metagenomeSeq analysis...")

  # Convert metadata to data.frame and ensure sample names are correct
  metadata <- as.data.frame(metadata)

  # Ensure abundance_mat and metadata have consistent sample order
  abundance_mat <- abundance_mat[, rownames(metadata)]

  # Create phenoData
  phenoData <- new("AnnotatedDataFrame",
                   data = metadata,
                   varMetadata = data.frame(
                     labelDescription = c("Sample ID", "Group"),
                     row.names = colnames(metadata)
                   ))

  # Ensure data is an integer matrix
  counts <- round(as.matrix(abundance_mat))

  # Create MRexperiment object
  obj <- try({
    metagenomeSeq::newMRexperiment(
      counts = counts,
      phenoData = phenoData,
      featureData = NULL,
      libSize = NULL,
      normFactors = NULL
    )
  }, silent = TRUE)

  if (inherits(obj, "try-error")) {
    stop("Failed to create MRexperiment object: ", attr(obj, "condition")$message)
  }

  # Normalize
  obj <- metagenomeSeq::cumNorm(obj)

  # Create model matrix
  mod <- stats::model.matrix(as.formula(paste0("~", group)), data = metadata)

  # Fit model
  fit <- metagenomeSeq::fitFeatureModel(obj, mod)

  # Extract results
  results <- data.frame(
    feature = rownames(abundance_mat),
    method = "metagenomeSeq",
    group1 = Level[1],
    group2 = Level[2],
    p_values = fit@pvalues,
    stringsAsFactors = FALSE
  )

  return(results)
}

# Helper function: Perform Maaslin2 analysis
.microfun_maaslin2_transform <- function(transform) {
  choices <- c("AST", "LOG", "LOGIT", "NONE")
  if (!is.character(transform) || length(transform) != 1L || is.na(transform) || !transform %in% choices)
    stop("maaslin2_transform must be one of: ", paste(choices, collapse = ", "), ".")
  transform
}

perform_maaslin2_analysis <- function(abundance_mat, metadata, group, confounder, reference, Level, length_Level,
                                     transform = "AST") {
  message("Running Maaslin2 analysis...")

  # Check if Maaslin2 is available
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    stop("Maaslin2 package is required but not installed")
  }

  # Transpose abundance matrix
  abundance_mat_t <- t(abundance_mat)
  # colnames(abundance_mat_t) <- rownames(abundance_mat)
  # print(colnames(abundance_mat_t))

  # Prepare metadata
  metadata <- as.data.frame(metadata)

  # Create temporary output directory
  output_dir <- tempfile("MicroFun_Maaslin2_")

  # Run Maaslin2 analysis
  fit_data <- tryCatch({
    Maaslin2::Maaslin2(
      input_data = abundance_mat_t,
      input_metadata = metadata,
      output = output_dir,
      transform = transform,
      fixed_effects = c(group, confounder),
      reference = paste0(group, ",", if (is.null(reference)) Level[1] else reference),
      normalization = "TSS",
      standardize = TRUE,
      min_prevalence = 0.1,
      min_variance = 0,
      cores = 1,
      plot_heatmap = FALSE,
      plot_scatter = FALSE
    )
  }, error = function(e) {
    message("  -> Maaslin2 failed: ", e$message)
    return(NULL)
  })

  if (is.null(fit_data)) return(NULL)

  # Read all results instead of significant results
  results_file <- file.path(output_dir, "all_results.tsv")

  if (file.exists(results_file)) {
    maaslin2_results <- utils::read.table(results_file,
                                          header = TRUE,
                                          sep = "\t",
                                          stringsAsFactors = FALSE)

    # Keep the requested primary group term and all its observed contrasts.
    mr <- maaslin2_results[maaslin2_results$metadata == group, , drop = FALSE]
    ref_level <- if (is.null(reference)) Level[1] else reference
    mr <- mr[as.character(mr$value) %in% setdiff(Level, ref_level), , drop = FALSE]
    results <- data.frame(
      feature = as.character(mr$feature), method = "Maaslin2",
      group1 = rep(ref_level, nrow(mr)), group2 = as.character(mr$value),
      p_values = mr$pval, logFC = mr$coef,
      effect_type = rep(paste("Maaslin2", transform, "coefficient"), nrow(mr)),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Maaslin2 analysis failed to produce results file")
  }

  return(results)
}

# Helper function: Perform Lefser analysis
perform_lefser_analysis <- function(abundance_mat, metadata, group, Level) {
  message("Running Lefser analysis...")

  # Check if lefser package is available
  if (!requireNamespace("lefser", quietly = TRUE)) {
    stop("The 'lefser' package is required for Lefser analysis. Please install it using BiocManager::install('lefser')")
  }

  metadata[[group]] <- factor(metadata[[group]])

  # Create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = abundance_mat),
    colData = metadata
  )

  # 3) Convert counts → relative abundances
  se <- lefser::relativeAb(se)
  # colData(se)[[group]] <- factor(colData(se)[[group]])

  # Perform Lefser analysis
  lefser_results <- lefser::lefser(se, classCol = group)  # Using classCol instead of deprecated groupCol

  # Check if results are empty
  if (length(lefser_results$Names) == 0 || is.null(lefser_results$Names)) {
    # If no significant features found, return an empty data frame with correct columns
    message("No significant features found by Lefser analysis.")
    return(data.frame(
      feature = character(0),
      method = character(0),
      group1 = character(0),
      group2 = character(0),
      effect_scores = numeric(0)
    ))
  } else {
    # If significant features found, return normal results
    results <- data.frame(
      feature = lefser_results$Names,
      method = "Lefser",
      group1 = Level[1],
      group2 = Level[2],
      effect_scores = lefser_results$scores
    )
    return(results)
  }
}

# Helper function: Perform LinDA analysis
perform_linda_analysis <- function(abundance, metadata, group, reference, Level, length_Level) {
  message("Running LinDA analysis...")

  # Ensure necessary packages are loaded
  if (!requireNamespace("MicrobiomeStat", quietly = TRUE)) {
    stop("MicrobiomeStat package is required but not installed")
  }

  # Prepare data
  feature.dat <- as.matrix(abundance)
  meta.dat <- metadata

  # Build formula
  formula <- paste0("~ ", group)

  # Run LinDA analysis
  linda_obj <- MicrobiomeStat::linda(
    feature.dat = feature.dat,
    meta.dat = meta.dat,
    formula = formula,
    feature.dat.type = "count",
    prev.filter = 0,
    mean.abund.filter = 0,
    adaptive = TRUE,
    zero.handling = "pseudo-count",
    pseudo.cnt = 0.5,
    p.adj.method = "BH",
    alpha = 0.05,
    n.cores = 1,
    verbose = FALSE
  )

  # Extract results
  if (length(linda_obj$output) > 0) {
    # Create an empty list to store results for each comparison
    results_list <- list()

    # Get all non-reference groups
    non_ref_groups <- Level[Level != reference]

    # Process each output element (each corresponds to a group comparison)
    for (i in seq_along(linda_obj$output)) {
      if (i <= length(non_ref_groups)) {
        comparison_df <- linda_obj$output[[i]]

        # Create a results data frame for this comparison
        results_list[[i]] <- data.frame(
          feature = rownames(comparison_df),
          method = "LinDA",
          group1 = reference,
          group2 = non_ref_groups[i],
          p_values = comparison_df$pvalue,
          log2FoldChange = comparison_df$log2FoldChange,
          stringsAsFactors = FALSE
        )
      }
    }

    # Combine all results
    results <- do.call(rbind, results_list)
  } else {
    stop("LinDA analysis failed to produce results")
  }

  return(results)
}
