suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(jsonlite)
  library(susieR)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(file.path(getwd()), mustWork = TRUE)
}

table_dir <- file.path(root, "results", "tables")
variant_path <- file.path(table_dir, "prkcb_ld_analysis_variants.csv")
ld_path <- file.path(table_dir, "prkcb_ld_matrix.csv")

variants <- fread(variant_path)
ld <- as.matrix(fread(ld_path), rownames = 1)
variants <- as.data.frame(variants)
rownames(variants) <- variants$rsid
if (!identical(colnames(ld), variants$rsid)) {
  stop("LD matrix columns are not aligned with PRKCB variant rows.")
}

symmetrize <- function(matrix) {
  matrix <- (matrix + t(matrix)) / 2
  eigen_values <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  minimum <- min(eigen_values)
  if (minimum <= 1e-6) {
    matrix <- matrix + diag(abs(minimum) + 1e-6, nrow(matrix))
  }
  matrix
}

ld <- symmetrize(ld)
colnames(ld) <- rownames(ld) <- variants$rsid

gwas <- list(
  beta = variants$gwas_beta,
  varbeta = variants$gwas_se^2,
  snp = variants$rsid,
  position = variants$position,
  N = 398668,
  type = "cc",
  s = 45106 / 398668,
  LD = ld
)

eqtl <- list(
  beta = variants$eql_beta,
  varbeta = variants$eql_se^2,
  snp = variants$rsid,
  position = variants$position,
  MAF = variants$maf_eur,
  N = 191,
  type = "quant",
  LD = ld
)

credible_set_rows <- list()
coloc_rows <- list()
fit_rows <- list()

for (L in c(1, 3, 5)) {
  message("Running SuSiE with L=", L)
  gwas_susie <- runsusie(
    gwas,
    maxit = 500,
    repeat_until_convergence = TRUE,
    L = L,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    check_prior = FALSE,
    max_iter = 500
  )
  eqtl_susie <- runsusie(
    eqtl,
    maxit = 500,
    repeat_until_convergence = TRUE,
    L = L,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    check_prior = FALSE,
    max_iter = 500
  )

  fit_rows[[length(fit_rows) + 1]] <- data.frame(
    L = L,
    gwas_converged = isTRUE(gwas_susie$converged),
    eqtl_converged = isTRUE(eqtl_susie$converged),
    gwas_credible_sets = length(gwas_susie$sets$cs),
    eqtl_credible_sets = length(eqtl_susie$sets$cs),
    gwas_min_abs_credible_pip = if (length(gwas_susie$sets$cs)) {
      min(gwas_susie$sets$purity$min.abs.corr)
    } else {
      NA_real_
    },
    eqtl_min_abs_credible_pip = if (length(eqtl_susie$sets$cs)) {
      min(eqtl_susie$sets$purity$min.abs.corr)
    } else {
      NA_real_
    }
  )

  for (trait_name in c("GWAS", "eQTL")) {
    fit <- if (trait_name == "GWAS") gwas_susie else eqtl_susie
    if (is.null(fit$sets$cs) || !length(fit$sets$cs)) {
      next
    }
    for (cs_index in seq_along(fit$sets$cs)) {
      members <- fit$sets$cs[[cs_index]]
      member_snps <- variants$rsid[members]
      member_pip <- fit$pip[members]
      credible_set_rows[[length(credible_set_rows) + 1]] <- data.frame(
        L = L,
        trait = trait_name,
        credible_set = cs_index,
        n_snps = length(members),
        top_snp = variants$rsid[members[which.max(member_pip)]],
        top_pip = max(member_pip),
        credible_set_coverage = fit$sets$coverage[cs_index],
        snps = paste(member_snps, collapse = ";")
      )
    }
  }

  if (length(gwas_susie$sets$cs) && length(eqtl_susie$sets$cs)) {
    for (p12 in c(1e-6, 5e-6, 1e-5)) {
      coloc_result <- coloc.susie(
        gwas_susie,
        eqtl_susie,
        p1 = 1e-4,
        p2 = 1e-4,
        p12 = p12,
        overlap.min = 0.5
      )
      result_summary <- as.data.frame(coloc_result$summary)
      if (!nrow(result_summary)) {
        next
      }
      result_summary$L <- L
      result_summary$p12 <- p12
      coloc_rows[[length(coloc_rows) + 1]] <- result_summary
    }
  }
}

fit_table <- rbindlist(fit_rows, fill = TRUE)
credible_table <- if (length(credible_set_rows)) {
  rbindlist(credible_set_rows, fill = TRUE)
} else {
  data.table()
}
coloc_table <- if (length(coloc_rows)) {
  rbindlist(coloc_rows, fill = TRUE)
} else {
  data.table()
}

fwrite(
  fit_table,
  file.path(table_dir, "prkcb_susie_fit_quality.csv")
)
fwrite(
  credible_table,
  file.path(table_dir, "prkcb_susie_credible_sets.csv")
)
fwrite(
  coloc_table,
  file.path(table_dir, "prkcb_susie_coloc_sensitivity.csv")
)

primary_rows <- coloc_table[
  L == 5 & p12 == 5e-6,
]
primary_h4 <- if (nrow(primary_rows)) {
  suppressWarnings(max(primary_rows$PP.H4.abf, na.rm = TRUE))
} else {
  NA_real_
}
if (!is.finite(primary_h4)) {
  primary_h4 <- NA_real_
}

summary <- list(
  method = "SuSiE-coloc with 1000 Genomes EUR LD",
  variant_region = "chr16:23.70-23.95 Mb",
  variants_used = nrow(variants),
  gwas_samples = 398668,
  eqtl_samples = 191,
  credible_sets_gwas_primary = if (nrow(credible_table)) {
    sum(credible_table$L == 5 & credible_table$trait == "GWAS")
  } else {
    0
  },
  credible_sets_eqtl_primary = if (nrow(credible_table)) {
    sum(credible_table$L == 5 & credible_table$trait == "eQTL")
  } else {
    0
  },
  primary_L = 5,
  primary_p12 = 5e-6,
  primary_max_H4 = primary_h4,
  primary_shared_credible_pair = if (
    nrow(primary_rows) && is.finite(primary_h4)
  ) {
    primary_rows[which.max(primary_rows$PP.H4.abf), c("hit1", "hit2")]
  } else {
    c(NA_character_, NA_character_)
  },
  supported_H4_ge_0_8 = is.finite(primary_h4) && primary_h4 >= 0.8,
  sensitivity_max_H4 = if (nrow(coloc_table)) {
    suppressWarnings(max(coloc_table$PP.H4.abf, na.rm = TRUE))
  } else {
    NA_real_
  },
  note = paste(
    "This is fine-mapped colocalization with an external LD reference.",
    "Shared credible sets support a shared causal variant hypothesis but",
    "do not establish the molecular mechanism experimentally."
  )
)

write_json(
  summary,
  file.path(table_dir, "prkcb_susie_coloc_summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)
print(summary)
print(head(credible_table, 20))
print(head(coloc_table, 20))
