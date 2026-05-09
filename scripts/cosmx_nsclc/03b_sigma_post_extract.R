#!/usr/bin/env Rscript
# Post-process the sigma sweep: re-extract per-sigma weights and cosines
# from the saved RDS. Standalone so we don't need to rerun the kernels/CCA
# if the cosine logic changes.
#
# Usage:
#   Rscript scripts/cosmx_nsclc/03b_sigma_post_extract.R [slide]

suppressPackageStartupMessages({
  library(CoPro)
})

args <- commandArgs(trailingOnly = TRUE)
SLIDE <- if (length(args) >= 1L) args[[1]] else "Lung5_Rep1"

DATA_ROOT <- "/home/zmiao2/spatial/cosmx_nsclc"
OUT_DIR   <- file.path(DATA_ROOT, "copro_results", paste0(SLIDE, "_sigma_sweep"))
RDS_PATH  <- file.path(OUT_DIR, sprintf("%s_sigma_sweep.rds", SLIDE))

obj <- readRDS(RDS_PATH)

# Keys are "geneScores|sigma<value>|<celltype>" (no underscore between
# 'sigma' and the value).
parse_key <- function(k) {
  parts <- strsplit(k, "\\|", perl = TRUE)[[1]]
  list(sigma = as.numeric(sub("^sigma", "", parts[2])),
       cell_type = parts[3])
}

keys <- names(obj@geneScores)
parsed <- lapply(keys, parse_key)
SIGMAS <- sort(unique(vapply(parsed, function(p) p$sigma, numeric(1))))
N_CC   <- ncol(obj@geneScores[[keys[1]]])

cat(sprintf("sigmas: %s; nCC: %d\n",
            paste(SIGMAS, collapse = ", "), N_CC))

gs_long <- do.call(rbind, lapply(seq_along(keys), function(i) {
  pk <- parsed[[i]]
  m  <- obj@geneScores[[keys[i]]]
  do.call(rbind, lapply(seq_len(N_CC), function(cc) {
    data.frame(
      sigma = pk$sigma,
      cell_type = pk$cell_type,
      cc = cc,
      gene = rownames(m),
      weight = m[, paste0("CC_", cc)],
      row.names = NULL
    )
  }))
}))

stopifnot(all(!is.na(gs_long$sigma)))

# Top-30 per (sigma, cell_type, cc)
top_n <- 30
gs_top <- do.call(rbind, lapply(split(gs_long,
                                       list(gs_long$sigma, gs_long$cell_type, gs_long$cc),
                                       drop = TRUE),
  function(df) head(df[order(-abs(df$weight)), , drop = FALSE], top_n)))
write.csv(gs_top,  file.path(OUT_DIR, sprintf("%s_topgene_weights_per_sigma.csv", SLIDE)),
          row.names = FALSE)
write.csv(gs_long, file.path(OUT_DIR, sprintf("%s_full_weights_per_sigma.csv", SLIDE)),
          row.names = FALSE)

# All-pairs cosine of CC weight vectors per (cell_type, cc).
# Sign of the eigenvector is arbitrary, so use abs(cosine).
cell_types <- sort(unique(gs_long$cell_type))
cosine_rows <- list()
for (ct in cell_types) {
  for (cc in seq_len(N_CC)) {
    sub <- gs_long[gs_long$cell_type == ct & gs_long$cc == cc, ]
    W <- reshape(sub[, c("gene", "sigma", "weight")],
                 idvar = "gene", timevar = "sigma", direction = "wide")
    rownames(W) <- W$gene; W$gene <- NULL
    sigma_cols <- as.numeric(sub("^weight\\.", "", colnames(W)))
    ord <- order(sigma_cols)
    W <- as.matrix(W[, ord, drop = FALSE]); sigma_cols <- sigma_cols[ord]

    norms <- sqrt(colSums(W^2))
    M <- abs((t(W) %*% W) / outer(norms, norms))

    # Spearman-style cosine on |weight| ranks (rank-stability)
    R <- apply(W, 2, function(v) rank(-abs(v)))
    R <- scale(R, center = TRUE, scale = FALSE)
    rnorms <- sqrt(colSums(R^2))
    M_rank <- (t(R) %*% R) / outer(rnorms, rnorms)

    for (i in seq_along(sigma_cols)) {
      for (j in seq_along(sigma_cols)) {
        cosine_rows[[length(cosine_rows) + 1L]] <- data.frame(
          cell_type = ct, cc = cc,
          sigma_i = sigma_cols[i], sigma_j = sigma_cols[j],
          cosine_weights = M[i, j],
          cosine_ranks   = M_rank[i, j],
          row.names = NULL
        )
      }
    }
  }
}
cosine_df <- do.call(rbind, cosine_rows)
write.csv(cosine_df, file.path(OUT_DIR, sprintf("%s_sigma_cosine.csv", SLIDE)),
          row.names = FALSE)

# Adjacent-sigma cosines (the acceptance check)
ord_sigma <- sort(SIGMAS)
adj_rows <- list()
for (ct in cell_types) {
  for (cc in seq_len(N_CC)) {
    for (k in seq_len(length(ord_sigma) - 1L)) {
      mask <- cosine_df$cell_type == ct & cosine_df$cc == cc &
        cosine_df$sigma_i == ord_sigma[k] &
        cosine_df$sigma_j == ord_sigma[k + 1L]
      adj_rows[[length(adj_rows) + 1L]] <-
        cosine_df[mask, , drop = FALSE]
    }
  }
}
adjacent_df <- do.call(rbind, adj_rows)
adjacent_df$pass_0p9 <- adjacent_df$cosine_weights >= 0.9
write.csv(adjacent_df, file.path(OUT_DIR, sprintf("%s_sigma_adjacent.csv", SLIDE)),
          row.names = FALSE)

cat("\n=== Adjacent-sigma cosine (acceptance: cosine_weights >= 0.9) ===\n")
print(adjacent_df, row.names = FALSE, digits = 4)
cat(sprintf("\nPasses: %d / %d\n",
            sum(adjacent_df$pass_0p9, na.rm = TRUE), nrow(adjacent_df)))

# Per (cell_type, cc) min-cosine across adjacent pairs
agg <- aggregate(cosine_weights ~ cell_type + cc, data = adjacent_df,
                 FUN = function(x) min(x, na.rm = TRUE))
names(agg)[names(agg) == "cosine_weights"] <- "min_adjacent_cosine"
cat("\n=== Per (cell_type, CC) minimum adjacent-pair cosine ===\n")
print(agg, row.names = FALSE, digits = 4)

# Pretty-print full pairwise matrix per (cell_type, cc)
cat("\n=== Full pairwise cosine matrices ===\n")
for (ct in cell_types) {
  for (cc in seq_len(N_CC)) {
    cat(sprintf("\n%s, CC%d:\n", ct, cc))
    sub <- cosine_df[cosine_df$cell_type == ct & cosine_df$cc == cc, ]
    M <- reshape(sub[, c("sigma_i", "sigma_j", "cosine_weights")],
                 idvar = "sigma_i", timevar = "sigma_j", direction = "wide")
    rownames(M) <- M$sigma_i; M$sigma_i <- NULL
    colnames(M) <- sub("^cosine_weights\\.", "sigma=", colnames(M))
    print(round(as.matrix(M), 3))
  }
}
