# Gene-Space Average Per-Slide Canonical Correlation Optimization

These functions implement the Proposal 1b algorithm: batch-robust CCA in
gene space using average per-slide canonical correlation. Each slide's
contribution is normalized by its own score variance, preventing batch
axes from inflating the objective.
