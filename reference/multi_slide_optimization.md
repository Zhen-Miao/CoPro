# Multi-slide optimization functions for CoPro

These functions extend the CoPro formulation to handle multiple slides
where weight vectors are shared across slides.

## Details

The objective is: \$\$Maximize\_{\\w_i,w_j\\} \sum_q w_i^T X\_{i,q}^T
K\_{ij,q} X\_{j,q} w_j\$\$ subject to \\\|\|w_i\|\| \leq 1\\ for all
\\i\\, where \\q\\ is the slide index.
