# CoPro Analysis Scripts

Reproducibility scripts for the analyses in Miao et al.

These scripts reproduce the main results and figures from the manuscript.
They are organized by dataset. Each script has a header documenting its
inputs, outputs, and dependencies.

**Note:** Before running these scripts, update the `DATA_DIR` variable at
the top of each script to point to your local data directory. The full
datasets are available on Zenodo (DOI: TBD).

## Directory structure

```
scripts/
  colon/
    01_run_copro_D3.R             Colon Day 3 (cross-type, orthogonal axes)
    02_run_copro_D9.R             Colon Day 9 (multi-slide pipeline)
    03_D9_cross_slide_transfer.R  Score transfer between D9 slides
    04_D9_LMM_analysis.R          Linear mixed model analysis (D9)
    05_gene_weights_comparison.R  PCA vs regression gene weights (D0)
  kidney/
    01_run_copro_supervised.R     Supervised nephron-vascular axis
    02_scrnaseq_transfer.R        Transfer scores to scRNA-seq
    03_downstream_GO.R            GO enrichment and MA plots
  liver/
    01_run_copro.R                CoPro on Visium liver data
    02_zonation_transfer.R        Zonation score transfer
  brain/
    01_run_copro_D1_D2.R          D1/D2 neuron co-progression (MERFISH)
  simulation/
    01_generate_data.R            Generate simulated spatial data
    02_run_copro.R                Run CoPro on simulated data
    03_run_alternative_methods.R  DIALOGUE comparison
    04_evaluate_and_plot.R        Method comparison figures
```
