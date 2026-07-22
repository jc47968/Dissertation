# ASD Brain Connectivity Analysis Pipeline

`Complete Code_Final.R` runs the full analysis pipeline top to bottom:
phenotype prep → connectivity correlation → merge → PCA → dataset creation →
exploratory data analysis → brain-network visuals → MANCOVA /
permutation-MANCOVA → SVM classification → SVM-based MANCOVA / post-hoc
tests. It is organized into 27 numbered stages, each printed to the console
as it starts (`STAGE N: ...`).

## 0. Generate the connectivity data (only if you don't already have it)

1. Run `ABIDE_Data_Extract.R` (in `Desktop\NCU\Dissertation Dataset\`). It
   reads the download URLs listed in `URL_ABIDE_Preprocessing_TS_Data.xlsx`
   (one worksheet per site) and downloads each subject's raw time-series
   file, saving them as individual CSVs into the `CC400_CPAC/` folder.
2. The main script's own Stage 2 then combines those per-subject CSVs into
   a single `CC400_combined.csv` automatically the first time it runs (it
   checks whether the file already exists and only does this combine step
   if it doesn't).

This only needs to be done once — after that, `CC400_combined.csv` can be
reused for every future run.

## 1. Required data files

These are the only files you need to supply — everything else the pipeline
needs is generated automatically as it runs. All three are already collected
in the `Data/` folder here:

| File | Where it must be | Used in |
|---|---|---|
| `Phenotypic_V1_0b_v1.csv` | `Data/` | Stage 1 |
| `CC400_combined.csv` (~840 MB) | `Data/` | Stage 3 (see Step 0 above if missing) |
| `CC400_ROI_labels.csv` | `Data/` | Stages 9–12 (brain-network visuals) |

## 2. Configure the paths

Near the top of the script is a CONFIG block:

```r
setwd(r"[...]")             # working directory — all pipeline files are read/written here
RAW_CC400_DIR <- r"[...]"    # folder of raw per-subject CC400 CSVs (only needed for Step 0)
CC400_COMBINED_FILE <- r"[...]"  # combined CC400 file
ROI_LABELS_FILE <- r"[...]"  # CC400 ROI label lookup table
```

**Check these before running** — as currently written, these four paths
still point at an older folder layout, not the `Data/` folder in this
project. Update them to match where your files actually live now, e.g.:

```r
setwd(r"[C:\Users\jason\OneDrive\NCU\Dissertation Dataset\Disseration Analysis\Full Stack_Final]")
CC400_COMBINED_FILE <- r"[C:\Users\jason\OneDrive\NCU\Dissertation Dataset\Disseration Analysis\Full Stack_Final\Data\CC400_combined.csv]"
ROI_LABELS_FILE <- r"[C:\Users\jason\OneDrive\NCU\Dissertation Dataset\Disseration Analysis\Full Stack_Final\Data\CC400_ROI_labels.csv]"
```

`RAW_CC400_DIR` only matters if `CC400_COMBINED_FILE` doesn't exist yet (see
Step 0).

Once these are correct, everything else the script reads or writes uses
relative filenames within the working directory set here.

## 3. Required R packages

Install these before running (all via `install.packages()`):

```
dplyr, data.table, readr, reshape2, ggplot2, tidyr, corrplot, GPArotation,
psych, igraph, ggraph, tidygraph, ggimage, tidyverse, MVN, car, biotools,
vegan, e1071, caret, pROC, multcomp, rstatix
```

## 4. Running the script

```r
source("Complete Code_Final.R")
```

or from the command line:

```
Rscript "Complete Code_Final.R"
```

**Expect a long runtime** if Step 0 or Stage 3 has to run from scratch — the
correlation analysis computes pairwise correlations across roughly 76,600
connectivity variables for over 1,000 subjects, the PCA runs on the same
scale, four separate stages each run a 10,000-permutation MANCOVA, and the
SVM stage performs a 121-combination grid search with 5-fold
cross-validation. A full run can take several hours.

If a stage fails (for example, an input file isn't where it's expected), the
script prints `STAGE N FAILED: <error message>` and continues on to the next
stage rather than stopping the whole run. Check the console output at the
end for any such messages, and re-run after fixing the underlying issue —
you don't need to re-run the whole script, just make sure the stage's inputs
are present before continuing.

## 5. Things to check before trusting the output

- **Component loading files (Stages 9–12):** these stages read
  `component_<N>_loadings.csv` files out of the `Varimax Factor Importance/`
  and `Promax Factor Importance/` subfolders that Stage 5 creates. Confirm
  the component number in each stage matches the component you intend to
  visualize (PC4, PC17, PC15, etc.).
- **Stages 22 and 24** (SVM MANCOVA, "Model" versions) build their formula
  around `DX_GROUP` (actual diagnosis) rather than the SVM's predicted class
  column (`pred`), even though `pred` is created earlier in the stage. If you
  intend to test against the model's *predicted* labels rather than actual
  diagnosis, update the formula in these two stages before running.
- **Existing output files:** if `varimax_alldata.csv` and
  `promax_alldata.csv` already exist in this folder from a previous run,
  Stage 6 will overwrite them. Back them up first if you want to keep the
  existing versions.

## 6. Outputs

Each stage writes its own CSVs and/or PNG plots directly into the working
directory, except the per-component loading CSVs from Stage 5, which are
written into `Varimax Factor Importance/` and `Promax Factor Importance/`
subfolders (kept separate because both rotations produce identically-named
`component_<N>_loadings.csv` files).
