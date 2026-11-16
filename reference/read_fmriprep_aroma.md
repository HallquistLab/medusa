# Read ICA-AROMA mixing matrices and component classifications for MEDuSA denoising

Ingests `fMRIPost-AROMA` (`*_desc-melodic_mixing.tsv` +
`*_desc-aroma_metrics.tsv`), legacy fMRIPrep `<= 23.0`
(`*_desc-MELODIC_mixing.tsv` + `*_AROMAnoiseICs.csv`), or in-memory
MELODIC mixing matrices and component classifications into a run-aligned
`"medusa_aroma"` bundle for
[`apply_aroma`](https://hallquistlab.github.io/medusa/reference/apply_aroma.md),
[`read_fmriprep_confounds`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_confounds.md),
and
[`read_fmriprep_subject`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_subject.md).

## Usage

``` r
read_fmriprep_aroma(
  mixing,
  classification = NULL,
  bold = NULL,
  run_names = NULL,
  na_strings = c("n/a", "N/A", "NA", ".", "NaN", "")
)
```

## Arguments

- mixing:

  Character vector of MELODIC mixing file paths
  (`*_desc-melodic_mixing.tsv`, `*_desc-MELODIC_mixing.tsv`, `.tsv.gz`,
  `.csv`), a single numeric `T x K` matrix or `data.frame` (single run),
  or a list of per-run `T x K` matrices / `data.frame`s.

- classification:

  Component classification specification identifying rejected motion
  components:

  - `NULL` (default when `mixing` is file paths): auto-discovers sibling
    `*_desc-aroma_metrics.tsv(.gz)` or `*_AROMAnoiseICs.csv` files in
    the same directory as each mixing file.

  - Character vector of paths to `*_desc-aroma_metrics.tsv`
    (`classification == "rejected"`) or `*_AROMAnoiseICs.csv`
    (comma-separated 1-based component indices; empty file indicates
    zero rejected components).

  - Single 1-based integer vector (e.g., `c(2L, 5L)`), logical mask of
    length `K` (`TRUE` = rejected), or metrics `data.frame` with a
    `classification` column (for single-run inputs), or a per-run list
    of any of these.

- bold:

  Optional `"medusa_nifti_ts"` object, named list of run BOLD matrices,
  single BOLD matrix, or character vector of run names used to validate
  acquisition entities, align run ordering, and verify volume counts.

- run_names:

  Optional explicit character vector of unique run names when `bold` is
  `NULL`.

- na_strings:

  Character vector of missing-value tokens passed to
  [`data.table::fread`](https://rdrr.io/pkg/data.table/man/fread.html).
  Default `c("n/a", "N/A", "NA", ".", "NaN", "")`.

## Value

An S3 list of class `"medusa_aroma"` containing:

- `mixing`: Run-named list of finite `T x K` numeric MELODIC mixing
  matrices.

- `noise_ics`: Run-named list of 1-based integer vectors of rejected
  motion component indices.

- `signal_ics`: Run-named list of 1-based integer vectors of accepted
  signal component indices.

- `metrics`: Run-named list of per-component classification
  `data.frame`s.

- `summary`: Per-run `data.frame` summarizing volume counts and
  accepted/rejected component counts.

## Details

Both headerless MELODIC mixing TSVs (standard FSL/`fMRIPost-AROMA`
format) and headered TSV/CSV tables are detected automatically. BIDS
file inputs are matched to `bold` runs using acquisition entities
(`sub-`, `ses-`, `task-`, `acq-`, `dir-`, `run-`) without requiring
template `space-` or `res-` tags to match, since `fMRIPost-AROMA`
derives run-level MELODIC timecourses in `MNI152NLin6Asym` (`res-2`)
regardless of the analysis output space.

## See also

[`apply_aroma`](https://hallquistlab.github.io/medusa/reference/apply_aroma.md),
[`read_fmriprep_subject`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_subject.md),
[`read_fmriprep_confounds`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_confounds.md)

## Examples

``` r
f_mix <- tempfile("sub-01_task-mid_run-01_", fileext = "_desc-melodic_mixing.tsv")
f_met <- sub("_desc-melodic_mixing\\.tsv$", "_desc-aroma_metrics.tsv", f_mix)
mix_mat <- matrix(rnorm(30 * 4), nrow = 30, ncol = 4)
utils::write.table(mix_mat, f_mix, sep = "\t", row.names = FALSE, col.names = FALSE)
met_df <- data.frame(
  IC = 1:4,
  classification = c("accepted", "rejected", "accepted", "rejected")
)
utils::write.table(met_df, f_met, sep = "\t", row.names = FALSE, quote = FALSE)

ar <- read_fmriprep_aroma(f_mix)
ar
#> <medusa_aroma> 1 run | 2/4 total component(s) classified as rejected
#>   run n_volumes n_components n_accepted n_rejected pct_rejected
#>  run1        30            4          2          2           50
unlink(c(f_mix, f_met))
```
