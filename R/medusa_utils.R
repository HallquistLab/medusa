# Internal spatial geometry, filename, duration, and batch-execution helpers

roi_image_geometry <- function(image) {
  header <- RNifti::niftiHeader(image)
  list(
    dim = as.integer(header$dim[2:4]),
    volumes = if (header$dim[1L] > 3L) prod(header$dim[5L:(header$dim[1L] + 1L)]) else 1,
    spacing = as.numeric(header$pixdim[2:4]),
    units = bitwAnd(as.integer(header$xyzt_units), 7L),
    affine = matrix(as.numeric(RNifti::xform(header, useQuaternionFirst = FALSE)), 4L, 4L)
  )
}

validate_roi_geometry <- function(mask, image, mask_file, image_file, tolerance = 1e-4) {
  close <- function(x, y) {
    length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
      all(abs(x - y) <= tolerance)
  }
  compatible <- identical(mask$dim, image$dim) && image$volumes == 1 &&
    identical(mask$units, image$units) && close(mask$spacing, image$spacing) &&
    close(mask$affine, image$affine)
  if (!isTRUE(compatible)) {
    stop("ROI mask and statistic image have incompatible spatial geometry (dimensions, voxel sizes, ",
      "units, or voxel-to-world transform). Resample/register the mask to the statistic image grid ",
      "before extraction. Mask: ", mask_file, "; image: ", image_file, call. = FALSE)
  }
  invisible(TRUE)
}

file_sans_ext <- function(file, withdot = TRUE) {
  ext_pat <- "^(.*)\\.(csv|dat|hdr|img|brik|head|nii|txt|tsv|yaml|json|1d|rds|parquet|pq)(\\.gz|\\.bz2|\\.zip|\\.xz)*$"
  matches <- grepl(ext_pat, file, ignore.case = TRUE)
  fout <- rep(NA_character_, length = length(file))
  fout[matches] <- sub(ext_pat, "\\1", file[matches], ignore.case = TRUE)
  fout
}

validate_dhms <- function(str) {
  checkmate::assert_string(str)
  if (grepl("^\\d+:\\d+?$", str, perl = TRUE)) {
    return(paste0("00:", str))
  } else if (grepl("^(\\d+-)?\\d+:\\d+:\\d+?$", str, perl = TRUE)) {
    return(str)
  } else if (grepl("^(\\d+-)?\\d+:\\d+?$", str, perl = TRUE)) {
    return(paste0(str, ":00"))
  } else if (grepl("^\\d+-\\d+$", str, perl = TRUE)) {
    return(paste0(str, ":00:00"))
  } else if (grepl("^\\d+$", str, perl = TRUE)) {
    return(paste0("00:", str, ":00"))
  } else {
    stop("Invalid duration string: ", str, call. = FALSE)
  }
}

named_list <- function(...) {
  vnames <- as.character(match.call())[-1L]
  stats::setNames(list(...), vnames)
}

create_decon_batch_job <- function(
  job_name,
  n_cpus = 1L,
  mem_per_cpu = "4g",
  wall_time = "01:00:00",
  scheduler = "local",
  input_objects = list(),
  r_packages = "medusa",
  r_code = character()
) {
  if (identical(scheduler, "local")) {
    return(list(
      submit = function() {
        eval_env <- list2env(input_objects, parent = asNamespace("medusa"))
        for (expr_txt in r_code) {
          eval(parse(text = expr_txt), envir = eval_env)
        }
        invisible(NULL)
      }
    ))
  }

  if (requireNamespace("fmri.pipeline", quietly = TRUE) &&
      exists("R_batch_job", envir = asNamespace("fmri.pipeline"), inherits = FALSE)) {
    rbj <- get("R_batch_job", envir = asNamespace("fmri.pipeline"), inherits = FALSE)
    return(rbj$new(
      job_name = job_name,
      n_cpus = n_cpus,
      mem_per_cpu = mem_per_cpu,
      wall_time = wall_time,
      scheduler = scheduler,
      input_objects = input_objects,
      r_packages = r_packages,
      r_code = r_code
    ))
  }

  list(
    submit = function() {
      work_dir <- if (!is.null(input_objects$out_file)) dirname(input_objects$out_file) else tempdir()
      dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
      safe_name <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(job_name))
      rds_file <- file.path(work_dir, paste0(safe_name, "_inputs.rds"))
      r_file <- file.path(work_dir, paste0(safe_name, ".R"))
      saveRDS(input_objects, rds_file)
      pkg_lines <- sprintf("library(%s)", r_packages)
      script_lines <- c(
        pkg_lines,
        sprintf("list2env(readRDS(%s), envir = environment())", deparse(rds_file)),
        r_code
      )
      writeLines(script_lines, r_file)
      if (identical(scheduler, "slurm")) {
        sh_file <- file.path(work_dir, paste0(safe_name, ".sbatch"))
        writeLines(c(
          "#!/bin/bash",
          sprintf("#SBATCH --job-name=%s", safe_name),
          sprintf("#SBATCH --cpus-per-task=%d", as.integer(n_cpus)),
          sprintf("#SBATCH --mem-per-cpu=%s", mem_per_cpu),
          sprintf("#SBATCH --time=%s", wall_time),
          sprintf("Rscript %s", shQuote(r_file))
        ), sh_file)
        system2("sbatch", sh_file)
      } else if (identical(scheduler, "torque")) {
        sh_file <- file.path(work_dir, paste0(safe_name, ".pbs"))
        writeLines(c(
          "#!/bin/bash",
          sprintf("#PBS -N %s", safe_name),
          sprintf("#PBS -l nodes=1:ppn=%d,pmem=%s,walltime=%s", as.integer(n_cpus), mem_per_cpu, wall_time),
          sprintf("Rscript %s", shQuote(r_file))
        ), sh_file)
        system2("qsub", sh_file)
      } else {
        stop("Unsupported scheduler: ", scheduler, call. = FALSE)
      }
      invisible(NULL)
    }
  )
}
