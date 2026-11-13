#define RNIFTI_NIFTILIB_VERSION 2

#include <Rcpp.h>
#include "RNifti.h"
#include "RNiftiAPI.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr size_t kTileSize = 64;

template <typename T>
inline double elem_to_double(T val) {
  return static_cast<double>(val);
}

template <>
inline double elem_to_double<int32_t>(int32_t val) {
  return (val == NA_INTEGER) ? NA_REAL : static_cast<double>(val);
}

template <>
inline double elem_to_double<double>(double val) {
  return val;
}

template <typename T, bool Scaled>
inline double read_scaled_elem(T raw, double slope, double inter) {
  const double val = elem_to_double<T>(raw);
  if (!Scaled) {
    return val;
  }
  return R_IsNA(val) ? NA_REAL : (val * slope + inter);
}

template <typename T, bool Scaled>
void extract_typed_loop(
    const T* base,
    size_t n_spatial,
    size_t n_time,
    const std::vector<size_t>& indices,
    double slope,
    double inter,
    bool time_first,
    double* out_ptr
) {
  const size_t n_vox = indices.size();
  if (n_vox == 0 || n_time == 0) {
    return;
  }

  if (!time_first) {
    for (size_t t = 0; t < n_time; ++t) {
      if ((t & 63) == 0) Rcpp::checkUserInterrupt();
      const T* vol_ptr = base + t * n_spatial;
      double* col_ptr = out_ptr + t * n_vox;
      for (size_t v = 0; v < n_vox; ++v) {
        col_ptr[v] = read_scaled_elem<T, Scaled>(vol_ptr[indices[v]], slope, inter);
      }
    }
    return;
  }

  double tile[kTileSize * kTileSize];
  for (size_t t0 = 0; t0 < n_time; t0 += kTileSize) {
    Rcpp::checkUserInterrupt();
    const size_t tb = std::min(kTileSize, n_time - t0);
    for (size_t v0 = 0; v0 < n_vox; v0 += kTileSize) {
      const size_t vb = std::min(kTileSize, n_vox - v0);
      for (size_t dt = 0; dt < tb; ++dt) {
        const T* vol_ptr = base + (t0 + dt) * n_spatial;
        for (size_t dv = 0; dv < vb; ++dv) {
          tile[dv * tb + dt] = read_scaled_elem<T, Scaled>(vol_ptr[indices[v0 + dv]], slope, inter);
        }
      }
      for (size_t dv = 0; dv < vb; ++dv) {
        std::memcpy(out_ptr + (v0 + dv) * n_time + t0, &tile[dv * tb], tb * sizeof(double));
      }
    }
  }
}

template <typename T>
void extract_typed_dispatch_scale(
    const void* blob,
    size_t n_spatial,
    size_t n_time,
    const std::vector<size_t>& indices,
    bool scaled,
    double slope,
    double inter,
    bool time_first,
    double* out_ptr
) {
  const T* base = static_cast<const T*>(blob);
  if (scaled) {
    extract_typed_loop<T, true>(base, n_spatial, n_time, indices, slope, inter, time_first, out_ptr);
  } else {
    extract_typed_loop<T, false>(base, n_spatial, n_time, indices, 1.0, 0.0, time_first, out_ptr);
  }
}

RNifti::NiftiImage load_reference_header(SEXP reference) {
  if (Rf_isString(reference)) {
    if (Rf_length(reference) != 1) {
      Rcpp::stop("reference file path must be a single string");
    }
    const std::string ref_path = Rcpp::as<std::string>(reference);
    return RNifti::NiftiImage(ref_path, false);
  }
  return RNifti::NiftiImage(reference, false, true);
}

int resolve_output_datatype(const std::string& datatype) {
  if (datatype == "float" || datatype == "float32") {
    return DT_FLOAT32;
  }
  if (datatype == "double" || datatype == "float64") {
    return DT_FLOAT64;
  }
  Rcpp::stop("Unsupported output NIfTI datatype '%s'; use 'float' or 'double'", datatype.c_str());
  return DT_FLOAT32;
}

RNifti::NiftiImage prepare_output_image(
    const RNifti::NiftiImage& ref,
    size_t& n_spatial_out,
    size_t n_time,
    double tr,
    double time_offset,
    int dt_code,
    int& filetype_out
) {
  if (ref.isNull()) {
    Rcpp::stop("Reference NIfTI image is null or unreadable");
  }
  const auto dims = ref.dim();
  if (dims.size() < 2) {
    Rcpp::stop("Reference NIfTI image must have at least two spatial dimensions");
  }
  const size_t dim_x = static_cast<size_t>(dims[0]);
  const size_t dim_y = static_cast<size_t>(dims[1]);
  const size_t dim_z = (dims.size() >= 3) ? static_cast<size_t>(dims[2]) : 1;
  n_spatial_out = dim_x * dim_y * dim_z;

  if (n_time == 0) {
    Rcpp::stop("Cannot write a NIfTI image with zero volumes");
  }

  nifti_image* copied = nifti2_copy_nim_info(ref);
  if (copied == nullptr) {
    Rcpp::stop("Failed to copy NIfTI reference header");
  }
  RNifti::NiftiImage out(copied, false);
  nifti_image* nim = out;

  nim->dim[0] = (n_time > 1) ? 4 : 3;
  nim->dim[1] = static_cast<RNifti::NiftiImage::dim_t>(dim_x);
  nim->dim[2] = static_cast<RNifti::NiftiImage::dim_t>(dim_y);
  nim->dim[3] = static_cast<RNifti::NiftiImage::dim_t>(dim_z);
  nim->dim[4] = static_cast<RNifti::NiftiImage::dim_t>(n_time);
  for (int d = 5; d < 8; ++d) {
    nim->dim[d] = 1;
    nim->pixdim[d] = 0.0;
  }

  if (std::isfinite(tr) && tr > 0.0) {
    nim->pixdim[4] = static_cast<RNifti::NiftiImage::pixdim_t>(tr);
  } else if (n_time <= 1) {
    nim->pixdim[4] = 0.0;
  }

  if (std::isfinite(time_offset)) {
    nim->toffset = static_cast<RNifti::NiftiImage::scale_t>(time_offset);
  }

  nifti2_update_dims_from_array(nim);

  nim->datatype = dt_code;
  nifti_datatype_sizes(dt_code, &nim->nbyper, &nim->swapsize);
  nim->scl_slope = 0.0;
  nim->scl_inter = 0.0;

  const size_t total_elems = n_spatial_out * n_time;
  void* raw_buf = std::calloc(total_elems, static_cast<size_t>(nim->nbyper));
  if (raw_buf == nullptr && total_elems > 0) {
    Rcpp::stop("Failed to allocate NIfTI output buffer (%zu elements)", total_elems);
  }
  nim->data = raw_buf;

  const size_t max_dim = std::max(std::max(dim_x, dim_y), std::max(dim_z, n_time));
  if (max_dim > 32767) {
    filetype_out = NIFTI_FTYPE_NIFTI2_1;
  } else if (nim->nifti_type <= 0 || nim->nifti_type == NIFTI_FTYPE_NIFTI1_2) {
    filetype_out = NIFTI_FTYPE_NIFTI1_1;
  } else {
    filetype_out = -1;
  }

  return out;
}

template <typename OutT>
void fill_indexed_buffer(
    OutT* buf,
    size_t n_spatial,
    size_t n_time,
    const std::vector<size_t>& indices,
    const double* ts_ptr,
    bool time_first,
    double fill_value,
    double& cal_min_out,
    double& cal_max_out
) {
  const size_t n_vox = indices.size();
  const size_t total_elems = n_spatial * n_time;
  const OutT fill_cast = static_cast<OutT>(fill_value);

  // Buffer comes from calloc (all zero bytes). Only fill explicitly if nonzero/NaN.
  if (!(fill_value == 0.0 && !std::signbit(fill_value))) {
    std::fill_n(buf, total_elems, fill_cast);
  }

  bool has_finite = false;
  double cmin = std::numeric_limits<double>::infinity();
  double cmax = -std::numeric_limits<double>::infinity();

  if (n_vox < n_spatial) {
    const double fill_rt = static_cast<double>(fill_cast);
    if (std::isfinite(fill_rt)) {
      has_finite = true;
      cmin = fill_rt;
      cmax = fill_rt;
    }
  }

  if (!time_first) {
    for (size_t t = 0; t < n_time; ++t) {
      if ((t & 63) == 0) Rcpp::checkUserInterrupt();
      OutT* vol_ptr = buf + t * n_spatial;
      const double* col_ptr = ts_ptr + t * n_vox;
      for (size_t v = 0; v < n_vox; ++v) {
        const OutT val = static_cast<OutT>(col_ptr[v]);
        vol_ptr[indices[v]] = val;
        const double dval = static_cast<double>(val);
        if (std::isfinite(dval)) {
          has_finite = true;
          if (dval < cmin) cmin = dval;
          if (dval > cmax) cmax = dval;
        }
      }
    }
  } else {
    OutT tile[kTileSize * kTileSize];
    for (size_t t0 = 0; t0 < n_time; t0 += kTileSize) {
      Rcpp::checkUserInterrupt();
      const size_t tb = std::min(kTileSize, n_time - t0);
      for (size_t v0 = 0; v0 < n_vox; v0 += kTileSize) {
        const size_t vb = std::min(kTileSize, n_vox - v0);
        for (size_t dv = 0; dv < vb; ++dv) {
          const double* col_ptr = ts_ptr + (v0 + dv) * n_time + t0;
          for (size_t dt = 0; dt < tb; ++dt) {
            const OutT val = static_cast<OutT>(col_ptr[dt]);
            tile[dt * vb + dv] = val;
            const double dval = static_cast<double>(val);
            if (std::isfinite(dval)) {
              has_finite = true;
              if (dval < cmin) cmin = dval;
              if (dval > cmax) cmax = dval;
            }
          }
        }
        for (size_t dt = 0; dt < tb; ++dt) {
          OutT* vol_ptr = buf + (t0 + dt) * n_spatial;
          const OutT* row_tile = &tile[dt * vb];
          for (size_t dv = 0; dv < vb; ++dv) {
            vol_ptr[indices[v0 + dv]] = row_tile[dv];
          }
        }
      }
    }
  }

  cal_min_out = has_finite ? cmin : 0.0;
  cal_max_out = has_finite ? cmax : 0.0;
}

template <typename OutT>
void fill_parcel_buffer(
    OutT* buf,
    size_t n_spatial,
    size_t n_time,
    const int* slots,
    size_t n_parcels,
    const std::vector<uint8_t>& slot_present,
    bool has_background,
    const double* parcel_ptr,
    double fill_value,
    double& cal_min_out,
    double& cal_max_out
) {
  const OutT fill_cast = static_cast<OutT>(fill_value);
  std::vector<OutT> lut(n_parcels + 1);
  lut[0] = fill_cast;

  bool has_finite = false;
  double cmin = std::numeric_limits<double>::infinity();
  double cmax = -std::numeric_limits<double>::infinity();

  if (has_background) {
    const double fill_rt = static_cast<double>(fill_cast);
    if (std::isfinite(fill_rt)) {
      has_finite = true;
      cmin = fill_rt;
      cmax = fill_rt;
    }
  }

  for (size_t t = 0; t < n_time; ++t) {
    if ((t & 31) == 0) Rcpp::checkUserInterrupt();
    const double* col_ptr = parcel_ptr + t * n_parcels;
    for (size_t m = 0; m < n_parcels; ++m) {
      const OutT val = static_cast<OutT>(col_ptr[m]);
      lut[m + 1] = val;
      if (slot_present[m + 1]) {
        const double dval = static_cast<double>(val);
        if (std::isfinite(dval)) {
          has_finite = true;
          if (dval < cmin) cmin = dval;
          if (dval > cmax) cmax = dval;
        }
      }
    }

    OutT* vol_ptr = buf + t * n_spatial;
    for (size_t i = 0; i < n_spatial; ++i) {
      vol_ptr[i] = lut[static_cast<size_t>(slots[i])];
    }
  }

  cal_min_out = has_finite ? cmin : 0.0;
  cal_max_out = has_finite ? cmax : 0.0;
}

}  // namespace

//' Extract selected voxel time series in a specified order
//'
//' @param infile Path to a scalar NIfTI image, optionally gzipped.
//' @param voxel_indices One-based linear indices within the spatial volume.
//' @param time_first Logical; if \code{TRUE}, returns a \code{volumes x voxels}
//'   matrix directly via cache-blocked extraction without an intermediate
//'   \code{voxels x volumes} allocation. Defaults to \code{FALSE} for backward
//'   compatibility.
//' @return A numeric matrix (\code{voxels x volumes} when \code{time_first = FALSE},
//'   or \code{volumes x voxels} when \code{time_first = TRUE}) in the supplied
//'   voxel order.
//' @details Uses the RNifti C++ API with single-dispatch typed pointer access,
//'   including NIfTI slope/intercept scaling, without creating a full
//'   double-precision R array or per-element virtual dispatch.
//'   Spatial alignment with the atlas/mask must be checked by the caller.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix extract_nifti_voxel_ts(
    std::string infile,
    Rcpp::NumericVector voxel_indices,
    bool time_first = false
) {
  RNifti::NiftiImage image(infile, true);
  if (image.isNull()) Rcpp::stop("Could not read NIfTI image: %s", infile.c_str());
  const auto dims = image.dim();
  if (dims.size() < 2) Rcpp::stop("Input image must have at least two spatial dimensions");
  for (size_t d = 4; d < dims.size(); ++d) {
    if (dims[d] != 1) Rcpp::stop("Input must be a scalar 3D or 4D NIfTI image");
  }
  size_t n_spatial = static_cast<size_t>(dims[0]) * dims[1];
  if (dims.size() >= 3) n_spatial *= dims[2];
  const size_t n_time = image.nVoxels() / n_spatial;
  if (n_time > static_cast<size_t>(std::numeric_limits<int>::max()) ||
      voxel_indices.size() > std::numeric_limits<int>::max()) {
    Rcpp::stop("Selected time series exceed R matrix dimensions");
  }

  std::vector<size_t> indices(voxel_indices.size());
  for (R_xlen_t v = 0; v < voxel_indices.size(); ++v) {
    const double index = voxel_indices[v];
    if (!std::isfinite(index) || index != std::floor(index) || index < 1 || index > n_spatial) {
      Rcpp::stop("voxel_indices must be finite, whole, one-based indices within the spatial volume");
    }
    indices[v] = static_cast<size_t>(index - 1);
  }

  const RNifti::NiftiImageData data = image.data();
  if (data.isComplex() || data.isRgb()) Rcpp::stop("Input NIfTI data must be real-valued scalar data");
  if (data.length() != image.nVoxels()) Rcpp::stop("NIfTI image data are unavailable");

  const R_xlen_t out_nrow = time_first ? static_cast<R_xlen_t>(n_time) : voxel_indices.size();
  const R_xlen_t out_ncol = time_first ? voxel_indices.size() : static_cast<R_xlen_t>(n_time);
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(out_nrow, out_ncol);
  if (indices.empty() || n_time == 0) {
    return out;
  }

  const void* blob = data.blob();
  const bool scaled = data.isScaled();
  const double slope = data.slope;
  const double inter = data.intercept;
  double* out_ptr = REAL(out);

  if (blob != nullptr) {
    switch (data.datatype()) {
      case DT_FLOAT32:
        extract_typed_dispatch_scale<float>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      case DT_FLOAT64:
        extract_typed_dispatch_scale<double>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      case DT_INT16:
        extract_typed_dispatch_scale<int16_t>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      case DT_INT32:
        extract_typed_dispatch_scale<int32_t>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      case DT_UINT8:
        extract_typed_dispatch_scale<uint8_t>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      case DT_INT8:
        extract_typed_dispatch_scale<int8_t>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      case DT_UINT16:
        extract_typed_dispatch_scale<uint16_t>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      case DT_UINT32:
        extract_typed_dispatch_scale<uint32_t>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      case DT_INT64:
        extract_typed_dispatch_scale<int64_t>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      case DT_UINT64:
        extract_typed_dispatch_scale<uint64_t>(blob, n_spatial, n_time, indices, scaled, slope, inter, time_first, out_ptr);
        return out;
      default:
        break;
    }
  }

  // Fallback for any non-standard scalar type
  for (size_t t = 0; t < n_time; ++t) {
    Rcpp::checkUserInterrupt();
    const size_t offset = t * n_spatial;
    if (!time_first) {
      for (size_t v = 0; v < indices.size(); ++v) {
        out(v, t) = static_cast<double>(data[offset + indices[v]]);
      }
    } else {
      for (size_t v = 0; v < indices.size(); ++v) {
        out(t, v) = static_cast<double>(data[offset + indices[v]]);
      }
    }
  }
  return out;
}

//' Write indexed voxel time series or static maps directly to a 3D/4D NIfTI file
//'
//' @param reference File path to a reference NIfTI image or an in-memory NIfTI object.
//' @param out_file Output file path (\code{.nii} or \code{.nii.gz}).
//' @param ts_mat Numeric matrix (\code{time x voxels} if \code{time_first = TRUE},
//'   else \code{voxels x time}).
//' @param voxel_indices Unique one-based spatial linear indices matching the
//'   voxel dimension of \code{ts_mat}.
//' @param tr Repetition time in seconds stored in \code{pixdim[5]} for 4D outputs.
//' @param time_offset Acquisition start offset stored in \code{toffset}.
//' @param fill_value Value written to unindexed background voxels.
//' @param datatype Output NIfTI datatype (\code{"float"} or \code{"double"}).
//' @param compression Gzip compression level (\code{0}--\code{9}) for \code{.nii.gz}.
//' @param time_first Logical; whether \code{ts_mat} is \code{time x voxels}.
//' @keywords internal
// [[Rcpp::export]]
std::string write_nifti_indexed_cpp(
    SEXP reference,
    std::string out_file,
    Rcpp::NumericMatrix ts_mat,
    Rcpp::NumericVector voxel_indices,
    double tr = 1.0,
    double time_offset = 0.0,
    double fill_value = 0.0,
    std::string datatype = "float",
    int compression = 2,
    bool time_first = true
) {
  if (compression < 0 || compression > 9) {
    Rcpp::stop("compression must be an integer between 0 and 9");
  }
  const int dt_code = resolve_output_datatype(datatype);
  const size_t n_time = static_cast<size_t>(time_first ? ts_mat.nrow() : ts_mat.ncol());
  const size_t n_vox = static_cast<size_t>(time_first ? ts_mat.ncol() : ts_mat.nrow());
  if (static_cast<size_t>(voxel_indices.size()) != n_vox) {
    Rcpp::stop("Length of voxel_indices (%lld) must match the voxel dimension of ts_mat (%zu)",
               static_cast<long long>(voxel_indices.size()), n_vox);
  }

  RNifti::NiftiImage ref = load_reference_header(reference);
  size_t n_spatial = 0;
  int filetype = -1;
  RNifti::NiftiImage out = prepare_output_image(
      ref, n_spatial, n_time, tr, time_offset, dt_code, filetype);

  std::vector<size_t> indices(n_vox);
  std::vector<uint8_t> seen(n_spatial, 0);
  for (size_t v = 0; v < n_vox; ++v) {
    const double idx = voxel_indices[v];
    if (!std::isfinite(idx) || idx != std::floor(idx) || idx < 1 || idx > n_spatial) {
      Rcpp::stop("voxel_indices must be finite, whole, one-based indices within the spatial volume");
    }
    const size_t zero_idx = static_cast<size_t>(idx - 1);
    if (seen[zero_idx]) {
      Rcpp::stop("voxel_indices must not contain duplicate voxel positions");
    }
    seen[zero_idx] = 1;
    indices[v] = zero_idx;
  }

  nifti_image* nim = out;
  double cmin = 0.0;
  double cmax = 0.0;
  if (dt_code == DT_FLOAT32) {
    fill_indexed_buffer<float>(
        static_cast<float*>(nim->data), n_spatial, n_time, indices,
        REAL(ts_mat), time_first, fill_value, cmin, cmax);
  } else {
    fill_indexed_buffer<double>(
        static_cast<double*>(nim->data), n_spatial, n_time, indices,
        REAL(ts_mat), time_first, fill_value, cmin, cmax);
  }
  nim->cal_min = static_cast<RNifti::NiftiImage::scale_t>(cmin);
  nim->cal_max = static_cast<RNifti::NiftiImage::scale_t>(cmax);

  const auto paths = out.toFile(out_file, dt_code, filetype, compression);
  return paths.second;
}

//' Write parcel-level time series or statistics directly to a 3D/4D NIfTI file
//'
//' @param reference File path to a reference 3D NIfTI mask or an in-memory NIfTI object.
//' @param out_file Output file path (\code{.nii} or \code{.nii.gz}).
//' @param parcel_slot Integer vector of length \code{n_spatial} with \code{0} for
//'   background voxels and \code{1..M} indexing rows of \code{parcel_values}.
//' @param parcel_values Numeric matrix of dimensions \code{M x T} (\code{parcels x time}).
//' @param tr Repetition time in seconds stored in \code{pixdim[5]} for 4D outputs.
//' @param fill_value Value written to background voxels (\code{parcel_slot == 0}).
//' @param datatype Output NIfTI datatype (\code{"float"} or \code{"double"}).
//' @param compression Gzip compression level (\code{0}--\code{9}) for \code{.nii.gz}.
//' @keywords internal
// [[Rcpp::export]]
std::string write_nifti_parcel_cpp(
    SEXP reference,
    std::string out_file,
    Rcpp::IntegerVector parcel_slot,
    Rcpp::NumericMatrix parcel_values,
    double tr = 1.0,
    double fill_value = 0.0,
    std::string datatype = "float",
    int compression = 2
) {
  if (compression < 0 || compression > 9) {
    Rcpp::stop("compression must be an integer between 0 and 9");
  }
  const int dt_code = resolve_output_datatype(datatype);
  const size_t n_parcels = static_cast<size_t>(parcel_values.nrow());
  const size_t n_time = static_cast<size_t>(parcel_values.ncol());

  RNifti::NiftiImage ref = load_reference_header(reference);
  size_t n_spatial = 0;
  int filetype = -1;
  RNifti::NiftiImage out = prepare_output_image(
      ref, n_spatial, n_time, tr, 0.0, dt_code, filetype);

  if (static_cast<size_t>(parcel_slot.size()) != n_spatial) {
    Rcpp::stop("parcel_slot length (%lld) must match reference spatial voxel count (%zu)",
               static_cast<long long>(parcel_slot.size()), n_spatial);
  }

  const int* slots = INTEGER(parcel_slot);
  std::vector<uint8_t> slot_present(n_parcels + 1, 0);
  bool has_background = false;
  for (size_t i = 0; i < n_spatial; ++i) {
    const int s = slots[i];
    if (s == NA_INTEGER || s < 0 || static_cast<size_t>(s) > n_parcels) {
      Rcpp::stop("parcel_slot values must be non-missing integers in 0..%zu", n_parcels);
    }
    if (s == 0) {
      has_background = true;
    } else {
      slot_present[static_cast<size_t>(s)] = 1;
    }
  }

  nifti_image* nim = out;
  double cmin = 0.0;
  double cmax = 0.0;
  if (dt_code == DT_FLOAT32) {
    fill_parcel_buffer<float>(
        static_cast<float*>(nim->data), n_spatial, n_time, slots, n_parcels,
        slot_present, has_background, REAL(parcel_values), fill_value, cmin, cmax);
  } else {
    fill_parcel_buffer<double>(
        static_cast<double*>(nim->data), n_spatial, n_time, slots, n_parcels,
        slot_present, has_background, REAL(parcel_values), fill_value, cmin, cmax);
  }
  nim->cal_min = static_cast<RNifti::NiftiImage::scale_t>(cmin);
  nim->cal_max = static_cast<RNifti::NiftiImage::scale_t>(cmax);

  const auto paths = out.toFile(out_file, dt_code, filetype, compression);
  return paths.second;
}
