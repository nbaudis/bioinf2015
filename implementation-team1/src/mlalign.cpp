#include "mlalign.h"

#include "simdwrapper.h"

// Please leave include order as is.
// This is to keep my IDE from complaining about redundant imports.
#include <cmath>
#include <cstdio>
#include <memory>
#include <algorithm>
#include <vector>
#include <array>

#include "indexers.h"

/**
 * Memoization unnecessary, only used in the initialization step and not even multiple times.
 */
inline double gamma2(size_t n, double lambda, double mu);

/**
 * Memoization over n: Only value used is 1.
 */
inline double p_desc(size_t n, double beta, double lambda, double mu, double tau);

/**
 * Memoization over n is possible and more than viable given that there are only 2 values used (0 and 1).
 */
inline double p_desc_inv(size_t n, double beta, double lambda, double mu, double tau);

/**
 * Memoization over the 4*4 possible values of i and j sensible.
 */
inline double p_trans(mlalign_Nucleotide i, mlalign_Nucleotide j, double delta, double pi[mlalign_NUMBER_OF_BASES],
                      double tau);

/**
 * Only used in the initialization step. No memoization necessary.
 */
inline double zeta(size_t n, double beta, double lambda);

template<typename T>
size_t max_idx(const std::initializer_list<T> &values);

void mlalign_TKF91(mlalign_Nucleotide *as, size_t na, mlalign_Nucleotide *bs,
                   size_t nb, double lambda, double mu, double tau, double pi[mlalign_NUMBER_OF_BASES],
                   mlalign_Site *alignment, size_t *n, double *score) {

  /***********************/
  /* initialize matrices */
  /***********************/

  size_t rows = na + 1;
  size_t cols = nb + 1;
  auto idxM = diagonalizedIndexer(rows, cols, VECSIZE);
  auto idxA = identityIndexer(na);
  auto idxB = identityIndexer(nb);
  auto idxP1 = identityIndexer(mlalign_NUMBER_OF_BASES);
  auto idxP2 = rowMajorIndexer(mlalign_NUMBER_OF_BASES, mlalign_NUMBER_OF_BASES);
  auto idx2 = identityIndexer(2);


  // Allocate enough bytes for the matrices and properly align that.
  // We won't use vector because
  //  1) it initializes every element to 0
  //  2) has some extra payload we won't need (e.g. capacity and size for bounds checks)
  //  3) we don't resize
  const size_t doubles_per_matrix = to_next_multiple(idxM(rows - 1, cols - 1) + 1, VECSIZE);
  std::unique_ptr<double[]> buffer(new double[3 * doubles_per_matrix + VECSIZE]);
  auto M0 = reinterpret_cast<double*>(
          to_next_multiple(
                  reinterpret_cast<size_t>(buffer.get()),
                  VECSIZE * sizeof(double)));
  auto M1 = M0 + doubles_per_matrix;
  auto M2 = M1 + doubles_per_matrix;
  assert(M0 - buffer.get() < static_cast<std::ptrdiff_t>(VECSIZE));

  //
  // Memoization
  //
  const double beta = (1 - std::exp((lambda - mu) * tau)) / (mu - lambda * std::exp((lambda - mu) * tau));
  const double delta = 1 / (1 - pi[0] * pi[0] - pi[1] * pi[1] - pi[2] * pi[2] - pi[3] * pi[3]);
  const double lambda_div_mu = lambda / mu;
  const double lambda_mul_beta = lambda * beta;
  double p_desc_1 = p_desc(1, beta, lambda, mu, tau);
  double M0_factor[mlalign_NUMBER_OF_BASES];
  double M1_factor[mlalign_NUMBER_OF_BASES * mlalign_NUMBER_OF_BASES];
  double M2_factor[mlalign_NUMBER_OF_BASES];
  double p_desc_inv_n[2];

  for (size_t i = 0; i < 2; ++i) {
    p_desc_inv_n[idx2(i)] = p_desc_inv(i, beta, lambda, mu, tau);
  }

  for (size_t i = 0; i < mlalign_NUMBER_OF_BASES; ++i) {
    const auto a = static_cast<mlalign_Nucleotide>(i);
    M0_factor[idxP1(a)] = std::log(lambda_div_mu * pi[idxP1(a)] * p_desc_inv_n[idx2(0)]);
  }

  for (size_t i = 0; i < mlalign_NUMBER_OF_BASES; ++i) {
    for (size_t j = 0; j < mlalign_NUMBER_OF_BASES; ++j) {
      const auto a = static_cast<mlalign_Nucleotide>(i);
      const auto b = static_cast<mlalign_Nucleotide>(j);
      M1_factor[idxP2(a, b)] = std::log(lambda_div_mu * pi[idxP1(a)]
                                 * std::max(p_trans(a, b, delta, pi, tau) * p_desc_1, pi[idxP1(b)] * p_desc_inv_n[idx2(1)]));
    }
  }

  for (size_t i = 0; i < mlalign_NUMBER_OF_BASES; ++i) {
    const auto b = static_cast<mlalign_Nucleotide>(i);
    M2_factor[idxP1(b)] = std::log(pi[idxP1(b)] * lambda_mul_beta);
  }

  size_t ix = idxM(0, 0); // We'll reuse ix throughout this function for caching idx calculations
  M0[ix] = -std::numeric_limits<double>::max();
  M1[ix] = std::log(gamma2(0, lambda, mu) * zeta(1, beta, lambda));
  M2[ix] = -std::numeric_limits<double>::max();

  double accum = 0.0;
  for (size_t i = 1; i < rows; i++) {
    ix = idxM(i, 0);
    accum += std::log(pi[idxP1(as[idxA(i-1)])] * p_desc_inv_n[idx2(0)]);

    M0[ix] = std::log(gamma2(i, lambda, mu) * zeta(i, beta, lambda)) + accum;
    M1[ix] = -std::numeric_limits<double>::max();
    M2[ix] = -std::numeric_limits<double>::max();
  }

  accum = 0.0;
  for (size_t j = 1; j < cols; j++) {
    ix = idxM(0, j);
    accum += std::log(pi[idxP1(bs[idxB(j-1)])]);

    M0[ix] = -std::numeric_limits<double>::max();
    M1[ix] = -std::numeric_limits<double>::max();
    M2[ix] = std::log(gamma2(0, lambda, mu) * zeta(j + 1, beta, lambda)) + accum;
  }

  /****************/
  /* DP algorithm */
  /****************/

  size_t top[3] = { 0, 0, SIZE_MAX };
  size_t bottom[3] { 1, 0, SIZE_MAX };
  // This is tricky. We store the corrected offsets for diagonals offset -1 and -2 (7, 3 resp.) and the
  // uncorrected (e.g. not aligned) offset for the current diagonal (which is 7 + 2 = 9).
  // We inititialize with the values for diagonal 2, 1 and 0, we start at diagonal = 2 anyway.
  size_t diagonal_offsets[3] = { 2*VECSIZE + 1, 2*VECSIZE - 1, VECSIZE - 1 };
  const size_t number_of_diagonals = rows + cols - 1;
  for (size_t diagonal = 2; diagonal < number_of_diagonals; ++diagonal) {

    // Topmost and bottommost row index of the diagonal,
    // skipping the already initialized first row and column (diagonal - 1 for that reason).
    // It's helpful for intuition to visualize this (rows=7, cols=5):
    //
    //   3  7 11 15 19
    //   8 12 16 20 24
    //   5 17 21 25 32
    //  18 22 26 33 40
    //  23 27 34 41 44
    //  28 35 42 45 48
    //  36 43 46 49 52
    //
    // diagonal             | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
    // bottom_without_first | 1 | 2 | 3 | 4 | 5 | 6 | 6 | 6 |
    // top_without_first    | 1 | 1 | 1 | 1 | 2 | 3 | 4 | 5 |
    //
    // We will then iterate over the diagonal through the indices
    // [top_without_first, bottom_without_first] inclusive. Column
    // indices are easily computed from the diagonal and row index.
    // We start with diagonal = 2, e.g. the third diagonal, because
    // it's the first containing uninitialized values.
    //
    const size_t bottom_without_first = std::min(diagonal - 1, rows - 1);
    const size_t top_without_first = diagonal - std::min(diagonal - 1, cols - 1);

    std::rotate(top, top + 2, top + 3); // This pushes all elements on further back
    std::rotate(bottom, bottom + 2, bottom + 3);
    bottom[0] = std::min(diagonal, rows - 1);
    top[0] = diagonal - std::min(diagonal, cols - 1);
    const size_t length_of_diagonal = bottom[0] + 1 - top[0];

    // We have to correct the diagonal offset for the alignment rules layed out in diagonalizedIndexer.
    //    1. If the diagonal starts in the first row (top[0]=0), then the offset must have alignment -1.
    //    2. If the diagonal doesn't start in the first row, then the offset must have alignment 0.
    diagonal_offsets[0] += top[0] == 0
                  ? VECSIZE - 1 - diagonal_offsets[0] % VECSIZE
                  : (VECSIZE - diagonal_offsets[0] % VECSIZE) % VECSIZE;

    assert((top[0] == 0 && diagonal_offsets[0] % VECSIZE == VECSIZE - 1)
           || (top[0] != 0 && diagonal_offsets[0] % VECSIZE == 0));

    // INDEX EXPLANATIONS
    // row_offset   offset of the row where the current SIMD vector begins
    // i            index into the SIMD vector
    // col          index of the column where the currently processed vector element is located in the DP matrix
    // row          index of the row where the currently processed vector element is located in the DP matrix
    
    // Step in vector size over diagonal entries
    for (size_t row_offset = top_without_first; row_offset <= bottom_without_first; row_offset += VECSIZE) {


      // Copy the factors values from the lookuptable into the vector
      // Apparently it is not better to copy into temp memory and the use loaddvec, but it is faster to just use the [] operator.
      dvec factors[3];
      for (size_t i = 0; i < VECSIZE; ++i) {
        const size_t row = row_offset + i;
        const size_t col = diagonal - row;

        if (row > bottom_without_first) {
          // We would load from somewhere out of the matrix, so we just leave
          // the
          // other vector elements untouched
          break;
        }

        const mlalign_Nucleotide a = as[idxA(row - 1)];
        const mlalign_Nucleotide b = bs[idxB(col - 1)];

        factors[0][i] = M0_factor[idxP1(a)];
        factors[1][i] = M1_factor[idxP2(a, b)];
        factors[2][i] = M2_factor[idxP1(b)];
      }

      // Compute the indexes of the top, left and top left cells and load them.

      const size_t ix_top = diagonal_offsets[1] + row_offset - 1 - top[1];
      assert(ix_top == idxM(row_offset - 1, diagonal - row_offset));
      const dvec max_top[3] = {
              loaddvec(M0 + ix_top),
              loaddvec(M1 + ix_top),
              loaddvec(M2 + ix_top),
      };

      const size_t ix_topleft = diagonal_offsets[2] + row_offset - 1 - top[2];
      assert(ix_topleft == idxM(row_offset - 1, diagonal - row_offset - 1));
      const dvec max_topleft[3] = {
              loaddvec(M0 + ix_topleft),
              loaddvec(M1 + ix_topleft),
              loaddvec(M2 + ix_topleft),
      };

      // Against all symmetry M0 isn't part of the formula.
      const size_t ix_left = ix_top + 1;
      assert(ix_left == idxM(row_offset, diagonal - row_offset - 1));
      const dvec max_left[2] = {
              loaddvec(M1 + ix_left),
              loaddvec(M2 + ix_left),
      };

      // Maxima for each cells on top/topleft/left of the calculated cell
      const dvec mt = maxdvec(max_top[0], maxdvec(max_top[1], max_top[2]));
      const dvec mtl = maxdvec(max_topleft[0], maxdvec(max_topleft[1], max_topleft[2]));
      const dvec ml = maxdvec(max_left[0], max_left[1]);

      // Vectorised add
      const dvec m0v = factors[0] + mt;
      const dvec m1v = factors[1] + mtl;
      const dvec m2v = factors[2] + ml;

      // Write back to matrix
      const size_t occupied_vector_elements = std::min(bottom_without_first + 1 - row_offset, VECSIZE);
      ix = diagonal_offsets[0] + row_offset - top[0];
      assert(ix % VECSIZE == 0);
      assert(ix == idxM(row_offset, diagonal - row_offset));
      storedvec(M0 + ix, m0v, occupied_vector_elements);
      storedvec(M1 + ix, m1v, occupied_vector_elements);
      storedvec(M2 + ix, m2v, occupied_vector_elements);
    }

    std::rotate(diagonal_offsets, diagonal_offsets + 2, diagonal_offsets + 3);
    diagonal_offsets[0] = diagonal_offsets[1] + length_of_diagonal;
  }

  /****************/
  /* Backtracking */
  /****************/

  ix = idxM(na, nb);
  *n = 0;
  *score = std::exp(std::max({M0[ix], M1[ix], M2[ix]}));
  size_t cur = max_idx({M0[ix], M1[ix], M2[ix]});

  for (size_t i = na, j = nb; i > 0 || j > 0;) {
    mlalign_Site site = {mlalign_Gap, mlalign_Gap};
    switch (cur) {
      case 0: // insert
        assert(i > 0);
        site.a = as[--i];
        site.b = mlalign_Gap;

        ix = idxM(i, j);
        cur = max_idx({M0[ix], M1[ix], M2[ix]});
        break;
      case 1: // substitution
        assert(i > 0);
        assert(j > 0);
        site.a = as[--i];
        site.b = bs[--j];

        ix = idxM(i, j);
        cur = max_idx({M0[ix], M1[ix], M2[ix]});
        break;
      case 2: // deletion
        assert(j > 0);
        site.a = mlalign_Gap;
        site.b = bs[--j];

        ix = idxM(i, j);
        cur = max_idx({M1[ix], M2[ix]}) + 1;
        break;
      default:
        assert(0);
    }

    alignment[(*n)++] = site; // Doesn't get any scarier
  }

  // The alignment is reversed, so we have to undo that.
  for (size_t i = 0; i < (*n) / 2; ++i) {
    std::swap(alignment[i], alignment[*n - i - 1]);
  }
}

double p_trans(mlalign_Nucleotide i, mlalign_Nucleotide j, double delta,
               double pi[mlalign_NUMBER_OF_BASES], double tau) {
  auto idxP = identityIndexer(mlalign_NUMBER_OF_BASES);

  double e_delta_t = std::exp(-delta * tau);

  if (i == j)
    return e_delta_t + pi[idxP(j)] * (1 - e_delta_t);
  else
    return pi[idxP(j)] * (1 - e_delta_t);
}

double p_desc(size_t n, double beta, double lambda, double mu, double tau) {
  return std::exp(-mu * tau) * (1 - lambda * beta) * std::pow(lambda * beta, n - 1);
}

double p_desc_inv(size_t n, double beta, double lambda, double mu, double tau) {
  if (n == 0)
    return mu * beta;
  else
    return (1 - std::exp(-mu * tau) - mu * beta) * (1 - lambda * beta)
           * std::pow(lambda * beta, n - 1);
}

double gamma2(size_t n, double lambda, double mu) {
  return (1 - lambda / mu) * std::pow(lambda / mu, n);
}

double zeta(size_t n, double beta, double lambda) {
  return (1 - lambda * beta) * std::pow(lambda * beta, n - 1);
}

template<typename T>
size_t max_idx(const std::initializer_list<T> &values) {
  size_t i = 0, mi = 0;
  const T *max = nullptr;
  for (const auto &v : values) {
    if (max == nullptr || v > *max) {
      max = &v;
      mi = i;
    }
    ++i;
  }
  return mi;
}

