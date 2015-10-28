#ifndef MLALIGN_ARRAY_VIEW_H
#define MLALIGN_ARRAY_VIEW_H

#include <cstddef>
#include <cassert>
#include <algorithm>
#include <tuple>
#include <functional>

auto identityIndexer = [] (size_t size) {
  return [size](size_t i) {
    assert(i < size);
    return i;
  };
};

auto rowMajorIndexer = [] (size_t rows, size_t cols) {
  return [rows, cols](size_t row, size_t col) {
    assert(row < rows);
    assert(col < cols);
    return row * cols + col;
  };
};

size_t aligned_stair(size_t n, size_t align) {
  size_t repeats = n / align;
  size_t rest = n % align;
  size_t regularPart = repeats*(repeats+1)/2*align*align;
  size_t restPart = rest*(repeats + 1)*align;
  return regularPart + restPart;
}

size_t to_next_multiple(size_t n, size_t align) {
  const size_t m = ((n + align - 1) / align) * align;
  assert(m % align == 0);
  assert(m - n < align);
  return m;
}

auto diagonalizedIndexer = [] (size_t rows, size_t cols, size_t align) {
  return [=](size_t row, size_t col) {
    assert(row < rows);
    assert(col < cols);
    size_t smaller_rank, larger_rank;
    std::tie(smaller_rank, larger_rank) = std::minmax(rows, cols);

    const size_t diagonal = row + col; // The indexed diagonal

    //
    // We will compute the indices based on the Cantor pairing function, adapted
    // to work on a matrix rather than on the typical upper left triangle of the N^2 space.
    // We also insert padding elements at the end of each diagonal to achieve a certain alignment:
    //    1. Elements of the first row (row=0) will have alignment -1 (so that row=1 is aligned).
    //       Remember that the first row contains initialized values that we don't want to change.
    //    2. The beginning of all other diagonals which don't start in the first row should be aligned.
    //
    // For a 7x5 matrix and an alignment of 4, the indices should be like this:
    //
    //   3  7 11 15 19 <- until this diagonal incl, we have the opening part
    //   8 12 16 20 24
    //   5 17 21 25 32 <- until this diagonal incl, we have the middle part
    //  18 22 26 33 40
    //  23 27 34 41 44
    //  28 35 42 45 48
    //  36 43 46 49 52 <- and until this diagonal incl, we have the closing part
    //
    // Note that I identified the parts which will have different indexing formulas by
    // the opening, middle and closing part.
    //
    // Computing the offset of the ith diagonal will need to compute the offsets of the
    // respective parts and then adding them together. When we have the offset of a
    // particular diagonal, we find the offset into that diagonal by
    //
    //     std::min(row, cols - 1 - col)
    //
    // (Just check for yourself)
    //

    // We first compute the size of the closing, middle and opening part
    size_t d = diagonal;

    size_t closing = std::max(d, larger_rank) - larger_rank; // essentially std::min(diagonal-larger_rank, 0)
    d -= closing;
    size_t middle = std::max(d, smaller_rank) - smaller_rank;
    d -= middle;
    size_t opening = d;

    // Now for computing the actual offsets.
    // The opening and closing parts are modeled with the aligned_staircase function,
    // which aligns the 'staircase' series (sum of [1, ..., n]) to multiples of align.
    size_t offset = 0;

    offset += aligned_stair(opening, align);

    // The following line ensures that we align the second row, e.g. row 0 is on alignment -1:
    offset += diagonal < cols ? align - 1 : smaller_rank % align == 1 ? 0 : align;

    offset += middle * to_next_multiple(smaller_rank, align);

    // It's smaller_rank - 1 because we don't count the last (full) smaller_rank-sized
    // diagonal to the closing part, so we start with the first non-full diagonal of size
    // smaller_rank - 1.
    offset += aligned_stair(smaller_rank - 1, align);
    offset -= aligned_stair(smaller_rank - 1 - closing, align);

    // What remains is indexing into the diagonal:
    offset += std::min(row, cols - 1 - col);

    return offset;
  };
};


#endif //MLALIGN_ARRAY_VIEW_H
