#ifndef MLALIGN_INIT_SIMDWRAPPER_H
#define MLALIGN_INIT_SIMDWRAPPER_H

#include <assert.h>


#ifdef SIMD_ARCH_X86_SSE3

#include <xmmintrin.h>
const size_t VECSIZE = 2;
typedef __m128d dvec;

#elif defined(SIMD_ARCH_X86_AVX)

#include <immintrin.h>
const size_t VECSIZE = 4;
typedef __m256d dvec;

#else

#error("Specify a simd architecture (SIMD_ARCH_X86_SSE3 or SIMD_ARCH_X86_AVX)")

#endif


static inline dvec maxdvec(dvec a, dvec b)
{
  #ifdef SIMD_ARCH_X86_SSE3
  return _mm_max_pd(a, b);
  #elif defined(SIMD_ARCH_X86_AVX)
  return _mm256_max_pd(a, b);
  #endif
}


static inline dvec loaddvec(double* src)
{
#ifdef SIMD_ARCH_X86_SSE3
  return _mm_loadu_pd(src);
#elif defined(SIMD_ARCH_X86_AVX)
  return _mm256_loadu_pd(src);
#endif 
}


#ifdef SIMD_ARCH_X86_AVX

const static __m256i STORE_MASKS[] = {
        _mm256_set_epi64x( 0,  0,  0,  0), // This should not be accessed anyway, just for padding
        _mm256_set_epi64x( 0,  0,  0, ~0),
        _mm256_set_epi64x( 0,  0, ~0, ~0),
        _mm256_set_epi64x( 0, ~0, ~0, ~0),
};

#endif


static inline void storedvec(double* dest, dvec src, size_t n)
{
  assert(n > 0);
  assert(n <= VECSIZE); // This is forbidden as part of the contract.

#ifdef SIMD_ARCH_X86_SSE3

  if (n == 2) {
      _mm_storeu_pd(dest, src); // Most likely case first
  } else if (n == 1) {
      _mm_storel_pd(dest, src);
  }

#elif defined(SIMD_ARCH_X86_AVX)

  if (n == VECSIZE) {
    _mm256_storeu_pd(dest, src);
  } else {
    _mm256_maskstore_pd(dest, STORE_MASKS[n], src);
  }

#endif

}

#endif /* end of include guard: MLALIGN_INIT_SIMDWRAPPER_H */
