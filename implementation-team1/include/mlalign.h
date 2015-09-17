#ifndef MLALIGN_MLALIGN_H
#define MLALIGN_MLALIGN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif


enum mlalign_Nucleotide {
  mlalign_A = 0,
  mlalign_C = 1,
  mlalign_G = 2,
  mlalign_T = 3,
  mlalign_Gap = 4,
  mlalign_NUMBER_OF_BASES = 4,
};

struct mlalign_Site {
  mlalign_Nucleotide a : 4;
  mlalign_Nucleotide b : 4;
};

void mlalign_TKF91(mlalign_Nucleotide *as, size_t na,
                   mlalign_Nucleotide *bs, size_t nb,
                   double lambda,
                   double mu,
                   double tau,
                   double pi[mlalign_NUMBER_OF_BASES],
                   mlalign_Site *alignment,
                   size_t *n, // alignment must be long enough to fit in the whole alignment, e.g. na + nb
                   double *score);

#ifdef __cplusplus

}; // extern "C"

#endif

#endif // MLALIGN_MLALIGN_H
