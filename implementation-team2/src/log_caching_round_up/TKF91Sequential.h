/*
 * TKF91Sequential.h
 *
 *  Created on: May 7, 2015
 *      Authors: Sarah Lutteropp, Pierre Barbera, Marin Vlastelica Pogančić
 */

#ifndef TKF91SEQUENTIAL_H_
#define TKF91SEQUENTIAL_H_

#include <string>
#include <string.h>
#include <math.h>
#include <cassert>
#include <map>
#include <iostream>
#include <vector>
#include <limits>
#include <utility> // for std::pair

typedef struct MatrixEntry_data {
	double m_0;
	double m_1;
	double m_2;
	
	short getMaxIdx012() {
		if (m_0 >= m_1 && m_0 >= m_2) {
			return 0;
		} else if (m_1 >= m_0 && m_1 >= m_2) {
			return 1;
		} else {
			return 2;
		}
	}
	short getMaxIdx12() {
		if (m_1 >= m_2) {
			return 1;
		} else {
			return 2;
		}
	}
} MatrixEntry;

class TKF91Sequential {
public:
	TKF91Sequential(std::string seq1, std::string seq2, const double &birthRate,
			const double &deathRate, const double &time,
			const double equilibriumProb[4], const bool outputMatrices = false, 
			const bool outputAlignment = false);
	~TKF91Sequential();
private:
	std::string seq1;
	std::string seq2;
	std::vector<size_t> seq1AsInt;
	std::vector<size_t> seq2AsInt;
	unsigned long numRows;
	unsigned long numCols;
	double mu; // death rate
	double logLambda; // log birth rate
	double logMu; // log death rate
	std::vector<double> pi; // equilibrium probabilities
	std::vector<double> logPi; // log equilibrium probabilities
	double delta; // substitution rate parameter for the F81 model
	double t; // time
	double logBeta; // log beta from paper

	double logOneMinusLambdaByMu;
	double logOneMinusLambdaTimesBeta;
	double longLog; // log(1-exp(-deathRate*t)-deathRate*beta)

	MatrixEntry *m; // matrices
	std::string alignedSeq1;
	std::string alignedSeq2;
	double alignmentScore;

	void initMatrices();
	void computeMatrices();
	void retrieveOptimalAlignment();
	void backtrackAlignment(short matrixNum);
	void printAlignment();
	void printMatrices();

	size_t CO(size_t row, size_t col);

	static inline size_t nucToIntFunc(char N) {
		size_t ret = 4;
		switch (N) {
			case 'A':
				ret = 0;
				break;
			case 'C':
				ret = 1;
				break;
			case 'G':
				ret = 2;
				break;
			case 'T':
				ret = 3;
				break;
			default:
				std::cout << "nuc to int failed: no such nucleotide\n";		
		}
		assert(ret != 4);
		return ret;
	}

	static std::vector<size_t> seqToIntVec(const std::string seq) {
		std::vector<size_t> ret(seq.length());
		for (size_t i = 0; i < seq.length(); ++i)
		{
			ret[i] = nucToIntFunc(seq.at(i));
		}
		return ret;
	}

	
	static void* malloc_aligned(size_t size, size_t align) 
	{
	  void *ptr = (void *)NULL;  
	  int res;
	  
	  ptr = malloc(size);
	  
	  if(ptr == (void*)NULL) 
	   assert(0);

	  res = posix_memalign( &ptr, align, size );

	  if(res != 0) 
	    assert(0);
	   
	  return ptr;
	}

};

#endif /* TKF91SEQUENTIAL_H_ */
