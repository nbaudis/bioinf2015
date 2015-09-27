/*
 * TKF91Sequential.cpp
 *
 *  Created on: May 7, 2015
 *      Authors: Sarah Lutteropp, Pierre Barbera, Marin Vlastelica Pogančić
 */

#include "TKF91Sequential.h"
#include <crlibm.h>

using namespace std;

const double MINUS_INFINITY = -std::numeric_limits<double>::infinity();

TKF91Sequential::TKF91Sequential(std::string pseq1, std::string pseq2,
		const double &birthRate, const double &deathRate, const double &time,
		const double equilibriumProb[4], const bool outputMatrices, const bool outputAlignment) :
		seq1(pseq1), seq2(pseq2), mu(deathRate), logLambda(log_ru(birthRate)), logMu(
				log_ru(deathRate)), t(time) {
	crlibm_init();
	assert(birthRate >= 0.0);
	assert(deathRate >= 0.0);
	assert(time > 0);
	assert(birthRate < deathRate);

	TKF91Sequential::numRows = seq1.size() + 1;
	TKF91Sequential::numCols = seq2.size() + 1;

	TKF91Sequential::alignedSeq1 = seq1;
	TKF91Sequential::alignedSeq2 = seq2;

	TKF91Sequential::seq1AsInt = seqToIntVec(seq1);
	TKF91Sequential::seq2AsInt = seqToIntVec(seq2);

	// equilibrium probabilities
	//const char* nucleotides = "ACGT";
	for (int i = 0; i < 4; ++i) {
		TKF91Sequential::pi.push_back(equilibriumProb[i]);
		TKF91Sequential::logPi.push_back(log_ru(equilibriumProb[i]));
	}
	// log beta
	double exponential = exp((birthRate - deathRate) * time);
	TKF91Sequential::logBeta = log_ru(1 - exponential)- log_ru(deathRate - birthRate * exponential);
	//TKF91Sequential::logBeta = log_ru((1- exponential)/(deathRate - birthRate * exponential));
	TKF91Sequential::logOneMinusLambdaByMu = log_ru(1 - birthRate / deathRate);
	TKF91Sequential::logOneMinusLambdaTimesBeta = log_ru(
			1 - birthRate * exp(logBeta));
	TKF91Sequential::longLog = log_ru(
			1 - exp(-deathRate * t) - deathRate * exp(logBeta));

	// compute substitution rate parameter for the F81 model
	TKF91Sequential::delta = 1
			/ (1 - equilibriumProb[0] * equilibriumProb[0]
					- equilibriumProb[1] * equilibriumProb[1]
					- equilibriumProb[2] * equilibriumProb[2]
					- equilibriumProb[3] * equilibriumProb[3]);

	initMatrices();
	computeMatrices();
	retrieveOptimalAlignment();
	
	if (outputMatrices)
		printMatrices();
	if (outputAlignment)
		printAlignment();
	
}

TKF91Sequential::~TKF91Sequential() {
	free(m);
}

size_t TKF91Sequential::CO(size_t row, size_t col) {
	return row * numCols + col;
}

void TKF91Sequential::initMatrices() {
	size_t size = numRows * numCols * sizeof(MatrixEntry);
	m = (MatrixEntry*) malloc(size);
}

void TKF91Sequential::computeMatrices() {
	std::vector<double> bestLogProb(16);
	double tempLog = log_ru(1 - exp(-delta * t));
	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			double way1;
			if (i == j) {
				way1 = log_ru(exp(-delta * t) + pi[j] * (1 - exp(-delta * t)));
			} else {
				way1 = logPi[j] + tempLog;
			}
			way1 -= mu * t;
			double way2 = logPi[j] + longLog;
			bestLogProb[i * 4 + j] = std::max(way1, way2);
		}
	}

	double preStuff_m0 = logLambda + logBeta;
	double preStuff_m1 = logLambda - logMu + logOneMinusLambdaTimesBeta;
	double preStuff_m2 = logLambda + logBeta;

	m[CO(0, 0)].m_0 = MINUS_INFINITY;
	m[CO(0, 0)].m_1 = logOneMinusLambdaByMu + logOneMinusLambdaTimesBeta;
	m[CO(0, 0)].m_2 = MINUS_INFINITY;
	m[CO(0, 1)].m_0 = MINUS_INFINITY;
	m[CO(0, 1)].m_2 = logOneMinusLambdaByMu + logOneMinusLambdaTimesBeta + logLambda
			+ logBeta + logPi[seq2AsInt[0]];
	m[CO(0, 1)].m_1 = MINUS_INFINITY;

	double oldVal_m2 = m[CO(0, 1)].m_2;
	for (size_t i = 2; i < numCols; ++i) {
		m[CO(0, i)].m_0 = MINUS_INFINITY;
		m[CO(0, i)].m_1 = MINUS_INFINITY;
		m[CO(0, i)].m_2 = oldVal_m2 + preStuff_m2 + logPi[seq2AsInt[i - 1]];
		oldVal_m2 = m[CO(0, i)].m_2;
	}
	

	m[CO(1, 0)].m_0 = logOneMinusLambdaByMu + logLambda + logOneMinusLambdaTimesBeta
		+ logBeta + logPi[seq1AsInt[0]];
	m[CO(1, 0)].m_1 = MINUS_INFINITY;
	m[CO(1, 0)].m_2 = MINUS_INFINITY;

	double oldVal_m0 = m[CO(1, 0)].m_0;
	for (size_t i = 2; i < numRows; ++i) {
		m[CO(i, 0)].m_0 = oldVal_m0 + 2 * logLambda + 2 * logBeta
				+ logPi[seq1AsInt[i - 1]];
		oldVal_m0 = m[CO(i, 0)].m_0;
		m[CO(i, 0)].m_1 = MINUS_INFINITY;
		m[CO(i, 0)].m_2 = MINUS_INFINITY;
	}

	for (size_t i = 1; i < numRows; ++i) {
		size_t nuc1 = seq1AsInt[i - 1];
		double curLogPi = logPi[nuc1];
		nuc1 *= 4;
		for (size_t j = 1; j < numCols; ++j) {
			size_t coord = CO(i, j);
			size_t up = CO(i - 1, j);
			size_t left = CO(i, j - 1);
			size_t diag = CO(i - 1, j - 1);
			m[coord].m_0 = preStuff_m0 + curLogPi;
			m[coord].m_0 += std::max(m[up].m_0, std::max(m[up].m_1, m[up].m_2));
			m[coord].m_1 = preStuff_m1 + curLogPi + bestLogProb[nuc1 + seq2AsInt[j - 1]];
			m[coord].m_1 += std::max(m[diag].m_0, std::max(m[diag].m_1, m[diag].m_2));
			m[coord].m_2 = preStuff_m2 + logPi[seq2AsInt[j - 1]];
			m[coord].m_2 += std::max(m[left].m_1, m[left].m_2);
		}
	}	
}

void TKF91Sequential::retrieveOptimalAlignment() {
	size_t coord = CO(numRows - 1, numCols - 1);
	short maxPos = m[coord].getMaxIdx012();
	backtrackAlignment(maxPos);
}

void TKF91Sequential::backtrackAlignment(short matrixNum) {
	size_t i = numRows - 1;
	size_t j = numCols - 1;
	size_t coord = CO(i, j);

	switch (matrixNum) {
	case 0:
		TKF91Sequential::alignmentScore = m[coord].m_0;
		break;
	case 1:
		TKF91Sequential::alignmentScore = m[coord].m_1;
		break;
	case 2:
		TKF91Sequential::alignmentScore = m[coord].m_2;
		break;
	default:
		assert(0);
	}

	while (coord != CO(0, 0)) {
		size_t up = CO(i - 1, j);
		size_t diag = CO(i - 1, j - 1);
		size_t left = CO(i, j - 1);

		short maxIndex = 0;

		switch (matrixNum) {
		case 0: // up
			maxIndex = m[up].getMaxIdx012();
			alignedSeq2.insert(j, 1, '-');
			i--;
			break;
		case 1: // diag
			maxIndex = m[diag].getMaxIdx012();
			i--;
			j--;
			break;
		case 2:
			// left
			maxIndex = m[left].getMaxIdx12();
			alignedSeq1.insert(i, 1, '-');
			j--;
			break;
		default:
			assert(0);
		}

		matrixNum = maxIndex;

		coord = CO(i, j);
	}
}

void TKF91Sequential::printAlignment() {
	std::cout << alignedSeq1 << "\n";
	std::cout << alignedSeq2 << "\n";
	cout.precision(100);
	std::cout << fixed << alignmentScore << "\n";
}

void TKF91Sequential::printMatrices() {
	cout << "Freqs: " << pi[0] << " " << pi[1] << " " << pi[2] << " "
			<< pi[3] << "\n";
	cout << "Beta: " << exp(logBeta) << "\n";
	cout << "\n";

	cout << "log m0:\n";
	for (size_t i = 0; i < numRows; ++i) {
		for (size_t j = 0; j < numCols; ++j) {
			cout << (m[CO(i, j)].m_0) << " ";
		}
		cout << "\n";
	}
	cout << "\n";

	cout << "log m1:\n";
	for (size_t i = 0; i < numRows; ++i) {
		for (size_t j = 0; j < numCols; ++j) {
			cout << (m[CO(i, j)].m_1) << " ";
		}
		cout << "\n";
	}
	cout << "\n";

	cout << "log m2:\n";
	for (size_t i = 0; i < numRows; ++i) {
		for (size_t j = 0; j < numCols; ++j) {
			cout << (m[CO(i, j)].m_2) << " ";
		}
		cout << "\n";
	}
	cout << "\n";
}


