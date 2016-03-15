#ifndef SPECMATH_H_
#define SPECMATH_H_

#include "config.h"
#include "ServiceFunc.h"

using namespace std;

	double RMS(const raw_t* source, uint32_t start, uint32_t finish);
	double Entropy(const double* source, uint32_t start, uint32_t finish, uint8_t binsCount, double minRaw, double maxRaw);
	double* calcMFCC(const double* source, uint32_t start, uint32_t finish, uint8_t mfccSize, uint32_t frequency, uint32_t freqMin, uint32_t freqMax);
	void fourierTransformFastRecursion(valarray<complex<double>>& data);
	double* fourierTransformFast(const double* source, uint32_t length, bool useWindow);
	double** getMelFilters(uint8_t mfccSize, uint32_t filterLength, uint32_t frequency, uint32_t freqMin, uint32_t freqMax);
	double* calcPower(const double* fourierRaw, uint32_t fourierLength, double** melFilters, uint8_t mfccSize);
	double* dctTransform(const double* data, uint32_t length);
	double calcDistance(double* seq1, uint32_t seq1size, double* seq2, uint32_t seq2size);
    double findDistance(uint32_t seq1size, uint32_t seq2size, double** diffM);
    double calcDistanceFor2Matrices(const vector<double*> matrix1, const vector<double*> matrix2, const uint16_t matricesWidth);

#endif /* SPECMATH_H_ */
