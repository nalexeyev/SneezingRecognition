#include "SpecMath.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 

using namespace std;

/************************************************************
 *  Audio Splitting and MFCC Calculation
 ************************************************************/

double RMS(const raw_t* source, uint32_t start, uint32_t finish) {
	double value = 0;

	for (uint32_t i = start; i <= finish; i++) {
		value += source[i] * source[i];
	}
	value /= (finish - start + 1);

	return sqrt(value);
}

double Entropy(const double* source, uint32_t start, uint32_t finish,
		uint8_t binsCount, double minRaw, double maxRaw) {
	double entropy = 0;

	double binSize = abs(maxRaw - minRaw) / static_cast<double>(binsCount);
	if (fabs(binSize) < numeric_limits<double>::epsilon()) {
		return 0;
	}

	double* p = new double[binsCount];
	for (uint8_t i = 0; i < binsCount; i++) {
		p[i] = 0.;
	}

	// Calculate probabilities
	uint8_t index;
	for (uint32_t i = start; i <= finish; i++) {
		double value = source[i];
		index = floor((value - minRaw) / binSize);

		if (index >= binsCount) {
			index = binsCount - 1;
		}

		p[index] += 1.;
	}

	// Normalize probabilities
	uint32_t size = finish - start + 1;
	for (uint8_t i = 0; i < binsCount; i++) {
		p[i] /= size;
	}

	// Calculate entropy
	for (uint8_t i = 0; i < binsCount; i++) {
		if (p[i] > numeric_limits<double>::epsilon()) {
			entropy += p[i] * log2(p[i]);
		}
	}

	entropy = -entropy;

	if (entropy<0){
		fprintf(stderr, "NEGATIVE ENTROPY\n");
		//getchar();
	}

	delete [] p;

	return entropy;
}

double* calcMFCC(const double* source, uint32_t start, uint32_t finish, uint8_t mfccSize, uint32_t frequency, uint32_t freqMin, uint32_t freqMax) {
	uint32_t sampleLength = finish - start + 1;
	uint32_t p2length = pow(2, floor(log2(sampleLength)));

	// Calc
	double* fourierRaw = fourierTransformFast(source + start, p2length, true);

	double** melFilters = getMelFilters(mfccSize, p2length, frequency, freqMin, freqMax);
	double* logPower = calcPower(fourierRaw, p2length, melFilters, mfccSize);
	double* dctRaw = dctTransform(logPower, mfccSize);

	// Clean up
	delete [] logPower;
	delete [] fourierRaw;

	for (unsigned short m = 0; m < mfccSize; m++) {
		delete [] melFilters[m];
	}
	delete [] melFilters;

	return dctRaw;
}

double* fourierTransformFast(const double* source, uint32_t length, bool useWindow) {

	// Extend source length to the power of 2
	uint32_t p2length = length;

	bool powerOfTwo = (length > 0) && !(length & (length - 1));
	assert("FFT input data size must have 2^n size" && powerOfTwo);

	// Move to complex calculations
	double* fourierRaw = new double[length];
	valarray<complex<double>> fourierRawTmp(p2length);

	for (uint32_t i = 0; i < p2length; i++) {

		// Even element is the real part of complex number
		if (i < length) {
			fourierRawTmp[i] = complex<double>(source[i], 0.);

			if (useWindow) {
				fourierRawTmp[i] *= (0.54 - 0.46 * cos(2 * M_PI * i / (length - 1)));
			}

		} else {
			fourierRawTmp[i] = complex<double>(0, 0);
		}
	}

	// Perform recursive calculations
	fourierTransformFastRecursion(fourierRawTmp);

	// As for magnitude, let's use Euclid's distance for its calculation
	for (uint32_t i = 0; i < length; i++) {
		fourierRaw[i] = sqrt(norm(fourierRawTmp[i]));
	}

	return fourierRaw;
}

void fourierTransformFastRecursion(valarray<complex<double>>& data) {

	// Exit from recursion
	const size_t n = data.size();
	if (n <= 1) {
		return;
	}

	// Divide into Even/Odd
	valarray<complex<double>> even = data[std::slice(0, n/2, 2)];
	valarray<complex<double>> odd = data[std::slice(1, n/2, 2)];

	// Compute recursion
	fourierTransformFastRecursion(even);
	fourierTransformFastRecursion(odd);

	// Combine
	for (size_t i = 0; i < n / 2; i++) {
		complex<double> t = polar(1.0, -2 * M_PI * i / n) * odd[i];
		data[i]       = even[i] + t;
		data[i + n/2] = even[i] - t;
	}
}

double convertToMel(double f) { return 1125. * log(1. + f/700.); }

double convertFromMel(double m) { return 700. * (exp(m/1125.) - 1); }


/**
 * Create triangular filters spaced on mel scale
 */
double** getMelFilters(uint8_t mfccSize, uint32_t filterLength, uint32_t frequency,	uint32_t freqMin, uint32_t freqMax) {

	// Create points for filter banks
	double* fb = new double[mfccSize + 2];
	fb[0] = convertToMel(freqMin);
	fb[mfccSize + 1] = convertToMel(freqMax);

	// Create mel bin
	for (uint8_t m = 1; m < mfccSize + 1; m++) {
		fb[m] = fb[0] + m * (fb[mfccSize + 1] - fb[0]) / (mfccSize + 1);
	}

	//frequency = 0.5 * frequency;
	for (uint8_t m = 0; m < mfccSize + 2; m++) {

		// Convert them from mel to frequency
		fb[m] = convertFromMel(fb[m]);

		// Map those frequencies to the nearest FT bin
		fb[m] = floor((filterLength + 1) * fb[m] / (double) frequency);

		assert("FT bin too small" &&
				!(m > 0 && (fb[m] - fb[m-1]) < numeric_limits<double>::epsilon()));
	}

	// Calc filter banks
	double** filterBanks = new double*[mfccSize];
	for (uint8_t m = 0; m < mfccSize; m++) {
		filterBanks[m] =  new double[filterLength];
	}

	for (uint8_t m = 1; m < mfccSize + 1; m++) {
		for (uint32_t k = 0; k < filterLength; k++) {

			if (fb[m - 1] <= k && k <= fb[m]) {
				filterBanks[m - 1][k] = (k - fb[m - 1]) / (fb[m] - fb[m - 1]);

			} else if (fb[m] < k && k <= fb[m + 1]) {
				filterBanks[m - 1][k] = (fb[m + 1] - k) / (fb[m + 1] - fb[m]);

			} else {
				filterBanks[m - 1][k] = 0;
			}
		}
	}

	delete [] fb;

	return filterBanks;
}

/**
 * Apply mel filters to spectrum's magnitudes, take the logs of the powers
 */
double* calcPower(const double* fourierRaw, uint32_t fourierLength, double** melFilters, uint8_t mfccCount) {

	double* logPower = new double[mfccCount];

	for (uint8_t m = 0; m < mfccCount; m++) {
		logPower[m] = 0.;

		for (uint32_t k = 0; k < fourierLength; k++) {
			logPower[m] += melFilters[m][k] * pow(fourierRaw[k], 2);
		}

		assert("Spectrum power is less than zero" &&
				!(logPower[m] < numeric_limits<double>::epsilon()));

		logPower[m] = log(logPower[m]);
	}

	return logPower;
}

/**
 * Take the discrete cosine transform of the list of mel log powers
 */
double* dctTransform(const double* data, uint32_t length) {

	double* dctTransform = new double[length];

	for (uint8_t n = 0; n < length; n++) {
		dctTransform[n] = 0;

		for (uint8_t m = 0; m < length; m++) {
			dctTransform[n] += data[m] * cos(M_PI * n * (m + 1./2.) / length);
		}
	}

	return dctTransform;
}

/************************************************************
 *  DTW Implementation
 ************************************************************/

double calcDistance(double* seq1, uint32_t seq1size, double* seq2, uint32_t seq2size) {

	// Create diff matrix
	double** diffM = new double*[seq1size];
	for (uint32_t i = 0; i < seq1size; i++) {
		diffM[i] = new  double[seq2size];
	}

	for (uint32_t i = 0; i < seq1size; i++) {
		for (uint32_t j = 0; j < seq2size; j++) {
			diffM[i][j] = fabs(seq1[i] - seq2[j]);
		}
	}

	// Compute distance
	double distance = findDistance(seq1size, seq2size, diffM);

	// Clean up
	for (uint32_t i = 0; i < seq1size; i++) {
		delete [] diffM[i];
	}
	delete [] diffM;

	return distance;
}

double calcDistanceFor2Matrices(const vector<double*> matrix1, const vector<double*> matrix2, const uint16_t matricesWidth)
{
	extern bool option_WriteLog;
	extern ofstream log_file;

	double d1 = 0., distance = 0.;

	if (option_WriteLog)
	{
		log_file << "have matrix sizes " << matrix1.size() << "x" << matricesWidth <<" and " << setw (2) << matrix2.size() << "x" << matricesWidth << " ... " << flush;
	}

	for (uint16_t i = 0; i < matricesWidth; i++)
	{
		d1 = calcDistance(matrixSlicer(matrix1, i), matrix1.size(), matrixSlicer(matrix2, i), matrix2.size());
		distance += d1;
	}
	return distance;
}

double findDistance(uint32_t seq1size, uint32_t seq2size, double** diffM) {

	// Create distance matrix (forward direction)
	double** pathM = new double*[seq1size];
	for (uint32_t i = 0; i < seq1size; i++) {
		pathM[i] = new double[seq2size];
	}

	pathM[0][0] = diffM[0][0];
	for (uint32_t i = 1; i < seq1size; i++) {
		pathM[i][0] = diffM[i][0] + pathM[i - 1][0];
	}
	for (uint32_t j = 1; j < seq2size; j++) {
		pathM[0][j] = diffM[0][j] + pathM[0][j - 1];
	}

	for (uint32_t i = 1; i < seq1size; i++) {
		for (uint32_t j = 1; j < seq2size; j++) {
			if (pathM[i - 1][j - 1] < pathM[i - 1][j]) {
				if (pathM[i - 1][j - 1] < pathM[i][j - 1]) {
					pathM[i][j] = diffM[i][j] + pathM[i - 1][j - 1];
				} else {
					pathM[i][j] = diffM[i][j] + pathM[i][j - 1];
				}
			} else {
				if (pathM[i - 1][j] < pathM[i][j - 1]) {
					pathM[i][j] = diffM[i][j] + pathM[i - 1][j];
				} else {
					pathM[i][j] = diffM[i][j] + pathM[i][j - 1];
				}
			}
		}
	}

	// Find the warping path (backward direction)
	uint32_t warpSize = seq1size * seq2size;
	double* warpPath = new double[warpSize];

	uint32_t warpPathIndex = 0;
	uint32_t i = seq1size - 1, j = seq2size - 1;

	warpPath[warpPathIndex] = pathM[i][j];

	do {
		if (i > 0 && j > 0) {

			if (pathM[i - 1][j - 1] < pathM[i - 1][j]) {
				if (pathM[i - 1][j - 1] < pathM[i][j - 1]) {
					i--;
					j--;
				} else {
					j--;
				}

			} else {
				if (pathM[i - 1][j] < pathM[i][j - 1]) {
					i--;
				} else {
					j--;
				}
			}

		} else {
			if (0 == i) {
				j--;
			} else {
				i--;
			}
		}

		warpPath[++warpPathIndex] = pathM[i][j];

	} while (i > 0 || j > 0);

	// Calculate path measure
	double distance = 0.;
	for (uint32_t k = 0; k < warpPathIndex + 1; k++) {
		distance += warpPath[k];
	}
	distance = distance / (warpPathIndex + 1);

	// Clean up
	delete[] warpPath;
	for (uint32_t i = 0; i < seq1size; i++) {
		delete[] pathM[i];
	}
	delete[] pathM;

	return distance;
}
