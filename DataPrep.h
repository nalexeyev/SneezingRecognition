/**
 * Represents WAV file data
 *
 * Currently supports only PCM format.
 *
 * @see http://en.wikipedia.org/wiki/WAV
 * @see http://en.wikipedia.org/wiki/Linear_pulse-code_modulation
 * @see https://ccrma.stanford.edu/courses/422/projects/WaveFormat/
 */
#ifndef DATA_PREP_H_
#define DATA_PREP_H_

#include "config.h"

using namespace std;

typedef int16_t raw_t;

/***************************************************************************************************
 * WAV header
 ***************************************************************************************************/
struct WavHeader {
    char     riff[4];        // RIFF Header
    uint32_t chunkSize;      // RIFF Chunk Size
    char     wave[4];        // WAVE Header

    char     fmt[4];         // FMT header
    uint32_t subchunk1Size;  // Size of the fmt chunk
    uint16_t audioFormat;    // Audio format 1=PCM (Other formats are unsupported)
    uint16_t numOfChan;      // Number of channels 1=Mono, 2=Stereo
    uint32_t samplesPerSec;  // Sampling Frequency in Hz
    uint32_t bytesPerSec;    // bytes per second
    uint16_t blockAlign;     // 2=16-bit mono, 4=16-bit stereo
    uint16_t bitsPerSample;  // Number of bits per sample

    // The data below depends on audioFormat, but we work only with PCM cases
    char     data[4];        // DATA header
    uint32_t subchunk2Size;  // Sampled data length
};

/***************************************************************************************************
 * WAV data
 ***************************************************************************************************/
class WavData {
public:

	~WavData() {
		if (NULL != this->rawData) {
			delete [] this->rawData;
		}
		if (NULL != this->normalizedData) {
			delete [] this->normalizedData;
		}
	}

	static WavData* readFromFile(const std::string& file);

	uint32_t getNumberOfSamples() const { return numberOfSamples; }
	void setNumberOfSamples(uint32_t numberOfSamples) { this->numberOfSamples = numberOfSamples; }

	raw_t getMaxVal() const { return maxVal; }
	void setMaxVal(raw_t maxVal) { this->maxVal = maxVal; }

	raw_t getMinVal() const { return minVal; }
	void setMinVal(raw_t minVal) { this->minVal = minVal; }

	const WavHeader& getHeader() const { return header; }
	const raw_t* getRawData() const { return rawData; }
	const double* getNormalizedData() const { return normalizedData; }

private:
	WavHeader		header;
	raw_t*			rawData;
	double*			normalizedData;

	raw_t			maxVal;
	raw_t			minVal;
	uint32_t		numberOfSamples;

	WavData(WavHeader header) {
		this->header = header;
		this->rawData = NULL;
		this->normalizedData = NULL;

		this->maxVal = 0;
		this->minVal = 0;
		this->numberOfSamples = 0;
	}

	static bool checkHeader(const WavHeader& wavHeader);
	static void readData(std::fstream& fs, const WavHeader& wavHeader, WavData& wavFile);
};

typedef struct
{
	WavHeader		header;
	raw_t*			rawData;
	double*			normalizedData;

	raw_t			maxVal;
	raw_t			minVal;
	uint32_t		numberOfSamples;
} sWave;

typedef struct
{
	uint32_t id;
	uint32_t firstSample;
	uint32_t lastSample;

	double rms;
	double entropy;

	double* mfcc;
} sFrame;

typedef struct
{
	uint32_t id;
	uint32_t firstFrame;
	uint32_t lastFrame;
	std::vector<double*> vMFCC;
	void reset()
		{
		id = 0;
		firstFrame = 0;
		lastFrame = 0;
		vMFCC.clear();
		}
} sSound;

typedef struct
{
	std::string waveFilename;
	uint32_t numberOfFrames;
	std::vector<double*> vMFCC;
	void reset()
	{
		waveFilename.clear();
		numberOfFrames = 0;
		vMFCC.clear();
	}
} svMFCC;

bool separateSamplesToFrames(WavData * data, std::vector<sFrame> *frames);
bool separateFramesToSounds(std::vector<sFrame> frames, std::vector<sSound> * sounds);
bool saveSoundAsAudio(const std::string& file, const std::vector<sFrame> frames, const sSound sound, WavData * wavdata);
bool makeTrainDataFile(const std::string& trainInputFolder, const std::string& trainDataFile);
bool readTrainDataFile(const std::string& trainDataFilePath, std::vector<svMFCC> *trainMatrices);
double* matrixSlicer(std::vector<double*> matrix, uint16_t sliceNo);
int makeDecision (std::string inputFile, std::vector<svMFCC> etalonframes);

#endif /* DATA_PREP_H_ */
