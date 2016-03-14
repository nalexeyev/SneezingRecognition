#ifdef _MSC_VER

#define NOMINMAX
#include <algorithm>

#endif

#include "DataPrep.h"
#include "ServiceFunc.h"

using namespace std;

/**
 * Read Wav data from a file
 */
WavData* WavData::readFromFile(const std::string& file) {
	WavHeader wavHeader;

	// Open file
	std::fstream fs;
	fs.open(file.c_str(), std::ios::in | std::ios::binary);

	if (!fs.good()) {
		fprintf(stderr, "Input file not found: %s\n", file.c_str());
		return NULL;
	}

	// Read header
	fs.read((char*)(&wavHeader), sizeof(WavHeader));
	if (!checkHeader(wavHeader)) {
		return NULL;
	}

	// Read raw data
	WavData* wavData = new WavData(wavHeader);
	readData(fs, wavHeader, *wavData);
	fs.close();

	return wavData;
}

/**
 * Checks a set of restrictions
 */
bool WavData::checkHeader(const WavHeader& wavHeader) {

	if ((0 != strncmp(wavHeader.riff, "RIFF", sizeof(wavHeader.riff))) || (0 != strncmp(wavHeader.wave, "WAVE", sizeof(wavHeader.wave)))) {
		fprintf(stderr, "Invalid RIFF/WAVE format\n");
		return false;
	}

	if (1 != wavHeader.audioFormat) {
		fprintf(stderr, "Invalid WAV format: only PCM audio format is supported\n");
		return false;
	}

	if (wavHeader.numOfChan > 2) {
		fprintf(stderr, "Invalid WAV format: only 1 or 2 channels audio is supported\n");
		return false;
	}

	uint32_t bitsPerChannel = wavHeader.bitsPerSample / wavHeader.numOfChan;
	if (16 != bitsPerChannel) {
		fprintf(stderr, "Invalid WAV format: only 16-bit per channel is supported\n");
		return false;
	}

	assert(wavHeader.subchunk2Size > 0);
	if (wavHeader.subchunk2Size > LONG_MAX) {
		fprintf(stderr, "File too big\n");
		return false;
	}

	return true;
}

void WavData::readData(std::fstream& fs, const WavHeader& wavHeader, WavData& wavFile) {
	raw_t value, minValue = 0, maxValue = 0; // in dB
	int16_t value16, valueLeft16, valueRight16;

	uint32_t bytesPerSample = static_cast<uint32_t>(wavHeader.bitsPerSample / 8);
	uint32_t numberOfSamplesXChannels = wavHeader.subchunk2Size /
			(wavHeader.numOfChan * bytesPerSample);

	wavFile.rawData = new raw_t[numberOfSamplesXChannels];

	uint32_t sampleNumber;
	for (sampleNumber=0; sampleNumber < numberOfSamplesXChannels && !fs.eof(); sampleNumber++) {

		if (1 == wavHeader.numOfChan) {
			fs.read((char*)(&value16), sizeof(int16_t));
			value = static_cast<raw_t>(value16);

		} else {
			fs.read((char*)(&valueLeft16), sizeof(int16_t));
			fs.read((char*)(&valueRight16), sizeof(int16_t));
			value = static_cast<raw_t>((abs(valueLeft16) + abs(valueRight16)) / 2);
		}

		if (maxValue < value) {
			maxValue = value;
		}

		if (minValue > value) {
			minValue = value;
		}

		wavFile.rawData[sampleNumber] = value;
	}
	assert(sampleNumber > 0);

	// Normalization
	wavFile.normalizedData = new double[numberOfSamplesXChannels];
	double maxAbs = max(fabs(minValue), fabs(maxValue));

	for (sampleNumber = 0; sampleNumber < numberOfSamplesXChannels; sampleNumber++) {
		wavFile.normalizedData[sampleNumber] = static_cast<double>(wavFile.rawData[sampleNumber]) / maxAbs;
		//	cout << wavFile.normalizedData[sampleNumber] << " : " << wavFile.rawData[sampleNumber] << endl;
	}

	// Update values
	wavFile.setMinVal(minValue);
	wavFile.setMaxVal(maxValue);
	wavFile.setNumberOfSamples(numberOfSamplesXChannels);
}

bool separateSamplesToFrames(WavData* data, std::vector<sFrame> *frames)
{
	extern bool option_SilentMode;

	sFrame sfr;

	uint32_t bytesPerFrame = static_cast<uint32_t>((data->getHeader()).bytesPerSec * FRAME_LENGTH / 1000.0);
	uint32_t bytesPerSample = static_cast<uint32_t>((data->getHeader()).bitsPerSample / 8);
	uint32_t samplesPerFrame = static_cast<uint32_t>(bytesPerFrame / bytesPerSample);
	//assert("Number of samples per frame cannot be less or equal than 0" && samplesPerFrame > 0);

	if (samplesPerFrame <= 0)
	{
		if (option_SilentMode)
			cout << "wrong header, file will be skipped. ";
		else
			cout << "wrong header, file will be skipped (number of samples per frame = " << samplesPerFrame << ", but cannot be less or equal than 0. ";

		return false;
	}

	uint32_t samplesPerNonOverlap =	static_cast<uint32_t>(samplesPerFrame * (1 - FRAME_OVERLAP));
	uint32_t framesCount =	((data->getHeader()).subchunk2Size / bytesPerSample) / samplesPerNonOverlap;
	//assert("File header is corrupted: subchunk2Size expected not less than 22040" && (data->getHeader()).subchunk2Size > 22040);

	if ((data->getHeader()).subchunk2Size < 10000)
	{
		if (option_SilentMode)
			cout << "Wrong header, file will be skipped.";
		else
			cout << "Wrong header, file will be skipped (subchunk 2 size = " << (data->getHeader()).subchunk2Size << ", but expected not less than 22040). ";

		return false;
	}


	uint32_t indexBegin = 0, indexEnd = 0;
	uint32_t size = data->getNumberOfSamples();
	uint32_t frameId = 0;

	for (frameId = 0; frameId < framesCount; ++frameId) {

		indexBegin = frameId * samplesPerNonOverlap;
		indexEnd = indexBegin + samplesPerFrame;
		if (indexEnd <= size) {
			sfr.id = frameId;
			sfr.firstSample = frameId * samplesPerNonOverlap;
			sfr.lastSample = sfr.firstSample + samplesPerFrame;
			sfr.rms = RMS(data->getRawData(), sfr.firstSample, sfr.lastSample);
			sfr.entropy = Entropy(data->getNormalizedData(), sfr.firstSample, sfr.lastSample, ENTROPY_BINS, -1, 1);

			sfr.mfcc = calcMFCC (data->getNormalizedData(), sfr.firstSample, sfr.lastSample, MFCC_SIZE, data->getHeader().samplesPerSec, MFCC_FREQ_MIN, MFCC_FREQ_MAX);

			frames->insert(frames->begin() + frameId, sfr);
		} else {
			break;
		}
	}
	return true;
}

bool findSilenceThreshold(std::vector<sFrame> frames, bool *hasSilence, double *rmsMax, double *soundThreshold) {

	*rmsMax = 0;
	double rms, rmsSilence = 0., entropyMax=-100., entropyMin=100., entropyAvg = 0.;
	rms = *rmsMax = frames.at(0).rms;

	*hasSilence = false;
	uint32_t cnt = 0;
//	cout << endl;

	for (vector<sFrame>::const_iterator frmIt = frames.begin(); frmIt != frames.end(); ++frmIt) {
		entropyMax = std::max(entropyMax, frmIt->entropy);
		entropyMin = std::min(entropyMin, frmIt->entropy);
		entropyAvg = (entropyMax + entropyMin) / 2;
	}

	for (vector<sFrame>::const_iterator frmIt = frames.begin(); frmIt != frames.end(); ++frmIt) {

		rms = frmIt->rms;
		*rmsMax = std::max(*rmsMax, rms);

		if (frmIt->entropy < entropyAvg*0.9) {
			*hasSilence = true;
			rmsSilence += frmIt->rms;
			cnt++;
		}


//		cout << frmIt->id << " " <<frmIt->entropy << " " << frmIt->rms << endl;
	}
	rmsSilence /= cnt;

	*soundThreshold = rmsSilence * SOUND_THRESHOLD_COEF;

	return true;
}

bool separateFramesToSounds(std::vector<sFrame> frames, std::vector<sSound> * sounds) {
	std::vector<sSound> ressounds;

	//assert(frames.size() > 10);

	bool hasSilence;
	double rmsMax = 0., soundThreshold = 0.;

	findSilenceThreshold(frames, &hasSilence, &rmsMax, &soundThreshold);
	//cout << hasSilence << " " << rmsMax << " " << soundThreshold << endl;

	int32_t SoundId = -1;
	int32_t firstFrameInCurrentSoundNumber = -1;
	uint32_t lastFrameInCurrentSoundNumber = 0;
	int32_t silenceCnt = 0;

	sSound tmpsSound;


	if (hasSilence) {
		for (vector<sFrame>::const_iterator frmIt = frames.begin(); frmIt != frames.end(); ++frmIt) {

			// Got a sound at 1st time or after a silence
			if ((frmIt->rms > soundThreshold) && (firstFrameInCurrentSoundNumber == -1)) {
				firstFrameInCurrentSoundNumber = frmIt->id;
				silenceCnt = 0;
			}
			// Sound continues
			if ((frmIt->rms > soundThreshold) && (firstFrameInCurrentSoundNumber != -1)) {
				lastFrameInCurrentSoundNumber = frmIt->id;
				silenceCnt = 0;
			}
			// Silence started after a sound
			if ((frmIt->rms <= soundThreshold) && (firstFrameInCurrentSoundNumber >= 0)) {
				++silenceCnt;

				if (silenceCnt > SOUNDS_MIN_DISTANCE) {

					if ((lastFrameInCurrentSoundNumber - firstFrameInCurrentSoundNumber) >= SOUND_MIN_SIZE) { ///////////////////////// Вот тут проверку на слишком длинные звуки!?
						++SoundId;

						tmpsSound.reset();
						tmpsSound.id = SoundId;
						tmpsSound.firstFrame = firstFrameInCurrentSoundNumber;
						tmpsSound.lastFrame = lastFrameInCurrentSoundNumber;
						for (uint32_t i = firstFrameInCurrentSoundNumber; i <= lastFrameInCurrentSoundNumber; i++)
						{
							tmpsSound.vMFCC.push_back(frames.at(i).mfcc);
						}
						ressounds.push_back(tmpsSound);
						sounds->push_back(tmpsSound);
					}

					firstFrameInCurrentSoundNumber = -1;
				}
			}
			// If there's not enough silence after the last sound and before the end of file
			if ((lastFrameInCurrentSoundNumber >= (frames.back().id - SOUNDS_MIN_DISTANCE)) && ((lastFrameInCurrentSoundNumber - firstFrameInCurrentSoundNumber) >= SOUND_MIN_SIZE)) {
				++SoundId;

				tmpsSound.reset();
				tmpsSound.id = SoundId;
				tmpsSound.firstFrame = firstFrameInCurrentSoundNumber;
				tmpsSound.lastFrame = lastFrameInCurrentSoundNumber;
				for (uint32_t i = firstFrameInCurrentSoundNumber; i <= lastFrameInCurrentSoundNumber; i++)
				{
					tmpsSound.vMFCC.push_back(frames.at(i).mfcc);
				}
				ressounds.push_back(tmpsSound);
				sounds->push_back(tmpsSound);

				break;
			}
		}
	}

	if (!hasSilence || (SoundId==-1))// There's no silence, whole file is one sound
	{
		double soundRMS1 = 0., soundRMS = 0.;
		tmpsSound.reset();
		firstFrameInCurrentSoundNumber = 0;
		lastFrameInCurrentSoundNumber = frames.back().id;

		for (auto frmIt = frames.begin(); frmIt != frames.end(); frmIt++)
		{
			soundRMS1 += frmIt->rms;
		}

		for (sFrame &frmIt : frames)

		{
			soundRMS += frmIt.rms;
		}



		soundRMS /= lastFrameInCurrentSoundNumber - firstFrameInCurrentSoundNumber;

		if (soundRMS < soundThreshold / SOUND_THRESHOLD_COEF) return false;

		tmpsSound.id = 0;
		tmpsSound.firstFrame = firstFrameInCurrentSoundNumber;
		tmpsSound.lastFrame = lastFrameInCurrentSoundNumber;
		for (uint32_t i = firstFrameInCurrentSoundNumber; i <= lastFrameInCurrentSoundNumber; i++)
			{
				tmpsSound.vMFCC.push_back(frames.at(i).mfcc);
			}

		ressounds.push_back(tmpsSound);
		sounds->push_back(tmpsSound);

		SoundId++;
	}

	if (SoundId==-1) return false;
			else	return true;

}

/*
bool enhanceSounds(std::vector<sSound> sounds, std::vector<sSound> * enhSounds) {


	for (sSound &item : sounds)  // range based
		{

		//	resfile << std::setw(2) << item.id << " ; "<< std::setw(5) << item.firstFrame << " ; " << std::setw(5) << item.lastFrame <<  endl;
		}

	return true;
}
*/


bool saveSoundAsAudio(const std::string& file, const std::vector<sFrame> frames, const sSound sound, WavData *wavdata) {

	uint32_t sampleStart = frames.at(sound.firstFrame).firstSample;
	uint32_t sampleFinish = frames.at(sound.lastFrame).lastSample;
	uint32_t waveSize = (sampleFinish - sampleStart) * sizeof(raw_t);

	// prepare a new header and write it to file stream
	WavHeader headerNew;
	strncpy(headerNew.riff, wavdata->getHeader().riff, 4);
	headerNew.chunkSize = waveSize + sizeof(WavHeader);
	strncpy(headerNew.wave, wavdata->getHeader().wave, 4);
	strncpy(headerNew.fmt, wavdata->getHeader().fmt, 4);
	headerNew.subchunk1Size = wavdata->getHeader().subchunk1Size;
	headerNew.audioFormat = wavdata->getHeader().audioFormat;
	headerNew.numOfChan = 1;
	headerNew.samplesPerSec = wavdata->getHeader().samplesPerSec;
	headerNew.bytesPerSec = wavdata->getHeader().samplesPerSec * sizeof(raw_t);
	headerNew.blockAlign = sizeof(raw_t);
	headerNew.bitsPerSample = sizeof(raw_t) * 8;
	strncpy(headerNew.data, wavdata->getHeader().data, 4);
	headerNew.subchunk2Size = waveSize;

	std::ofstream fs;
	fs.open(file.c_str(), std::ios::out | std::ios::binary);
	if (!fs.good()) {
		fprintf(stderr, "Output file is not created: %s\n", file.c_str());
		return false;
	}
	fs.write((char*)&headerNew, sizeof(WavHeader));

	raw_t* data = new raw_t[waveSize / sizeof(raw_t)];

	uint32_t i = 0;


	for (uint32_t currentSample = sampleStart; currentSample <= sampleFinish; currentSample++) {

		data[i] = wavdata->getRawData()[currentSample];

		++i;
	}

	fs.write((char*)data, waveSize);
	fs.close();
	delete [] data;

	return true;
}


bool makeTrainDataFile(const std::string& trainInputFolder, const std::string& trainDataFilePath)
{

	const char *  trainWaveFilesFolder = trainInputFolder.c_str();
	WavData* trainWave;
	std::string filename, path ,fullpath;

	std::vector<sFrame> trainFrames;

	std::fstream trainfile;
	trainfile.open(trainDataFilePath.c_str(), std::ios::out | std::ios::binary);

	cout<< "Writing of MFCC from wave files to a train data file " << trainDataFilePath << " is started..."<< endl;

	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (trainWaveFilesFolder)) != NULL)
	{
		while ((ent = readdir (dir)) != NULL)
		{
			if (ent->d_name[0]!='.')
			{
				cout<< ent->d_name << " - " << flush;

				fullpath = trainInputFolder + filename.assign(ent->d_name);

				trainFrames.clear();

				trainWave = WavData::readFromFile(fullpath);
				if (!trainWave)
				{
					return false;
				}

				separateSamplesToFrames(trainWave, &trainFrames);
		//		path = filename + "_trainframes_formatted.txt";
		//		printframes(trainFrames,path);

				trainfile << ent->d_name << std::endl;
				trainfile << trainFrames.size() << std::endl;
				cout << trainFrames.size() << " frames" << endl;

				for (std::vector<sFrame>::const_iterator item = trainFrames.begin(); item!=trainFrames.end(); ++item)
				{
					for (int i = 0; i < MFCC_SIZE; ++i)
					{
						trainfile << std::setprecision(12) << item->mfcc[i] << std::endl;
					}
				}
			}
		}
		closedir (dir);
		trainfile.close();
	}
	else
	{
		cerr << "Cannot open directory " << trainDataFilePath << endl;
		return false;
	}
	return true;
}

bool readTrainDataFile(const std::string& trainDataFilePath, std::vector<svMFCC> *trainMatrices)
{
	extern bool option_SilentMode;

	svMFCC item;
	std::fstream trainfile;
	trainfile.open(trainDataFilePath, std::ios::in | std::ios::binary);

	if (!trainfile.good())
	{
		cerr << "Train data file not found: " << trainDataFilePath << endl;
		return false;
	}

	if (!option_SilentMode) {cout<< "Reading of MFCC from a train data file " << trainDataFilePath << " is started..."<< flush;}

	std::string waveFileName, strFramesNumber, strMFCC;

	uint32_t FramesNumber = 0;
	double * MFCC1 = new double[100000];

	while (!trainfile.eof())
	{
		item.reset();

		getline(trainfile, waveFileName);
		if (waveFileName == "")
		{
			//cout<< "Empty string" << endl;
			break;
		}

		getline(trainfile, strFramesNumber);

		FramesNumber = atoi(strFramesNumber.c_str());

		item.waveFilename = waveFileName;
		item.numberOfFrames = FramesNumber;

		for (uint32_t i = 0; i < FramesNumber; i++)
		{
			item.vMFCC.push_back(MFCC1);

			for (uint32_t j = 0; j < MFCC_SIZE; j++, MFCC1++)
			{
				getline(trainfile, strMFCC);
				*MFCC1 = stod(strMFCC);
			}
		}


		trainMatrices->push_back(item);

	}
	trainfile.close();
	if (!option_SilentMode) {cout << "Completed. " << endl << endl ;}
	return true;
}

double* matrixSlicer(std::vector<double*> matrix, uint16_t sliceNo)
{
	double* arr = new double[matrix.size()];

	uint16_t i = 0;
	for (auto it = matrix.begin(); it!=matrix.end(); ++it)
	{
		arr[i] = (*it)[sliceNo];
		++i;
	}

	return arr;
}


int makeDecision (std::string inputFile, std::vector<svMFCC> etalonframes)
{
	extern bool option_WriteLog;
	extern bool option_SilentMode;
	extern bool option_WriteSplittedSounds;
	extern bool option_Recognize;

	extern std::string setting_InputDataFolder;
	extern std::string setting_ResultsFileName;

	extern ofstream log_file;

	bool isOk;
	std::vector<sFrame> frames;
	std::vector<sSound> sounds;


	WavData* wavData = WavData::readFromFile(setting_InputDataFolder + inputFile);

	isOk = separateSamplesToFrames(wavData, &frames);
	if (!isOk) {cout << "Something went wrong during separation of samples to frames. " << flush; return -1;}

	//cout<< endl; printframes(frames, "frames.txt");

	isOk = separateFramesToSounds(frames, &sounds);
	if (!isOk) {cout << "Something went wrong during separation of frames to sounds. " << endl; return -1;}


	double distance = 0.;
	double absoluteMinimum = 3000.;

	if (option_WriteLog)
		{
		log_file << "================================================================================================" << endl;
		log_file << " file: " << inputFile << endl;
		log_file << "================================================================================================" << endl;
		}


	if (option_Recognize)
	{
		for (std::vector<sSound>::const_iterator sitem = sounds.begin(); sitem != sounds.end(); ++sitem)
		{
			double distMin = 3000.;
			if (!option_SilentMode) {cout << "Sound #" << sitem->id << " " << flush;}

			for (std::vector<svMFCC>::const_iterator item = etalonframes.begin(); item != etalonframes.end(); ++item)
			{
				if (option_WriteLog) ////////////////////////// ???????????????????????????????????
					{
						log_file << "Sound #" << sitem->id << " and etalon file " << setw (8) << item->waveFilename << " ... " << flush;
					}

				distance = calcDistanceFor2Matrices(sitem->vMFCC, item->vMFCC, MFCC_SIZE);

					if (option_WriteLog)
					{
						log_file << "Distance is " << distance << endl;
					}

				distMin = std::min(distMin, distance);
			}

			if (!option_SilentMode) {cout << "min. distance = " << fixed << setprecision(2) << distMin << "; " << flush; }
			if (option_WriteLog)
				{
				log_file <<  "=================================" << endl;
				log_file <<  "  Minimal distance = " << distMin << "; " << endl;
				log_file <<  "=================================" << endl << endl;}

			absoluteMinimum = std::min(absoluteMinimum, distMin);
		}
		if (!option_SilentMode) {cout << " WAVE min. distance = " << fixed << setprecision(2) << absoluteMinimum << " ==> " << flush; }
		//cout << ":" << flush;
	}
	/** (Optional) Writing of Sounds extracted from input files *************************************/

	if (option_WriteSplittedSounds)
	{
		if (!option_SilentMode) { cout << endl << "Writing of Sounds extracted from input files is started..." <<flush;}

		system(("mkdir -p " + setting_InputDataFolder + "splitted" + separator()).c_str());

		//for (uint i = 0; i <= sounds.back().id; ++i)
		int i = 0;
		for (std::vector<sSound>::const_iterator it = sounds.begin(); it != sounds.end(); ++it)
		{
			std::string path(setting_InputDataFolder + "splitted" + separator()  + inputFile + "_Sound_" + to_string(i++) + ".wav");
			isOk = saveSoundAsAudio(path, frames, *it, wavData);
			if (!isOk) {cout << "Something went wrong during writing of Sounds extracted from input files" << endl; return -2;}
		}
		if (!option_SilentMode) {cout << "Done." << endl << endl;}
	}


	if (absoluteMinimum <= DECISION_DISTANCE_THRESHOLD)
		return 1;
	else
		return 0;
}


