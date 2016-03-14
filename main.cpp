#ifdef _MSC_VER
	#define NOMINMAX
	#include <algorithm>
#endif

#include "config.h"
#include "DataPrep.h"
#include "SpecMath.h"
#include "ServiceFunc.h"

using namespace std;

/** Options *************************************************************************************/
bool option_CreateTrainDataFile = false;
bool option_Recognize = true;
bool option_SilentMode = false;
bool option_WriteLog = true;
bool option_ResultsToFile = true;
bool option_WriteSplittedSounds = false;

/** Settings*************************************************************************************/
	std::string setting_TrainDataFileFolder ("training/");
	std::string setting_InputDataFolder("input");
	std::string setting_ResultsFileName("results.delmtd.csv");
	std::string setting_LogFileName("log.txt");

	ofstream res_file;
	ofstream log_file;


int main(int argc, char** argv) {

	bool isOk;
	int decision;
	std::vector<svMFCC> etalonframes;

/** (Optional) Creation of train data file ******************************************************/

	if (option_CreateTrainDataFile)
	{
		std::string trainInputFolder(setting_TrainDataFileFolder);
		std::string trainDataFile("trainData.dat");
		isOk = makeTrainDataFile(trainInputFolder, trainDataFile);
		if (!isOk) {cerr << "Something went wrong during creation of Train Data File" << endl; return EXIT_FAILURE;}
	}

/** Reading of train data file ******************************************************************/

	isOk = readTrainDataFile("trainData.dat", &etalonframes);
	if (!isOk) {cerr << "Something went wrong during reading of Train Data File" << endl; return EXIT_FAILURE;}

/** Initialize output files *********************************************************************/

	if (option_ResultsToFile)
			{
				res_file.open (setting_ResultsFileName);
			}

	if (option_WriteLog)
			{
				log_file.open (setting_LogFileName);
			}

/** Reading of input data files in folder and recognize it **************************************/

	DIR *inputDIR;
	struct dirent *DIRent;
	std::string inputFile ("");
	int numberOfFilesTotal = 0, numberOfFilesPositive = 0, numberOfFilesNegative = 0;

	setting_InputDataFolder = setting_InputDataFolder + SEPARATOR;

	if (!option_SilentMode) {cout << "Reading and analysis of files in folder " << setting_InputDataFolder << " is started:" << endl;}

	if ((inputDIR = opendir (setting_InputDataFolder.c_str())) != NULL)
		{
			while ((DIRent = readdir (inputDIR)) != NULL)
			{
				if (DIRent->d_name[0]!='.' && strcmp(DIRent->d_name, "splitted")!=0)
				{
					inputFile.assign(DIRent->d_name);

					cout << inputFile << ": " << flush;

					decision = makeDecision (inputFile, etalonframes);

					++numberOfFilesTotal;

					if (option_Recognize)
					{
					if (decision==1)
						{ cout << "Detected" << endl; ++numberOfFilesPositive; }
					else if (decision==0)
						{ cout << "Not detected" << endl; ++numberOfFilesNegative; }
					else if (decision==-2)
						{ cout << "Something went wrong during writing splitted sounds."<< endl << "Check the availability of folder named splitted" << endl; }
					else
						{ cout << "Something went wrong during making decision" << endl;}


						if (option_ResultsToFile)
						{
							res_file << inputFile << "," << decision << endl;
						}

					}
				}
			}
		}
	else
	{
		cerr << "Something went wrong while opening folder with input data files" << endl;
		return EXIT_FAILURE;
	}

	if (!option_SilentMode) {cout << "Completed. " << endl;}
	cout << endl << "Summary:" << endl;
	cout << "Total number of files processed: "<< numberOfFilesTotal << endl;
	cout << "Sneezing   sound   detected: "<< numberOfFilesPositive << endl;
	cout << "Sneezing sound NOT detected: "<< numberOfFilesNegative << endl << endl;

	if (option_ResultsToFile)
		{
		cout << "Results were exported to CSV file " << setting_ResultsFileName << endl;
		res_file.close();
		}

	if (option_WriteLog)
		{
		cout << "Log of analysis process was saved to file" << setting_LogFileName << endl;
		log_file.close();
		}


	return EXIT_SUCCESS;
}

