/*
 * funclib.cpp
 *
 *  Created on: Feb 19, 2016
 *      Author: nick
 */

#include "ServiceFunc.h"


using namespace std;

int printframes(std::vector<sFrame> framestoprint){

 	cout << "Id ; Start ; Finish; Entropy;  RMS    ;  MFCC1  ;  MFCC2 ;  MFCC3 ;  MFCC4 ;  MFCC5 ;  MFCC6 ;  MFCC7 ;  MFCC8 ;  MFCC9 ;  MFCC10;  MFCC11;  MFCC12;" << endl;
    cout << "---+-------+-------+--------+---------+-------" << endl;
	for (auto it = framestoprint.begin(); it != framestoprint.end(); ++it)
	{

		cout << std::setw(2) << it->id << " ; "<< std::setw(5) << it->firstSample << " ; " << std::setw(5) << it->lastSample << " ; "
			 << std::right   << std::fixed << std::setprecision(2)
		     << std::setw(6) << it->entropy  << " ; "
		     << std::setw(7) << it->rms      << " ; ["
		     << std::setw(6) << it->mfcc[0]  << " ; "
		     << std::setw(6) << it->mfcc[1]  << " ; "
		     << std::setw(6) << it->mfcc[2]  << " ; "
		     << std::setw(6) << it->mfcc[3]  << " ; "
		     << std::setw(6) << it->mfcc[4]  << " ; "
		     << std::setw(6) << it->mfcc[5]  << " ; "
		     << std::setw(6) << it->mfcc[6]  << " ; "
		     << std::setw(6) << it->mfcc[7]  << " ; "
		     << std::setw(6) << it->mfcc[8]  << " ; "
		     << std::setw(6) << it->mfcc[9]  << " ; "
		     << std::setw(6) << it->mfcc[10]  << " ; "
		     << std::setw(6) << it->mfcc[11]  << " ] " << endl;
	}

	return 0;
}

int printframes(std::vector<sFrame> framestoprint, const std::string& file){

	ofstream resfile;
	resfile.open (file);

	resfile << "Id ; Start ; Finish; Entropy;  RMS    ;  MFCC1  ;  MFCC2 ;  MFCC3 ;  MFCC4 ;  MFCC5 ;  MFCC6 ;  MFCC7 ;  MFCC8 ;  MFCC9 ;  MFCC10;  MFCC11;  MFCC12;" << endl;

	for (auto it = framestoprint.begin(); it != framestoprint.end(); ++it)
	{

     resfile << std::setw(2) << it->id << " ; "<< std::setw(5) << it->firstSample << " ; " << std::setw(5) << it->lastSample << " ; "
			 << std::right   << std::fixed << std::setprecision(2)
		     << std::setw(6) << it->entropy  << " ; "
		     << std::setw(7) << it->rms      << " ; ["
		     << std::setw(6) << it->mfcc[0]  << " ; "
		     << std::setw(6) << it->mfcc[1]  << " ; "
		     << std::setw(6) << it->mfcc[2]  << " ; "
		     << std::setw(6) << it->mfcc[3]  << " ; "
		     << std::setw(6) << it->mfcc[4]  << " ; "
		     << std::setw(6) << it->mfcc[5]  << " ; "
		     << std::setw(6) << it->mfcc[6]  << " ; "
		     << std::setw(6) << it->mfcc[7]  << " ; "
		     << std::setw(6) << it->mfcc[8]  << " ; "
		     << std::setw(6) << it->mfcc[9]  << " ; "
		     << std::setw(6) << it->mfcc[10]  << " ; "
		     << std::setw(6) << it->mfcc[11]  << " ] " << endl;
	}
    resfile.close();
	return 0;
}

	int printsounds(std::vector<sSound> soundsToPrint)
	{

		for (sSound &item : soundsToPrint)  // range based
		{
		 	cout << "Id ; 1stFrm ; LastFrm; " << endl;
		    cout << "---+-------+---------" << endl;
			cout << std::setw(2) << item.id << " ; "<< std::setw(5) << item.firstFrame << " ; " << std::setw(5) << item.lastFrame <<  endl;
		    cout << "---+-------+---------" << endl;

			cout << "  MFCC1  ;  MFCC2 ;  MFCC3 ;  MFCC4 ;  MFCC5 ;  MFCC6 ;  MFCC7 ;  MFCC8 ;  MFCC9 ;  MFCC10;  MFCC11;  MFCC12;"  <<  endl;
			cout << "+-----------------------------------------------------------------------------------------------------------+"  <<  endl;

//			for (double *MFCC : item.vMFCC)
//			cout << std::setw(6) << MFCC[0]  << " ; "

			for (std::vector<double*>::const_iterator MFCC = item.vMFCC.begin(); MFCC != item.vMFCC.end(); ++MFCC)
			{
				//double * temp = *MFCC;  !!! iterator eto object
				cout << "[ ";
				cout << std::right   << std::fixed << std::setprecision(2);
				cout << std::setw(6) << (*MFCC)[0]  << " ; "
				     << std::setw(6) << (*MFCC)[1]   << " ; "
				     << std::setw(6) << (*MFCC)[2]  << " ; "
				     << std::setw(6) << (*MFCC)[3]  << " ; "
				     << std::setw(6) << (*MFCC)[4]  << " ; "
				     << std::setw(6) << (*MFCC)[5]  << " ; "
				     << std::setw(6) << (*MFCC)[6]  << " ; "
				     << std::setw(6) << (*MFCC)[7]  << " ; "
				     << std::setw(6) << (*MFCC)[8]  << " ; "
				     << std::setw(6) << (*MFCC)[9]  << " ; "
				     << std::setw(6) << (*MFCC)[10]  << " ; "
				     << std::setw(6) << (*MFCC)[11]  << " ] " << endl;
			}
			cout << endl << endl;
		}
		return 0;
	}


	int printsounds(std::vector<sSound> soundsToPrint, const std::string& file)
	{

		ofstream resfile;
		resfile.open (file);

		for (sSound &item : soundsToPrint)  // range based
		{
			resfile << "Id ; 1stFrm ; LastFrm; " << endl;
			resfile << "---+-------+---------" << endl;
			resfile << std::setw(2) << item.id << " ; "<< std::setw(5) << item.firstFrame << " ; " << std::setw(5) << item.lastFrame <<  endl;
			resfile << "---+-------+---------" << endl;

			resfile << "  MFCC1  ;  MFCC2 ;  MFCC3 ;  MFCC4 ;  MFCC5 ;  MFCC6 ;  MFCC7 ;  MFCC8 ;  MFCC9 ;  MFCC10;  MFCC11;  MFCC12;"  <<  endl;
			resfile << "+-----------------------------------------------------------------------------------------------------------+"  <<  endl;

//			for (double *MFCC : item.vMFCC)
//			cout << std::setw(6) << MFCC[0]  << " ; "

			for (std::vector<double*>::const_iterator MFCC = item.vMFCC.begin(); MFCC != item.vMFCC.end(); ++MFCC)
			{
				//double * temp = *MFCC;  !!! iterator eto object
				resfile << "[ ";
				resfile << std::right   << std::fixed << std::setprecision(2);
				resfile << std::setw(6) << (*MFCC)[0]  << " ; "
				     << std::setw(6) << (*MFCC)[1]   << " ; "
				     << std::setw(6) << (*MFCC)[2]  << " ; "
				     << std::setw(6) << (*MFCC)[3]  << " ; "
				     << std::setw(6) << (*MFCC)[4]  << " ; "
				     << std::setw(6) << (*MFCC)[5]  << " ; "
				     << std::setw(6) << (*MFCC)[6]  << " ; "
				     << std::setw(6) << (*MFCC)[7]  << " ; "
				     << std::setw(6) << (*MFCC)[8]  << " ; "
				     << std::setw(6) << (*MFCC)[9]  << " ; "
				     << std::setw(6) << (*MFCC)[10]  << " ; "
				     << std::setw(6) << (*MFCC)[11]  << " ] " << endl;
			}
			resfile << endl << endl;
		}
	    resfile.close();
		return 0;
	}

	int printMFCCmatrix(std::vector<double*> vmfcc, uint32_t itemLength)
	{
		for(auto it = vmfcc.begin(); it!= vmfcc.end(); ++it)
		{
			printMFCCarr((*it), itemLength);
		}
		return 0;
	}

	int printMFCCarr(double* mfcc, uint32_t itemLength)
	{
		for (uint32_t i=0; i< itemLength; ++i)
				{
					cout << mfcc[i] << "; " << flush;
				}
				cout << endl;
		return 0;
	}
