#ifndef SERVICEFUNC_H_
#define SERVICEFUNC_H_

#include "config.h"
#include "DataPrep.h"
#include "SpecMath.h"

int printframes(std::vector<sFrame> framestoprint);
int printframes(std::vector<sFrame> framestoprint, const std::string& file);
int printsounds(std::vector<sSound> soundstoprint);
int printsounds(std::vector<sSound> soundsToPrint, const std::string& file);
int printMFCCmatrix(std::vector<double*> vmfcc, uint32_t itemLength);
int printMFCCarr(double* mfcc, uint32_t itemLength);

#endif /* SERVICEFUNC_H_ */
