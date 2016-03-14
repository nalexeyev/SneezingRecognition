/*
 * config.h
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */


#ifndef CONFIG_H_
#define CONFIG_H_


#define __GXX_EXPERIMENTAL_CXX0X__ 1


#ifdef _WIN32
#define SEPARATOR "\\"
#else
#define SEPARATOR "/"
#endif

inline char separator()
{
#if defined _WIN32 || defined __CYGWIN__
    return '\\';
#else
    return '/';
#endif
}

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <fstream>
#include <cmath>
#include <complex>
#include <valarray>
#include <cstring>
#include <limits.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <iomanip>


	/**
	 * Raw audio data type
	 */
	typedef int16_t raw_t;

	/**
	 * Length of frame (ms)
	 */
	const uint16_t FRAME_LENGTH = 50;

	/**
	 * Percentage of overlap for frames (0 <= x < 1)
	 */
	const double FRAME_OVERLAP = 0.5;

	/**
	 * Minimal size of sound (in frames)
	 */
	const uint16_t SOUND_MIN_SIZE = (100 / FRAME_LENGTH) / (1 - FRAME_OVERLAP); // 4 frames, 100 ms

	/**
	 * Maximal length of sound (in frames)
	 */
	const uint16_t SOUND_MAX_SIZE = (1000 / FRAME_LENGTH) / (1 - FRAME_OVERLAP); // 40 frames, 1000 ms


	/**
	 * Minimal amount of framer between two sounds
	 */
	const uint16_t SOUNDS_MIN_DISTANCE = SOUND_MIN_SIZE * 0.25;

	/**
	 * Amount of MFCC coefficients
	 */
	const uint16_t MFCC_SIZE = 12;

	/**
	 * Frequency bounds
	 */
	const uint16_t MFCC_FREQ_MIN = 3000; //300
	const uint16_t MFCC_FREQ_MAX = 4000;

	/**
	 * Entropy parameters
	 */
	const uint16_t ENTROPY_BINS = 75;
	const double ENTROPY_THRESHOLD = 4.75; //2

	/**
	 * Sound threshold coefficient
	 */
	const uint16_t SOUND_THRESHOLD_COEF = 4;

	/**
	 * DTW distance decision threshold
	 */
	const uint16_t DECISION_DISTANCE_THRESHOLD = 273; // 260

#endif /* CONFIG_H_ */
