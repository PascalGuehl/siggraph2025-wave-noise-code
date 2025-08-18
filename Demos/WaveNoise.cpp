/*
 * Publication: Multi-Dimensional Procedural Wave Noise
 *
 * ACM Transactions on Graphics (TOG), Vol. 44, No. 4, Article 37, August 2025
 * SIGGRAPH 2025, August 2025, Vancouver, Canada
 *
 * Authors: Pascal Guehl (1), Rémi Allègre (2), Guillaume Gilet (3), Basile Sauvage (2), Marie-Paule Cani (1), Jean-Michel Dischler (2)
 *
 * (1) LIX, Ecole Polytechnique, CNRS, Institut Polytechnique de Paris, France
 * (2) ICube, Université de Strasbourg, CNRS, France
 * (3) Université de Sherbrooke, Canada
 * 
 * Website: https://pascalguehl.github.io/siggraph2025-wave-noise/
 */

#include "WaveNoise.h"

/******************************************************************************
 ******************************* INCLUDE SECTION ******************************
 ******************************************************************************/

 // System
#include <climits>
#include <cmath>
#include <ctime>

// STL
#include <iostream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iomanip>
#include <vector>
#include <fstream>
#include <chrono>

// Project
#include "hvPair.h"

/******************************************************************************
 ****************************** NAMESPACE SECTION *****************************
 ******************************************************************************/

// Project
using namespace Wn;

// STL
using namespace std;

/******************************************************************************
 ************************* DEFINE AND CONSTANT SECTION ************************
 ******************************************************************************/

/******************************************************************************
 ***************************** TYPE DEFINITION ********************************
 ******************************************************************************/

/******************************************************************************
 ***************************** METHOD DEFINITION ******************************
 ******************************************************************************/

/******************************************************************************
 * Constructor
 ******************************************************************************/
WaveNoise::WaveNoise()
{
}

/******************************************************************************
 * Destructor
 ******************************************************************************/
WaveNoise::~WaveNoise()
{
	//finalize();
}

/******************************************************************************
 * Initialize
 ******************************************************************************/
void WaveNoise::initialize()
{
	// precision of the pre-computed wave (sampling rate)
	NARRAY = 512;
	MAX_FREQ = NARRAY / 2; // half of half because of FFT symetry and Nyquist
	NDIR = 20;			 // 40+50;

	FREQ_LOW = 1;
	sFREQ_LOW = 1;
	FREQ_HIGH = 32;
	sFREQ_HIGH = 32;
	Ffreq_low = 1.0 / 64.0;
	Ffreq_high = 32.0 / 64.0;

	// Isotropic noise pre-computed arrays
	fs.resize( NARRAY );
	for ( auto& row : fs )
	{
		row.resize( 3, 0.0 );
	}
	vmax = 0.0;
	vmin = 0.0;
	fs_cr.resize( 3 * NARRAY, 0u );
	fsd_cr.resize( 3 * NARRAY, 0u );

	// the user defined discrete spectral energy distribution used for noise control
	spectralEnergyDistribution.resize( NDIR );
	for ( auto& row : spectralEnergyDistribution )
	{
		row.resize( MAX_FREQ, 0.0 );
	}
}

/******************************************************************************
 * Create a (procedural) energy distribution for amplitudes for the isotropic case.
 * Hence, it will use only one amplitudes distribution for all sampled directions.
 ******************************************************************************/
void WaveNoise::createIsotropicProceduralEnergyDistri()
{
	printf( "create isotropic energy distri\n" );

	const int FF = 4;
	
	// Iterate over the frequency domain
	for ( int ii = 0; ii < MAX_FREQ; ii++ )
	{
		double ampli = 0.0;
		// Check current allowed frequency bandwith
		if ( ii >= FREQ_LOW && ii <= FREQ_HIGH )
		{
			double freq = (double)ii / (double)(NARRAY / 2);
			switch ( item_current )
			{
				case 0:
					// ampli decrease as gaussian -> gaussian noise 
					freq = (double)(ii - FREQ_LOW) / (double)(FREQ_HIGH - FREQ_LOW);
					ampli = 1.0 / 16.0 * exp(-freq * freq * 3.0 * 3.0);
					break;

				case 1:
					// ampli constant -> white noise  
					ampli = 1.0 / 8.0;
					break;

				case 2:
					// ampli increases -> blue noise  
					ampli = pow((double)(ii - FREQ_LOW) / (double)(FREQ_HIGH - FREQ_LOW), 2.0) / 16.0;
					break;

				case 3:
					// ampli decrease -> brown noise
					ampli = 1.0 / 16.0 * (1.0 - pow((double)(ii - FREQ_LOW) / (double)(FREQ_HIGH - FREQ_LOW), 0.045));
					break;

				case 10:
					// ampli constant on two-steps -> two-level white noise
					if (ii < 64 / FF)
						ampli = 0.0;
					else if (ii >= 64 / FF && ii < 80 / FF)
						ampli = 0.8 / 16.0;
					else if (ii >= 80 / FF && ii < 256 / FF)
						ampli = 0.275 / 16.0;
					else
						ampli = 0.0;
					break;

				default:
					// ampli decrease as gaussian -> gaussian noise 
					freq = (double)(ii - FREQ_LOW) / (double)(FREQ_HIGH - FREQ_LOW);
					ampli = 1.0 / 16.0 * exp(-freq * freq * 3.0 * 3.0);
					break;
				}
		}

		// Fill only one amplitudes distribution for all sampled directions (slot direction #0)
		spectralEnergyDistribution[ 0 ][ ii ] = ampli;
	}
}

/******************************************************************************
 * ...
 ******************************************************************************/
void WaveNoise::precomputePlanarWave(float scale)
{
	int ii, k;

	int Nfreq = MAX_FREQ;
	int finit = 0;
	vmax = 0.0;
	//double rphases[MAX_FREQ];
	std::vector< double > rphases( MAX_FREQ, 0.0 );
	srand(10);
	for (k = 0; k < MAX_FREQ; k++)
		rphases[k] = 2.0 * (double)rand() / (double)RAND_MAX - 1.0;

	for (ii = 0; ii < NARRAY; ii++)
	{
		fs[ii][0] = 0.0;
		fs[ii][1] = 0.0;
		fs[ii][2] = 0.0;
		for (k = finit; k < finit + Nfreq; k++)
		{
			// random phase
			double phase = 2.0 * M_PI * rphases[k]; // inoise(2 * k, 0, 0);
			double freq = 1.0 / (float)NARRAY;
			// get user defined amplitude
			double ampli = spectralEnergyDistribution[0][k];
			// compute sum real and imaginary parts
			double vcos = ampli * cos(scale * 2.0 * M_PI * (double)ii * freq * (double)k + phase + 2.0 * M_PI / 10.0);
			fs[ii][0] += vcos;
			double vsin = ampli * sin(scale * 2.0 * M_PI * (double)ii * freq * (double)k + phase + 2.0 * M_PI / 10.0);
			fs[ii][1] += vsin;
		}
		if (abs(fs[ii][0]) > vmax)
			vmax = abs(fs[ii][0]);
		if (abs(fs[ii][1]) > vmax)
			vmax = abs(fs[ii][1]);
	}

	// Normalize
	for (ii = 0; ii < NARRAY; ii++)
	{
		fs[ii][0] /= vmax;
		fs_cr[3 * ii + 0] = (unsigned char)(fs[ii][0] * 127.0 + 128.0);
		fs[ii][1] /= vmax;
		fs_cr[3 * ii + 1] = (unsigned char)(fs[ii][1] * 127.0 + 128.0);
		fs_cr[3 * ii + 2] = 0;
	}
	// Gradient
	for (ii = 0; ii < NARRAY; ii++)
	{
		fsd_cr[3 * ii + 0] = (unsigned char)((fs[(ii + NARRAY - 1) % NARRAY][0] - fs[(ii + 1) % NARRAY][0]) * 64.0 + 128.0);
		fsd_cr[3 * ii + 1] = (unsigned char)((fs[(ii + NARRAY - 1) % NARRAY][1] - fs[(ii + 1) % NARRAY][1]) * 64.0 + 128.0);
		fsd_cr[3 * ii + 2] = 0;
	}
}

/******************************************************************************
 * ...
 ******************************************************************************/
double* WaveNoise::MakeSpatialWaveProfile(int pow_2)
{
	//double rphases[MAX_FREQ];
	std::vector< double > rphases( MAX_FREQ, 0.0 );
	srand(10);
	for (int k = 0; k < MAX_FREQ; k++)
		rphases[k] = 2.0 * (double)rand() / (double)RAND_MAX - 1.0;

	int Nfreq = MAX_FREQ;
	int finit = 0;
	int N = 1 << pow_2;
	int i;
	const float scale = 1.0;
	double* A = (double*)malloc(sizeof(double) * N);
	for (i = 0; i < N; i++)
	{
		// DIFFERENT WAVE PROFILES (None Gaussian)
		///////////////////////
		//
		double ff = (double)i / (double)N * 8.0;
		ff = ff - (int)ff;
		int ii = (int)(ff * (double)N);
		double x = (double)i / (double)(N / 2);
		x *= 8.0;
		int y = (int)(x);
		x = x - (int)x;
		x *= 2.0;
		double step = 1.0, aa = 4.0;
		float sumr = 0.0, sumi = 0.0;
		switch (item_current)
		{
		// periodic peaks
		case 4:
			if (x <= 1.0)
				A[i] = 0.8 * (2.0 * pow(x, Power) - 1.0); // Power=50.0
			else
				A[i] = 0.8 * (2.0 * pow(2.0 - x, Power) - 1.0);
			if (y % 2 == 0)
				A[i] = -A[i];
			break;
		// other peaks
		case 5:
			if (ii < 20)
				A[i] = (ii < 10 ? (double)ii / 10.0 : 1.0 - (double)(ii - 10) / 10.0);
			else if (ii < 60)
				A[i] = (ii < 40 ? (double)(ii - 20) / 20.0 : 1.0 - (double)(ii - 40) / 20.0);
			else if (ii < 120)
				A[i] = (ii < 90 ? -(double)(ii - 60) / 30.0 : -1.0 + (double)(ii - 90) / 30.0);
			else if (ii < 170)
				A[i] = (ii < 145 ? (double)(ii - 120) / 25.0 : 1.0 - (double)(ii - 145) / 25.0);
			else if (ii < 250)
				A[i] = (ii < 210 ? (double)(ii - 170) / 40.0 : 1.0 - (double)(ii - 210) / 40.0);
			else if (ii < 300)
				A[i] = (ii < 275 ? -(double)(ii - 250) / 25.0 : -1.0 + (double)(ii - 275) / 25.0);
			else if (ii < 400)
				A[i] = (ii < 350 ? (double)(ii - 300) / 50.0 : 1.0 - (double)(ii - 350) / 50.0);
			else if (ii < 512)
				A[i] = (ii < 456 ? (double)(ii - 400) / 56.0 : 1.0 - (double)(ii - 456) / 56.0);
			else
				A[i] = (ii < 506 ? (double)(ii - 500) / 6.0 : 1.0 - (double)(ii - 506) / 6.0);
			if (A[i] < 0.0)
				A[i] = -2.0 * pow(-A[i], Power); // Power=25.0
			else
				A[i] = 2.0 * pow(A[i], Power);
			break;
		// triangular functions
		case 6:
			A[i] = 0.0;
			for (int k = 0; k < 5; k++)
			{
				double val = 0.0;
				ff = (double)i / (double)N * step;
				ff = ff - (int)ff;
				ii = (int)(ff * (double)N);
				if (ii < 20)
					val = (ii < 10 ? (double)ii / 10.0 : 1.0 - (double)(ii - 10) / 10.0);
				else if (ii < 60)
					val = (ii < 40 ? (double)(ii - 20) / 20.0 : 1.0 - (double)(ii - 40) / 20.0);
				else if (ii < 120)
					val = (ii < 90 ? -(double)(ii - 60) / 30.0 : -1.0 + (double)(ii - 90) / 30.0);
				else if (ii < 170)
					val = (ii < 145 ? (double)(ii - 120) / 25.0 : 1.0 - (double)(ii - 145) / 25.0);
				else if (ii < 250)
					val = (ii < 210 ? (double)(ii - 170) / 40.0 : 1.0 - (double)(ii - 210) / 40.0);
				else if (ii < 300)
					val = (ii < 275 ? -(double)(ii - 250) / 25.0 : -1.0 + (double)(ii - 275) / 25.0);
				else if (ii < 400)
					val = (ii < 350 ? (double)(ii - 300) / 50.0 : 1.0 - (double)(ii - 350) / 50.0);
				else if (ii < 512)
					val = (ii < 456 ? (double)(ii - 400) / 56.0 : 1.0 - (double)(ii - 456) / 56.0);
				else
					val = (ii < 506 ? (double)(ii - 500) / 6.0 : 1.0 - (double)(ii - 506) / 6.0);
				if (val < 0.0)
					val = -pow(-val, Power); // Power=8.0
				else
					val = pow(val, Power);
				A[i] += val * aa;
				aa /= 1.25;
				step *= 1.5;
			}
			break;
		// multi step function
		case 7:
			if (ii < 60)
				A[i] = 0.0;
			else if (ii < 120)
				A[i] = 0.2;
			else if (ii < 170)
				A[i] = -0.2;
			else if (ii < 250)
				A[i] = 0.0;
			else if (ii < 400)
				A[i] = 0.2;
			else
				A[i] = 0.0;
			break;
		// threshold from spectrum
		case 8:
			A[i] = 0.0;
			for (int k = finit; k < finit + Nfreq; k++)
			{
				double phase = 2.0 * M_PI * rphases[k]; // inoise(2 * k, 0, 0);
				double freq = 1.0 / (float)NARRAY;
				double ampli;
				if (item_current == 8)
					ampli = pow((double)(ii - FREQ_LOW) / (double)(FREQ_HIGH - FREQ_LOW), 2.0) / 16.0 / (double)NARRAY;
				double vcos = ampli * cos(2.0 * M_PI * (double)i * freq * (double)k + phase + 2.0 * M_PI / 10.0);
				A[i] += vcos;
			}
			// apply threshold
			if (A[i] > 1.0)
				A[i] = 1.0;
			else if (A[i] < -1.0)
				A[i] = -1.0;
			if (A[i] >= 0.0)
				A[i] = 5.0 * pow(A[i], item_current == 8 ? Power * 0.1 : Power); // Power=2.0
			else
				A[i] = -5.0 * pow(-A[i], item_current == 8 ? Power * 0.1 : Power);
			break;
		// contrast augmented spectrum
		case 9:
			for (int k = FREQ_LOW; k < FREQ_HIGH; k++)
			{
				// random phase
				double phase = 2.0 * M_PI * rphases[k]; // inoise(2 * k, 0, 0);
				double freq = 1.0 / (float)NARRAY;
				// get defined amplitude
				double ampli = 1.0 / 16.0 * (1.0 - pow((double)(k - FREQ_LOW) / (double)(FREQ_HIGH - FREQ_LOW), 0.05));
				// compute sum real and imaginary parts
				double vcos =
					ampli * cos(scale * 2.0 * M_PI * (double)ii * freq * (double)k + phase + 2.0 * M_PI / 10.0);
				double vsin =
					ampli * sin(scale * 2.0 * M_PI * (double)ii * freq * (double)k + phase + 2.0 * M_PI / 10.0);
				sumr += vcos;
				sumi += vsin;
			}
			A[i] = pow(10.0 * (sumr * sumr + sumi * sumi), 0.005);
			if (A[i] > 1.0)
				A[i] = 1.0;
			break;
		}
	}
	if (item_current == 9) // normalize wave
	{
		float vmin = A[0], vmax = A[0];
		for (i = 0; i < N; i++)
		{
			if (vmin > A[i])
				vmin = A[i];
			if (vmax < A[i])
				vmax = A[i];
		}
		for (i = 0; i < N; i++)
		{
			A[i] = -1.0 + 2.0 * (A[i] - vmin) / (vmax - vmin);
		}
	}
	return A;
}

/******************************************************************************
 * ...
 ******************************************************************************/
void WaveNoise::precomputePlanarWaveFromFFT1D(double* A, int N, int pow_2)
{
	int ii, k;
	vmax = 0.0;

	hvArray1<hvPair<double, double>>* F = new hvArray1<hvPair<double, double>>(N, hvPair<double, double>(0.0, 0.0));
	for (ii = 0; ii < N; ii++)
		F->update(ii, hvPair<double, double>(A[ii], 0.0));
	printf("compute FFT of profile: N=%d, pow2=%d...\n", N, pow_2);
	hvArray1<double>::fft(*F, pow_2, 1, 0, false);
	printf("done.\n");
	printf("compute complex wave front...\n");
	int Nfreq = MAX_FREQ;
	int finit = 0;
	for (ii = 0; ii < NARRAY; ii++)
	{
		fs[ii][0] = 0.0;
		fs[ii][1] = 0.0;
		fs[ii][2] = 0.0;
		for (k = finit; k < finit + 2 * Nfreq; k++)
		{
			// random phase
			double phase = F->get(k).phase();
			double freq = 1.0 / (float)(NARRAY);

			// get user defined amplitude
			double ampli = F->get(k).mod();
			// compute sum real and imaginary parts
			double vcos = ampli * cos(2.0 * M_PI * (double)ii * freq * (double)k + phase);
			fs[ii][0] += vcos;
			double vsin = ampli * sin(2.0 * M_PI * (double)ii * freq * (double)k + phase);
			fs[ii][1] += vsin;
		}
		if (abs(fs[ii][0]) > vmax)
			vmax = abs(fs[ii][0]);
		if (abs(fs[ii][1]) > vmax)
			vmax = abs(fs[ii][1]);
	}
	if (vmax < 1.0)
		vmax = 1.0;
	// Normalize
	for (ii = 0; ii < NARRAY; ii++)
	{
		fs[ii][0] /= vmax;
		// printf("%d=%g\n", ii, fs[ii][0]);
		fs_cr[3 * ii + 0] = (unsigned char)(fs[ii][0] * 127.0 + 128.0);
		fs[ii][1] /= vmax;
		fs_cr[3 * ii + 1] = (unsigned char)(fs[ii][1] * 127.0 + 128.0);
		fs_cr[3 * ii + 2] = 0;
	}
	// Gradient
	for (ii = 0; ii < NARRAY; ii++)
	{
		fsd_cr[3 * ii + 0] = (unsigned char)((fs[(ii + NARRAY - 1) % NARRAY][0] - fs[(ii + 1) % NARRAY][0]) * 64.0 + 128.0);
		fsd_cr[3 * ii + 1] = (unsigned char)((fs[(ii + NARRAY - 1) % NARRAY][1] - fs[(ii + 1) % NARRAY][1]) * 64.0 + 128.0);
		fsd_cr[3 * ii + 2] = 0;
	}
}
