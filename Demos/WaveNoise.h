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

#ifndef _WN_WAVE_NOISE_H_
#define _WN_WAVE_NOISE_H_

// STL
#include <vector>

/******************************************************************************
 ******************************** CLASS USED **********************************
 ******************************************************************************/

namespace Wn
{

class WaveNoise
{
    /**************************************************************************
     ***************************** PUBLIC SECTION *****************************
     **************************************************************************/

public:

    /****************************** INNER TYPES *******************************/

    /******************************* ATTRIBUTES *******************************/

    // precision of the pre-computed wave (sampling rate)
    int NARRAY;
    int MAX_FREQ; // half of half because of FFT symetry and Nyquist
    int NDIR;			 // 40+50;

    int FREQ_LOW;
    int sFREQ_LOW;
    int FREQ_HIGH;
    int sFREQ_HIGH;
    float Ffreq_low;
    float Ffreq_high;

    int Ndir, Nc, Na, Oper, NRec;
    float tX, tY, tZ, tV;
    float Orient, Period, Zoom, Time, Ratio, iRatio, Power, old_power, contrast, Proba, Anisodd;
    int item_current, old_item;
    int complex_current;

    // Isotropic noise pre-computed arrays
    std::vector< std::vector< double > > fs;
    double vmax;
    double vmin;
    std::vector< unsigned char > fs_cr;
    std::vector< unsigned char > fsd_cr;

    // The user defined discrete spectral energy distribution used for noise control
    std::vector< std::vector< double > > spectralEnergyDistribution;
    
    /******************************** METHODS *********************************/
        
    /**
     * Constructor
     */
    WaveNoise();

    /**
     * Destructor
     */
    ~WaveNoise();

    /**
     * Initialize
     */
    void initialize();

    /**
     * Compute
     */
    int execute();

    // void iso3dangles();
    void createIsotropicProceduralEnergyDistri();
    void precomputePlanarWave( float scale );
    double* MakeSpatialWaveProfile( int pow_2 );
    void precomputePlanarWaveFromFFT1D( double* A, int N, int pow_2 );

    /**************************************************************************
     **************************** PROTECTED SECTION ***************************
     **************************************************************************/

protected:

    /****************************** INNER TYPES *******************************/

    /******************************* ATTRIBUTES *******************************/

    /******************************** METHODS *********************************/

    /**************************************************************************
     ***************************** PRIVATE SECTION ****************************
     **************************************************************************/
  
private:

    /****************************** INNER TYPES *******************************/

    /******************************* ATTRIBUTES *******************************/

    /******************************** METHODS *********************************/

};

}

#endif // _WN_WAVE_NOISE_H_
