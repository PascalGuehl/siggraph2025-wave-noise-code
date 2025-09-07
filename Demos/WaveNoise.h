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

/******************************************************************************
 ******************************* INCLUDE SECTION ******************************
 ******************************************************************************/

 // EasyCppOGL
#include <gl_eigen.h>
#include <texture1d.h>

// STL
#include <vector>
#include <string>

/******************************************************************************
 ******************************** CLASS USED **********************************
 ******************************************************************************/

namespace Wn
{

 /**
  * @class WaveNoise
  *
  * @brief The WaveNoise class provides interface to multi-dimensional procedural wave noise.
  *
  * Multi-dimensional procedural wave noise.
  */
class WaveNoise
{
    /**************************************************************************
     ***************************** PUBLIC SECTION *****************************
     **************************************************************************/

public:

    /****************************** INNER TYPES *******************************/

    /**
     * Value type
     */
    enum class EValueType
    {
        eReal = 0,
        eImaginary,
        eModulus,
        eArgument,
        eNbValueTypes
    };
    EValueType mValueType;
    static const std::vector< std::string > mValueTypeNames;

    /**
     * Wave type
     */
    enum class EWaveType
    {
        eNoiseGaussian = 0,
        eNoiseWhite,
        eNoiseBlue,
        eNoiseBrown,
        eNonGaussianCrystal1,
        eNonGaussianWeb,
        eNonGaussianMarble,
        eNonGaussianCrystal2,
        eNonGaussianScratches,
        eNonGaussianSmoothCells,
        eNoiseTwoAmpliLevels,
        eNbWaveTypes
    };
    EWaveType mWaveType;
    static const std::vector< std::string > mWaveTypeNames;

    /******************************* ATTRIBUTES *******************************/

    // precision of the pre-computed wave (sampling rate)
    int NARRAY;
    int NDIR;

    int Oper;
    float Orient, contrast, Anisodd;

    // Isotropic noise pre-computed arrays
    std::vector< std::vector< double > > fs;
    double vmax;
    double vmin;
    std::vector< unsigned char > fs_cr;
    std::vector< unsigned char > fsd_cr;

    /**
     * 1D profile texture
     */
	EZCOGL::Texture1D::SP tex;

    /**
     * 1D profile's gradient texture
     */
    EZCOGL::Texture1D::SP texd;
   
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

    int getNbOrientations() const;
    void setNbOrientations( int pValue );

    int getMinFrequency() const;
    void setMinFrequency( int pValue );

    int getMaxFrequency() const;
    void setMaxFrequency( int pValue );

    int getRecursionLevel() const;
    void setRecursionLevel( int pValue );
    
    float getRecursionProbability() const;
    void setRecursionProbability( float pValue );

    float getVelocity() const;
    void setVelocity( float pValue );
    
    float getTime() const;
    void setTime( float pValue );

    int getMaximumFrequency() const; // half of half because of FFT symetry and Nyquist

    const EZCOGL::GLVec3& getTranslation() const;
    void setTranslation( const EZCOGL::GLVec3& pT );

    float getZoom() const;
    void setZoom( float pValue );

    float getRatio() const;
    void setRatio( float pValue );

    float getPower() const;
    void setPower( float pValue );

    // The user defined discrete spectral energy distribution used for noise control
    const std::vector< std::vector< double > >& getSpectralEnergyDistribution() const;
    std::vector< std::vector< double > >& editSpectralEnergyDistribution();
    void setSpectralEnergyDistribution( const std::vector< std::vector< double > >& pData );

    /**************************************************************************
     **************************** PROTECTED SECTION ***************************
     **************************************************************************/

protected:

    /****************************** INNER TYPES *******************************/

    /******************************* ATTRIBUTES *******************************/

    /**
     * Number of orientations
     */
    int Ndir;

    /**
     * Minimum frequency
     */
    int FREQ_LOW;

    /**
     * Maximum frequency
     */
    int FREQ_HIGH;

    /**
     * Maximum recursion level for cellular noise
     */
    int NRec;

    /**
     * Recursion probability for cellular noise
     */
    float Proba;

    /**
     * Wave maximum velocity
     */
    float tV;

    /**
     * Time
     */
    float Time;

    /**
     * Maximum 1D profile frequency
     */
    int MAX_FREQ; // half of half because of FFT symetry and Nyquist

    /**
     * Translation
     */
    EZCOGL::GLVec3 mTranslation;

    /**
     * Zoom
     */
    float Zoom;

    /**
     * ...
     */
    float Ratio;

    /**
     * [non-gaussian] wave sharpness
     */
    float Power;

    // The user defined discrete spectral energy distribution used for noise control
    /**
     * ...
     */
    std::vector< std::vector< double > > spectralEnergyDistribution;

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
