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

/**
 * @version 1.0
 */

/******************************************************************************
 ******************************* INCLUDE SECTION ******************************
 ******************************************************************************/

// Project
#include "Viewer.h"

// STL
#include <string>

/******************************************************************************
 ****************************** NAMESPACE SECTION *****************************
 ******************************************************************************/

/******************************************************************************
 ************************* DEFINE AND CONSTANT SECTION ************************
 ******************************************************************************/

// NVIDIA Optimus
// A key feature of Optimus configurations is to support rendering applications using NVIDIA High Performance Graphics while displaying on monitors connected to the Integrated Graphics.
// https://docs.nvidia.com/gameworks/content/technologies/desktop/optimus.htm
extern "C"
{
	__declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001;
}

/******************************************************************************
 ***************************** TYPE DEFINITION ********************************
 ******************************************************************************/

/******************************************************************************
 ***************************** METHOD DEFINITION ******************************
 ******************************************************************************/


/******************************************************************************
 * Main entry program
 *
 * @param pArgc Number of arguments
 * @param pArgv List of arguments
 *
 * @return flag telling whether or not it succeeds
 ******************************************************************************/
int main( int, char** )
{
	Viewer v;

	const std::string windowTitle = "Multi-Dimensional Procedural Wave Noise";
	v.setTitle( windowTitle.c_str() );
	v.set_size( 1800, 1000 );

	return v.launch3d();
}
