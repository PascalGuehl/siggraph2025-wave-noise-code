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

// EasyCppOGL
#include <gl_viewer.h>
#include <mesh.h>
#include <shader_program.h>
#include <texture1d.h>

/******************************************************************************
 ****************************** NAMESPACE SECTION *****************************
 ******************************************************************************/

namespace Wn
{
	class WaveNoise;
}

/******************************************************************************
 ************************* DEFINE AND CONSTANT SECTION ************************
 ******************************************************************************/

/******************************************************************************
 ***************************** TYPE DEFINITION ********************************
 ******************************************************************************/

/******************************************************************************
 ***************************** METHOD DEFINITION ******************************
 ******************************************************************************/

//////////////////////////////////////////////////////////////////////////////////////
//  MAIN  VIEWER
//////////////////////////////////////////////////////////////////////////////////////

// Creation du VIEWER
class Viewer : public EZCOGL::GLViewer
{
	/**************************************************************************
	 ***************************** PUBLIC SECTION *****************************
	 **************************************************************************/

public:

	/****************************** INNER TYPES *******************************/

	/******************************* ATTRIBUTES *******************************/

	/******************************** METHODS *********************************/

	Viewer();
	~Viewer();

	void init_ogl() override;
	void draw_ogl() override;
	void interface_ogl() override;
	void close_ogl() override;

	void initializeNoise();

	bool display1DProfileWidget();
	void displayPeformance();

	void setTitle( const char* pText );

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

	EZCOGL::ShaderProgram::UP prg_p;

	std::vector< EZCOGL::MeshRenderer::UP > renderer_p;
	int nbMeshParts;

	// GPU timer
	GLuint mQueryTimeElapsed;
	GLuint64 mGPUTimeElapsed;

	/**
	 * Wave noise
	 */
	Wn::WaveNoise* waveNoise;
	bool mUseContinuousAnimation;

	float iRatio;

	/******************************** METHODS *********************************/

};
