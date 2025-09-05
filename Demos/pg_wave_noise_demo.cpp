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
#include <fbo.h>
#include <gl_viewer.h>
#include <mesh.h>
#include <shader_program.h>
#include <texture1d.h>
#include <texture3d.h>

// STL
#include <iostream>
#include <iomanip>

// GLFW
#include <GLFW/glfw3.h>

// Project
#include "WaveNoise.h"

/******************************************************************************
 ****************************** NAMESPACE SECTION *****************************
 ******************************************************************************/

using namespace EZCOGL;

/******************************************************************************
 ************************* DEFINE AND CONSTANT SECTION ************************
 ******************************************************************************/

#define macro_str(s) #s
#define macro_xstr(s) macro_str(s)
#define DATA_PATH std::string(macro_xstr(MYAPP_DATA_PATH))
#define SHADERS_PATH std::string(macro_xstr(MYAPP_SHADERS_PATH))

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

//////////////////////////////////////////////////////////////////////////////////////
//  SHADERS
//////////////////////////////////////////////////////////////////////////////////////

// VERTEX shader
static const std::string p_vert = R"(
////////////////////////////////////////////////////////////////////////////////
// VERSION
////////////////////////////////////////////////////////////////////////////////

#version 460

/******************************************************************************
 * INPUTS
 ******************************************************************************/

// Mesh attributes
layout( location = 1 ) in vec3 position_in;
layout( location = 2 ) in vec3 normal_in;
layout( location = 3 ) in vec3 text_in;
layout( location = 4 ) in vec3 tangents_in;

/******************************************************************************
 * OUTPUTS
 ******************************************************************************/

out vec3 Po;
out vec3 No;
out vec3 Co;
out vec3 NCo, TCo;
out mat3 NoMat;

/******************************************************************************
 * UNIFORMS
 ******************************************************************************/

// Matrix transforms
layout( location = 1 ) uniform mat4 uProjectionMatrix;
layout( location = 2 ) uniform mat4 uViewMatrix;
layout( location = 3 ) uniform mat3 uNormalMatrix;

/******************************************************************************
 * Main Function
 ******************************************************************************/
void main()
{
	// Output data
	Co = vec3( 0.5 ) + position_in * 0.5;
	NCo = normal_in;
	TCo = tangents_in;
	No = uNormalMatrix * normal_in;
	NoMat = uNormalMatrix;
	vec4 Po4 = uViewMatrix * vec4( position_in, 1.0 );
	Po = Po4.xyz;

	// Send position to clip space
	gl_Position = uProjectionMatrix * Po4;
}
)";

// FRAGMENT shader
static const std::string p_frag = R"(
////////////////////////////////////////////////////////////////////////////////
// VERSION
////////////////////////////////////////////////////////////////////////////////

#version 460
precision highp float;

////////////////////////////////////////////////////////////////////////////////
// INPUTS
////////////////////////////////////////////////////////////////////////////////

in vec3 Po;
in vec3 No;
in vec3 Co;
in vec3 NCo, TCo;
in mat3 NoMat;

////////////////////////////////////////////////////////////////////////////////
// OUTPUTS
////////////////////////////////////////////////////////////////////////////////
out vec3 frag_out;

////////////////////////////////////////////////////////////////////////////////
// UNIFORMS
////////////////////////////////////////////////////////////////////////////////

layout( location = 10 ) uniform vec3 color_diff;
layout( location = 11 ) uniform vec3 color_amb;
layout( location = 12 ) uniform vec3 color_spec;
layout( location = 13 ) uniform float Ns;
layout( location = 14 ) uniform vec3 light_pos;
layout( location = 15 ) uniform int NDIR;
layout( location = 16 ) uniform float anisodd;

layout( location = 18 ) uniform int NN;
layout( location = 19 ) uniform float zoom;
layout( location = 20 ) uniform float tV;
layout( location = 21 ) uniform float time;
layout( location = 22 ) uniform float RATIO;
layout( location = 23 ) uniform int QQ;
layout( location = 24 ) uniform float contrast;
layout( location = 25 ) uniform float tX;
layout( location = 26 ) uniform float tY;
layout( location = 27 ) uniform float tZ;

layout( location = 30 ) uniform int Operator;
layout( location = 31 ) uniform int NRec;
layout( location = 32 ) uniform float Proba;

layout( binding = 0 ) uniform sampler1D wave;
layout( binding = 1 ) uniform sampler1D waved;

////////////////////////////////////////////////////////////////////////////////
// LOCAL VARIABLES
////////////////////////////////////////////////////////////////////////////////

const int NDIRMAX = 100;
const float orient = 1.0;
const int NNmin = 0;
const float M_PI = 3.14159265358979;

// PRNG (Pseudo-Random Number Generation)
// - per-thread seed
uint prng_seed = 1;

////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

// PRNG (Pseudo-Random Number Generation)
uint prng_hash( uint x );
uint prng_hash( uvec2 v );
void prng_seeding( uint x, uint y, uint z );
float prng_random();
float prng_next();

// Wave Noise
vec2 wavenoise_3D_isotropic( float x, float y, float z );
vec2 wavenoise_3D_anisotropic( float x, float y, float z, int nd );
vec2 wavenoise_3D_cellular( float x, float y, float z );
vec2 wavenoise_3D_cellular_recursive( float x, float y, float z, uint rec );

/******************************************************************************
 * 1D hash function
 ******************************************************************************/
uint prng_hash( uint x )
{
    x += ( x << 10u );
    x ^= ( x >>  6u );
    x += ( x <<  3u );
    x ^= ( x >> 11u );
    x += ( x << 15u );
    return x;
}

/******************************************************************************
 * 2D hash function
 ******************************************************************************/
uint prng_hash( uvec2 v )
{
	return prng_hash( v.x ^ prng_hash( v.y ) );
}

/******************************************************************************
 * Set the seed for a given integer triplet (x,y,z)
* - all next drawn random numbers will be based on this current seed (for replicability of PRNG)
 ******************************************************************************/
void prng_seeding( uint x, uint y, uint z )
{
	prng_seed = prng_hash( x + prng_hash( y + prng_hash( z ) ) );
}

/******************************************************************************
 * Random value between 0 and 1.0
 ******************************************************************************/
float prng_random()
{
	prng_seed = prng_seed * 22695477u + 1u;
	uint v = prng_seed & 0xffffu;
	return float( v ) / float( 0xffffu );
}

/******************************************************************************
 * Random value between -1.0 and 1.0
 ******************************************************************************/
float prng_next()
{
	return 2.0 * prng_random() - 1.0;
}

)"
R"(

/******************************************************************************
 * Computes a sum of random waves in 3D at spatial location (x,y,z) using a single precomputed 1D table
 * - the sum is isotropic (same energy distribution for all directions)
 ******************************************************************************/
vec2 wavenoise_3D_isotropic( float x, float y, float z )
{
	vec2 sum = vec2( 0.0 ); // result is complex number

	const float poids = 1.0;

	// Iterate through directions (3D stratified sampling)
	const int strat = 4;
	const int NF = NDIR / strat;  // very simple stratified sampling
	int dir = 0;
	for ( int j = 0; j <= NF; j++ ) // stratified alpha
	{
		for ( int i = 0; i < strat; i++ ) // stratified beta
		{
			// Configure per "direction" probabilities
			// strat angular sector
			prng_seeding( uint( dir ) * 5u, 2u, 3u );
			// - configure animation
			float tspeed = tV * sign( prng_next() ) * time;
			// - configure orientation stratified sampling
			float anbeta = acos( ( float( i ) + prng_random() ) / float( strat ) ); // in [0;Pi/2]
			float alpha1 = 2.0 * M_PI * float( j ) / float( NF + 1 );
			float alpha2 = 2.0 * M_PI * ( float( j ) + 1.0 ) / float( NF + 1 );
			float analpha = ( alpha1 + alpha2 ) / 2.0;
			float analphawidth = alpha2 - alpha1;

			// Slicing
			// - draw a random orientation (stratified sampling)
			float aa =  M_PI / 4.0 * orient + analpha + prng_next() * analphawidth * 0.5;
			float beta = anbeta + M_PI / 4.0 * orient;
			// - project the (x,y,z) position on the 1D line
			float id = sin( beta ) * cos( aa ) * x + sin( beta ) * sin( aa ) * y + cos( beta ) * z;
			// - interger part (slice id)
			int iid = ( id > 0.0 ? int( id ) : int( id ) - 1 ); //FASTFLOOR( id );
			// - fractional part (position inside the slice)
			float oid = id - float( iid );

			// Configure per "slice id" along 1D line probabilities (according to direction)
			prng_seeding( uint( 7 * dir ), uint( 4 * iid ), 4u );
			// - jittering along the 1D line
			float pl = 0.5 + 0.2 * prng_next();
			// - blending weight between current slice and next slice waves
			if ( oid < pl - 0.3 ) oid = 0.0;
			else if ( oid > pl + 0.3 ) oid = 1.0;
			else oid = ( oid - pl + 0.3 ) / 0.6;

			// Compute CURRENT (slice) wave (along 1D line, according to direction)
			// - configure per "wave" probabilities
			prng_seeding( uint( 4 * dir ), uint( 4 * iid ), 0 );
			// - draw a random orientation
			float alpha = ( analpha + prng_next() * analphawidth ) * orient + aa * ( 1.0 - orient );
			float bb = acos( ( float( i ) + ( 0.5 * prng_next() + 0.5 ) ) / float( strat ) ) * orient + beta * ( 1.0 - orient );
			// - rotation, scale, then translation of the wave
			float dd = ( 1.0 / RATIO * ( sin( bb ) * cos( alpha ) * x + sin( bb ) * sin( alpha ) * y + cos( bb ) * z ) + 2.0 * prng_next() - tspeed );
			// - add wave contribution
			if ( dir >= NNmin && dir < NN && dir < NDIR ) {
				sum += poids * ( 1.0 - oid ) * ( 2.0 * texture( wave, fract( dd ) ).xy - vec2( 1.0 ) );
			}

			// Compute NEXT (slice) wave (along 1D line, according to direction)
			// - configure per "wave" probabilities
			prng_seeding( uint( 4 * dir ), uint( 4 * ( iid + 1 ) ), 0 );
			// - draw a random orientation
			alpha = ( analpha + prng_next() * analphawidth ) * orient + aa * ( 1.0 - orient );
			bb = acos( ( float( i ) + ( 0.5 * prng_next() + 0.5 ) ) / float( strat ) ) * orient + beta * ( 1.0 - orient );
			// - rotation, scale, then random translation of the wave
			dd = ( 1.0 / RATIO * ( sin( bb ) * cos( alpha ) * x + sin( bb ) * sin( alpha ) * y + cos( bb ) * z ) + 2.0 * prng_next() - tspeed );
			// - add wave contribution
			if ( dir >= NNmin && dir < NN && dir < NDIR ) {
				sum += poids * oid * ( 2.0 * texture( wave, fract( dd ) ).xy - vec2( 1.0 ) );
			}

			// Update next direction ID
			dir += 1;
		}
	}

	// Result (with normalization)
    return sum / vec2( 2.0 + 0.2 * float( NN - NNmin ) );
}

)"
R"(

/******************************************************************************
 * Computes a sum of random waves in 3D at spatial location (x,y,z) using a single precomputed 1D table
 * - the sum is anisotropic (different energy distribution for some directions)
 ******************************************************************************/
vec2 wavenoise_3D_anisotropic( float x, float y, float z, int nd )
{
	vec2 sum = vec2( 0.0 ); // result is complex number

	// Iterate through directions (3D stratified sampling)
	const int strat = 2;
	const int NF = NDIR / strat / nd;  // very simple stratified sampling
	int dir = 0;
	for ( int j = 0; j < nd * NF; j++ ) // stratified alpha
	{
		for ( int i = 0; i < strat; i++ ) // stratified beta
		{
			// Configure per "direction" probabilities
			prng_seeding( uint( dir ) * 5u, 2u, 3u );
			// - configure animation
			float tspeed = tV * sign( prng_next() ) * time;
			// - configure orientation stratified sampling
			float anbeta = acos( 0.25 + anisodd * 0.5 * ( float( i ) + ( 0.5 * prng_next() + 0.5 ) ) / float( strat ) );
			float alpha1 = ( j < NF ? M_PI / 4.0: 3.0 * M_PI / 4.0 ) + anisodd * M_PI / 8.0 * float( j ) / float( NF + 1 );
			float alpha2 = ( j < NF ? M_PI / 4.0: 3.0 * M_PI / 4.0 ) + anisodd * M_PI / 8.0 * ( float( j ) + 1.0) / float( NF + 1 );
			float analpha = ( alpha1 + alpha2 ) / 2.0;
			float analphawidth = alpha2 - alpha1;

			// Slicing
			// - draw a random orientation (stratified sampling)
			float aa =  M_PI / 4.0 + analpha + prng_next() * analphawidth * 0.5;
			float beta = anbeta + M_PI / 4.0; 
			// - project the (x,y,z) position on the 1D line
			float id = sin( beta ) * cos( aa ) * x + sin( beta ) * sin( aa ) * y + cos( beta ) * z;
			// - interger part (slice id)
			int iid = ( id > 0.0 ? int( id ) : int( id ) - 1 ); //FASTFLOOR( id );
			// - fractional part (position inside the slice)
			float oid = id - float( iid );

			// Configure per "slice id" along 1D line probabilities (according to direction)
			prng_seeding( uint( 7 * dir ), uint( 4 * iid ), 4u );
			// - jittering along the 1D line
			float pl = 0.5 + 0.2 * prng_next();
			// - blending weight between current slice and next slice waves
			if ( oid < pl - 0.3 ) oid = 0.0;
			else if ( oid > pl + 0.3 ) oid = 1.0;
			else oid = ( oid - pl + 0.3 ) / 0.6;

			// Compute CURRENT (slice) wave (along 1D line, according to direction)
			// - configure per "wave" probabilities
			prng_seeding( uint( 4 * dir ), uint( 4 * iid ), 0 );
			// - draw a random orientation
			float alpha = ( analpha + prng_next() * analphawidth );
			float bb = acos( 0.5 + 0.1 * ( float( i ) + ( 0.5 * prng_next() + 0.5 ) ) / float( strat ) );
			// - rotation, scale, then translation of the wave
			float dd = ( 1.0 / RATIO * ( sin( bb ) * cos( alpha ) * x + sin( bb ) * sin( alpha ) * y + cos( bb ) * z ) + 2.0 * prng_next() - tspeed );
			// - add wave contribution
			if ( dir >= NNmin && dir < NN && dir < NDIR )
			{
				sum += ( 1.0 - oid ) * ( 2.0 * texture( wave, fract( dd ) ).xy - vec2( 1.0 ) );
			}

			// Compute NEXT (slice) wave (along 1D line, according to direction)
			// - configure per "wave" probabilities
			prng_seeding( uint( 4 * dir ), uint( 4 * ( iid + 1 ) ), 0 );
			// - draw a random orientation
			alpha = ( analpha + prng_next() * analphawidth );
			bb = acos( 0.5 + 0.1 * ( float( i ) + ( 0.5 * prng_next() + 0.5 ) ) / float( strat ) );
			// - rotation, scale, then random translation of the wave
			dd = ( 1.0 / RATIO * ( sin( bb ) * cos( alpha ) * x + sin( bb ) * sin( alpha ) * y + cos( bb ) * z ) + 2.0 * prng_next() - tspeed );
			// - add wave contribution
			if ( dir >= NNmin && dir < NN && dir < NDIR )
			{
				sum += oid * ( 2.0 * texture( wave, fract( dd ) ).xy - vec2( 1.0 ) );
			}

			// Update next direction ID
			dir += 1;
		}
	}

	// Result (with normalization)
	return sum / vec2( 2.0 + 0.2 * float( NN - NNmin ) );
}

/******************************************************************************
 * Computes cellular wave noise at a given recursion level
 ******************************************************************************/
vec2 wavenoise_3D_cellular_recursive( float x, float y, float z, uint rec )
{
	uint resv = 0;
	float rvalue = 2.0;

	// Iterate through directions (3D stratified sampling)
	const int strat = 2;
	const int NF = NDIR / strat;  // very simple stratified sampling
	int dir = 0;
	for ( int j = 0; j <= NF; j++ ) // stratified alpha
	{
		for ( int i = 0; i < strat; i++ ) // stratified beta
		{
			// Configure per "direction" probabilities
			prng_seeding( uint( dir ) * 5u, 2u, 3u + rec );
			// - configure animation
			float tspeed = tV * sign( prng_next() ) * time;
			// - configure orientation stratified sampling
			float anbeta = acos( ( float( i ) + ( 0.5 * prng_next() + 0.5 ) ) / float( strat ) );
			float alpha1 = M_PI * float( j ) / float( NF + 1 );
			float alpha2 = M_PI * ( float( j ) + 1.0 ) / float( NF + 1 );
			float analpha = ( alpha1 + alpha2 ) / 2.0;
			float analphawidth = alpha2 - alpha1;

			// Slicing
			// - draw a random orientation (stratified sampling)
			float aa =  M_PI / 4.0 * orient + analpha + prng_next() * analphawidth * 0.5;
			float beta = anbeta + M_PI / 4.0 * orient;
			// - project the (x,y,z) position on the 1D line
			float id = sin( beta ) * cos( aa ) * x + sin( beta ) * sin( aa ) * y + cos( beta ) * z  + 5.0 + tspeed;
			// - interger part (slice id)
			int iid = ( id > 0.0 ? int( id ) : int( id ) - 1 ); //FASTFLOOR( id );
			// - fractional part (position inside the slice)
			float oid = id - float( iid );

			// Configure per "slice id" along 1D line probabilities (according to direction)
			prng_seeding( uint( 7 * dir ), uint( 4 * iid ), 4u );
			// - jittering along the 1D line
			float pl = 0.5 + 0.3 * prng_next();
			// - determine closest cell border between current, previous and next cells, along the 1D line (dist min)
			if ( oid <= pl )
			{
				// current cell
				float dist1 = pl - oid;
				rvalue = min( dist1, rvalue );
				// previous cell
				prng_seeding( uint( 7 * dir ), uint( 4 * ( iid - 1 ) ), 4u );
				float pl2 = 0.5 + 0.3 * prng_next();
				float dist2 = oid + 1.0 - pl2;
				rvalue = min( dist2, rvalue );
			}
			else
			{
				// current cell
				float dist1 = oid - pl;
				rvalue = min( dist1, rvalue );
				// next cell
				prng_seeding( uint( 7 * dir ), uint( 4 * ( iid + 1 ) ), 4u );
				float pl2 = 0.5 + 0.3 * prng_next();
				float dist2 = 1.0 - oid + pl2;
				rvalue = min( dist2, rvalue );
			}
			// cell index
			int iidabs = abs( iid );
			if ( ( oid <= pl && iidabs % 2 == 0 ) || ( oid > pl && iidabs % 2 == 1 ) ) resv += 1;
			resv = resv * 2u;

			// Orthogonal tesselation (to current 1D line) to create regular cellular patterns
			if ( oid <= pl ) prng_seeding( uint( 4 * dir ), uint( 4 * iid ), 0 );
			else prng_seeding( uint( 4 * dir ), uint( 4 * ( iid + 1 ) ), 0 );
			// - draw a random orthogonal orientation (stratified sampling in a cone)
			float alpha = ( analpha + prng_next() * analphawidth ) * orient + aa * ( 1.0 - orient );
			float bb = acos( ( float( i ) + ( 0.5 * prng_next() + 0.5 ) ) / float( strat ) ) * orient + beta * ( 1.0 - orient );
			// - project the (x,y) position on the 1D line
			float dd = sin( bb ) * cos( alpha ) * x + sin( bb ) * sin( alpha ) * y + cos( bb ) * z + 5.0 + tspeed;
			// - interger part (slice id)
			int p = ( dd > 0.0 ? int( dd ) : int( dd ) - 1 ); //FASTFLOOR( dd );
			// - fractional part (position inside the slice)
			dd = dd - float( p);
			prng_seeding( uint( 7 * dir ), uint( 4 * p ), 4u );
			// - jittering along the 1D line
			pl = 0.5 + 0.3 * prng_next();
			// - determine in which we are, current or previous cell, along the 1D line (dist min)
			if ( dd <= pl )
			{
				// current cell
				float dist1 = pl - dd;
				rvalue = min( dist1, rvalue );
				// previous cell
				prng_seeding( uint( 7 * dir ), uint( 4 * ( p - 1 ) ), 4u );
				float pl2 = 0.5 + 0.3 * prng_next();
				float dist2 = dd + 1.0 - pl2;
				rvalue = min( dist2, rvalue );
			}
			else
			{
				// current cell
				float dist1 = dd - pl;
				rvalue = min( dist1, rvalue );
				// next cell
				prng_seeding( uint( 7 * dir ), uint( 4 * ( p + 1 ) ), 4u );
				float pl2 = 0.5 + 0.3 * prng_next();
				float dist2 = 1.0 - dd + pl2;
				rvalue = min( dist2, rvalue );
			}
			// cell index
			iidabs = abs( p );
			if ( ( dd <= pl && iidabs % 2 == 0 ) || ( dd > pl && iidabs % 2 == 1 ) ) resv += 1;
			resv = resv * 2u;

			// Update next direction ID
			dir += 1;
		}
	}

	// Result
	uint cellident = ( resv * 1453u ) % 255u;
	return vec2( float( cellident ) / 255.0 * 2.0 - 1.0, rvalue * 20.0 - 1.0 );
}

/******************************************************************************
 * Computes cellular wave noise
 ******************************************************************************/
vec2 wavenoise_3D_cellular( float x, float y, float z )
{
	// Computer cellular noise
	vec2 res = wavenoise_3D_cellular_recursive( x, y, z, 0 );

	// Recursive call (simulate STIT patterns)
	bool cont = true;
	int rr = 0;
	for ( rr = 1; rr < NRec && cont; rr++ )
	{
		uint cellident = uint( 255.0 * ( res.x + 1.0 ) * 0.5 );
		prng_seeding( cellident, uint( rr ), 3 );
		if ( 0.5 * ( prng_next() + 1.0 ) < Proba )
		{
			vec2 res2 = wavenoise_3D_cellular_recursive( x, y, z, uint( rr ) );

			res = vec2( res2.x, min( res.y, res2.y ) );
		}
		else cont = false;
	}

	return res;
}

/******************************************************************************
 * Computes multi-dimensional procedural wave noise
 ******************************************************************************/
vec2 wavenoise_3D( vec3 p )
{
	vec2 noise = vec2( 0.0 );
    
	if ( Operator == 0 )
	{
		noise = wavenoise_3D_isotropic( p.x, p.y, p.z );
	}
	else if ( Operator == 1 )
	{
		noise = wavenoise_3D_anisotropic( p.x, p.y, p.z, 1 );
	}
	else if ( Operator == 2 )
	{
		noise = wavenoise_3D_anisotropic( p.x, p.y, p.z, 2 );
	}
	else
	{
		noise = wavenoise_3D_cellular( 0.1 * p.x, 0.1 * p.y, 0.1 * p.z );
	}

	return noise;
}

////////////////////////////////////////////////////////////////////////////////
// MAIN PROGRAM
////////////////////////////////////////////////////////////////////////////////
void main()
{
	// Noise computation
	vec3 Xpos = vec3( zoom * Co.x + tX, zoom * Co.y + tY, zoom * Co.z + tZ );
	vec2 noise = wavenoise_3D( Xpos );

	// Noise post-processing
	float val = 0.0;
	if ( Operator <= 2 ) { 
		if ( QQ == 0 ) val = clamp( ( contrast * noise.x + 0.5 ), 0.0, 1.0 );
		else if ( QQ == 1 ) val = clamp( ( contrast * noise.y + 0.5 ), 0.0, 1.0 );
		else if ( QQ == 2 ) val = clamp( length( noise ) * 2.0 * contrast, 0.0, 1.0 );
		else val = atan( noise.y, noise.x ) / M_PI + 0.5; //smoothstep( -1.0, 1.0, cos( atan( noise.y, noise.x ) ) ); 
	}
	else {
		if ( Operator == 3 ) val = clamp( ( 0.5 * noise.x + 0.5 ), 0.0, 1.0 );
		else if ( Operator == 4 ) val = 0.8 * pow( clamp( ( 0.5 * noise.y + 0.5 ), 0.0, 1.0 ), contrast ) + 0.2;
		else if ( Operator == 5 ) val = step( 0.05, 0.5 * noise.y + 0.5 );
		else val = 1.0 - pow( clamp( ( 0.5 * noise.y + 0.5 ), 0.0, 1.0 ), contrast );
	}

	// Shading (ADS: ambient, diffuse, specular)
	vec3 Nco = normalize( NCo );
	vec3 Nnn = NCo;
	if ( gl_FrontFacing == false ) Nnn = -Nnn;
	vec3 L = normalize( light_pos - Po );
	//vec3 Npp = normalize( TCo );
	//vec3 Ncc = cross( Nnn, Npp );
	vec3 N = normalize( NoMat * normalize( Nnn ) );
	float lamb = abs( dot( N, L ) );
	vec3 E = normalize( -Po );
	vec3 R = reflect( -L, N );
	float spec = Ns == 0.0 ? 0.0 : pow( max( dot( R, E ), 0.0 ), Ns );

	// Write output color
	frag_out = min( color_amb * val + color_diff * lamb * val + color_spec * spec, vec3( 1.0 ) );
}
)";

//////////////////////////////////////////////////////////////////////////////////////
//  MAIN  VIEWER
//////////////////////////////////////////////////////////////////////////////////////

// Creation du VIEWER
class Viewer : public GLViewer
{
	ShaderProgram::UP prg_p;
	
	std::vector<MeshRenderer::UP> renderer_p;
	int nbMeshParts;

	Texture1D::SP tex, texd;
	
	// GPU timer
	GLuint mQueryTimeElapsed;
	GLuint64 mGPUTimeElapsed;

	Wn::WaveNoise* waveNoise;
	bool mUseContinuousAnimation;

public:

	Viewer();
	~Viewer();

	void init_ogl() override;
	void draw_ogl() override;
	void interface_ogl() override;

	bool display1DProfileWidget();
	void displayPeformance();
};

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

	GLFWwindow* window = v.window();
	glfwSetWindowTitle( window, "Multi-Dimensional Procedural Wave Noise" );
	v.set_size(1800, 1000);

	return v.launch3d();
}

/******************************************************************************
 * Constructor
 ******************************************************************************/
Viewer::Viewer()
: nbMeshParts(0)
, waveNoise( nullptr )
, mQueryTimeElapsed( 0 )
, mGPUTimeElapsed( 0 )
, mUseContinuousAnimation( false )
{
	// Create and initialize noise
	waveNoise = new Wn::WaveNoise();
	waveNoise->initialize();

	// Dafault parameters value
	waveNoise->tX = 0.0;
	waveNoise->tY = 0.0;
	waveNoise->tZ = 0.0;
	waveNoise->setVelocity( 0.0 );
	waveNoise->Ndir = 40;
	waveNoise->Orient = 1.0f;
	waveNoise->Period = 0.0;
	waveNoise->Zoom = 0.2f;
	waveNoise->Time = 0.0f;
	waveNoise->item_current = 0; // Gaussian
	waveNoise->Oper = 0;		 // Isowave
	waveNoise->old_item = 0;
	waveNoise->Ratio = 64.0f;
	waveNoise->iRatio = 6;
	waveNoise->complex_current = 0;
	waveNoise->Power = 25.0;
	waveNoise->old_power = 1.0;
	waveNoise->contrast = 0.5;
	waveNoise->setRecursionProbability( 0.5 );
	waveNoise->setRecursionLevel( 3 );
	waveNoise->Anisodd = 0.5;
	
	// Precompute 1D profile wave(s)
	// iso3dangles();
	waveNoise->createIsotropicProceduralEnergyDistri();
	waveNoise->precomputePlanarWave( 4.0 );
}

/******************************************************************************
 * Destructor
 ******************************************************************************/
Viewer::~Viewer()
{
	// Device timer
	glDeleteQueries( 1, &mQueryTimeElapsed );
}

/******************************************************************************
 * Callback for custom initialization function
 ******************************************************************************/
void Viewer::init_ogl()
{
	prg_p = ShaderProgram::create({ {GL_VERTEX_SHADER, p_vert}, {GL_FRAGMENT_SHADER, p_frag} }, "prog");

	// Load OBJ file mesh
	auto mesh = Mesh::load(DATA_PATH + "/models/manual_cutaway_cube.obj");

	nbMeshParts = mesh.size();
	// set the renderer and the materials for all the meshes parts
	for (int i = 0; i < nbMeshParts; ++i)
	{
		renderer_p.push_back(mesh[i]->renderer(1, 2, 3, 4, -1));
	}
	set_scene_center(mesh[0]->BB()->center());
	set_scene_radius(3.0 * mesh[0]->BB()->radius());
	// set_scene_radius(50.f);

	tex = Texture1D::create();
	// tex->update(0, N, wprofile);
	tex->alloc( waveNoise->NARRAY, GL_RGB8, waveNoise->fs_cr.data() );

	texd = Texture1D::create();
	// tex->update(0, N, wprofile);
	texd->alloc( waveNoise->NARRAY, GL_RGB8, waveNoise->fsd_cr.data() );

	// Device timer
	glCreateQueries( GL_TIME_ELAPSED, 1, &mQueryTimeElapsed );
}

/******************************************************************************
 * Callback for custom rendering/compute function
 ******************************************************************************/
void Viewer::draw_ogl()
{
	GLMat4 sc = Transfo::scale(2.5);
	GLMat4 rotx = Transfo::rotateX(-60.0);
	GLMat4 rotz = Transfo::rotateZ(25.0);
	const GLMat4& proj = this->get_projection_matrix();
	const GLMat4& mv = this->get_view_matrix() * sc * rotx * rotz;

	glEnable(GL_DEPTH_TEST);
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	waveNoise->FREQ_LOW = (int)(waveNoise->Ffreq_low * 64.0);
	waveNoise->FREQ_HIGH = (int)(waveNoise->Ffreq_high * 64.0);
	prg_p->bind();

	// Bind textures
	// - precomputed 1D wave profile
	tex->bind( 0 );
	// - precomputed 1D wave gradient
	texd->bind( 1 );

	if (waveNoise->FREQ_LOW != waveNoise->sFREQ_LOW || waveNoise->FREQ_HIGH != waveNoise->sFREQ_HIGH || waveNoise->item_current != waveNoise->old_item ||
		(waveNoise->item_current >= 4 && waveNoise->old_power != waveNoise->Power))
	{
		waveNoise->sFREQ_LOW = waveNoise->FREQ_LOW;
		waveNoise->sFREQ_HIGH = waveNoise->FREQ_HIGH;
		waveNoise->old_item = waveNoise->item_current;
		waveNoise->old_power = waveNoise->Power;
		if (waveNoise->item_current >= 4 && waveNoise->item_current <= 9)
		{
			int pow_2 = 1, N = 1;
			while (N < waveNoise->NARRAY)
			{
				pow_2++;
				N *= 2;
			}
			pow_2--;
			printf("NARRAY=%d, N=%d, pow2=%d\n", waveNoise->NARRAY, N, pow_2);
			double* A = waveNoise->MakeSpatialWaveProfile( pow_2 );
			waveNoise->precomputePlanarWaveFromFFT1D( A, N, pow_2 );
			delete A;
		}
		else
		{
			waveNoise->createIsotropicProceduralEnergyDistri();
			waveNoise->precomputePlanarWave( 4.0 );
		}
		tex->alloc(waveNoise->NARRAY, GL_RGB8, waveNoise->fs_cr.data());
		texd->alloc(waveNoise->NARRAY, GL_RGB8, waveNoise->fsd_cr.data());
	}
	set_uniform_value(1, proj);
	set_uniform_value(2, mv);
	set_uniform_value(3, Transfo::inverse_transpose(mv));

	set_uniform_value(14, GLVec3(0, 2.0, 3.0));

	set_uniform_value(25, waveNoise->tX+5.0);
	set_uniform_value(26, waveNoise->tY+5.0);
	set_uniform_value(27, waveNoise->tZ+5.0);
	set_uniform_value(15, waveNoise->Ndir);
	set_uniform_value(16, waveNoise->Anisodd);
	waveNoise->Period = 1.0 - waveNoise->Orient;
	// set_uniform_value(17, Period);
	set_uniform_value(18, waveNoise->Ndir);
	set_uniform_value( 20, waveNoise->getVelocity() );
	if ( ! mUseContinuousAnimation )
	{
		set_uniform_value( 21, waveNoise->Time );
	}
	else
	{
		set_uniform_value( 21, current_time() );
	}
	waveNoise->Ratio = (float)pow(2.0, waveNoise->iRatio );
	set_uniform_value( 22, waveNoise->Ratio );
	set_uniform_value( 19, waveNoise->Zoom * waveNoise->Ratio );
	set_uniform_value( 23, waveNoise->complex_current );
	set_uniform_value( 24, waveNoise->contrast );
	set_uniform_value( 30, waveNoise->Oper );
	set_uniform_value( 31, waveNoise->getRecursionLevel() );
	set_uniform_value( 32, waveNoise->getRecursionProbability() );

	// GPU timer
	GLuint64 result = 0;
	mGPUTimeElapsed = 0;

	for (int i = 0; i < renderer_p.size(); i++)
	{
		auto mat = renderer_p[i]->material();

		set_uniform_value(10, mat->Kd);
		set_uniform_value(11, mat->Ka);
		set_uniform_value(12, mat->Ks);
		set_uniform_value(13, mat->Ns);
		set_uniform_value("light_pos", GLVec3(0, 2, 2.5));

		glBeginQuery( GL_TIME_ELAPSED, mQueryTimeElapsed );

		renderer_p[ i ]->draw( GL_TRIANGLES );

		// GPU timer
		glEndQuery( GL_TIME_ELAPSED );
		glGetQueryObjectui64v( mQueryTimeElapsed, GL_QUERY_RESULT, &result );
		mGPUTimeElapsed += result;
	}
}

/******************************************************************************
 * Callback for custom GUI (graphical user interface)
 ******************************************************************************/
void Viewer::interface_ogl()
{
		ImGui::Begin("Wave Noise Gui", nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize(ImVec2(560, 525));
		ImGui::SetWindowPos(ImVec2(10, 10));
		// ImGui::Text("FPS: %2.2lf", fps_);

		const float sliderWidth = 100.0f; // largeur personnalis�e pour chaque slider
		if (ImGui::CollapsingHeader("Spatial Transformation"))
		{
			ImGui::PushItemWidth(sliderWidth);
			ImGui::SliderFloat("X", &waveNoise->tX, 0.0, 10.0);
			ImGui::SameLine();
			ImGui::SliderFloat("Y", &waveNoise->tY, 0.0, 10.0);
			ImGui::SameLine();
			ImGui::SliderFloat("Z", &waveNoise->tZ, 0.0, 10.0);
			ImGui::PopItemWidth(); // remet la largeur par d�faut

			ImGui::SliderFloat("Zoom", &waveNoise->Zoom, 0.01, 2.0);
		}
	
		// ImGui::Text("MEM:  %2.2lf %%", 100.0 * mem_);

		if (ImGui::CollapsingHeader("Spectral Configuration")) // PG
		{
			//ImGui::Text("Quality vs.Framerate");
			if (ImGui::TreeNode("Quality vs. Framerate"))
			{
				ImGui::SliderInt("NDir", &waveNoise->Ndir, 1, 100);
				ImGui::SliderFloat("Slice size", &waveNoise->iRatio, 0.0, 8.0);

				ImGui::TreePop();
			}

			//ImGui::Text("Frequency Ranges");
			if (ImGui::TreeNode("Frequency Ranges"))
			{
				ImGui::SliderFloat("aniso wave dir width", &waveNoise->Anisodd, 0.0, 1.0);

				ImGui::PushItemWidth(sliderWidth);
				ImGui::SliderFloat("FreqMin", &waveNoise->Ffreq_low, 1.0 / 64.0, 1.0);
				ImGui::SameLine();
				ImGui::SliderFloat("FreqMax", &waveNoise->Ffreq_high, 1.0 / 64.0, 1.0);
				ImGui::PopItemWidth(); // remet la largeur par d�faut

				ImGui::TreePop();
			}
		}

		if (ImGui::CollapsingHeader("Waveforms & Noise Types")) // PG
		{
			const char* items[] = { "[noise] gaussian",		"[noise] white (0dB)",		"[noise] blue (+3dB)",
								   "[noise] brown (-6dB)",	"[non-gaussian] crystal1",		"[non-gaussian] web",
								   "[non-gaussian] marble",	"[non-gaussian] crystal2",	"[non-gaussian] scratches",
								   "[non-gaussian] smooth cells", "[noise] two ampli levels" };
			ImGui::Combo("Wave type", &waveNoise->item_current, items, IM_ARRAYSIZE(items));

			const char* itemsop[] = {
				"[isotropic] Sum", "[anisotropic] Sum - one direction", "[anisotropic] Sum - two directions", "[cellular] Random polytopes", "[cellular] Cellular",
				"[cellular] Hyperplan",	 "[cellular] Reversed Cellular" };
			ImGui::Combo("Operator", &waveNoise->Oper, itemsop, IM_ARRAYSIZE(itemsop));
			ImGui::SliderFloat("nongauss wave sharpness", &waveNoise->Power, 0.2, 50.0);
		}

		if (ImGui::CollapsingHeader("Cellular Noise Settings (STIT)")) // PG
		{
			ImGui::PushItemWidth(sliderWidth);
			float recursionProbability = waveNoise->getRecursionProbability();
			if ( ImGui::SliderFloat( "Subdivision Probability", &recursionProbability, 0.0, 1.0 ) )
			{
				waveNoise->setRecursionProbability( recursionProbability );
			}
			ImGui::SameLine();
			int recursionLevel = waveNoise->getRecursionLevel();
			if ( ImGui::SliderInt( "Nb Recursions", &recursionLevel, 1, 5 ) )
			{
				waveNoise->setRecursionLevel( recursionLevel );
			}
			ImGui::PopItemWidth(); // remet la largeur par d�faut
		}

		if (ImGui::CollapsingHeader("Post-Processing & Display")) // PG
		{
			const char* itemscplx[] = { "real", "imaginary", "modulus", "phasor" };
			ImGui::Combo("Value", &waveNoise->complex_current, itemscplx, IM_ARRAYSIZE(itemscplx));
			ImGui::SliderFloat("Contrast", &waveNoise->contrast, 0.1, 1.0);
		}

		if ( ImGui::CollapsingHeader( "Temporal & Animation Controls" ) )
		{
			ImGui::PushItemWidth( sliderWidth );
			float velocity = waveNoise->getVelocity();
			if ( ImGui::SliderFloat( "Speed", &velocity, 0.0, 0.1 ) )
			{
				waveNoise->setVelocity( velocity );
			}
			ImGui::SameLine();
			ImGui::SliderFloat( "Time", &waveNoise->Time, 0.0, 5.0 );
			ImGui::SameLine();
			ImGui::Checkbox( "Continuous", &mUseContinuousAnimation );
			ImGui::PopItemWidth(); // remet la largeur par d�faut
		}

		// ImGui::SetWindowSize({ 0, 0 });
		ImGui::End();

	// Display wave noise's' 1D profile
	display1DProfileWidget();
	displayPeformance();
}

/******************************************************************************
 * ...
 ******************************************************************************/
bool Viewer::display1DProfileWidget()
{
	bool updateRequested = false;

	// Force UI on top right
	ImGuiIO& io = ImGui::GetIO();
	const float PAD = 10.0f; // petite marge du bord
	ImGui::SetNextWindowPos(ImVec2( io.DisplaySize.x - PAD, PAD ), ImGuiCond_Always, ImVec2( 1.0f, 0.0f ) );
	ImGui::SetNextWindowBgAlpha( 0.9f );
	ImGuiWindowFlags flags = ImGuiWindowFlags_NoMove
		| ImGuiWindowFlags_NoResize
		| ImGuiWindowFlags_NoCollapse
		| ImGuiWindowFlags_NoSavedSettings;

	ImGui::Begin( "PLANE WAVE##hpn", nullptr, flags );
	{
		ImGui::TextColored( ImVec4( 0.f, 1.f, 0.f, 1.f ), "Amplitude (energy distribution)" );

		const char* items[] = { "noise-gaussian",	"noise-white",	 "noise-blue",	  "noise-brown",
		"nongauss-crystal1", "nongauss-web", "nongauss-marble", "nongauss-crystal2", "nongauss-scratches", "nongauss-smooth cells", "noise-two ampli levels" };
		const std::string wave_type = "[" + std::string( items[ waveNoise->item_current ] ) + "]";
		ImGui::TextColored( ImVec4( 0.f, 1.f, 0.f, 1.f ), wave_type.c_str() );

		// Display amplitude's energy distribution
		{
			const unsigned int MAX_FREQ = waveNoise->getMaximumFrequency();
			const unsigned int directionIndex = 0;
			const std::vector< std::vector< double > >& spectralEnergyDistribution = waveNoise->spectralEnergyDistribution;
			const std::vector< double >& mfs = spectralEnergyDistribution[ directionIndex ];
			std::vector< float > samples( MAX_FREQ );
			for ( int n = 0; n < MAX_FREQ; n++ )
			{
				samples[ n ] = mfs[ n ];
			}
			const int values_offset = 0;
			const char* overlay_text = NULL;
			const float scale_min = FLT_MAX;
			const float scale_max = FLT_MAX;
			ImVec2 graph_size = ImVec2( 256, 128 ); // hard-coded
			int stride = sizeof( float );
			ImGui::PlotLines( "", samples.data(), MAX_FREQ, values_offset, overlay_text, scale_min, scale_max, graph_size);
		}

		ImGui::TextColored( ImVec4( 0.f, 1.f, 0.f, 1.f ), "1D Profile (real part)" );

		// Display pre-computed 1D wave
		{
			const unsigned int mNARRAY = waveNoise->NARRAY;
			const std::vector< std::vector< double > >& mfs = waveNoise->fs;
			std::vector< float > samples( mNARRAY );
			for ( int n = 0; n < mNARRAY; n++ )
			{
				samples[ n ] = mfs[ n ][ 0 ]; // real part
			}
			const int values_offset = 0;
			const char* overlay_text = NULL;
			const float scale_min = FLT_MAX;
			const float scale_max = FLT_MAX;
			ImVec2 graph_size = ImVec2( 512, 128 ); // hard-coded
			int stride = sizeof( float );
			ImGui::PlotLines( "", samples.data(), mNARRAY, values_offset, overlay_text, scale_min, scale_max, graph_size );
		}
	}
	ImGui::End();

	return updateRequested;
}

/******************************************************************************
 * ...
 ******************************************************************************/
void Viewer::displayPeformance()
{
	const float DISTANCE = 10.0f;
	ImGuiIO& io = ImGui::GetIO();
	ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - DISTANCE, io.DisplaySize.y - DISTANCE),
		ImGuiCond_Always, ImVec2(1.0f, 1.0f));

	ImGuiWindowFlags flags = ImGuiWindowFlags_NoTitleBar
		| ImGuiWindowFlags_NoResize
		| ImGuiWindowFlags_AlwaysAutoResize
		| ImGuiWindowFlags_NoMove
		| ImGuiWindowFlags_NoSavedSettings
		| ImGuiWindowFlags_NoFocusOnAppearing
		| ImGuiWindowFlags_NoNav;

	if ( ImGui::Begin( "Example: Fixed Overlay", nullptr, flags ) )
	{
		ImGui::Text( "Performance" );
		ImGui::Separator();
		ImGui::Text( "FPS          %.1f", ImGui::GetIO().Framerate );
		ImGui::Text( "Frame        %.3f ms", 1000.0f / ImGui::GetIO().Framerate );
		ImGui::Text( "Window       %dx%d", this->width(), this->height() );
		ImGui::Separator();

		// LOG info
		/*std::cout << "\tTOTAL: ";
		std::cout << "\t" << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << std::setfill( ' ' ) << ( mGPUTimeElapsed * 1.e-6 ) << " ms\n";
		const std::string time = 
		ImGui::TextColored( ImVec4( 0.f, 1.f, 0.f, 1.f ),  );*/
		//ImGui::Text("%9.3f ms", mGPUTimeElapsed * 1e-6);
		ImGui::Text("GPU time");
		ImGui::SameLine(0); // colonne
		ImGui::TextColored( ImVec4( 0.f, 1.f, 0.f, 1.f ), "%9.3f ms", mGPUTimeElapsed * 1e-6 );

	}
	ImGui::End();
}
