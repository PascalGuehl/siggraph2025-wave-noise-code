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

#include "Viewer.h"

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
layout( location = 18 ) uniform int NN;

layout( location = 16 ) uniform float anisodd;

layout( location = 19 ) uniform float zoom;
layout( location = 25 ) uniform vec3 tX;

layout( location = 20 ) uniform float tV;
layout( location = 21 ) uniform float time;

layout( location = 22 ) uniform float RATIO;

layout( location = 23 ) uniform int QQ;
layout( location = 24 ) uniform float contrast;

layout( location = 30 ) uniform int Operator;
layout( location = 31 ) uniform int NRec;
layout( location = 32 ) uniform float Proba;

// 1D wave profile
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
float wavenoise_3D_postprocess( vec2 noise );
vec2 wavenoise_3D( vec3 p );

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
			float bb = acos( ( float( i ) + prng_random() ) / float( strat ) ) * orient + beta * ( 1.0 - orient );
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
			bb = acos( ( float( i ) + prng_random() ) / float( strat ) ) * orient + beta * ( 1.0 - orient );
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
			float anbeta = acos( 0.25 + anisodd * 0.5 * ( float( i ) + prng_random() ) / float( strat ) );
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
			float bb = acos( 0.5 + 0.1 * ( float( i ) + prng_random() ) / float( strat ) );
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
			bb = acos( 0.5 + 0.1 * ( float( i ) + prng_random() ) / float( strat ) );
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
			float anbeta = acos( ( float( i ) + prng_random() ) / float( strat ) );
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
			float bb = acos( ( float( i ) + prng_random() ) / float( strat ) ) * orient + beta * ( 1.0 - orient );
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
		if ( prng_random() < Proba )
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

/******************************************************************************
 * Post-process multi-dimensional procedural wave noise
 ******************************************************************************/
float wavenoise_3D_postprocess( vec2 noise )
{
	float noise_value = 0.0;
    
	if ( Operator <= 2 ) // non-cellular noise
	{ 
		if ( QQ == 0 ) // real part
		{
			noise_value = clamp( ( contrast * noise.x + 0.5 ), 0.0, 1.0 );
		}
		else if ( QQ == 1 ) // imaginary part
		{
			noise_value = clamp( ( contrast * noise.y + 0.5 ), 0.0, 1.0 );
		}
		else if ( QQ == 2 ) // modulus
		{
			noise_value = clamp( length( noise ) * 2.0 * contrast, 0.0, 1.0 );
		}
		else // phase
		{
			noise_value = atan( noise.y, noise.x ) / M_PI + 0.5; //smoothstep( -1.0, 1.0, cos( atan( noise.y, noise.x ) ) );
		}
	}
	else // cellular noise
	{
		if ( Operator == 3 )
		{
			noise_value = clamp( ( 0.5 * noise.x + 0.5 ), 0.0, 1.0 );
		}
		else if ( Operator == 4 )
		{
			noise_value = 0.8 * pow( clamp( ( 0.5 * noise.y + 0.5 ), 0.0, 1.0 ), contrast ) + 0.2;
		}
		else if ( Operator == 5 )
		{
			noise_value = step( 0.05, 0.5 * noise.y + 0.5 );
		}
		else
		{
			noise_value = 1.0 - pow( clamp( ( 0.5 * noise.y + 0.5 ), 0.0, 1.0 ), contrast );
		}
	}

	return noise_value;
}

////////////////////////////////////////////////////////////////////////////////
// MAIN PROGRAM
////////////////////////////////////////////////////////////////////////////////
void main()
{
	// Noise computation
	vec3 Xpos = vec3( zoom * Co.x + tX.x, zoom * Co.y + tX.y, zoom * Co.z + tX.z );
	vec2 noise = wavenoise_3D( Xpos );

	// Noise post-processing
	float noise_value = wavenoise_3D_postprocess( noise );

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
	frag_out = min( color_amb * noise_value + color_diff * lamb * noise_value + color_spec * spec, vec3( 1.0 ) );
}
)";

/******************************************************************************
 * Constructor
 ******************************************************************************/
Viewer::Viewer()
: nbMeshParts(0)
, waveNoise( nullptr )
, mQueryTimeElapsed( 0 )
, mGPUTimeElapsed( 0 )
, mUseContinuousAnimation( false )
, iRatio( 0.f )
, sFREQ_LOW( 0 )
, sFREQ_HIGH( 0 )
, Ffreq_low( 0.f )
, Ffreq_high( 0.f )
, old_power( 0.f )
{
	// Create and initialize noise
	waveNoise = new Wn::WaveNoise();
}

/******************************************************************************
 * Destructor
 ******************************************************************************/
Viewer::~Viewer()
{
	delete waveNoise;
	waveNoise = nullptr;
}

/******************************************************************************
 * Set title
 ******************************************************************************/
void Viewer::setTitle( const char* pText )
{
	GLFWwindow* w = window();
	glfwSetWindowTitle( w, "Multi-Dimensional Procedural Wave Noise" );
}

/******************************************************************************
 * ...
 ******************************************************************************/
void Viewer::initializeNoise()
{
	// Initialize noise
	waveNoise->initialize();

	// Dafault parameters value
	waveNoise->setTranslation( GLVec3( 0.f, 0.f, 0.f ) );
	waveNoise->setVelocity( 0.0 );
	waveNoise->setNbOrientations( 40 );
	waveNoise->Orient = 1.0f;
	waveNoise->setZoom( 0.2f );
	waveNoise->setTime( 0.0f );
	waveNoise->item_current = 0; // Gaussian
	waveNoise->Oper = 0;		 // Isowave
	waveNoise->old_item = 0;
	waveNoise->setRatio( 64.0f );
	waveNoise->mValueType = Wn::WaveNoise::EValueType::eReal;
	waveNoise->setPower( 25.f );
	waveNoise->contrast = 0.5;
	waveNoise->setRecursionProbability( 0.5 );
	waveNoise->setRecursionLevel( 3 );
	waveNoise->Anisodd = 0.5;
	
	// Precompute 1D profile wave(s)
	waveNoise->createIsotropicProceduralEnergyDistri();
	waveNoise->precomputePlanarWave( 4.0 );

	// Wave 1D profile
	waveNoise->tex = Texture1D::create();
	waveNoise->tex->alloc( waveNoise->NARRAY, GL_RGB8, waveNoise->fs_cr.data() );
	// Wave 1D profile's gradient
	waveNoise->texd = Texture1D::create();
	waveNoise->texd->alloc( waveNoise->NARRAY, GL_RGB8, waveNoise->fsd_cr.data() );
}

/******************************************************************************
 * Callback for custom initialization function
 ******************************************************************************/
void Viewer::init_ogl()
{
	prg_p = ShaderProgram::create({ {GL_VERTEX_SHADER, p_vert}, {GL_FRAGMENT_SHADER, p_frag} }, "prog");

	// Load OBJ file mesh
	auto mesh = Mesh::load(DATA_PATH + "/models/manual_cutaway_cube.obj");

	nbMeshParts = static_cast< int >( mesh.size() );
	// set the renderer and the materials for all the meshes parts
	for (int i = 0; i < nbMeshParts; ++i)
	{
		renderer_p.push_back(mesh[i]->renderer(1, 2, 3, 4, -1));
	}
	set_scene_center(mesh[0]->BB()->center());
	set_scene_radius(3.0 * mesh[0]->BB()->radius());
	// set_scene_radius(50.f);

	// Initialize the noise
	initializeNoise();

	// Initialize wave noise's helper parameters
	iRatio = std::log2f( waveNoise->getRatio() );
	sFREQ_LOW = waveNoise->getMinFrequency(); // 1
	sFREQ_HIGH = waveNoise->getMaxFrequency(); // 32
	Ffreq_low = 1.0 / 64.0; // TODO: explain this value!
	Ffreq_high = 32.0 / 64.0; // TODO: explain this value!
	old_power = 1.f;

	// Device timer
	glCreateQueries( GL_TIME_ELAPSED, 1, &mQueryTimeElapsed );

	// Set global GL states
	glEnable( GL_DEPTH_TEST );
	glClearColor( 1.0, 1.0, 1.0, 0.0 );
}

/******************************************************************************
 * ...
 ******************************************************************************/
void Viewer::close_ogl()
{
	// Device timer
	glDeleteQueries( 1, &mQueryTimeElapsed );
}

/******************************************************************************
 * Update noise if requested (not clean and not optimized...)
 ******************************************************************************/
void Viewer::updateNoise()
{
	waveNoise->setMinFrequency( (int)( Ffreq_low * 64.0 ) );
	waveNoise->setMaxFrequency( (int)( Ffreq_high * 64.0 ) );

	if ( waveNoise->getMinFrequency() != sFREQ_LOW ||
		 waveNoise->getMaxFrequency() != sFREQ_HIGH ||
		 waveNoise->item_current != waveNoise->old_item ||
		( waveNoise->item_current >= 4 && old_power != waveNoise->getPower() ) )
	{
		sFREQ_LOW = waveNoise->getMinFrequency();
		sFREQ_HIGH = waveNoise->getMaxFrequency();
		waveNoise->old_item = waveNoise->item_current;
		old_power = waveNoise->getPower();
		if ( waveNoise->item_current >= 4 && waveNoise->item_current <= 9 )
		{
			int pow_2 = 1, N = 1;
			while ( N < waveNoise->NARRAY )
			{
				pow_2++;
				N *= 2;
			}
			pow_2--;
			printf( "NARRAY=%d, N=%d, pow2=%d\n", waveNoise->NARRAY, N, pow_2 );
			double* A = waveNoise->MakeSpatialWaveProfile( pow_2 );
			waveNoise->precomputePlanarWaveFromFFT1D( A, N, pow_2 );
			delete A;
		}
		else
		{
			waveNoise->createIsotropicProceduralEnergyDistri();
			waveNoise->precomputePlanarWave( 4.0 );
		}
		waveNoise->tex->alloc( waveNoise->NARRAY, GL_RGB8, waveNoise->fs_cr.data() );
		waveNoise->texd->alloc( waveNoise->NARRAY, GL_RGB8, waveNoise->fsd_cr.data() );
	}
}

/******************************************************************************
 * Callback for custom rendering/compute function
 ******************************************************************************/
void Viewer::draw_ogl()
{
	glClear( GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );

	// Update noise if requested (not clean and not optimized...)
	updateNoise();

	// Bind shader program
	prg_p->bind();

	// Bind textures
	// - precomputed 1D wave profile
	waveNoise->tex->bind( 0 );
	// - precomputed 1D wave gradient
	waveNoise->texd->bind( 1 );

	// Model-View-Projection transforms
	GLMat4 sc = Transfo::scale( 2.5 );
	GLMat4 rotx = Transfo::rotateX( -60.0 );
	GLMat4 rotz = Transfo::rotateZ( 25.0 );
	const GLMat4& proj = this->get_projection_matrix();
	const GLMat4& mv = this->get_view_matrix() * sc * rotx * rotz;
	set_uniform_value( 1, proj );
	set_uniform_value( 2, mv );
	set_uniform_value( 3, Transfo::inverse_transpose( mv ) );

	set_uniform_value( 14, GLVec3( 0, 2.0, 3.0 ) ); // light position

	GLVec3 translation = waveNoise->getTranslation() + GLVec3( 5.f, 5.f, 5.f );
	set_uniform_value( 25, translation );
	
	set_uniform_value( 15, waveNoise->getNbOrientations() ); // problem? uniforms 15 and 18 are the same! why?
	set_uniform_value( 16, waveNoise->Anisodd );
	
	set_uniform_value( 18, waveNoise->getNbOrientations() ); // problem? uniforms 15 and 18 are the same! why?

	set_uniform_value( 20, waveNoise->getVelocity() );
	if ( ! mUseContinuousAnimation )
	{
		set_uniform_value( 21, waveNoise->getTime() );
	}
	else
	{
		set_uniform_value( 21, current_time() );
	}

	set_uniform_value( 22, waveNoise->getRatio() );
	set_uniform_value( 19, waveNoise->getZoom() * waveNoise->getRatio() );
	set_uniform_value( 23, static_cast< int >( waveNoise->mValueType ) );
	set_uniform_value( 24, waveNoise->contrast );
	set_uniform_value( 30, waveNoise->Oper );

	set_uniform_value( 31, waveNoise->getRecursionLevel() );
	set_uniform_value( 32, waveNoise->getRecursionProbability() );

	// GPU timer
	GLuint64 result = 0;
	mGPUTimeElapsed = 0;

	for ( int i = 0; i < renderer_p.size(); i++ )
	{
		auto mat = renderer_p[ i ]->material();

		set_uniform_value( 10, mat->Kd );
		set_uniform_value( 11, mat->Ka );
		set_uniform_value( 12, mat->Ks );
		set_uniform_value( 13, mat->Ns );
		set_uniform_value( "light_pos", GLVec3( 0, 2, 2.5 ) );

		// GPU timer
		glBeginQuery( GL_TIME_ELAPSED, mQueryTimeElapsed );

		// Render
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

		const float sliderWidth = 100.0f; // largeur personnalisee pour chaque slider
		if (ImGui::CollapsingHeader("Spatial Transformation"))
		{
			GLVec3 translation = waveNoise->getTranslation();
			bool translationModified = false;
			ImGui::PushItemWidth(sliderWidth);
			translationModified |= ImGui::SliderFloat( "X", &translation.x(), 0.0, 10.0 );
			ImGui::SameLine();
			translationModified |= ImGui::SliderFloat( "Y", &translation.y(), 0.0, 10.0 );
			ImGui::SameLine();
			translationModified |= ImGui::SliderFloat( "Z", &translation.z(), 0.0, 10.0 );
			ImGui::PopItemWidth(); // remet la largeur par defaut
			if ( translationModified )
			{
				waveNoise->setTranslation( translation );
			}

			float zoom = waveNoise->getZoom();
			if ( ImGui::SliderFloat( "Zoom", &zoom, 0.01f, 2.0f ) )
			{
				waveNoise->setZoom( zoom );
			}
		}
	
		// ImGui::Text("MEM:  %2.2lf %%", 100.0 * mem_);

		if (ImGui::CollapsingHeader("Spectral Configuration"))
		{
			//ImGui::Text("Quality vs.Framerate");
			if (ImGui::TreeNode("Quality vs. Framerate"))
			{
				int nbOrientations = waveNoise->getNbOrientations();
				if ( ImGui::SliderInt( "NDir", &nbOrientations, 1, 100 ) )
				{
					waveNoise->setNbOrientations( nbOrientations );
				}

				if ( ImGui::SliderFloat( "Slice size", &iRatio, 0.0, 8.0 ) )
				{
					waveNoise->setRatio( (float)pow( 2.0, iRatio ) );
				}

				ImGui::TreePop();
			}

			//ImGui::Text("Frequency Ranges");
			if (ImGui::TreeNode("Frequency Ranges"))
			{
				ImGui::SliderFloat("aniso wave dir width", &waveNoise->Anisodd, 0.0, 1.0);

				ImGui::PushItemWidth(sliderWidth);
				ImGui::SliderFloat("FreqMin", &Ffreq_low, 1.0 / 64.0, 1.0);
				ImGui::SameLine();
				ImGui::SliderFloat("FreqMax", &Ffreq_high, 1.0 / 64.0, 1.0);
				ImGui::PopItemWidth(); // remet la largeur par d�faut

				ImGui::TreePop();
			}
		}

		if (ImGui::CollapsingHeader("Waveforms & Noise Types"))
		{
			const char* items[] = { "[noise] gaussian",			   "[noise] white (0dB)",	  "[noise] blue (+3dB)",
								    "[noise] brown (-6dB)",		   "[non-gaussian] crystal1", "[non-gaussian] web",
								    "[non-gaussian] marble",	   "[non-gaussian] crystal2", "[non-gaussian] scratches",
								    "[non-gaussian] smooth cells", "[noise] two ampli levels" };
			ImGui::Combo("Wave type", &waveNoise->item_current, items, IM_ARRAYSIZE(items));

			const char* itemsop[] = {
				"[isotropic] Sum", "[anisotropic] Sum - one direction", "[anisotropic] Sum - two directions", "[cellular] Random polytopes", "[cellular] Cellular",
				"[cellular] Hyperplan",	 "[cellular] Reversed Cellular" };
			ImGui::Combo("Operator", &waveNoise->Oper, itemsop, IM_ARRAYSIZE(itemsop));
			float power = waveNoise->getPower();
			if ( ImGui::SliderFloat( "[non-gaussian] wave sharpness", &power, 0.2f, 50.0f ) )
			{
				waveNoise->setPower( power );
			}
		}

		if (ImGui::CollapsingHeader("Cellular Noise Settings (STIT)"))
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

		if ( ImGui::CollapsingHeader( "Post-Processing & Display" ) )
		{
			auto& valueTypeNames = Wn::WaveNoise::mValueTypeNames;
			static std::vector< const char* > valueTypeItems; // built only once!
			if ( valueTypeItems.size() != valueTypeNames.size() )
			{
				valueTypeItems.clear();
				valueTypeItems.reserve( valueTypeNames.size() );
				for ( auto& s : valueTypeNames )
				{
					valueTypeItems.push_back( s.c_str() );
				}
			}
			int valueType = static_cast< int >( waveNoise->mValueType );
			if ( ImGui::Combo( "Value", &valueType, valueTypeItems.data(), (int)valueTypeItems.size() ) )
			{
				waveNoise->mValueType = static_cast< Wn::WaveNoise::EValueType >( valueType );
			}

			ImGui::SliderFloat("Contrast", &waveNoise->contrast, 0.1f, 1.0f);
		}

		if ( ImGui::CollapsingHeader( "Temporal & Animation Controls" ) )
		{
			ImGui::PushItemWidth( sliderWidth );
			float velocity = waveNoise->getVelocity();
			if ( ImGui::SliderFloat( "Speed", &velocity, 0.0f, 0.1f ) )
			{
				waveNoise->setVelocity( velocity );
			}
			ImGui::SameLine();
			float time = waveNoise->getTime();
			if( ImGui::SliderFloat( "Time", &time, 0.0, 5.0 ) )
			{
				waveNoise->setTime( time );
			}
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
		
		const char* items[] = { "[noise] gaussian",			   "[noise] white (0dB)",	  "[noise] blue (+3dB)",
								    "[noise] brown (-6dB)",		   "[non-gaussian] crystal1", "[non-gaussian] web",
								    "[non-gaussian] marble",	   "[non-gaussian] crystal2", "[non-gaussian] scratches",
								    "[non-gaussian] smooth cells", "[noise] two ampli levels" };
		const std::string wave_type = std::string( items[ waveNoise->item_current ] );
		ImGui::TextColored( ImVec4( 0.f, 1.f, 0.f, 1.f ), wave_type.c_str() );

		// Display amplitude's energy distribution
		{
			const unsigned int MAX_FREQ = waveNoise->getMaximumFrequency();
			const unsigned int directionIndex = 0;
			const std::vector< std::vector< double > >& spectralEnergyDistribution = waveNoise->getSpectralEnergyDistribution();
			const std::vector< double >& mfs = spectralEnergyDistribution[ directionIndex ];
			std::vector< float > samples( MAX_FREQ );
			for ( unsigned int n = 0; n < MAX_FREQ; n++ )
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
			for ( unsigned int n = 0; n < mNARRAY; n++ )
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
