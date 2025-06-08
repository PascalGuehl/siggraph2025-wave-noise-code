//
#include "fbo.h"
#include "gl_viewer.h"
#include "mesh.h"
#include "shader_program.h"
#include "texture1d.h"
#include "texture3d.h"

#include <iostream>

#include <GLFW/glfw3.h>

#define macro_str(s) #s
#define macro_xstr(s) macro_str(s)
#define DATA_PATH std::string(macro_xstr(MYAPP_DATA_PATH))
#define SHADERS_PATH std::string(macro_xstr(MYAPP_SHADERS_PATH))

extern "C"
{
	__declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001;
}

using namespace EZCOGL;

//////////////////////////////////////////////////////////////////////////////////////
//  FOR USING FFT 1D
//////////////////////////////////////////////////////////////////////////////////////

template <class T, class U>
class hvPair
{
protected:
	T left;
	U right;

public:
	hvPair()
	{
		left = T(0);
		right = U(0);
	}
	hvPair(const T& x, const U& y)
	{
		left = x;
		right = y;
	}

	void setLeft(const T& x)
	{
		left = x;
	}
	void setRight(const U& x)
	{
		right = x;
	}
	T getLeft()
	{
		return left;
	}
	U getRight()
	{
		return right;
	}
	bool operator==(const hvPair<T, U>& pp) const
	{
		return left == pp.left && right == pp.right;
	}

	// complex numbers
	double mod()
	{
		return sqrt((double)left * (double)left + (double)right * (double)right);
	}
	double energy()
	{
		return (double)left * (double)left + (double)right * (double)right;
	}
	double phase()
	{
		double rr, r = (double)left, i = (double)right;

		if (r == 0.0 && i == 0.0)
			return (0.0);
		if ((r > 0.0 ? r : -r) > (i > 0.0 ? i : -i))
		{
			rr = i / r;
			if (r < 0.0)
				rr = M_PI + atan(rr);
			else
				rr = atan(rr);
		}
		else
		{
			rr = r / i;
			if (i > 0.0)
				rr = M_PI / 2.0 - atan(rr);
			else
				rr = 3.0 * M_PI / 2.0 - atan(rr);
		}
		if (rr > M_PI)
			return (rr - 2.0 * M_PI);
		else if (rr < -M_PI)
			return (rr + 2.0 * M_PI);
		return (rr);
	}
};

template <class T>
class hvArray1
{
protected:
	T* t;
	int sx;

public:
	hvArray1()
	{
		t = 0;
		sx = 0;
	}
	hvArray1(int x, T nil)
	{
		t = new T[x];
		if (t == 0)
		{
			sx = -1;
			return;
		}
		for (int i = 0; i < x; i++)
			t[i] = nil;
		sx = x;
	}
	T* data() const
	{
		return t;
	}

	// copy
	hvArray1(const hvArray1& a)
	{
		hvFatal("No temporary creation of hvArray1!");
	}

	// affectation
	hvArray1& operator=(const hvArray1& a)
	{
		if (this != &a)
		{
			if (t != 0)
				delete[] t;
			if (a.isInvalid())
			{
				t = 0;
				sx = -1;
				return *this;
			}
			sx = a.sx;
			t = new T[sx];
			if (t == 0)
			{
				sx = -1;
				return *this;
			}
			for (int i = 0; i < sx; i++)
				t[i] = a.t[i];
		}
		return *this;
	}

	void reset(int x, T nil)
	{
		if (t != 0)
			delete[] t;
		t = new T[x];
		if (t == 0)
		{
			sx = -1;
			return;
		}
		for (int i = 0; i < x; i++)
			t[i] = nil;
		sx = x;
	}
	void reset(int x)
	{
		if (t != 0)
			delete[] t;
		t = new T[x];
		if (t == 0)
		{
			sx = -1;
			return;
		}
		sx = x;
	}
	void reset()
	{
		if (t != 0)
			delete[] t;
		t = 0;
		sx = 0;
	}

	// isInvalid
	bool isInvalid() const
	{
		return sx == -1;
	}
	// isVoid
	bool isVoid() const
	{
		return t == 0;
	}

	// operations
	void clear(T nil)
	{
		for (int i = 0; i < sx; i++)
			t[i] = nil;
	}

	// selectors
	int size() const
	{
		return sx;
	}

	T& operator[](int x)
	{
		if (x < 0 || x >= sx)
		{
			fprintf(stderr, "out of range!");
		}
		if (t == 0)
		{
			fprintf(stderr, "hvArray1 is void!");
		}
		return t[x];
	}
	T get(int x) const
	{
		if (x < 0 || x >= sx)
		{
			fprintf(stderr, "out of range!");
		}
		if (t == 0)
		{
			fprintf(stderr, "hvArray1 is void!");
		}
		return t[x];
	}
	T* getPointer(int x)
	{
		if (x < 0 || x >= sx)
		{
			fprintf(stderr, "out of range!");
		}
		if (t == 0)
		{
			fprintf(stderr, "hvArray1 is void!");
		}
		return t + x;
	}
	void update(int x, T val)
	{
		if (x < 0 || x >= sx)
		{
			fprintf(stderr, "out of range!");
		}
		if (t == 0)
		{
			fprintf(stderr, "hvArray1 is void!");
		}
		t[x] = val;
	}
	template <class U>
	void updateRange(int start, int end, U* val)
	{
		if (start < 0 || start >= sx)
		{
			fprintf(stderr, "start out of range!");
		}
		if (end < 0 || end >= sx)
		{
			fprintf(stderr, "end out of range!");
		}
		if (t == 0)
		{
			fprintf(stderr, "hvArray1 is void!");
		}
		int i;
		for (i = start; i <= end; i++)
			t[i] = T(val[i - start]);
	}

	~hvArray1()
	{
		if (t != 0)
			delete[] t;
	}

	// FFT algorithm
	// F is pointer to pairs (integer, imaginary) part of complex numbers of type T
	// N is such that 2^N is the size of the array F
	static void fft(hvArray1<hvPair<T, T>>& F, int N, int s, int id, bool divn)
	{
		int n, i, j, k, l, nv2, ip;
		hvPair<T, T> u, w, tmp;

		n = 1 << N;	  /* n=2^N */
		nv2 = n >> 1; /* nv2=2^(N-1) */
		j = 1;
		for (i = 1; i <= n - 1; i++)
		{
			if (i < j)
			{
				tmp = F[(j - 1) * s + id];
				F[(j - 1) * s + id] = F[(i - 1) * s + id];
				F[(i - 1) * s + id] = tmp;
			}
			k = nv2;
			while (k < j)
			{
				j = j - k;
				k = k >> 1;
			}
			j = j + k;
		}
		for (l = 1; l <= N; l++)
		{
			int le = 1 << l;
			int le1 = le >> 1;
			u = hvPair<T, T>(T(1.0), T(0.0));
			double a = M_PI / (double)le1;
			w = hvPair<T, T>(T(cos(a)), T(-sin(a)));
			for (j = 1; j <= le1; j++)
			{
				for (i = j; i <= n; i += le)
				{
					ip = i + le1;
					tmp = hvPair<T, T>(F[((ip - 1) * s + id)].getLeft() * u.getLeft() -
						F[((ip - 1) * s + id)].getRight() * u.getRight(),
						F[((ip - 1) * s + id)].getLeft() * u.getRight() +
						F[((ip - 1) * s + id)].getRight() * u.getLeft());
					F[((ip - 1) * s + id)] = hvPair<T, T>(F[((i - 1) * s + id)].getLeft() - tmp.getLeft(),
						F[((i - 1) * s + id)].getRight() - tmp.getRight());
					F[((i - 1) * s + id)] = hvPair<T, T>(F[((i - 1) * s + id)].getLeft() + tmp.getLeft(),
						F[((i - 1) * s + id)].getRight() + tmp.getRight());
				}
				u = hvPair<T, T>(u.getLeft() * w.getLeft() - u.getRight() * w.getRight(),
					u.getLeft() * w.getRight() + u.getRight() * w.getLeft());
			}
		}
		if (divn)
		{
			double ff = 1.0 / (double)n;
			for (i = 1; i <= n; i++)
				F[((i - 1) * s + id)] =
				hvPair<T, T>(F[((i - 1) * s + id)].getLeft() * ff, F[((i - 1) * s + id)].getRight() * ff);
		}
	}
};

//////////////////////////////////////////////////////////////////////////////////////
//  SHADERS
//////////////////////////////////////////////////////////////////////////////////////

static const std::string p_vert = R"(
#version 430
layout(location=1) uniform mat4 projectionMatrix;
layout(location=2) uniform mat4 viewMatrix;
layout(location=3) uniform mat3 normalMatrix;

layout(location=1) in vec3 position_in;
layout(location=2) in vec3 normal_in;
layout(location=3) in vec3 text_in;
layout(location = 4) in vec3 tangents_in;

out vec3 Po;
out vec3 No;
out vec3 Co;
out vec3 NCo, TCo;
out mat3 NoMat;


void main()
{
	Co = vec3(0.5)+position_in*0.5;
	NCo = normal_in;
	TCo = tangents_in;
	No = normalMatrix * normal_in;
	NoMat = normalMatrix;
	vec4 Po4 = viewMatrix * vec4(position_in,1);
	Po = Po4.xyz;
	gl_Position = projectionMatrix * Po4;
}
)";

static const std::string p_frag = R"(
#version 430
precision highp float;
in vec3 Po;
in vec3 No;
in vec3 Co;
in vec3 NCo, TCo;
in mat3 NoMat;
out vec3 frag_out;

const int NDIRMAX = 100; 
layout(location=10) uniform vec3 color_diff;
layout(location=11) uniform vec3 color_amb;
layout(location=12) uniform vec3 color_spec;
layout(location=13) uniform float Ns;
layout(location=14) uniform vec3 light_pos;
layout(location=15) uniform int NDIR;
layout(location=16) uniform float anisodd;
const float orient=1.0;
const int NNmin=0;
layout(location=18) uniform int NN;
layout(location=19) uniform float zoom;
layout(location=20) uniform float tV;
layout(location=21) uniform float time;
layout(location=22) uniform float RATIO;
layout(location=23) uniform int QQ;
layout(location=24) uniform float contrast;
layout(location=25) uniform float tX;
layout(location=26) uniform float tY;
layout(location=27) uniform float tZ;

layout(location=30) uniform int Operator;
layout(location=31) uniform int NRec;
layout(location=32) uniform float Proba;

layout(binding=0) uniform sampler1D wave;
layout(binding=1) uniform sampler1D waved;

uint hash( uint x ) {
    x += ( x << 10u );
    x ^= ( x >>  6u );
    x += ( x <<  3u );
    x ^= ( x >> 11u );
    x += ( x << 15u );
    return x;
}
uint hash( uvec2 v ) { return hash( v.x ^ hash(v.y)); }

uint seed=1;
float random()
{
 seed = seed*22695477u+1u;
 uint v = seed & 0xffffu;
 return float(v)/float(0xffffu);
}
float next() { return 2.0*random()-1.0; }
void seeding(uint x, uint y, uint z)  { seed = hash(x + hash(y + hash(z))); }

//float RATIO=8.0;
const float M_PI = 3.14159265358979;
int NF = 1;
int NFS[3];
//const int NFST = 10; // 10;

//vec3 grad;

vec2 anisowavenoise3D(float x, float y, float z, int nd)
{
	const int strat=2;
	NF = NDIR/strat/nd;  // very simple stratified sampling

	vec2 sum = vec2(0.0); // result is complex number

	float analpha, analphawidth, anbeta;
	int i, j, k;
	int dir=0;
	for (j = 0; j < nd*NF; j++) //strat alpha
	{
		for (i = 0; i < strat; i++) // strat beta
		{
			seeding(uint(dir) * 5u, 2u, 3u);
			float tspeed = tV*sign(next())*time;

			anbeta = acos(0.25+anisodd*0.5*(float(i)+(0.5*next()+0.5))/float(strat)); 
			float alpha1 = (j<NF ? M_PI/4.0: 3.0*M_PI/4.0) +anisodd*M_PI/8.0 * float(j) / float(NF+1);
			float alpha2 = (j<NF ? M_PI/4.0: 3.0*M_PI/4.0) +anisodd*M_PI/8.0 *(float(j) + 1.0) / float(NF+1);
			analpha = (alpha1 + alpha2) / 2.0;
			analphawidth = alpha2 - alpha1;

			float aa =  M_PI / 4.0 + analpha + next() * analphawidth * 0.5 ;
			float beta = anbeta + M_PI / 4.0; 

			float id = sin(beta)*cos(aa) * x + sin(beta)*sin(aa) * y + cos(beta) * z;
			int iid = (id>0.0 ? int(id):int(id)-1); //FASTFLOOR(id);
			float oid = id - float(iid);

			seeding(uint(7 * dir), uint(4 * iid), 4u);
			float pl = 0.5 + 0.2 * next();
			if (oid < pl - 0.3) oid = 0.0;
			else if (oid > pl + 0.3) oid = 1.0;
			else oid = (oid - pl + 0.3) / 0.6;

			seeding(uint(4 * dir), uint(4 * iid), 0);
			float alpha = (analpha + next() * analphawidth);
			float bb = acos(0.5+0.1*(float(i)+(0.5*next()+0.5))/float(strat));

			float dd = (1.0 / RATIO * (sin(bb)*cos(alpha) * x + sin(bb)*sin(alpha) * y + cos(bb) * z) + 2.0*next() -tspeed);
			if (dir>=NNmin && dir<NN && dir<NDIR) sum += (1.0 - oid) * (2.0*texture(wave,fract(dd)).xy-vec2(1.0));

			seeding(uint(4 * dir), uint(4 * (iid + 1)), 0);
			alpha = (analpha + next() * analphawidth);
			bb = acos(0.5+0.1*(float(i)+(0.5*next()+0.5))/float(strat));

			dd = (1.0 / RATIO * (sin(bb) * cos(alpha) * x + sin(bb) * sin(alpha) * y + cos(bb) * z) + 2.0 * next() -tspeed);
			if (dir>=NNmin && dir<NN && dir<NDIR) sum += oid * (2.0*texture(wave,fract(dd)).xy-vec2(1.0));

			dir += 1;
		}
	}
	return sum /  vec2(2.0+0.2*float(NN-NNmin));
}

vec2 isowavenoise3D(float x, float y, float z)
{
	vec2 sum = vec2(0.0); // result is complex number
	float analpha, analphawidth, anbeta, anbetawidth;
	float poids = 1.0; 
	int i, j;
	const int strat=4;
	NF = NDIR/strat;  // very simple stratified sampling
	int dir=0;
	for (j = 0; j <= NF; j++) //strat alpha
	{
		for (i = 0; i < strat; i++) // strat beta
		{
			// strat angular sector
			seeding(uint(dir) * 5u, 2u, 3u);
			float tspeed = tV*sign(next())*time;
			anbeta = acos((float(i)+(0.5*next()+0.5))/float(strat)); 
			float alpha1 = 2.0 * M_PI * float(j) / float(NF+1);
			float alpha2 = 2.0 * M_PI * (float(j) + 1.0) / float(NF+1);
			analpha = (alpha1 + alpha2) / 2.0;
			analphawidth = alpha2 - alpha1;

			// slice
			float aa =  M_PI / 4.0*orient + analpha + next() * analphawidth * 0.5 ;
			float beta = anbeta + M_PI / 4.0*orient; 
			float id = sin(beta)*cos(aa) * x + sin(beta)*sin(aa) * y + cos(beta) * z;
			int iid = (id>0.0 ? int(id):int(id)-1); //FASTFLOOR(id);
			float oid = id - float(iid);
			seeding(uint(7 * dir), uint(4 * iid), 4u);
			float pl = 0.5 + 0.2 * next();
			if (oid < pl - 0.3) oid = 0.0;
			else if (oid > pl + 0.3) oid = 1.0;
			else oid = (oid - pl + 0.3) / 0.6;

			// wave left
			seeding(uint(4 * dir), uint(4 * iid), 0);
			float alpha = (analpha + next() * analphawidth)*orient+aa*(1.0-orient);
			float bb = acos((float(i)+(0.5*next()+0.5))/float(strat))*orient+beta*(1.0-orient);
			float dd = (1.0 / RATIO * (sin(bb)*cos(alpha) * x + sin(bb)*sin(alpha) * y + cos(bb) * z) + 2.0 * next()-tspeed);
			if (dir>=NNmin && dir<NN && dir<NDIR) {
				sum += poids*(1.0 - oid) * (2.0*texture(wave,fract(dd)).xy-vec2(1.0));
			}
			// wave right
			seeding(uint(4 * dir), uint(4 * (iid + 1)), 0);
			alpha = (analpha + next() * analphawidth)*orient+aa*(1.0-orient);
			bb = acos((float(i)+(0.5*next()+0.5))/float(strat))*orient+beta*(1.0-orient);
			dd = (1.0 / RATIO * (sin(bb) * cos(alpha) * x + sin(bb) * sin(alpha) * y + cos(bb) * z) + 2.0 * next() -tspeed);
			if (dir>=NNmin && dir<NN && dir<NDIR) {
				sum += poids* oid * (2.0*texture(wave,fract(dd)).xy-vec2(1.0));
			}
			dir += 1;
		}
	}
    return sum/vec2(2.0+0.2*float(NN-NNmin));
}

vec2 isowavenoise3Drec(float x, float y, float z, uint rec)
{
	uint resv = 0;
	float rvalue = 2.0;

	float analpha, analphawidth, anbeta, anbetawidth;
	int i, j;
	const int strat=2;
	NF = NDIR/strat;  // very simple stratified sampling
	int dir=0;
	for (j = 0; j <= NF; j++) //strat alpha
	{
		for (i = 0; i < strat; i++) // strat beta
		{
			seeding(uint(dir) * 5u, 2u, 3u+rec);
			float tspeed = tV*sign(next())*time;
			anbeta = acos((float(i)+(0.5*next()+0.5))/float(strat)); 
			float alpha1 = M_PI * float(j) / float(NF+1);
			float alpha2 = M_PI * (float(j) + 1.0) / float(NF+1);
			analpha = (alpha1 + alpha2) / 2.0;
			analphawidth = alpha2 - alpha1;

			float aa =  M_PI / 4.0*orient + analpha + next() * analphawidth * 0.5 ;
			float beta = anbeta + M_PI / 4.0*orient; 

			float id = sin(beta)*cos(aa) * x + sin(beta)*sin(aa) * y + cos(beta) * z  + 5.0 + tspeed;
			int iid = (id>0.0 ? int(id):int(id)-1); //FASTFLOOR(id);
			float oid = id - float(iid);

			seeding(uint(7 * dir), uint(4 * iid), 4u);
			float pl = 0.5 + 0.3 * next();

			if (oid <= pl)
			{
				float dist1 = pl - oid;
				rvalue = min(dist1,rvalue);
				seeding(uint(7 * dir), uint(4 * (iid - 1)), 4u);
				float pl2 = 0.5 + 0.3 * next();
				float dist2 = oid + 1.0 - pl2;
				rvalue = min(dist2,rvalue);
			}
			else
			{
				float dist1 = oid - pl;
				rvalue = min(dist1,rvalue);
				seeding(uint(7 * dir), uint(4 * (iid + 1)), 4u);
				float pl2 = 0.5 + 0.3 * next();
				float dist2 = 1.0 - oid + pl2;
				rvalue = min(dist2,rvalue);
			}
			// cell index
			int iidabs = abs(iid);
			if ((oid <= pl && iidabs % 2 == 0) || (oid > pl && iidabs % 2 == 1)) resv += 1;
			resv = resv * 2u;

			if (oid <= pl) seeding(uint(4 * dir), uint(4 * iid), 0);
			else seeding(uint(4 * dir), uint(4 * (iid+1)), 0);
			float alpha = (analpha + next() * analphawidth)*orient+aa*(1.0-orient);
			float bb = acos((float(i)+(0.5*next()+0.5))/float(strat))*orient+beta*(1.0-orient);

			float dd = sin(bb)*cos(alpha) * x + sin(bb)*sin(alpha) * y + cos(bb) * z + 5.0 + tspeed;
			int p = (dd>0.0 ? int(dd):int(dd)-1); //FASTFLOOR(dd);
			dd = dd - float(p);
			seeding(uint(7 * dir), uint(4 * p), 4u);
			pl = 0.5 + 0.3 * next();
			if (dd <= pl)
			{
				float dist1 = pl - dd;
				rvalue = min(dist1,rvalue);
				seeding(uint(7 * dir), uint(4 * (p - 1)), 4u);
				float pl2 = 0.5 + 0.3 * next();
				float dist2 = dd + 1.0 - pl2;
				rvalue = min(dist2,rvalue);
			}
			else
			{
				float dist1 = dd - pl;
				rvalue = min(dist1,rvalue);
				seeding(uint(7 * dir), uint(4 * (p + 1)), 4u);
				float pl2 = 0.5 + 0.3 * next();
				float dist2 = 1.0 - dd + pl2;
				rvalue = min(dist2,rvalue);
			}
			// cell index
			iidabs = abs(p);
			if ((dd <= pl && iidabs % 2 == 0) || (dd > pl && iidabs % 2 == 1)) resv += 1;
			resv = resv * 2u;

			dir += 1;
		}
	}
	uint cellident = (resv * 1453u) % 255u;
	return vec2(float(cellident)/255.0*2.0-1.0,rvalue*20.0-1.0);
}

vec2 cellwavenoise3D(float x, float y, float z)
{
int rr=0;
vec2 res = isowavenoise3Drec(x,y,z,0);
bool cont=true;
for (rr=1; rr<NRec && cont; rr++)
	{
	uint cellident = uint(255.0*(res.x+1.0)*0.5);
	seeding(cellident, uint(rr), 3);
	if (0.5*(next()+1.0)<Proba)
		{
			vec2 res2 = isowavenoise3Drec(x,y,z,uint(rr));
			res = vec2(res2.x, min(res.y,res2.y));
		}
	else cont=false;
	}
return res;
}

void main()
{
	vec3 Nco = normalize(NCo);
	vec3 Nnn = NCo;
	if (gl_FrontFacing==false) Nnn = -Nnn;
	vec3 L = normalize(light_pos-Po);
	vec3 col;

	vec3 Npp = normalize(TCo);
	vec3 Ncc = cross(Nnn,Npp);
	vec2 waven=vec2(0.0);
    float val = 0.0;
	vec3 Xpos = vec3(zoom*Co.x+tX, zoom*Co.y+tY, zoom*Co.z+tZ);
	if (Operator==0) waven = isowavenoise3D(Xpos.x, Xpos.y, Xpos.z);
	else if (Operator==1) waven = anisowavenoise3D(Xpos.x, Xpos.y, Xpos.z, 1);
	else if (Operator==2) waven = anisowavenoise3D(Xpos.x, Xpos.y, Xpos.z, 2);
	else waven = cellwavenoise3D(0.1*Xpos.x, 0.1*Xpos.y, 0.1*Xpos.z);
	if (Operator<=2) { 
		if (QQ==0) val = clamp((contrast*waven.x+0.5), 0.0, 1.0);
		else if (QQ==1) val = clamp((contrast*waven.y+0.5), 0.0, 1.0);
		else if (QQ==2) val = clamp(length(waven)*2.0*contrast, 0.0, 1.0);
		else val = atan(waven.y, waven.x)/M_PI+0.5; //smoothstep(-1.0, 1.0, cos(atan(waven.y, waven.x))); 
		}
	else {
		if (Operator==3) val = clamp((0.5*waven.x+0.5), 0.0, 1.0);
		else if (Operator==4) val = 0.8*pow(clamp((0.5*waven.y+0.5), 0.0, 1.0),contrast)+0.2;
		else if (Operator==5) val = step(0.05, 0.5*waven.y+0.5);
		else val = 1.0-pow(clamp((0.5*waven.y+0.5), 0.0, 1.0),contrast);
	}

	vec3 N = normalize(NoMat*normalize(Nnn));
	float lamb = abs(dot(N,L));
	vec3 E = normalize(-Po);
	vec3 R = reflect(-L, N);
	float spec = Ns==0.0 ? 0.0 : pow( max(dot(R,E), 0.0), Ns);

	frag_out = min(color_amb*val+color_diff*lamb*val+color_spec*spec,vec3(1));
}
)";

//////////////////////////////////////////////////////////////////////////////////////
//  MAIN  VIEWER
//////////////////////////////////////////////////////////////////////////////////////

// precision of the pre-computed wave (sampling rate)
const int NARRAY = 512;
const int MAX_FREQ = NARRAY / 2; // half of half because of FFT symetry and Nyquist
const int NDIR = 20;			 // 40+50;

int FREQ_LOW = 1, sFREQ_LOW = 1;
int FREQ_HIGH = 32, sFREQ_HIGH = 32;
float Ffreq_low = 1.0 / 64.0, Ffreq_high = 32.0 / 64.0;

// Creation du VIEWER
class Viewer : public GLViewer
{
	ShaderProgram::UP prg_p;
	std::vector<MeshRenderer::UP> renderer_p;
	int nbMeshParts;
	int Ndir, Nc, Na, Oper, NRec;
	float tX, tY, tZ, tV;
	float Orient, Period, Zoom, Time, Ratio, iRatio, Power, old_power, contrast, Proba, Anisodd;
	int item_current, old_item;
	int complex_current;

	Texture1D::SP tex, texd;

	// Isotropic noise pre-computed arrays
	double fs[NARRAY][3], vmax, vmin;
	GLubyte fs_cr[3 * NARRAY];
	GLubyte fsd_cr[3 * NARRAY];

	// the user defined discrete spectral energy distribution used for noise control
	double spectralEnergyDistribution[NDIR][MAX_FREQ];

public:
	Viewer();
	void init_ogl() override;
	void draw_ogl() override;
	void interface_ogl() override;

	// void iso3dangles();
	void createIsotropicProceduralEnergyDistri();
	void precomputePlanarWave(float scale);
	double* MakeSpatialWaveProfile(int pow_2);
	void precomputePlanarWaveFromFFT1D(double* A, int N, int pow_2);
};

// MAIN c'est juste creer un viewer
int main(int, char**)
{
	Viewer v;

	GLFWwindow* window = v.window();
	glfwSetWindowTitle( window, "Multi-Dimensional Procedural Wave Noise" );
	v.set_size(1800, 1000);

	return v.launch3d();
}

Viewer::Viewer() : nbMeshParts(0)
{
	tX = 0.0;
	tY = 0.0;
	tZ = 0.0;
	tV = 0.0;
	Ndir = 40;
	Nc = 40;
	Na = 0;
	Orient = 1.0f;
	Period = 0.0;
	Zoom = 0.2f;
	Time = 0.0f;
	item_current = 0; // Gaussian
	Oper = 0;		  // Isowave
	old_item = 0;
	Ratio = 64.0f;
	iRatio = 6;
	complex_current = 0;
	Power = 25.0;
	old_power = 1.0;
	contrast = 0.5;
	Proba = 0.5;
	NRec = 3;
	Anisodd = 0.5;
	// iso3dangles();
	createIsotropicProceduralEnergyDistri();
	precomputePlanarWave(4.0);
}

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
	tex->alloc(NARRAY, GL_RGB8, fs_cr);

	texd = Texture1D::create();
	// tex->update(0, N, wprofile);
	texd->alloc(NARRAY, GL_RGB8, fsd_cr);
}

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

	FREQ_LOW = (int)(Ffreq_low * 64.0);
	FREQ_HIGH = (int)(Ffreq_high * 64.0);
	prg_p->bind();
	tex->bind(0);
	texd->bind(1);
	if (FREQ_LOW != sFREQ_LOW || FREQ_HIGH != sFREQ_HIGH || item_current != old_item ||
		(item_current >= 4 && old_power != Power))
	{
		sFREQ_LOW = FREQ_LOW;
		sFREQ_HIGH = FREQ_HIGH;
		old_item = item_current;
		old_power = Power;
		if (item_current >= 4 && item_current <= 9)
		{
			int pow_2 = 1, N = 1;
			while (N < NARRAY)
			{
				pow_2++;
				N *= 2;
			}
			pow_2--;
			printf("NARRAY=%d, N=%d, pow2=%d\n", NARRAY, N, pow_2);
			double* A = MakeSpatialWaveProfile(pow_2);
			precomputePlanarWaveFromFFT1D(A, N, pow_2);
			delete A;
		}
		else
		{
			createIsotropicProceduralEnergyDistri();
			precomputePlanarWave(4.0);
		}
		tex->alloc(NARRAY, GL_RGB8, fs_cr);
		texd->alloc(NARRAY, GL_RGB8, fsd_cr);
	}
	set_uniform_value(1, proj);
	set_uniform_value(2, mv);
	set_uniform_value(3, Transfo::inverse_transpose(mv));

	set_uniform_value(14, GLVec3(0, 2.0, 3.0));

	set_uniform_value(25, tX+5.0);
	set_uniform_value(26, tY+5.0);
	set_uniform_value(27, tZ+5.0);
	set_uniform_value(15, Ndir);
	set_uniform_value(16, Anisodd);
	Period = 1.0 - Orient;
	// set_uniform_value(17, Period);
	set_uniform_value(18, Ndir);
	// set_uniform_value(20, Na);
	set_uniform_value(20, tV);
	set_uniform_value(21, Time);
	Ratio = (float)pow(2.0, iRatio);
	set_uniform_value(22, Ratio);
	set_uniform_value(19, Zoom * Ratio);
	set_uniform_value(23, complex_current);
	set_uniform_value(24, contrast);
	set_uniform_value(30, Oper);
	set_uniform_value(31, NRec);
	set_uniform_value(32, Proba);

	for (int i = 0; i < renderer_p.size(); i++)
	{
		auto mat = renderer_p[i]->material();

		set_uniform_value(10, mat->Kd);
		set_uniform_value(11, mat->Ka);
		set_uniform_value(12, mat->Ks);
		set_uniform_value(13, mat->Ns);
		set_uniform_value("light_pos", GLVec3(0, 2, 2.5));
		renderer_p[i]->draw(GL_TRIANGLES);
	}
}

void Viewer::interface_ogl()
{
	ImGui::Begin("Gui", nullptr, ImGuiWindowFlags_NoSavedSettings);
	ImGui::SetWindowSize(ImVec2(600, 500));
	ImGui::SetWindowPos(ImVec2(10, 10));
	//ImGui::Text("FPS: %2.2lf", fps_);
	ImGui::SliderFloat("X", &tX, 0.0, 10.0);
	ImGui::SliderFloat("Y", &tY, 0.0, 10.0);
	ImGui::SliderFloat("Z", &tZ, 0.0, 10.0);
	ImGui::SliderFloat("Zoom", &Zoom, 0.01, 2.0);
	ImGui::SliderFloat("Time", &Time, 0.0, 5.0);
	ImGui::SliderFloat("Speed", &tV, 0.0, 0.1);
	ImGui::SliderFloat("Contrast", &contrast, 0.1, 1.0);

	// ImGui::Text("MEM:  %2.2lf %%", 100.0 * mem_);
	ImGui::SliderInt("NDir", &Ndir, 1, 100);
	ImGui::SliderFloat("Slice size", &iRatio, 0.0, 8.0);
	ImGui::SliderFloat("FreqMin", &Ffreq_low, 1.0/64.0, 1.0);
	ImGui::SliderFloat("FreqMax", &Ffreq_high, 1.0 / 64.0, 1.0);
	const char* items[] = { "noise-gaussian",	"noise-white",	 "noise-blue",	  "noise-brown",	   
		"nongauss-crystal1", "nongauss-web", "nongauss-marble", "nongauss-crystal2", "nongauss-scratches", "nongauss-smooth cells", "noise-two ampli levels" };
	ImGui::Combo("Wave type", &item_current, items, IM_ARRAYSIZE(items));
	const char* itemscplx[] = { "real", "imaginary", "modulus", "phasor" };
	ImGui::Combo("Value", &complex_current, itemscplx, IM_ARRAYSIZE(itemscplx));
	const char* itemsop[] = { "Isotropic Sum", "Aniso Sum - one direction", "Aniso Sum - two directions", "Random polytopes", "Cellular", "Hyperplan", "Reversed Cellular" };
	ImGui::Combo("Operator", &Oper, itemsop, IM_ARRAYSIZE(itemsop));
	ImGui::SliderFloat("aniso wave dir width", &Anisodd, 0.0, 1.0);
	ImGui::SliderFloat("nongauss wave sharpness", &Power, 0.2, 50.0);
	ImGui::SliderFloat("STIT Probability", &Proba, 0.0, 1.0);
	ImGui::SliderInt("STIT Recursions", &NRec, 1, 5);

	// ImGui::SetWindowSize({ 0, 0 });
	ImGui::End();
}

void Viewer::createIsotropicProceduralEnergyDistri()
{
	const int FF = 4;
	printf("create isotropic energy distri\n");
	int ii = 0, jj = 0;
	for (ii = 0; ii < MAX_FREQ; ii++)
	{
		double ampli = 0.0;
		if (ii >= FREQ_LOW && ii <= FREQ_HIGH)
		{
			double freq = (double)ii / (double)(NARRAY / 2);
			switch (item_current)
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
			}
		}
		spectralEnergyDistribution[0][ii] = ampli;
	}
}

void Viewer::precomputePlanarWave(float scale)
{
	int ii, k;

	int Nfreq = MAX_FREQ;
	int finit = 0;
	vmax = 0.0;
	double rphases[MAX_FREQ];
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
		fs_cr[3 * ii + 0] = (GLubyte)(fs[ii][0] * 127.0 + 128.0);
		fs[ii][1] /= vmax;
		fs_cr[3 * ii + 1] = (GLubyte)(fs[ii][1] * 127.0 + 128.0);
		fs_cr[3 * ii + 2] = 0;
	}
	// Gradient
	for (ii = 0; ii < NARRAY; ii++)
	{
		fsd_cr[3 * ii + 0] = (GLubyte)((fs[(ii + NARRAY - 1) % NARRAY][0] - fs[(ii + 1) % NARRAY][0]) * 64.0 + 128.0);
		fsd_cr[3 * ii + 1] = (GLubyte)((fs[(ii + NARRAY - 1) % NARRAY][1] - fs[(ii + 1) % NARRAY][1]) * 64.0 + 128.0);
		fsd_cr[3 * ii + 2] = 0;
	}
}

double* Viewer::MakeSpatialWaveProfile(int pow_2)
{
	double rphases[MAX_FREQ];
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

void Viewer::precomputePlanarWaveFromFFT1D(double* A, int N, int pow_2)
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
		fs_cr[3 * ii + 0] = (GLubyte)(fs[ii][0] * 127.0 + 128.0);
		fs[ii][1] /= vmax;
		fs_cr[3 * ii + 1] = (GLubyte)(fs[ii][1] * 127.0 + 128.0);
		fs_cr[3 * ii + 2] = 0;
	}
	// Gradient
	for (ii = 0; ii < NARRAY; ii++)
	{
		fsd_cr[3 * ii + 0] = (GLubyte)((fs[(ii + NARRAY - 1) % NARRAY][0] - fs[(ii + 1) % NARRAY][0]) * 64.0 + 128.0);
		fsd_cr[3 * ii + 1] = (GLubyte)((fs[(ii + NARRAY - 1) % NARRAY][1] - fs[(ii + 1) % NARRAY][1]) * 64.0 + 128.0);
		fsd_cr[3 * ii + 2] = 0;
	}
}
