#ifndef INPAINTINGH
#define INPAINTINGH

#include <math.h>
#include <vector>
#include <list>
#include "CImg-1.3.2/CImg.h"
using namespace std;
using namespace cimg_library;

//Default parameters See (23)
/*#define EPSILON 5
#define KAPPA   25
#define SIGMA   1.4
#define RHO     4*/


#define EPSILON 6.
#define KAPPA   25.
#define SIGMA   2.
#define RHO     4.


#define DELTA   1  // See (20)
#define SCALE   4  // filter size is SCALE*SIGMA

/*---------------------------------------------------------------------------*/
// Basic types

typedef unsigned char uchar;
typedef double FLOAT;

typedef struct coord
{
	int x,y;
	coord() {x = 0; y = 0;}
	coord(int _x,int _y) {x = _x; y = _y;}
} coord;

class ord_coord
{
	public:
	int x,y;
	FLOAT T;
	
	ord_coord(int _x, int _y, FLOAT _T) {x = _x; y = _y; T = _T;}};
inline bool operator < (const ord_coord &p, const ord_coord &q) {return p.T>q.T;}

class vec2
{
	public:
	FLOAT x,y;      // x,y coordinate
	
	vec2() {x = 0.; y = 0.;}
	vec2(FLOAT _x, FLOAT _y) {x = _x; y = _y;}
	vec2(const coord &c) {x = c.x; y = c.y;}
	
	inline FLOAT norm2() const {return x*x+y*y;}
	inline FLOAT norm()  const {return sqrt(norm2());}
	inline vec2  operator - (const vec2 &v) const {return vec2(x-v.x,y-v.y);} 
	inline vec2  operator / (const FLOAT &s) const {return vec2(x/s,y/s);} 
	inline FLOAT operator * (const vec2 &v) const {return x*v.x + y*v.y;} 
};
class mat2
{
	public:
	FLOAT a[2][2];      // x,y coordinate
	
	mat2() {init(0.);}
	mat2(FLOAT v) {init(v);}
	mat2(FLOAT a00, FLOAT a01, FLOAT a10, FLOAT a11) {a[0][0] = a00; a[0][1] = a01; a[1][0] = a10; a[1][1] = a11;}
	
	inline void init(FLOAT v) {a[0][0] = v; a[0][1] = v; a[1][0] = v; a[1][1] = v;}
	
	inline mat2  operator - (const mat2 &m) const {return mat2(a[0][0]-m.a[0][0], a[0][1]-m.a[0][1], a[1][0]-m.a[1][0], a[1][1]-m.a[1][1]);} 
	inline mat2  operator + (const mat2 &m) const {return mat2(a[0][0]+m.a[0][0], a[0][1]+m.a[0][1], a[1][0]+m.a[1][0], a[1][1]+m.a[1][1]);} 	
	inline void  operator += (const mat2 &m) {a[0][0]+=m.a[0][0]; a[0][1]+=m.a[0][1]; a[1][0]+=m.a[1][0]; a[1][1]+=m.a[1][1];} 	
	inline void  operator -= (const mat2 &m) {a[0][0]-=m.a[0][0]; a[0][1]-=m.a[0][1]; a[1][0]-=m.a[1][0]; a[1][1]-=m.a[1][1];} 		
	inline mat2  operator * (const FLOAT &s) const {return mat2(a[0][0]*s, a[0][1]*s, a[1][0]*s, a[1][1]*s);} 	
	inline mat2  operator / (const FLOAT &s) const {return mat2(a[0][0]/s, a[0][1]/s, a[1][0]/s, a[1][1]/s);}
	inline CImg<FLOAT> toImg() const {return CImg<FLOAT>::matrix(a[0][0], a[0][1], a[1][0], a[1][1]);}
};

/*---------------------------------------------------------------------------*/
// Image Tools

// Computes the gradient of I at point c for color channel d
template<typename T>
vec2 gradient(const CImg<T> &I, const coord &c, int d = 0, CImg<bool> *val_ok = NULL);

// Computes column vector (1D) or matrix (2D) for gaussian convolution
CImg<FLOAT> gaussian1D(FLOAT sigma, int size);
CImg<FLOAT> gaussian2D(FLOAT sigma, int size);

// Apply a convolution filter centered on p on I
template<typename T>
T IFilter(const CImg<T> &I, const CImg<FLOAT> &filter, const coord& p, int d, const CImg<bool> *val_ok = NULL);

// Apply gaussian blur to an image
template<typename T>
CImg<T> IGaussianBlur(const CImg<T> &I, FLOAT sigma, const CImg<bool> *val_ok = NULL, const CImg<int> *do_compute = NULL);

// Divide pixels of I1 by those of I2
template<typename T>
CImg<T> IGaussianBlur(const CImg<T> &I1, const CImg<T> &I2);

/*---------------------------------------------------------------------------*/
// Weight classes:
// A weight class should have following methods:
// - a constructor of type Weight(const CImg<FLOAT> &img_source, const Img<FLOAT> &distance)
// - FLOAT operator() (const coord& p, const coord& q);   	// computes the contribution w(p,q) of pixel q to the total weight centered in p
// - void post_painted() (const coord& p) ; // signal that pixel p has been painted to proceed to eventuel updates
class TeleaWeight
{
	private:
	const CImg<FLOAT> &dist;
	CImg<vec2>  grad;	
	CImg<bool>  grad_ok;		
	
	public:
	TeleaWeight(const CImg<FLOAT>&, const CImg<FLOAT> &distance);
	FLOAT operator() (const coord& p, const coord& q);
	void pre_painting(const coord& p);
	void post_painting(const coord& p);
};

class BornemannWeight
{
  public:
//	private:
	const CImg<FLOAT> &I;
	const CImg<FLOAT> &dist;
	
	CImg<bool> _0,_1;	CImg<FLOAT> _1sigma, _1rho;
	
	CImg<FLOAT> IBlur;
	CImg<FLOAT> Vsigma;

	CImg<mat2>  Tensor;
	CImg<mat2>  subJ;
	
	CImg<FLOAT> Grho;
	CImg<FLOAT> Gsigma;
	 
  vec2  c;
  FLOAT mu;
	
  static CImg<int> comp_nb(CImg<int> do_comp, int margin);
  	
	template<typename T>
	static void bluring_update(const CImg<FLOAT>& gaussian_filter, const T& add_val, CImg<T>& map, const coord& p, int z=0, CImg<bool> *do_compute = NULL);
	
	public:
	BornemannWeight(const CImg<FLOAT>& img_source, const CImg<FLOAT>& distance);
	FLOAT operator() (const coord& p, const coord& q);
	
	void pre_painting(const coord& p);
	void post_painting(const coord& p);
};
/*---------------------------------------------------------------------------*/
// Inpainting functions

// Return a mask indicating the region to inpaint:  0 if should not be inpainted, 1 otherwise
CImg<FLOAT> get_inpainting_mask(const CImg<uchar> &I, uchar R, uchar G, uchar B);  

// Modify the mask: pixel value = boundary distance,  0 if outside the inpainting zone.
// Return the pixels to inpaint in inpainting order
vector<coord> fast_march_ordering(CImg<FLOAT> &dist); 

// Return the list of neighbours of p at euclian distance max_dist
list<coord> get_neighbours(CImg<FLOAT> &dist, const coord &p, FLOAT max_dist);   

// Head function to call for inpainting, works with color images
// Use parameters R, G and B to indicate the inpainting color
template <class Weight>
CImg<FLOAT> inpaint(const CImg<uchar> &I, uchar R, uchar G, uchar B);  


#endif
