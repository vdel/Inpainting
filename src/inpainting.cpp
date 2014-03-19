#include <limits>
#include <queue>
#include <iostream>
#include <sys/time.h>
#include <time.h>
#include "inpainting.h"
#include "progressbar.h"

#define ForNeigbours(_p,_r,_x,_w,_y,_h) \
      	  for(int _x=(_p.x-(_r)<0?0:_p.x-(_r)); _x<=_p.x+(_r) && _x<_w; _x++) \
	  	      for(int _y=(_p.y-(_r)<0?0:_p.y-(_r)); _y<=_p.y+(_r) && _y<_h; _y++)
	  	      
double get_time()
{
  struct timeval t;
  gettimeofday(&t,NULL); 
  return double(t.tv_sec)+double(t.tv_usec)*0.000001;  
}

int main(int argc, char** argv)
{
  if(argc != 5)
  {
    printf("The usage is ./inpaint filename R V B\n");
    return 1;
  }
  char filename[200];
  double d,f;
  CImg<unsigned char> image(argv[1]);
  CImgDisplay main_disp(image,"Original Picture");
  
  d = get_time();
  CImg<unsigned char> inpainted_telea     = inpaint<TeleaWeight>(image, atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));   
  f = get_time();  
  sprintf(filename,"telea_%s",argv[1]);
  inpainted_telea.save(filename);
  CImgDisplay visu_disp1(inpainted_telea,"Picture inpainted by Telea algorithm");  
  printf("Telea's algorithm : Computing Time = %.2fs\n",f-d);

  d = get_time();
  CImg<unsigned char> inpainted_bornemann = inpaint<BornemannWeight>(image, atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));   
  f = get_time();
  sprintf(filename,"bornemann_%.1f_%.1f_%.1f_%s",EPSILON,SIGMA,RHO,argv[1]);
  inpainted_bornemann.save(filename);
  CImgDisplay visu_disp2(inpainted_bornemann,"Picture inpainted by Bornemann algorithm");
  printf("Bornemann's algorithm : Computing Time = %.2fs\n",f-d);
  
  while (!main_disp.is_closed() && !visu_disp1.is_closed() && !visu_disp2.is_closed()) 
    main_disp.wait();

  return 0;
}


/*---------------------------------------------------------------------------*/ 
template <class Weight>
CImg<FLOAT> inpaint(const CImg<uchar> &I, uchar R, uchar G, uchar B)  
{	
	CImg<FLOAT> dist = get_inpainting_mask(I, R, G, B);
	vector<coord> p = fast_march_ordering(dist);
	CImg<FLOAT> IR = +I;
	CImg<FLOAT> U = I.get_RGBtoYUV();
	CImg<FLOAT> YUVpix(1,1,1,IR.spectrum());
	CImg<FLOAT> RGBpix(1,1,1,IR.spectrum());
	Weight w(IR,dist);
  
  list<coord> NB;
	FLOAT sum_weights;
	FLOAT *sum_channels = new FLOAT[IR.spectrum()];
	
	printf("Inpainting...    ");
	ProgressBar pgbar(0,p.size());
  for(unsigned int k=0;k<p.size();k++)
  {
	  w.pre_painting(p[k]);
	  
	  NB = get_neighbours(dist, p[k], EPSILON);
	  
	  sum_weights = 0.;
	  for(int d=0; d<IR.spectrum(); d++)
		  sum_channels[d] = 0.;

	  for(list<coord>::iterator q=NB.begin(); q!=NB.end(); q++)
	  {
		  FLOAT weight = w(p[k],*q);   // see equation (2)
		  sum_weights += weight;
		  for(int d=0; d<IR.spectrum(); d++)
			  sum_channels[d] += weight*U(q->x,q->y,0,d);
	  }		

	  for(int d=0; d<IR.spectrum(); d++)
	  {
	    YUVpix(0,0,0,d) = sum_channels[d]/sum_weights;	
		  U(p[k].x,p[k].y,0,d) = YUVpix(0,0,0,d);
	  }
	  RGBpix = YUVpix.get_YUVtoRGB();
	  for(int d=0; d<IR.spectrum(); d++)
	    IR(p[k].x,p[k].y,0,d) = RGBpix(0,0,0,d);

	  w.post_painting(p[k]);
	
	  pgbar.update(k+1);
  }
  pgbar.done();
	
	delete[] sum_channels;
	return IR;
}
/*---------------------------------------------------------------------------*/
// Return a mask indicating the region to inpaint:  0 if should not be inpainted, 1 otherwise
CImg<FLOAT> get_inpainting_mask(const CImg<uchar> &I, uchar R, uchar G, uchar B)  
{
  printf("Computing mask...\n");
	CImg<FLOAT> mask(I.width(),I.height());
	for(int x=0; x<I.width(); x++)
		for(int y=0; y<I.height(); y++)
			mask(x,y) = (I(x,y,0,0)==R && I(x,y,0,1)==G && I(x,y,0,2)==B)?1.:0.;
	return mask;
}

/*---------------------------------------------------------------------------*/
	
// Modify the mask: pixel value = boundary distance,  0 if outside the inpainting zone.
// Return the pixels to inpaint in inpainting order
enum {KNOWN, BAND, INSIDE};
vector<coord> fast_march_ordering(CImg<FLOAT> &dist)
{	
  printf("Computing ordering...\n");
	const int w = dist.width();
	const int h = dist.height();
	const int dx[] = {-1,0,1,0};
	const int dy[] = {0,-1,0,1};
	
	FLOAT v1,v2;	
	vector<coord> order;
	CImg<int> status(dist.width(), dist.height());
	priority_queue<ord_coord> band;

	// Initialization
	for(int x=0; x<w; x++)
		for(int y=0; y<h; y++)
			if(dist(x,y) == 0.)
			{	
	      bool nb = false;
				     if(x>0 && dist(x-1,y) != 0.)   nb=true;
				else if(y>0 && dist(x,y-1) != 0.)   nb=true;
				else if(x<w-1 && dist(x+1,y) != 0.) nb=true;
				else if(y<h-1 && dist(x,y+1) != 0.) nb=true;
				if(nb)
					band.push(ord_coord(x,y,0.));
				status(x,y) = KNOWN;		
			}
			else
			{
					status(x,y) = INSIDE;
					dist(x,y) = numeric_limits<FLOAT>::infinity();			
			}
	
	while(!band.empty())
	{
		ord_coord p = band.top();
		band.pop();
		if(status(p.x,p.y) != KNOWN)
     	order.push_back(coord(p.x,p.y));
		status(p.x,p.y) = KNOWN;
		for(int k=0;k<4;k++)
		{
			const int u = p.x+dx[k];
			const int v = p.y+dy[k];
			if(u<0 || v<0 || u>=w || v>=h) continue;
			if(status(u,v) == INSIDE)
			{
   			// See http://people.math.jussieu.fr/~chalons/Propagation_Fronts/slides_cours_fast_marching.pdf
   			if(u == 0)
   			  v1 = dist(u+1,v);
   			else if(u == w-1)
   			  v1 = dist(u-1,v);
   			else   			
  			  v1 = dist(u-1,v)<dist(u+1,v)?dist(u-1,v):dist(u+1,v);
   			if(v == 0)
   			  v2 = dist(u,v+1);
   			else if(v == h-1)
   			  v2 = dist(u,v-1);   			
   			else  			
			    v2 = dist(u,v-1)<dist(u,v+1)?dist(u,v-1):dist(u,v+1);
			  if(fabs(v1-v2)<1)
			    dist(u,v) = (v1+v2+sqrt(2.-(v2-v1)*(v2-v1)))/2.;
        else
          dist(u,v) = (v1<v2?v1:v2)+1.;
  			band.push(ord_coord(u,v,dist(u,v)));				
				status(u,v) = BAND;
			}
		}
	}
	return order;
}

/*---------------------------------------------------------------------------*/
// Return the list of neighbours of p at euclian distance max_dist
list<coord> get_neighbours(CImg<FLOAT> &dist, const coord &p, FLOAT max_dist)
{
	const int w = dist.width();
	const int h = dist.height();
	const FLOAT max = dist(p.x,p.y);
  int d = max_dist;
	max_dist *= max_dist;  // square dist
	list<coord> Binf;
	int xi = p.x<=d?0:p.x-d;
	int xf = p.x+d>=w?w-1:p.x+d;
	int yi = p.y<=d?0:p.y-d;
	int yf = p.y+d>=h?h-1:p.y+d;
		
	for(int x=xi; x<=xf; x++)
		for(int y=yi; y<=yf; y++)
		{
			vec2 v(x-p.x,y-p.y);
			if(v.norm2()<= max_dist && dist(x,y)<max)
				Binf.push_back(coord(x,y));
		}
	return Binf;
}

/*---------------------------------------------------------------------------*/
// Image Tools
/*---------------------------------------------------------------------------*/

// Computes the gradient of I at point c for color channel d
template<typename T>
vec2 gradient(const CImg<T> &I, const coord &c, int d, CImg<bool> *val_ok)
{
	FLOAT sx = 2.;
	FLOAT sy = 2.;
	int xi = c.x-1;
	int xf = c.x+1;
	int yi = c.y-1;
	int yf = c.y+1;
	if(xi<0           || (val_ok && !(*val_ok)(xi,c.y))) {xi++; sx = 1.;}
	if(xf>=I.width()  || (val_ok && !(*val_ok)(xf,c.y))) {xf--; sx = 1.;}
	if(yi<0           || (val_ok && !(*val_ok)(c.x,yi))) {yi++; sy = 1.;}
	if(yf>=I.height() || (val_ok && !(*val_ok)(c.x,yf))) {yf--; sy = 1.;}	
	
	FLOAT dx = (I(xf,c.y,0,d)-I(xi,c.y,0,d))/sx;
	FLOAT dy = (I(c.x,yf,0,d)-I(c.x,yi,0,d))/sy;
		
	return vec2(dx,dy);
}

/*---------------------------------------------------------------------------*/
// Computes line vector (1D) or matrix (2D) for gaussian convolution
CImg<FLOAT> gaussian1D(FLOAT sigma, int size)
{
	if(size%2==0) size++;
	int mid = size/2;
	FLOAT norm = 1./(2.*M_PI*sigma*sigma);
	FLOAT t = -2.*sigma*sigma;
	CImg<FLOAT>	G1D(size,1);
	for(int x=0; x<size; x++)
		G1D(x) = norm * exp(pow(x-mid,2)/t);
	return G1D;
}
CImg<FLOAT> gaussian2D(FLOAT sigma, int size)
{
	if(size%2==0) size++;
	int mid = size/2;
	FLOAT norm = 1./(sigma*sqrt(2.*M_PI));
	norm *= norm;
	FLOAT t = -2.*sigma*sigma;
	CImg<FLOAT>	G2D(size,size);
	for(int x=0; x<size; x++)
		for(int y=0; y<size; y++)
			G2D(x,y) = norm * exp(pow(x-mid,2)/t+pow(y-mid,2)/t);		
	return G2D;
}
/*---------------------------------------------------------------------------*/
// Apply a convolution filter centered on p on I
template<typename T>
T IFilter(const CImg<T> &I, const CImg<FLOAT> &filter, const coord& p, int d, const CImg<bool> *val_ok)
{
	const int w = I.width();
	const int h = I.height();
	int dx = filter.width()/2;
	int dy = filter.height()/2;
	int cx = p.x-dx;
	int cy = p.y-dy;
	int xi = p.x<=dx?0:cx;
	int xf = p.x+dx>=w?w-1:p.x+dx;
	int yi = p.y<=dy?0:cy;
	int yf = p.y+dy>=h?h-1:p.y+dy;
	T sum(0);
	
	if(!val_ok)
	{
	  for(int x=xi; x<=xf; x++)
		  for(int y=yi; y<=yf; y++)   		  sum += I(x,y,0,d)*filter(x-cx,y-cy);
	}
	else
	{
	  for(int x=xi; x<=xf; x++)
		  for(int y=yi; y<=yf; y++)
			  if((*val_ok)(x,y))
    		  sum += I(x,y,0,d)*filter(x-cx,y-cy);
	}

	return sum;
}

/*---------------------------------------------------------------------------*/
// Apply gaussian blur to an image
template<typename T>
CImg<T> IGaussianBlur(const CImg<T> &I, FLOAT sigma, const CImg<bool> *val_ok, const CImg<int> *do_compute)
{
	const int w = I.width();
	const int h = I.height();
	const int d = I.spectrum();
	CImg<T> tmp(w, h);
	CImg<T> res(w, h, 1, d);
	CImg<FLOAT> G = gaussian1D(sigma,SCALE*sigma);
	if(!do_compute)
	{
	  for(int z=0; z<d; z++)
	  {
		  for(int x=0; x<w; x++)
			  for(int y=0; y<h; y++)
 				  tmp(x,y) = IFilter<T>(I,G,coord(x,y),z,val_ok);
		  for(int x=0; x<w; x++)
			  for(int y=0; y<h; y++)   				res(x,y,0,z) = IFilter<T>(tmp,G.get_transpose(),coord(x,y),0,NULL);
    }
	}
	else
	{
	  for(int z=0; z<d; z++)
	  {
		  for(int x=0; x<w; x++)
			  for(int y=0; y<h; y++)
			    if((*do_compute)(x,y))
   				  tmp(x,y) = IFilter<T>(I,G,coord(x,y),z,val_ok);
		  for(int x=0; x<w; x++)
			  for(int y=0; y<h; y++)
			    if((*do_compute)(x,y))     				res(x,y,0,z) = IFilter<T>(tmp,G.get_transpose(),coord(x,y),0,NULL);	
     			else
     			  res(x,y,0,z) = T(0);
   	}
	}
	return res;
}

/*---------------------------------------------------------------------------*/
// Telea Weight
/*---------------------------------------------------------------------------*/
TeleaWeight::TeleaWeight(const CImg<FLOAT>&, const CImg<FLOAT> &distance):dist(distance)
{
  grad = CImg<vec2>(dist.width(),dist.height());
  grad_ok = CImg<bool>(dist.width(),dist.height());  
  grad_ok = false;
}
/*---------------------------------------------------------------------------*/
FLOAT TeleaWeight::operator() (const coord& p, const coord& q)
{
  if(!grad_ok(p.x,p.y))
  {
    grad(p.x,p.y) = gradient<FLOAT>(dist,p);
    grad_ok(p.x,p.y) = true;
  }
  vec2 v = vec2(p)-vec2(q);
	return (fabs(grad(p.x,p.y)*v)+0.000001)/v.norm2();
}
/*---------------------------------------------------------------------------*/
void TeleaWeight::pre_painting(const coord&)
{
}
/*---------------------------------------------------------------------------*/
void TeleaWeight::post_painting(const coord&) 
{
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
// Bornemann Weight
/*---------------------------------------------------------------------------*/
BornemannWeight::BornemannWeight(const CImg<FLOAT> &img_source, const CImg<FLOAT> &distance):I(img_source),dist(distance)
{
  printf("Initializing tensors...\n");
	const int w = I.width();
	const int h = I.height();
	const int d = I.spectrum();

	Grho      = gaussian2D(RHO, SCALE*RHO);
	Gsigma    = gaussian2D(SIGMA, SCALE*SIGMA);
	
	_0 = CImg<bool>(w, h);
	_1 = CImg<bool>(w, h);	
	for(int x=0; x<w; x++)
		for(int y=0; y<h; y++)
		{
			_1(x,y) = dist(x,y) == 0.;
			_0(x,y) = dist(x,y) != 0.;
		}
		
  CImg<int> nb_rho = comp_nb(_0,Grho.width()-1); 
  CImg<int> nb_sigma = comp_nb(nb_rho,Gsigma.width()-1);
	
	_1sigma = _1;
  _1rho   = _1;
	_1sigma = IGaussianBlur<FLOAT>(_1sigma, SIGMA, &_1, &nb_sigma);
	_1rho   = IGaussianBlur<FLOAT>(_1rho,   RHO,   &_1, &nb_rho);
	IBlur   = IGaussianBlur<FLOAT>(I,  SIGMA, &_1, &nb_sigma);
	
	Vsigma  = CImg<FLOAT>(w, h, 1, d);
	for(int x=0; x<w; x++)
		for(int y=0; y<h; y++)
		  if(_1(x,y) && nb_rho(x,y))		
  			for(int z=0; z<d; z++)
          Vsigma(x,y,0,z) = IBlur(x,y,0,z) / _1sigma(x,y); 
    
	vec2 grad;
	Tensor  = CImg<mat2>(w, h, 1, d);
	for(int x=0; x<w; x++)
		for(int y=0; y<h; y++)
  	  if(_1(x,y) && nb_rho(x,y))
	  		for(int z=0; z<d; z++)
  			{
	  			grad = gradient<FLOAT>(Vsigma, coord(x,y), z, &_1);
	  			Tensor(x,y,0,z).a[0][0] = grad.x * grad.x;
	  			Tensor(x,y,0,z).a[1][0] = Tensor(x,y,0,z).a[0][1] = grad.x * grad.y;
	  			Tensor(x,y,0,z).a[1][1] = grad.y * grad.y;
	  		}
				
  printf("Bluring tensors...\n");		
	subJ = IGaussianBlur<mat2>(Tensor, RHO, &_1, &nb_rho);
}
/*---------------------------------------------------------------------------*/
CImg<int> BornemannWeight::comp_nb(CImg<int> do_comp, int margin)
{
	const int w = do_comp.width();
	const int h = do_comp.height();
	int *b = new int[2*margin+1];
	int d;
	int k;
	  
  CImg<int> tmp(w,h);
  for(int x=0; x<w; x++)
  {
    d = 0;
    for(k=0;k<margin && k<h;k++)
    {
      b[k] = do_comp(x,k)?1:0;
      d += b[k];      
    }
    for(; k<2*margin+1; k++)    
      b[k] = 0;
    k = margin;
    for(int y=0; y<h; y++, k=(k+1)%(2*margin+1))
    {
      d -= b[k];
      if(y+margin<h)
      {
        b[k] = do_comp(x,y+margin)?1:0;
        d += b[k];
      }
      tmp(x,y) = d;
    }
  }
  
  CImg<int> res(w,h);
  for(int y=0; y<h; y++)
  {
    d = 0;
    for(k=0;k<margin && k<w;k++)
    {
      b[k] = tmp(k,y);
      d += b[k];      
    }
    for(; k<2*margin+1; k++)    
      b[k] = 0;
    k = margin;
    for(int x=0; x<w; x++, k=(k+1)%(2*margin+1))
    {
      d -= b[k];
      if(x+margin<w)
      {
        b[k] = tmp(x+margin,y);
        d += b[k];
      }
      res(x,y) = d;
    }
  }
  
  delete[] b;
  return res;
}
/*---------------------------------------------------------------------------*/
FLOAT BornemannWeight::operator() (const coord& p, const coord& q)
{
	vec2 v(p.x-q.x, p.y-q.y);
	return sqrt(M_PI_2/v.norm2())*mu*exp(-pow(mu*(c*v)/(FLOAT)EPSILON,2.)/2.);
}
/*---------------------------------------------------------------------------*/
void BornemannWeight::pre_painting(const coord& p)
{
  mat2 J = subJ(p.x,p.y,0,0)*(0.299/_1rho(p.x,p.y))+
           subJ(p.x,p.y,0,1)*(0.587/_1rho(p.x,p.y))+
           subJ(p.x,p.y,0,2)*(0.114/_1rho(p.x,p.y));
				 
	CImg<FLOAT> val;
	CImg<FLOAT> vec;
	J.toImg().eigen(val,vec);

	FLOAT L1 = val(0,0);
	FLOAT L2 = val(0,1);  
	int i = L1<L2?0:1;
	
	mu = L1==L2?1.:(1 + KAPPA * exp( -DELTA*DELTA*DELTA*DELTA / pow(L1-L2,2) ));

	c = vec2(-vec(i,1), vec(i,0));  // take the orthogonal
	c = c / sqrt(c.x*c.x+c.y*c.y); 	
}
/*---------------------------------------------------------------------------*/
void BornemannWeight::post_painting(const coord& p) 
{
	const int w = I.width();
	const int h = I.height();
	const int d = I.spectrum();

	_1(p.x,p.y) = true;
	_0(p.x,p.y) = false;	

	bluring_update(Grho,   1., _1rho,   p);
	bluring_update(Gsigma, 1., _1sigma, p);
	for(int z=0; z<d; z++)	
	{
  	bluring_update(Gsigma, I(p.x,p.y,0,z), IBlur, p, z);	
    Vsigma(p.x,p.y,0,z) = IBlur(p.x,p.y,0,z) / _1sigma(p.x,p.y);
    if(0<=p.x-1) Vsigma(p.x-1,p.y,0,z) = IBlur(p.x-1,p.y,0,z) / _1sigma(p.x-1,p.y);
    if(p.x+1<w) Vsigma(p.x+1,p.y,0,z) = IBlur(p.x+1,p.y,0,z) / _1sigma(p.x+1,p.y);
    if(0<=p.y-1) Vsigma(p.x,p.y-1,0,z) = IBlur(p.x,p.y-1,0,z) / _1sigma(p.x,p.y-1);
    if(p.y+1<h) Vsigma(p.x,p.y+1,0,z) = IBlur(p.x,p.y+1,0,z) / _1sigma(p.x,p.y+1);
  }
        
  // Approximation with only the new tensor
	vec2 grad;
  for(int z=0; z<d; z++)
  {
	  grad = gradient<FLOAT>(Vsigma, p, z, &_1);			
	  Tensor(p.x,p.y,0,z).a[0][0] = grad.x * grad.x;
	  Tensor(p.x,p.y,0,z).a[0][1] = Tensor(p.x,p.y,0,z).a[1][0] = grad.x * grad.y;
	  Tensor(p.x,p.y,0,z).a[1][1] = grad.y * grad.y;				
    bluring_update(Grho, Tensor(p.x,p.y,0,z), subJ, p, z, &_0);
  }
}
/*---------------------------------------------------------------------------*/
template<typename T>
// Special function for applying local modification of a blured image.
// Use the symetries of the gaussian filter to avoid recomputing 4 times values to add
// (useful if values are matrix).
void BornemannWeight::bluring_update(const CImg<FLOAT>& gaussian_filter, const T& add_val, CImg<T>& map, const coord& p, int z, CImg<bool> *do_compute)
{
	const int w = map.width();
	const int h = map.height();
	const int r = gaussian_filter.width() / 2;
	const int dx = p.x-r;
	const int dy = p.y-r;	
	
	ForNeigbours(p,r,x,w,y,h) 
	  if(!do_compute || (*do_compute)(x,y))
      map(x,y,0,z) += add_val*gaussian_filter(x-dx,y-dy);
}
/*---------------------------------------------------------------------------*/
