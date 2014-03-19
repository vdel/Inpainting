/*
 #
 #  File        : gmic.cpp
 #                ( C++ source file )
 #
 #  Description : GREYC's Magic Image Converter - G'MIC Language interpreter
 #                ( http://gmic.sourceforge.net )
 #                This file is a part of the CImg Library project.
 #                ( http://cimg.sourceforge.net )
 #
 #  Copyright   : David Tschumperle
 #                ( http://www.greyc.ensicaen.fr/~dtschump/ )
 #
 #  License     : CeCILL v2.0
 #                ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
 #
 #  This software is governed by the CeCILL  license under French law and
 #  abiding by the rules of distribution of free software.  You can  use,
 #  modify and/ or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and  rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
 #
*/

// Add specific G'MIC methods to the CImg<T> class.
//-------------------------------------------------
#ifdef cimg_plugin

template<typename t>
CImg<T>& replace(CImg<t>& img) {
  return img.move_to(*this);
}

template<typename t>
CImg<T> get_replace(const CImg<t>& img) const {
  return +img;
}

CImg<T> get_gmic_set(const double value, const int x, const int y, const int z, const int v) const {
  return (+*this).gmic_set(value,x,y,z,v);
}

CImg<T>& gmic_set(const double value, const int x, const int y, const int z, const int v) {
  (*this).atXYZC(x,y,z,v,0) = (T)value;
  return *this;
}

CImg<T> get_draw_point(const int x, const int y, const int z, const T *const col, const float opacity) const {
  return (+*this).draw_point(x,y,z,col,opacity);
}

CImg<T> get_draw_line(const int x0, const int y0, const int x1, const int y1, const T *const col, const float opacity) const {
  return (+*this).draw_line(x0,y0,x1,y1,col,opacity);
}

template<typename t>
CImg<T> get_draw_polygon(const CImg<t>& pts, const T *const col, const float opacity) const {
  return (+*this).draw_polygon(pts,col,opacity);
}

CImg<T> get_draw_spline(const int x0, const int y0, const float u0, const float v0,
                        const int x1, const int y1, const float u1, const float v1,
                        const T *const color, const float opacity, const float precision) const {
  return (+*this).draw_spline(x0,y0,u0,v0,x1,y1,u1,v1,color,opacity,precision);
}

CImg<T> get_draw_ellipse(const int x, const int y, const float r0, const float r1,
                         const float angle, const T *const col, const float opacity) const {
  return (+*this).draw_ellipse(x,y,r0,r1,angle,col,opacity);
}

CImg<T> get_draw_text(const int x, const int y, const char *const text, const T *const col,
                      const int bg, const float opacity,const int siz) const {
  return (+*this).draw_text(x,y,text,col,bg,opacity,siz);
}

CImg<T> get_draw_image(const int x, const int y, const int z,
                       const CImg<T>& sprite, const CImg<T>& mask, const float opacity) const {
  return (+*this).draw_image(x,y,z,sprite,mask,opacity);
}

CImg<T> get_draw_image(const int x, const int y, const int z,
                       const CImg<T>& sprite, const float opacity) const {
  return (+*this).draw_image(x,y,z,sprite,opacity);
}

CImg<T> get_draw_plasma(const float alpha, const float beta, const float opacity) const {
  return (+*this).draw_plasma(alpha,beta,opacity);
}

CImg<T> get_draw_mandelbrot(const CImg<T>& color_palette, const float opacity,
                            const double z0r, const double z0i, const double z1r, const double z1i,
                            const unsigned int itermax, const bool normalized_iteration,
                            const bool julia_set, const double paramr, const double parami) const {
  return (+*this).draw_mandelbrot(color_palette,opacity,z0r,z0i,z1r,z1i,itermax,
                                  normalized_iteration,julia_set,paramr,parami);
}

template<typename t1, typename t2>
CImg<T> get_draw_quiver(const CImg<t1>& flow,
                        const t2 *const color, const float opacity=1,
                        const unsigned int sampling=25, const float factor=-20,
                        const bool arrows=true, const unsigned int pattern=~0U) {
  return (+*this).draw_quiver(flow,color,opacity,sampling,factor,arrows,pattern);
}

CImg<T> get_draw_fill(const int x, const int y, const int z,
                      const T *const col, const float opacity, const float tolerance) const {
  return (+*this).draw_fill(x,y,z,col,opacity,tolerance);
}

static bool is_almost(const T x, const T c) {
  return x>=c && x<c+1;
}

bool is_CImg3d() const {
  const bool is_header = (_width==1 && _height>=8 && _depth==1 && _spectrum==1 &&
                          is_almost((*this)[0],(T)'C') && is_almost((*this)[1],(T)'I') &&
                          is_almost((*this)[2],(T)'m') && is_almost((*this)[3],(T)'g') &&
                          is_almost((*this)[4],(T)'3') && is_almost((*this)[5],(T)'d'));
  if (!is_header) return false;
  const int nbv = (int)(*this)[6], nbp = (int)(*this)[7];
  if (nbv<0 || nbp<0) return false;
  if (!nbv && !nbp) return true;
  const T *ptrs = data() + 8 + 3*nbv, *const ptre = end();
  if (ptrs>=ptre) return false;
  for (int i = 0; i<nbp && ptrs<ptre; ++i) {
    const int N = (int)*(ptrs++);
    if (N<=0 || N>=8) return false;
    ptrs+=N;
  }
  ptrs+=4*nbp;
  if (ptrs>ptre) return false;
  return true;
}

template<typename tp, typename tf, typename tc, typename to>
CImg<T> get_draw_object3d(const float x0, const float y0, const float z0,
                          const CImg<tp>& vertices, const CImgList<tf>& primitives,
                          const CImgList<tc>& colors, const CImg<to>& opacities,
                          const unsigned int render_mode, const bool double_sided, const float focale,
                          const float light_x, const float light_y,const float light_z,
                          const float specular_light, const float specular_shine,
                          CImg<floatT>& zbuffer) const {
  return (+*this).draw_object3d(x0,y0,z0,vertices,primitives,colors,opacities,render_mode,double_sided,focale,
                                light_x,light_y,light_z,specular_light,specular_shine,zbuffer);
}

template<typename tp, typename tc, typename to>
CImg<T>& object3dtoCImg3d(CImgList<tp>& primitives, CImgList<tc>& colors, CImg<to>& opacities) {
  is_object3d(primitives,true,true,"object3dtoCImg3d");
  const unsigned int primitives_size = primitives.size();
  CImgList<floatT> res;
  (CImg<floatT>("CImg3d",1,6)+=0.5f).move_to(res);
  CImg<floatT>::vector((float)_width,(float)primitives.size()).move_to(res);
  if (!is_empty()) transpose().unroll('y').move_to(res);
  cimglist_for(primitives,p) {
    CImg<floatT>::vector((float)primitives[p].size()).move_to(res);
    primitives[p].unroll('y').move_to(res);
  }
  primitives.assign();
  const unsigned int defined_colors = colors.size();
  cimglist_for(colors,c) colors[c].resize(1,3,1,1,-1).move_to(res);
  colors.assign();
  if (defined_colors<primitives_size) CImg<floatT>(1,3*(primitives_size-defined_colors),1,1,200).move_to(res);
  const unsigned int defined_opacities = opacities.size();
  if (opacities) opacities.unroll('y').move_to(res);
  if (defined_opacities<primitives.size()) CImg<floatT>(1,primitives_size-defined_opacities,1,1,1).move_to(res);
  return (res>'y').move_to(*this);
}

template<typename tp, typename tc, typename to>
CImg<T>& CImg3dtoobject3d(CImgList<tp>& primitives, CImgList<tc>& colors, CImg<to>& opacities) {
  const T *ptrs = data() + 6;
  const unsigned int
    nbv = (unsigned int)*(ptrs++),
    nbp = (unsigned int)*(ptrs++);
  CImg<T> vertices(nbv,3);
  cimg_forX(vertices,x) { vertices(x,0) = (T)*(ptrs++); vertices(x,1) = (T)*(ptrs++); vertices(x,2) = (T)*(ptrs++); }
  primitives.assign(nbp);
  cimglist_for(primitives,p) {
    const unsigned int N = (unsigned int)*(ptrs++);
    primitives[p].assign(ptrs,1,N);
    ptrs+=N;
  }
  colors.assign(nbp,1,3,1,1);
  cimglist_for(colors,c) { colors(c,0) = (tc)*(ptrs++); colors(c,1) = (tc)*(ptrs++); colors(c,2) = (tc)*(ptrs++); }
  CImg<T>(ptrs,nbp).move_to(opacities);
  return vertices.move_to(*this);
}

CImg<T> get_appendCImg3d(const CImg<T>& img) const {
  CImg<T> res(1,img.size() + size() - 8);
  const T *ptrs = _data + 6, *ptrs0 = img._data + 6;
  T *ptrd = res._data;
  *(ptrd++) = (T)('C' + 0.5f); *(ptrd++) = (T)('I' + 0.5f);
  *(ptrd++) = (T)('m' + 0.5f); *(ptrd++) = (T)('g' + 0.5f);
  *(ptrd++) = (T)('3' + 0.5f); *(ptrd++) = (T)('d' + 0.5f);
  const unsigned int
    nbv = (unsigned int)*(ptrs++),
    nbv0 = (unsigned int)*(ptrs0++),
    nbp = (unsigned int)*(ptrs++),
    nbp0 = (unsigned int)*(ptrs0++);
  *(ptrd++) = (T)(nbv + nbv0);
  *(ptrd++) = (T)(nbp + nbp0);
  std::memcpy(ptrd,ptrs,sizeof(T)*nbv*3);
  ptrd+=3*nbv; ptrs+=3*nbv;
  std::memcpy(ptrd,ptrs0,sizeof(T)*nbv0*3);
  ptrd+=3*nbv0; ptrs0+=3*nbv0;
  for (unsigned int i = 0; i<nbp; ++i) {
    const unsigned int N = (unsigned int)*(ptrs++);
    *(ptrd++) = (T)N;
    std::memcpy(ptrd,ptrs,sizeof(T)*N);
    ptrd+=N; ptrs+=N;
  }
  for (unsigned int i = 0; i<nbp0; ++i) {
    const unsigned int N = (unsigned int)*(ptrs0++);
    *(ptrd++) = (T)N;
    for (unsigned int j = 0; j<N; ++j) *(ptrd++) = (T)(*(ptrs0++) + nbv);
  }
  std::memcpy(ptrd,ptrs,sizeof(T)*nbp*3);
  ptrd+=3*nbp; ptrs+=3*nbp;
  std::memcpy(ptrd,ptrs0,sizeof(T)*nbp0*3);
  ptrd+=3*nbp0; ptrs0+=3*nbp0;
  std::memcpy(ptrd,ptrs,sizeof(T)*nbp);
  ptrd+=nbp;
  std::memcpy(ptrd,ptrs0,sizeof(T)*nbp0);
  return res;
}

CImg<T>& appendCImg3d(const CImg<T>& img) {
  return get_appendCImg3d(img).move_to(*this);
}

CImg<T>& centerCImg3d() {
  const unsigned int nbv = (unsigned int)(*this)[6];
  const T *ptrs = data() + 8;
  float xm = cimg::type<float>::max(), ym = xm, zm = xm, xM = cimg::type<float>::min(), yM = xM, zM = xM;
  for (unsigned int i = 0; i<nbv; ++i) {
    const float x = (float)*(ptrs++), y = (float)*(ptrs++), z = (float)*(ptrs++);
    if (x<xm) xm = x; if (x>xM) xM = x;
    if (y<ym) ym = y; if (y>yM) yM = y;
    if (z<zm) zm = z; if (z>zM) zM = z;
  }
  const float xc = (xm + xM)/2, yc = (ym + yM)/2, zc = (zm + zM)/2;
  T *ptrd = data() + 8;
  for (unsigned int i = 0; i<nbv; ++i) { *(ptrd++)-=(T)xc; *(ptrd++)-=(T)yc; *(ptrd++)-=(T)zc; }
  return *this;
}

CImg<T> get_centerCImg3d() const {
  return (+*this).centerCImg3d();
}

CImg<T>& normalizeCImg3d() {
  const unsigned int nbv = (unsigned int)(*this)[6];
  const T *ptrs = data() + 8;
  float xm = cimg::type<float>::max(), ym = xm, zm = xm, xM = cimg::type<float>::min(), yM = xM, zM = xM;
  for (unsigned int i = 0; i<nbv; ++i) {
    const float x = (float)*(ptrs++), y = (float)*(ptrs++), z = (float)*(ptrs++);
    if (x<xm) xm = x; if (x>xM) xM = x;
    if (y<ym) ym = y; if (y>yM) yM = y;
    if (z<zm) zm = z; if (z>zM) zM = z;
  }
  const float delta = cimg::max(xM-xm,yM-ym,zM-zm);
  if (delta>0) {
    T *ptrd = data() + 8;
    for (unsigned int i = 0; i<3*nbv; ++i) *(ptrd++)/=(T)delta;
  }
  return *this;
}

CImg<T> get_normalizeCImg3d() const {
  return (+*this).normalizeCImg3d();
}

template<typename t>
CImg<T>& rotateCImg3d(const CImg<t>& rot) {
  const unsigned int nbv = (unsigned int)(*this)[6];
  const T *ptrs = data() + 8;
  const float
    a = (float)rot(0,0), b = (float)rot(1,0), c = (float)rot(2,0),
    d = (float)rot(0,1), e = (float)rot(1,1), f = (float)rot(2,1),
    g = (float)rot(0,2), h = (float)rot(1,2), i = (float)rot(2,2);
  T *ptrd = data() + 8;
  for (unsigned int j = 0; j<nbv; ++j) {
    const float x = (float)*(ptrs++), y = (float)*(ptrs++), z = (float)*(ptrs++);
    *(ptrd++) = (T)(a*x + b*y + c*z);
    *(ptrd++) = (T)(d*x + e*y + f*z);
    *(ptrd++) = (T)(g*x + h*y + i*z);
  }
  return *this;
}

template<typename t>
CImg<T> get_rotateCImg3d(const CImg<t>& rot) const {
  return (+*this).rotateCImg3d(rot);
}

CImg<T>& shiftCImg3d(const float tx, const float ty, const float tz) {
  const unsigned int nbv = (unsigned int)(*this)[6];
  T *ptrd = data() + 8;
  for (unsigned int j = 0; j<nbv; ++j) { *(ptrd++)+=(T)tx; *(ptrd++)+=(T)ty; *(ptrd++)+=(T)tz; }
  return *this;
}

CImg<T> get_shiftCImg3d(const float tx, const float ty, const float tz) const {
  return (+*this).shiftCImg3d(tx,ty,tz);
}

CImg<T>& coloropacityCImg3d(const float R, const float G, const float B,
                            const float opacity, const bool set_RGB, const bool set_opacity) {
  T *ptrd = data() + 6;
  const unsigned int
    nbv = (unsigned int)*(ptrd++),
    nbp = (unsigned int)*(ptrd++);
  ptrd+=3*nbv;
  for (unsigned int i = 0; i<nbp; ++i) { const unsigned int N = (unsigned int)*(ptrd++); ptrd+=N; }
  if (set_RGB) for (unsigned int c = 0; c<nbp; ++c) { *(ptrd++) = (T)R; *(ptrd++) = (T)G; *(ptrd++) = (T)B; } else ptrd+=3*nbp;
  if (set_opacity) for (unsigned int o = 0; o<nbp; ++o) *(ptrd++) = (T)opacity;
  return *this;
}

CImg<T> get_coloropacityCImg3d(const float R, const float G, const float B,
                               const float opacity, const bool set_RGB, const bool set_opacity) const {
  return (+*this).coloropacityCImg3d(R,G,B,opacity,set_RGB,set_opacity);
}

template<typename t>
CImg<T>& inpaint(CImg<t>& mask) {
  if (!is_sameXYZ(mask))
    throw CImgArgumentException("CImg<%s>::inpaint() : Invalid mask (%u,%u,%u,%u,%p) for instance image (%u,%u,%u,%u,%p).",
                                pixel_type(),mask._width,mask._height,mask._depth,mask._spectrum,mask._data,
                                _width,_height,_depth,_spectrum,_data);
  CImg<t> nmask(mask);
  CImg_3x3(M,t); Mpp = Mnp = Mpn = Mnn = 0;
  CImg_3x3(I,T); Ipp = Inp = Icc = Ipn = Inn = 0;
  bool is_pixel = false;
  do {
    is_pixel = false;
    cimg_forZ(mask,z) cimg_for3x3(mask,x,y,z,0,M,t) if (Mcc && (!Mpc || !Mnc || !Mcp || !Mcn)) {
      is_pixel = true;
      const float wcp = Mcp?0.0f:1.0f, wpc = Mpc?0.0f:1.0f, wnc = Mnc?0.0f:1.0f, wcn = Mcn?0.0f:1.0f, sumw = wcp + wpc + wnc + wcn;
      cimg_forC(*this,k) {
        cimg_get3x3(*this,x,y,z,k,I,T);
        (*this)(x,y,z,k) = (T)((wcp*Icp + wpc*Ipc + wnc*Inc + wcn*Icn)/sumw);
      }
      nmask(x,y,z) = 0;
    }
    mask = nmask;
  } while (is_pixel);
  return *this;
}

template<typename t>
CImg<T> get_inpaint(CImg<t>& mask) const {
  return (+*this).inpaint(mask);
}

CImgList<T> get_split_patch(const unsigned int px, const unsigned int py,
                            const unsigned int pz, const unsigned int pv, const bool borders) const {
  CImgList<T> res;
  const unsigned int
    px1 = px/2, px2 = px - px1 - 1,
    py1 = py/2, py2 = py - py1 - 1,
    pz1 = pz/2, pz2 = pz - pz1 - 1,
    pv1 = pv/2, pv2 = pv - pv1 - 1;
  unsigned int p = 0;
  if (pv) {
    res.assign(_width*_height*_depth*_spectrum);
    cimg_forXYZC(*this,x,y,z,v) get_crop(x-px1,y-py1,z-pz1,v-pv1,x+px2,y+py2,z+pz2,v+pv2,borders).move_to(res[p++]);
  } else if (pz) {
    res.assign(_width*_height*_depth);
    cimg_forXYZ(*this,x,y,z) get_crop(x-px1,y-py1,z-pz1,x+px2,y+py2,z+pz2,borders).move_to(res[p++]);
  } else if (py) {
    res.assign(_width*_height);
    cimg_forXY(*this,x,y) get_crop(x-px1,y-py1,x+px2,y+py2,borders).move_to(res[p++]);
  } else if (px) {
    res.assign(_width);
    cimg_forX(*this,x) get_crop(x-px1,x+px2,borders).move_to(res[p++]);
  }
  return res;
}

CImg<T>& mark_filename() {
  const unsigned int siz = size();
  if (siz<2 || (*this)[siz-2]!='*') {
    CImg<char> filename(siz+1);
    std::memcpy(filename.data(),data(),siz);
    filename[siz-1] = '*';
    filename[siz] = 0;
    filename.move_to(*this);
  }
  return *this;
}

CImg<T> get_mark_filename() const {
  return (+*this).mark_filename();
}

#else  // eq. to #ifndef cimg_plugin

#include "gmic.h"
using namespace cimg_library;
#undef min
#undef max

#if !defined(gmic_main) || !defined(gmic_separate_compilation)

// Define some useful variables and macros.
//------------------------------------------

// Character codes for evaluation braces.
const char brace1 = 26, brace2 = 28;

// End character, should be used as write only.
static char end = 0;

// Return current selection as a selection string.
#define gmic_selection selection2string(selection,filenames,true)

// Check validity of a selected indice.
#define gmic_check_indice(ind) { \
  const int indo = (int)ind; \
  if (ind<0) ind+=images.size(); \
  if (ind<0 || ind>=(int)images.size()) { \
    if (images.size()) error(images,"Command '%s' : Invalid indice '%d' (not in range -%u..%u).", \
                             command_name,indo,images.size(),images.size()-1); \
    else error(images,"Command '%s' : Invalid indice '%d' (no data available).",command_name,indo); \
  } \
}

// Code for having 'get' or 'non-get' versions of G'MIC commands.
#define gmic_apply(instance,function) { \
  unsigned int posi = 0; \
  const bool is_inlist = images.contains(instance,posi); \
  CImg<T> &src = instance; \
  if (patch_w) { \
    CImg<T> tmp, res(src.width(),src.height(),src.depth(),src.spectrum(),0); \
    const unsigned int \
      sx1 = patch_w/2, sx2 = patch_w - sx1 - 1, \
      sy1 = patch_h/2, sy2 = patch_h - sy1 - 1, \
      sz1 = patch_d/2, sz2 = patch_d - sz1 - 1, \
      sv1 = patch_c/2, sv2 = patch_c - sv1 - 1; \
    if (patch_c) { \
      cimg_forXYZC(src,x,y,z,v) { \
        src.get_crop(x-sx1,y-sy1,z-sz1,v-sv1,x+sx2,y+sy2,z+sz2,v+sv2,patch_borders?true:false).function.move_to(tmp); \
        res(x,y,z,v) = tmp(tmp.width()>>1,tmp.height()>>1,tmp.depth()>>1,tmp.spectrum()>>1); \
      } \
    } else if (patch_d) { \
      cimg_forXYZ(src,x,y,z) { \
        src.get_crop(x-sx1,y-sy1,z-sz1,x+sx2,y+sy2,z+sz2,patch_borders?true:false).function.move_to(tmp); \
        for (int vmin = cimg::min(src.spectrum(),tmp.spectrum()), v = 0; v<vmin; ++v) \
          res(x,y,z,v) = tmp(tmp.width()>>1,tmp.height()>>1,tmp.depth()>>1,v); \
      } \
    } else if (patch_h) { \
      cimg_forXY(src,x,y) { \
        src.get_crop(x-sx1,y-sy1,x+sx2,y+sy2,patch_borders?true:false).function.move_to(tmp); \
        for (int vmin = cimg::min(src.spectrum(),tmp.spectrum()), v = 0; v<vmin; ++v) \
          for (int zmin = cimg::min(src.depth(),tmp.depth()), z = 0; z<zmin; ++z) \
            res(x,y,z,v) = tmp(tmp.width()>>1,tmp.height()>>1,z,v); \
      } \
    } else { \
      cimg_forX(src,x) { \
        src.get_crop(x-sx1,x+sx2,patch_borders?true:false).function.move_to(tmp); \
        for (int vmin = cimg::min(src.spectrum(),tmp.spectrum()), v = 0; v<vmin; ++v) \
          for (int zmin = cimg::min(src.depth(),tmp.depth()), z = 0; z<zmin; ++z) \
            for (int ymin = cimg::min(src.height(),tmp.height()), y = 0; y<ymin; ++y) \
              res(x,y,z,v) = tmp(tmp.width()>>1,y,z,v); \
      } \
    } \
    if (get_version) { \
      res.move_to(images); \
      if (is_inlist) filenames.insert(filenames[posi].get_mark_filename()); \
      else CImg<char>("(unnamed)",10).move_to(filenames); \
    } else { res.move_to(instance); filenames[posi].mark_filename(); } \
  } else { \
    if (get_version) { \
      instance.get_##function.move_to(images); \
      if (is_inlist) filenames.insert(filenames[posi].get_mark_filename()); \
      else CImg<char>("(unnamed)",10).move_to(filenames); \
    } else { instance.function; filenames[posi].mark_filename(); } \
  } \
}

// Code for simple commands that has no parameters and act on images.
#define gmic_simple_item(option,function,description) \
  if (!std::strcmp(option,command_name)) { \
    print(images,description,gmic_selection); cimg_foroff(selection,l) gmic_apply(images[selection[l]],function()); \
    continue; \
}

// Code for the type cast command.
#define gmic_cast(pixel_type,st_type) \
  if (!std::strcmp(#pixel_type,argument)) { \
    print(images,"Set pixel data type to '%s'.",#pixel_type); ++position; \
    if (!std::strcmp(st_type,cimg::type<T>::string())) continue; \
    CImgList<pixel_type> casted_images; \
    while (images) { casted_images.insert(images[0]); images.remove(0); } \
    return parse_##pixel_type(command_line,position,casted_images,filenames,dowhiles,repeatdones,locals,initial_call); \
}

// Code for G'MIC arithmetic commands.
#define gmic_arithmetic_item(command1,command2,\
                             function1,description1,arg1_1,arg1_2,value_type1, \
                             function2,description2,arg2_1,arg2_2, \
                             description3,arg3_1,arg3_2, \
                             description4) \
 if (!std::strcmp(command1,command_name) || !std::strcmp(command2,command_name)) { \
   double value = 0; char inds[4096] = { 0 }, sep = 0; \
    if (std::sscanf(argument,"%lf%c",&value,&end)==1) { \
      print(images,description1 ".",arg1_1,arg1_2); \
      cimg_foroff(selection,l) \
       if (get_version) { \
         filenames.insert(filenames[selection[l]].get_mark_filename()); \
         images.insert(images[selection[l]]); images.back().function1((value_type1)value); \
       } else { \
         filenames[selection[l]].mark_filename(); \
         images[selection[l]].function1((value_type1)value); \
       } \
      ++position; \
    } else if (std::sscanf(argument,"[%4095[0-9.eE%+-]%c%c",inds,&sep,&end)==2 && sep==']') { \
      const CImg<unsigned int> ind = selection2cimg(inds,images.size(),command1,false); \
      if (ind.size()!=1) error(images,"Command '%s' : Multiple image%s specified in argument '%s', when only one was expected.", \
                               command1+1,selection2string(ind,filenames,true),argument_text); \
      print(images,description2 ".",arg2_1,arg2_2); \
      const CImg<T> img0 = images[ind[0]]; \
      cimg_foroff(selection,l) \
       if (get_version) { \
         filenames.insert(filenames[selection[l]].get_mark_filename()); \
         images.insert(images[selection[l]]); images.back().function2(img0); \
       } else { \
         filenames[selection[l]].mark_filename(); \
         images[selection[l]].function2(img0); \
       } \
      ++position; \
    } else if (std::sscanf(argument,"'%4095[^']%c%c",inds,&sep,&end)==2 && sep=='\'') { \
      print(images,description3 ".",arg3_1,arg3_2); \
      cimg_foroff(selection,l) \
        if (get_version) { \
          filenames.insert(filenames[selection[l]].get_mark_filename()); \
          images.insert(images[selection[l]]); images.back().function2((const char*)inds); \
        } else { \
          filenames[selection[l]].mark_filename(); \
          images[selection[l]].function2((const char*)inds); \
        } \
      ++position; \
    } else { \
      print(images,description4 ".",gmic_selection); \
      if (images && selection) { \
        if (get_version) { \
          CImg<T> img0 = images[selection[0]]; \
          for (unsigned int siz = selection.size(), l = 1; l<siz; ++l) img0.function2(images[selection[l]]); \
          filenames.insert(filenames[selection[0]].get_mark_filename()); \
          img0.move_to(images); \
        } else for (unsigned int siz = selection.size(), ind0 = selection[0], off = 0, l = 1; l<siz; ++l) { \
          const unsigned int ind = selection[l] - off; \
          filenames[ind0].mark_filename(); images[ind0].function2(images[ind]); \
          images.remove(ind); filenames.remove(ind); \
          ++off; \
        }}} continue; \
}

// Constructors / Destructors.
//----------------------------
#if defined(gmic_float) || !defined(gmic_separate_compilation)
#include "gmic_def.h"

gmic::gmic(const char *const command_line, const char *const custom_commands, const bool default_commands,
           float *const p_progress, int *const p_cancel) {
  CImgList<float> images;
  gmic(command_line,images,custom_commands,default_commands,p_progress,p_cancel);
}

gmic::~gmic() {
  if (tmpstr) delete[] tmpstr;
}

// Get current scope as a string.
//-------------------------------
CImg<char> gmic::scope2string() const {
  CImgList<char> res(scope);
  cimglist_for(res,l) res[l].back() = '/';
  CImg<char>::vector(0).move_to(res);
  return res>'x';
}

// Return command line items from a string.
//-----------------------------------------
CImgList<char> gmic::command_line_to_CImgList(const char *const command_line) const {
  if (!command_line) return CImgList<char>();
  const char *ptrs0 = command_line;
  while (*ptrs0 && *ptrs0==' ') ++ptrs0;
  bool is_dquoted = false, is_escaped = false;
  CImg<char> item(ptrs0,std::strlen(ptrs0)+1);
  CImgList<char> items;
  char *const ptrd0 = item.data(), *ptrd = ptrd0;
  for (const char *ptrs = ptrs0; *ptrs || ptrs0!=ptrs; ++ptrs) {
    switch (*ptrs) {
    case '\\' : is_escaped = !is_escaped; *(ptrd++) = *ptrs; break;
    case '\"' : if (is_escaped) { *(ptrd-1) = *ptrs; is_escaped = false; } else is_dquoted = !is_dquoted; break;
    case '{' : case '}' : if (is_escaped) { *(ptrd-1) = *ptrs; is_escaped = false; } else *(ptrd++) = *ptrs-'{'+26; break;
    case ',' : if (is_escaped) { *(ptrd-1) = 29; is_escaped = false; } else *(ptrd++) = *ptrs; break;
    case '#' : if (is_escaped) { *(ptrd-1) = *ptrs; is_escaped = false; } else *(ptrd++) = *ptrs; break;
    case ' ' :
      if (is_escaped) { *(ptrd-1) = *ptrs; is_escaped = false; }
      else if (is_dquoted) *(ptrd++) = *ptrs;
      else {
        *ptrd = 0; CImg<char>(ptrd0,ptrd-ptrd0+1).move_to(items); ptrd = ptrd0;
        ++ptrs; while (*ptrs && *ptrs==' ') ++ptrs; ptrs0 = (ptrs--);
      }
      break;
    case 0 : *ptrd = 0; CImg<char>(ptrd0,ptrd-ptrd0+1).move_to(items); ptrs0 = (ptrs--); break;
    default : *(ptrd++) = *ptrs; is_escaped = false;
    }
  }
  cimglist_for(items,l) cimg::strescape(items[l].data());
  if (is_debug) {
    debug("Decompose command line into %u items : ",items.size());
    cimglist_for(items,l) {
      std::strcpy(tmpstr,items[l].data());
      std::fprintf(cimg::output(),"%s[%u]='%s'%s%s ",
                   cimg::t_green,l,tmpstr,l<(int)items.size()-1?"":",",cimg::t_normal);
    }
  }
  return items;
}

// Set default values of G'MIC parameters and add custom commands.
//----------------------------------------------------------------
gmic& gmic::assign(const char *const custom_commands, const bool default_commands,
                   float *const p_progress, int *const p_cancel) {
  command_names.assign();
  command_definitions.assign();
  scope.assign(CImg<char>(".",2));
  position = 0;
  verbosity_level = 0;
  is_released = true;
  is_debug = false;
  is_begin = true;
  is_end = false;
  is_quit = false;
  check_elif = false;
  patch_w = patch_h = patch_d = patch_c = 0;
  patch_borders = false;
  background3d[0] = 100;
  background3d[1] = 100;
  background3d[2] = 110;
  render3d = 4;
  renderd3d = -1;
  is_double3d = true;
  focale3d = 500;
  light3d_x = 0;
  light3d_y = 0;
  light3d_z = -5000;
  specular_light3d = 0.15f;
  specular_shine3d = 0.8f;
  is_fullpath = false;
  tmpstr = new char[16384];
  std::memset(tmpstr,0,16384);
  if (p_progress) progress = p_progress; else { _progress = -1; progress = &_progress; }
  if (p_cancel) cancel = p_cancel; else { _cancel = 0; cancel = &_cancel; }
  if (default_commands) add_commands(data_gmic_def);
  add_commands(custom_commands);
  return *this;
}

// General error procedure.
//--------------------------
const gmic& gmic::error(const char *const format, ...) const {
  char message[4096] = { 0 };
  va_list ap;
  va_start(ap,format);
  std::vsprintf(message + std::sprintf(message,"*** Error in %s *** ",scope2string().data()),format,ap);
  va_end(ap);
  if (verbosity_level>=0) {
    std::fprintf(cimg::output(),
                 "\n%s[gmic]%s %s%s%s\n"
                 "[gmic]%s Abort G'MIC instance.\n",
                 cimg::t_red,scope2string().data(),cimg::t_bold,message,cimg::t_normal,
                 scope2string().data());
    std::fflush(cimg::output());
  }
  throw gmic_exception(0,message);
  return *this;
}

// Warning procedure.
//-------------------
const gmic& gmic::warning(const char *format, ...) const {
  if (verbosity_level<0) return *this;
  char message[4096] = { 0 };
  va_list ap;
  va_start(ap,format);
  std::vsprintf(message + std::sprintf(message,"*** Warning in %s *** ",scope2string().data()),format,ap);
  va_end(ap);
  std::fprintf(cimg::output(),"\n[gmic]%s %s%s%s%s",
               scope2string().data(),cimg::t_bold,cimg::t_red,message,cimg::t_normal);
  std::fflush(cimg::output());
  return *this;
}

// Print debug message.
//---------------------
const gmic& gmic::debug(const char *format, ...) const {
  if (!is_debug) return *this;
  char message[4096] = { 0 };
  va_list ap;
  va_start(ap,format);
  std::vsprintf(message,format,ap);
  va_end(ap);
  for (char *s = message; *s; ++s) *s = *s==brace1?'{':*s==brace2?'}':*s;
  std::fprintf(cimg::output(),"\n%s<gmic>%s %s%s",
               cimg::t_green,scope2string().data(),message,cimg::t_normal);
  std::fflush(cimg::output());
  return *this;
}

// Print status message.
//----------------------
const gmic& gmic::print(const char *format, ...) const {
  if (verbosity_level<0) return *this;
  char message[4096] = { 0 };
  va_list ap;
  va_start(ap,format);
  std::vsprintf(message,format,ap);
  va_end(ap);
  std::fprintf(cimg::output(),"\n[gmic]%s %s",
               scope2string().data(),message);
  std::fflush(cimg::output());
  return *this;
}

// Add custom commands from a char* buffer.
//------------------------------------------
gmic& gmic::add_commands(const char *const data_commands) {
  if (!data_commands || !*data_commands) return *this;
  char mac[4096] = { 0 }, com[256*1024] = { 0 }, line[256*1024] = { 0 }, sep = 0;
  unsigned int pos = 0, comsiz0 = command_names.size();
  bool is_continued = false, nis_continued = false;
  for (const char *data = data_commands; *data; is_continued = nis_continued) {

    // Read line and remove comments.
    if (std::sscanf(data,"%262143[^\n]",line)>0) {
      data+=std::strlen(line); if (*data) ++data;
      if (*line=='#') *line = 0; else { char *const ptrc = std::strchr(line,'#'); if (ptrc && *(ptrc-1)==' ') *(ptrc-1) = 0; }
      cimg::strpare(line,' ',false,true);
      nis_continued = false;
      if (*line) {
        const int l = std::strlen(line)-1;
        for (const char *s = line + l; s>=line && *s=='\\'; --s) nis_continued = !nis_continued;
        if (nis_continued) line[l] = 0;
      }
      if (!*line) continue;
    } else { ++data; nis_continued = false; continue; }

    // Check for command definition.
    mac[0] = com[0] = 0;
    if (!is_continued && std::sscanf(line,"%4095[ a-zA-Z0-9_]%c%262143[^\n]",mac,&sep,com)>=2 &&
        sep==':' && std::sscanf(mac,"%4095s",line)==1) {
      CImg<char>(line,std::strlen(line)+1).move_to(command_names,pos);
      cimg::strpare(com,' ',false,true); CImg<char>(com,std::strlen(com)+1).move_to(command_definitions,pos++);
    } else { // Line continuation.
      if (comsiz0==command_names.size()) error("Command 'command' : Syntax error in expression '%s'.");
      if (!is_continued) command_definitions[pos-1].back() = ' '; else --(command_definitions[pos-1]._width);
      command_definitions[pos-1].append(CImg<char>(line,std::strlen(line)+1),'x');
    }
  }
  return *this;
}

// Add commands from a file.
//---------------------------
gmic& gmic::add_commands(std::FILE *const file) {
  if (!file) return *this;
  unsigned int siz = 0;
  std::fseek(file,0,SEEK_END);
  siz = (unsigned int)std::ftell(file);
  std::rewind(file);
  if (siz) {
    CImg<char> buffer(siz+1);
    if (std::fread(buffer.data(),sizeof(char),siz,file)) add_commands(buffer.data());
    buffer[siz] = 0;
  }
  return *this;
}

// Return selection from a selection string.
//------------------------------------------
CImg<unsigned int> gmic::selection2cimg(const char *const string, const unsigned int indice_max,
                                        const char *const command, const bool is_selection) const {
  if (!string || !*string) {
    if (indice_max) return CImg<unsigned int>::sequence(indice_max,0,indice_max-1);
    else return CImg<unsigned int>();
  }
  const char *const stype = is_selection?"selection":"subset";
  CImgList<unsigned int> selection;
  const char *it = string;
  for (bool stopflag = false; !stopflag; ) {
    char sep = 0, item0[4096] = { 0 }, item1[4096] = { 0 };
    float ind0 = 0, ind1 = 0, step = 1;
    if (std::sscanf(it,"%4095[^,]%c",item0,&end)!=2) stopflag = true;
    else it+=1 + std::strlen(item0);
    const int err = std::sscanf(item0,"%4095[^:]%c%f%c",item1,&sep,&step,&end);
    if (err!=1 && err!=3) error("Command '%s' : Syntax error in %s [%s].",command,stype,string);
    if (std::sscanf(item1,"%f%%-%f%c%c",&ind0,&ind1,&sep,&end)==3 && sep=='%') {
      ind0 = (float)cimg::round(ind0*(indice_max-1)/100,1);
      ind1 = (float)cimg::round(ind1*(indice_max-1)/100,1);
    } else if (std::sscanf(item1,"%f%%-%f%c",&ind0,&ind1,&end)==2)
      ind0 = (float)cimg::round(ind0*(indice_max-1)/100,1);
    else if (std::sscanf(item1,"%f-%f%c%c",&ind0,&ind1,&sep,&end)==3 && sep=='%')
      ind1 = (float)cimg::round(ind1*(indice_max-1)/100,1);
    else if (std::sscanf(item1,"%f-%f%c",&ind0,&ind1,&end)==2) { }
    else if (std::sscanf(item1,"%f%c%c",&ind0,&sep,&end)==2 && sep=='%')
      ind1 = (ind0 = (float)cimg::round(ind0*(indice_max-1)/100,1));
    else if (std::sscanf(item1,"%f%c",&ind0,&end)==1)
      ind1 = ind0;
    else error("Command '%s' : Syntax error in %s [%s].",command,stype,string);
    if (ind0<0) ind0+=indice_max;
    if (ind1<0) ind1+=indice_max;
    const int
      iind0 = (int)ind0,
      _ind1 = (int)ind1, iind1 = (int)(_ind1 - cimg::mod((float)_ind1,step));
    if (ind0>ind1) cimg::swap(ind0,ind1);

    if (!indice_max) error("Command '%s' : Invalid %s [%s] (no data available).",command,stype,string);
    if (step<=0) error("Command '%s' : Invalid %s [%s] (defines step '%g', should be >0).",command,stype,string,step);
    if (iind0<0 || iind0>=(int)indice_max)
      error("Command '%s' : Invalid %s [%s] (contains indice '%d', not in range -%u..%u).",
            command,stype,string,iind0,indice_max,indice_max-1);
    if (iind1<0 || iind1>=(int)indice_max)
      error("Command '%s' : Invalid %s [%s] (contains indice '%d', not in range -%u..%u).",
            command,stype,string,iind1,indice_max,indice_max-1);

    if (iind0==iind1) CImg<unsigned int>::vector((unsigned int)iind0).move_to(selection);
    else (CImg<unsigned int>::sequence((unsigned int)(1+(iind1-iind0)/step),
                                       (unsigned int)iind0,
                                       (unsigned int)iind1)<'y').move_to(selection);
  }
  selection = (selection>'y').sort()<'y';  // Sort indices in increasing order.
  cimglist_for(selection,l)
    if (l!=(int)selection.size()-1 && selection(l,0)==selection(l+1,0)) selection.remove(l--); // Remove indices that appear multiple time.
  if (is_debug) {
    debug(stype);
    std::fprintf(cimg::output(),"%s",cimg::t_green);
    (selection>'y').print("");
    std::fprintf(cimg::output(),"%s",cimg::t_normal);
  }
  return (selection>'y').sort();
}

// Return selection or filename strings from image selection.
//-----------------------------------------------------------
char *gmic::selection2string(const CImg<unsigned int>& selection, const CImgList<char>& filenames, const bool display_indices) const {
  static char res0[4096] = { 0 }, res1[4096] = { 0 };
  const unsigned int siz = selection.size();
  if (display_indices) {
    switch (siz) {
    case 0: std::sprintf(res0," []"); break;
    case 1: std::sprintf(res0," [%u]",selection[0]); break;
    case 2: std::sprintf(res0,"s [%u,%u]",selection[0],selection[1]); break;
    case 3: std::sprintf(res0,"s [%u,%u,%u]",selection[0],selection[1],selection[2]); break;
    case 4: std::sprintf(res0,"s [%u,%u,%u,%u]",selection[0],selection[1],selection[2],selection[3]); break;
    default: std::sprintf(res0,"s [%u,..,%u]",selection[0],selection[siz-1]);
    }
    return res0;
  }
  switch (siz) {
  case 0: std::sprintf(res1," "); break;
  case 1: std::sprintf(res1,"%s",filenames[selection[0]].data()); break;
  case 2: std::sprintf(res1,"%s, %s",filenames[selection[0]].data(),filenames[selection[1]].data()); break;
  case 3: std::sprintf(res1,"%s, %s, %s",filenames[selection[0]].data(),filenames[selection[1]].data(),
                       filenames[selection[2]].data()); break;
  case 4: std::sprintf(res1,"%s, %s, %s, %s",filenames[selection[0]].data(),filenames[selection[1]].data(),
                       filenames[selection[2]].data(), filenames[selection[3]].data()); break;
  default: std::sprintf(res1,"%s, .., %s",filenames[selection[0]].data(),filenames[selection[siz-1]].data());
  }
  return res1;
}
#endif // #if defined(gmic_float) || !defined(gmic_separate_compilation)

// General error procedure.
//-------------------------
template<typename T>
const gmic& gmic::error(const CImgList<T>& list, const char *const format, ...) const {
  char message[4096] = { 0 };
  va_list ap;
  va_start(ap,format);
  std::vsprintf(message + std::sprintf(message,"*** Error in %s *** ",scope2string().data()),format,ap);
  va_end(ap);
  if (verbosity_level>=0) {
    std::fprintf(cimg::output(),
                 "\n%s[gmic]-%u%s %s%s%s\n"
                 "[gmic]-%u%s Abort G'MIC instance.\n",
                 cimg::t_red,list.size(),scope2string().data(),cimg::t_bold,message,cimg::t_normal,
                 list.size(),scope2string().data());
    std::fflush(cimg::output());
  }
  throw gmic_exception(0,message);
  return *this;
}

// Bad argument procedure.
//------------------------
#define arg_error(command) _arg_error(images,command,argument_text)
template<typename T>
const gmic& gmic::_arg_error(const CImgList<T>& list, const char *const command, const char *const argument) const {
  char message[4096] = { 0 };
  std::sprintf(message,"*** Error in %s *** Command '%s' : Invalid argument '%s'.",scope2string().data(),command,argument);
  if (verbosity_level>=0) {
    std::fprintf(cimg::output(),
                 "\n%s[gmic]-%u%s %s%s%s\n"
                 "[gmic]-%u%s Abort G'MIC instance.\n",
                 cimg::t_red,list.size(),scope2string().data(),cimg::t_bold,message,cimg::t_normal,
                 list.size(),scope2string().data());
    std::fflush(cimg::output());
  }
  throw gmic_exception(command,message);
  return *this;
}

// Warning procedure.
//-------------------
template<typename T>
const gmic& gmic::warning(const CImgList<T>& list, const char *format, ...) const {
  if (verbosity_level<0) return *this;
  char message[4096] = { 0 };
  va_list ap;
  va_start(ap,format);
  std::vsprintf(message + std::sprintf(message,"*** Warning in %s *** ",scope2string().data()),format,ap);
  va_end(ap);
  std::fprintf(cimg::output(),"\n[gmic]-%u%s %s%s%s%s",
               list.size(),scope2string().data(),cimg::t_bold,cimg::t_red,message,cimg::t_normal);
  std::fflush(cimg::output());
  return *this;
}

// Print debug message.
//---------------------
template<typename T>
const gmic& gmic::debug(const CImgList<T>& list, const char *format, ...) const {
  if (!is_debug) return *this;
  char message[4096] = { 0 };
  va_list ap;
  va_start(ap,format);
  std::vsprintf(message,format,ap);
  va_end(ap);
  for (char *s = message; *s; ++s) *s = *s==brace1?'{':*s==brace2?'}':*s;
  std::fprintf(cimg::output(),"\n%s<gmic>-%u%s %s%s",
               cimg::t_green,list.size(),scope2string().data(),message,cimg::t_normal);
  std::fflush(cimg::output());
  return *this;
}

// Print status message.
//----------------------
template<typename T>
const gmic& gmic::print(const CImgList<T>& list, const char *format, ...) const {
  if (verbosity_level<0) return *this;
  char message[4096] = { 0 };
  va_list ap;
  va_start(ap,format);
  std::vsprintf(message,format,ap);
  va_end(ap);
  std::fprintf(cimg::output(),"\n[gmic]-%u%s %s",
               list.size(),scope2string().data(),message);
  std::fflush(cimg::output());
  return *this;
}

// Template constructors.
//-----------------------
template<typename T>
gmic::gmic(const int argc, const char *const *const argv, CImgList<T>& images,
           const char *custom_commands, const bool default_commands,
           float *const p_progress, int *const p_cancel) {
  assign(custom_commands,default_commands,p_progress,p_cancel);
  CImgList<char> command_line;
  for (int pos = 1; pos<argc; ++pos) {
    CImg<char>(argv[pos],std::strlen(argv[pos])+1).move_to(command_line);
    if (pos<argc-1) command_line.back().back()=' ';
  }
  is_debug = (argc>=2 && !std::strcmp(argv[1],"-debug"));
  const CImg<char> command(command_line>'x');
  debug("Initial command line : %s",command.data());
  command_line = command_line_to_CImgList(command);
  is_released = is_debug = false;
  unsigned int position = 0;
  CImgList<char> filenames;
  CImgList<unsigned int> dowhiles, repeatdones, locals;
  parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,true);
}

template<typename T>
gmic::gmic(const char *const command_line, CImgList<T>& images,
           const char *custom_commands, const bool default_commands,
           float *const p_progress, int *const p_cancel) {
  assign(custom_commands,default_commands,p_progress,p_cancel);
  is_released = true;
  const CImgList<char> items = command_line_to_CImgList(command_line);
  unsigned int position = 0;
  CImgList<char> filenames;
  CImgList<unsigned int> dowhiles, repeatdones, locals;
  parse(items,position,images,filenames,dowhiles,repeatdones,locals,true);
}

// Display selected images.
//-------------------------
template<typename T>
bool gmic::display_images(const CImgList<T>& images, const CImgList<char>& filenames, const CImg<unsigned int>& selection,
                          const bool verbose) const {
  if (!images || !filenames || !selection) { if (verbose) print(images,"Display image []."); return false; }
#if cimg_display==0
  if (verbose) print(images,"Display image%s (skipped, no display available).",gmic_selection);
  return true;
#endif
  CImgList<unsigned int> inds = selection.get_unroll('x')<'x';
  CImgList<T> visu;
  int max_height = 0;
  cimglist_for(inds,l) {
    const CImg<T>& img = images[inds(l,0)];
    if (img.height()>max_height && !img.is_CImg3d()) max_height = img.height();
  }
  cimglist_for(inds,l) {
    const unsigned int ind = inds(l,0);
    const CImg<T> &img = images[ind];
    if (img) {
      if (!max_height || img.height()<max_height) visu.insert(img,~0U,true);
      else img.get_lines(0,max_height-1).move_to(visu);
    } else if (verbose) { warning(images,"Command 'display' : Image [%d] is empty.",ind); inds.remove(l--); }
  }
  const CImg<unsigned int> nselection = inds>'x';
  const char *const fnames = selection2string(nselection,filenames,false);
  if (verbose) print(images,"Display image%s = '%s'.\n\n",gmic_selection,fnames);
  if (visu.size()) {
    if (visu.size()!=1) visu.display(fnames,verbosity_level>=0,'x','p');
    else {
      const CImg<T> &img = visu[0];
      std::sprintf(tmpstr,"%s (%dx%dx%dx%d)",fnames,img.width(),img.height(),img.depth(),img.spectrum());
      img.display(tmpstr,verbosity_level>=0);
    }
  }
  return true;
}

// Display plots of selected images.
//----------------------------------
template<typename T>
bool gmic::display_plots(const CImgList<T>& images, const CImgList<char>& filenames, const CImg<unsigned int>& selection,
                         const unsigned int plot_type, const unsigned int vertex_type,
                         const double xmin, const double xmax,
                         const double ymin, const double ymax,
                         const bool verbose) const {
  if (!images || !filenames || !selection) { print(images,"Plot image []."); return false; }
#if cimg_display==0
  print(images,"Plot image%s (skipped, no display available).",gmic_selection);
  return true;
#endif
  CImgDisplay disp(cimg_fitscreen(640,480,1),0,0);
  cimg_forY(selection,l) {
    const unsigned int ind = selection[l];
    const CImg<T>& img = images[ind];
    if (img) {
      print(images,"Plot image%s = '%s'.\n",gmic_selection,selection2string(selection,filenames,false));
      if (verbosity_level>=0) { std::fputc('\n',cimg::output()); img.print(filenames[ind].data()); }
      std::sprintf(tmpstr,"%s (%dx%dx%dx%d)",
                   filenames[ind].data(),img.width(),img.height(),img.depth(),img.spectrum());
      img.display_graph(disp.set_title("%s",tmpstr),plot_type,vertex_type,0,xmin,xmax,0,ymin,ymax);
    } else if (verbose) warning(images,"Command 'plot' : Image [%d] is empty.",ind);
  }
  return true;
}

// Display selected 3D objects.
//-----------------------------
template<typename T>
bool gmic::display_objects3d(const CImgList<T>& images, const CImgList<char>& filenames, const CImg<unsigned int>& selection,
                             const bool verbose) const {
  if (!selection) { print(images,"Display 3D object []."); return false; }
#if cimg_display==0
  print(images,"Display 3D object%s (skipped, no display available).",gmic_selection);
  return true;
#endif
  CImg<unsigned char> background;
  bool exist3d = false;
  CImgDisplay disp;
  cimg_forY(selection,l) {
    const unsigned int ind = selection[l];
    const CImg<T> &img = images[ind];
    if (!img.is_CImg3d()) {
      if (verbose) warning(images,"Command 'display3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
    } else {
      exist3d = true;
      if (!background || !disp) {
        background.assign(cimg_fitscreen(640,480,1),1,3);
        cimg_forC(background,k) background.get_shared_channel(k).fill(background3d[k]);
        disp.assign(background);
      }
      CImgList<unsigned int> primitives;
      CImgList<unsigned char> colors;
      CImg<float> opacities, vertices(img);
      vertices.CImg3dtoobject3d(primitives,colors,opacities);
      print(images,"Display 3D object [%u] = '%s' (%d vertices, %u primitives).",
            ind,filenames[ind].data(),vertices.width(),primitives.size());
      disp.set_title("%s (%d vertices, %u primitives)",
                     filenames[ind].data(),vertices.width(),primitives.size());
      background.display_object3d(disp,vertices,primitives,colors,opacities,
                                  true,render3d,renderd3d,is_double3d,focale3d,
                                  light3d_x,light3d_y,light3d_z,specular_light3d,specular_shine3d);
      if (disp.is_closed()) break;
    }
  }
  return exist3d;
}

// Substitute '@' and '{}' expressions.
//--------------------------------------
template<typename T>
bool gmic::substitute_item(const char *const source, char *const destination, const CImgList<T>& images,
                           const CImgList<unsigned int>& repeatdones) const {
  if (!source || !destination) return false;
  bool substitution_done = false;
  CImgList<char> items;
  for (const char *nsource = source; *nsource; )
    if (*nsource!='@' && *nsource!=brace1) { // If not starting with '@' and '{'
      const char *nsource0 = nsource; unsigned int l = 0;
      for (l = 0; *nsource && *nsource!='@' && *nsource!=brace1; ++l) ++nsource;
      CImg<char>(nsource0,l).move_to(items);

    } else { // '@' or '{}' expression found.
      char argument[4096] = { 0 }, sep = 0;
      int ind = 0, larg = 0;
      bool no_braces = true;

      // Isolate arguments between '{}'.
      if (*nsource==brace1) {
        const char *const ptr_beg = nsource + 1, *ptr_end = ptr_beg; unsigned int p = 0;
        for (p = 1; p>0 && *ptr_end; ++ptr_end) { if (*ptr_end==brace1) ++p; if (*ptr_end==brace2) --p; }
        if (p) { CImg<char>(nsource++,1).move_to(items); continue; }
        larg = ptr_end - ptr_beg - 1;
        if (larg>0) { char s[4096] = { 0 }; std::memcpy(s,ptr_beg,larg); substitute_item(s,argument,images,repeatdones); }
        nsource+=larg+2;
        const CImg<T> empty, &img = images.size()?images.back():empty;
        try {
          std::sprintf(tmpstr,"%g",img.eval(argument));
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        }
        catch (CImgException&) {
          for (const char *s = argument; *s; ++s) {
            std::sprintf(tmpstr,"%u",(unsigned int)*s);
            CImg<char>(tmpstr,std::strlen(tmpstr)+1).move_to(items).back().back()=',';
          }
          --(items.back()._width);
        }
        substitution_done = true;
        continue;

        // Isolate arguments between '@{}'.
      } else if (nsource[1]==brace1) {
        const char *const ptr_beg = nsource + 2, *ptr_end = ptr_beg; unsigned int p = 0;
        for (p = 1; p>0 && *ptr_end; ++ptr_end) { if (*ptr_end==brace1) ++p; if (*ptr_end==brace2) --p; }
        if (p) { CImg<char>(nsource++,1).move_to(items); continue; }
        larg = ptr_end - ptr_beg - 1;
        if (larg>0) { char s[4096] = { 0 }; std::memcpy(s,ptr_beg,larg); substitute_item(s,argument,images,repeatdones); }
        no_braces = false;
      }

      if (nsource[1]=='#') {
        // Substitute '@#' -> number of images in the list.
        nsource+=2;
        std::sprintf(tmpstr,"%u",images.size());
        CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        substitution_done = true;

        // Substitute '@!' -> visibility state of the first instant display window.
      } else if (nsource[1]=='!') {
        nsource+=2;
#if cimg_display==0
        std::sprintf(tmpstr,"0");
#else
        std::sprintf(tmpstr,"%d",instant_window[0]?(instant_window[0].is_closed()?0:1):0);
#endif
        CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        substitution_done = true;

        // Substitute '@{!}', @{!1}, '@{!,subset}' and '@{!1,subset}' -> features of an instant display window.
      } else if (argument[0]=='!' && (argument[1]==0 ||
                                      (argument[1]>='0' && argument[1]<='9' && argument[2]==0) ||
                                      (argument[1]==',' && argument[2]) ||
                                      (argument[1]>='0' && argument[1]<='9' && argument[2]==',' && argument[3]))) {
        nsource+=3 + larg;
#if cimg_display==0
        std::sprintf(tmpstr,"0");
#else
        unsigned int wind = 0;
        const char *nargument = argument+1;
        if (*nargument>='0' && *nargument<='9') { wind = (unsigned int)(*nargument-'0'); ++nargument; }
        if (!*nargument) std::sprintf(tmpstr,"%d",instant_window[wind]?(instant_window[wind].is_closed()?0:1):0);
        else if (*nargument==',') switch(*(++nargument)) {
          case 'w' : std::sprintf(tmpstr,"%d",instant_window[wind].width()); break;
          case 'h' : std::sprintf(tmpstr,"%d",instant_window[wind].height()); break;
          case 'u' : std::sprintf(tmpstr,"%d",CImgDisplay::screen_width()); break;
          case 'v' : std::sprintf(tmpstr,"%d",CImgDisplay::screen_height()); break;
          case 'x' : std::sprintf(tmpstr,"%d",instant_window[wind].mouse_x()); break;
          case 'y' : std::sprintf(tmpstr,"%d",instant_window[wind].mouse_y()); break;
          case 'n' : std::sprintf(tmpstr,"%d",instant_window[wind].normalization()); break;
          case 'b' : std::sprintf(tmpstr,"%d",instant_window[wind].button()); break;
          case 'o' : std::sprintf(tmpstr,"%d",instant_window[wind].wheel()); break;
          case 'c' : std::sprintf(tmpstr,"%d",(int)instant_window[wind].is_closed()); break;
          case 'r' : std::sprintf(tmpstr,"%d",(int)instant_window[wind].is_resized()); break;
          case 'm' : std::sprintf(tmpstr,"%d",(int)instant_window[wind].is_moved()); break;
          default : std::sprintf(tmpstr,"%d",instant_window[wind].is_key(argument+2)); break;
          }
        else std::sprintf(tmpstr,"@{!%s}",argument);
#endif
        CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        substitution_done = true;

        // Substitute '@*' -> number of items in the global stack.
      } else if (nsource[1]=='*') {
        nsource+=2;
        std::sprintf(tmpstr,"%u",stack.size());
        CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        substitution_done = true;

        // Substitute '@{*}' and '@{*,subset}' -> content of the global stack.
      } else if (argument[0]=='*' && (argument[1]==0 || (argument[1]==',' && argument[2]))) {
        nsource+=3 + larg;
        const CImg<unsigned int> sub = selection2cimg(argument+2,stack.size(),"Item substitution",false);
        if (sub) {
          cimg_foroff(sub,i) items.insert(stack[sub[i]]).back().back() = ',';
          --(items.back()._width);
        }
        substitution_done = true;

        // Substitute '@>' and '@<' -> current number nested loops.
      } else if (nsource[1]=='>' || nsource[1]=='<') {
        nsource+=2;
        std::sprintf(tmpstr,"%u",repeatdones.size());
        CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        substitution_done = true;

        // Substitute '@{>}' and '@{>,subset}' -> forward values of loop indices.
      } else if (argument[0]=='>' && (argument[1]==0 || (argument[1]==',' && argument[2]))) {
        nsource+=3 + larg;
        const CImg<unsigned int> sub = selection2cimg(argument+2,repeatdones.size(),"Item substitution",false);
        if (sub) {
          cimg_foroff(sub,i) {
            std::sprintf(tmpstr,"%u",repeatdones(sub[i],2));
            CImg<char>(tmpstr,std::strlen(tmpstr)+1).move_to(items).back().back()=',';
          }
          --(items.back()._width);
        }
        substitution_done = true;

        // Substitute '@{<}' and '@{<,subset}' -> backward values of loop indices.
      } else if (argument[0]=='<' && (argument[1]==0 || (argument[1]==',' && argument[2]))) {
        nsource+=3 + larg;
        const CImg<unsigned int> sub = selection2cimg(argument+2,repeatdones.size(),"Item substitution",false);
        if (sub) {
          cimg_foroff(sub,i) {
            std::sprintf(tmpstr,"%u",repeatdones(sub[i],1)-1);
            CImg<char>(tmpstr,std::strlen(tmpstr)+1).move_to(items).back().back()=',';
          }
          --(items.back()._width);
        }
        substitution_done = true;

        // Substitute '@ind', '@{ind}' and '@{ind,argument}' -> image values or feature.
      } else if (std::sscanf(nsource+1,"%d",&ind)==1 ||
                 std::sscanf(argument,"%d%c",&ind,&(end=0))==1 ||
                 std::sscanf(argument,"%d,%c",&ind,&sep)==2) {
        const unsigned int lind = std::sprintf(tmpstr,"%d",ind);
        nsource+=no_braces?1 + lind:3 + larg;
        int nind = ind;
        if (nind<0) nind+=images.size();
        if (nind<0 || nind>=(int)images.size()) {
          if (images.size()) error(images,"Item substitution : Invalid indice '%d' in expression '%s' (not in range -%u..%u).",
                                   ind,!*argument?"@ind":end?"@{ind,...}":"@{ind}",images.size(),images.size()-1);
          else error(images,"Item substitution : Invalid indice '%d' in expression '@ind' (no data available).",ind);
        }
        const CImg<T>& img = images[nind];

        float x = 0, y = 0, z = 0, v = 0; char sepp = 0, sepx = 0, sepy = 0, sepz = 0, sepv = 0;
        char argx[256] = { 0 }, argy[256] = { 0 }, argz[256] = { 0 }, argv[256] = { 0 }; int bcond = 0;
        const char *subset = sep?argument+lind+1:&sep;
        const unsigned int l = std::strlen(subset);
        if (*subset=='w' && l==1) {  // Substitute by image width.
          std::sprintf(tmpstr,"%d",img.width());
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='h' && l==1) { // Substitute by image height.
          std::sprintf(tmpstr,"%d",img.height());
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='d' && l==1) { // Substitute by image depth.
          std::sprintf(tmpstr,"%d",img.depth());
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='s' && l==1) { // Substitute by image spectrum.
          std::sprintf(tmpstr,"%d",img.spectrum());
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='#' && l==1) { // Substitute by number of values.
          std::sprintf(tmpstr,"%u",img.size());
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='+' && l==1) { // Substitute by sum of values.
          double res = img?(double)img.front():0; for (const T *ptrs = img.data()+1, *ptre = img.end(); ptrs<ptre; res+=(double)*ptrs++) {}
          std::sprintf(tmpstr,"%g",res);
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='-' && l==1) { // Substitute by difference of values.
          double res = img?(double)img.front():0; for (const T *ptrs = img.data()+1, *ptre = img.end(); ptrs<ptre; res-=(double)*ptrs++) {}
          std::sprintf(tmpstr,"%g",res);
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='*' && l==1) { // Substitute by product of values.
          double res = img?(double)img.front():0; for (const T *ptrs = img.data()+1, *ptre = img.end(); ptrs<ptre; res*=(double)*ptrs++) {}
          std::sprintf(tmpstr,"%g",res);
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='/' && l==1) { // Substitute by division of values.
          double res = img?(double)img.front():0; for (const T *ptrs = img.data()+1, *ptre = img.end(); ptrs<ptre; res/=(double)*ptrs++) {}
          std::sprintf(tmpstr,"%g",res);
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='m' && l==1) { // Substitute by minimum value.
          std::sprintf(tmpstr,"%g",(double)img.min());
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='M' && l==1) { // Substitute by maximum value.
          std::sprintf(tmpstr,"%g",(double)img.max());
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='a' && l==1) { // Substitute by image average.
          std::sprintf(tmpstr,"%g",img.mean());
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='v' && l==1) { // Substitute by image variance.
          std::sprintf(tmpstr,"%g",img.variance());
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='t' && l==1) { // Substitute by text string from image values.
          const T *ptrs = img.data(); char *ptrd = tmpstr;
          cimg_foroff(img,l) *(ptrd++) = (char)*(ptrs++);
          *ptrd = 0;
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='c' && l==1) { // Substitute by coordinates of minimum value.
          const CImg<unsigned int> st = img.get_stats();
          std::sprintf(tmpstr,"%u,%u,%u,%u",st[4],st[5],st[6],st[7]);
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if (*subset=='C' && l==1) { // Substitute by coordinates of maximum value.
          const CImg<unsigned int> st = img.get_stats();
          std::sprintf(tmpstr,"%u,%u,%u,%u",st[8],st[9],st[10],st[11]);
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else if ((std::sscanf(subset,"(%255[0-9.eE%+-]%c%c",argx,&sepp,&end)==2 || // Substitute by value at specified coordinates.
                    std::sscanf(subset,"(%255[0-9.eE%+-],%255[0-9.eE%+-]%c%c",argx,argy,&sepp,&end)==3 ||
                    std::sscanf(subset,"(%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c%c",argx,argy,argz,&sepp,&end)==4 ||
                    std::sscanf(subset,"(%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c%c",argx,argy,argz,argv,&sepp,&end)==5 ||
                    std::sscanf(subset,"(%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%d%c%c",argx,argy,argz,argv,&bcond,&sepp,&end)==6) &&
                   sepp==')' &&
                   (!*argx || std::sscanf(argx,"%f%c",&x,&end)==1 || (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
                   (!*argy || std::sscanf(argy,"%f%c",&y,&end)==1 || (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
                   (!*argz || std::sscanf(argz,"%f%c",&z,&end)==1 || (std::sscanf(argz,"%f%c%c",&z,&sepz,&end)==2 && sepz=='%')) &&
                   (!*argv || std::sscanf(argv,"%f%c",&v,&end)==1 || (std::sscanf(argv,"%f%c%c",&v,&sepv,&end)==2 && sepv=='%'))) {
          const int
            nx = (int)cimg::round(sepx=='%'?x*(img.width()-1)/100:x,1),
            ny = (int)cimg::round(sepy=='%'?y*(img.height()-1)/100:y,1),
            nz = (int)cimg::round(sepz=='%'?z*(img.depth()-1)/100:z,1),
            nv = (int)cimg::round(sepv=='%'?v*(img.spectrum()-1)/100:v,1);
          std::sprintf(tmpstr,"%g",bcond?(double)img.atXYZC(nx,ny,nz,nv):(double)img.atXYZC(nx,ny,nz,nv,0));
          CImg<char>(tmpstr,std::strlen(tmpstr)).move_to(items);
        } else { // Substitute by value subset (default).
          CImg<T> values;
          if (!*subset) values = img.get_shared();
          else {
            const CImg<unsigned int> inds = selection2cimg(subset,img.size(),"Item substitution",false);
            values.assign(inds.size());
            cimg_foroff(inds,p) values[p] = img[inds(p)];
          }
          CImg<char> s_values = values.value_string();
          --(s_values._width); s_values.move_to(items);
        }
        substitution_done = true;

        // Substitute any other '@' expression.
      } else CImg<char>(nsource++,1).move_to(items);
    }
  CImg<char>::vector(0).move_to(items);
  const CImg<char> _items = items>'x';
  if (_items.size()>4095) {
    std::strcpy(tmpstr,source);
    for (char *s = tmpstr; *s; ++s) *s = *s==brace1?'{':*s==brace2?'}':*s;
    error(images,"Item substitution : Buffer overflow when substituting '%s'.",tmpstr);
  }
  char *d = destination;
  for (const char *s = _items.data(); *s; ++s) *(d++) = *s==brace1?'{':*s==brace2?'}':*s; *d = 0;
  if (substitution_done) debug(images,"Substitute item '%s' as '%s'.",source,destination);
  return substitution_done;
}

// Main parsing procedure.
//------------------------
template<typename T>
gmic& gmic::parse(const CImgList<char>& command_line, unsigned int& position, CImgList<T> &images, CImgList<char> &filenames,
                  CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones, CImgList<unsigned int>& locals,
                  const bool initial_call) {

  try {
    const int no_ind = (int)(~0U>>1);
    cimg::exception_mode() = 0;
    if (images.size()<filenames.size()) filenames.remove(images.size(),~0U);
    else if (images.size()>filenames.size()) filenames.insert(images.size() - filenames.size(),CImg<char>("(unnamed)",10));

    // Begin command line parsing.
    while (position<command_line.size() && !is_quit) {
      const char
        *const orig_item = command_line[position].data(),
        *const orig_argument = position+1<command_line.size()?command_line[position+1].data():"";

      // Get a constant reference to the last image.
      const CImg<T> _last_image, &last_image = images.size()?images.back():_last_image;

      // Check consistency between variables.
      if (filenames.size()!=images.size())
        error("Internal error : Images (%u) and filenames (%u) are in an inconsistent state.",filenames.size(),images.size());
      if (!scope.size())
        error("Internal error : Scope is empty.");
      if (scope.size()>64)
        error("Internal error : Scope overflow (possible infinite recursion).");
      if (stack.size()>65536)
        error("Internal error : Global stack overflow.");

      // Substitute '@' and '{}' expressions in 'item' and 'argument' if necessary.
      char _item[4096] = { 0 }, _argument[4096] = { 0 };
      bool sub_item = false, sub_argument = false;
      if (*orig_item=='-' || *orig_item=='@' || *orig_item==brace1 || *orig_item=='[' || *orig_item=='(') {
        if (std::strchr(orig_item,'@') || std::strchr(orig_item,brace1))
          sub_item = substitute_item(orig_item,_item,images,repeatdones);
        if (*orig_item=='-' &&
            (*orig_argument!='-' || orig_argument[1]=='.' || orig_argument[1]=='@' || orig_argument[1]==brace1 ||
             (orig_argument[1]>='0' && orig_argument[1]<='9')) &&
            (std::strchr(orig_argument,'@') || std::strchr(orig_argument,brace1)))
          sub_argument = substitute_item(orig_argument,_argument,images,repeatdones);
      }
      const char *item = sub_item?_item:orig_item, *argument = sub_argument?_argument:orig_argument;
      char argument_text[64] = { 0 };
      if (std::strlen(argument)>=64) {
        std::memcpy(argument_text,argument,60*sizeof(char));
        argument_text[60] = argument_text[61] = argument_text[62] = '.'; argument_text[63] = 0;
      } else std::strcpy(argument_text,argument);

      // Get current item/command from the command line.
      char command_name[4096] = { 0 }, command_restriction[4096] = { 0 };
      bool get_version = false, is_restriction = false;
      CImg<unsigned int> selection;
      if (item[0]=='-' && item[1] && item[1]!='.') {
        char sep0 = 0, sep1 = 0;
        if (item[1]=='-' && item[2] && item[2]!='[' && (item[2]!='3' || item[3]!='d')) { ++item; get_version = true; }
        const int err = std::sscanf(item,"%4095[^[]%c%4095[0-9.eE%,:+-]%c%c",command_name,&sep0,command_restriction,&sep1,&end);
        if (err==1) selection = CImg<unsigned int>::sequence(images.size(),0,images.size()-1);
        else if (err==2 && sep0=='[' && item[std::strlen(command_name)+1]==']') { selection.assign(); is_restriction = true; }
        else if (err==4 && sep1==']') {
          is_restriction = true;
          if (!std::strcmp("-push",command_name) || !std::strcmp("-p",command_name))
            selection = selection2cimg(command_restriction,stack.size()+1,command_name,true);
          else if (!std::strcmp("-pop",command_name) || !std::strcmp("-pp",command_name))
            selection = selection2cimg(command_restriction,stack.size(),command_name,true);
          else if (!std::strcmp("-wait",command_name))
            selection = selection2cimg(command_restriction,10,command_name,true);
          else if (!std::strcmp("-input",command_name) || !std::strcmp("-i",command_name))
            selection = selection2cimg(command_restriction,images.size()+1,command_name,true);
          else
            selection = selection2cimg(command_restriction,images.size(),command_name,true);
        } else { std::strcpy(command_name,item); *command_restriction = 0; }
      }
      ++position;

      // Check for verbosity commands.
      if (*item=='-') {
        if (!std::strcmp("-verbose+",item) || !std::strcmp("-v+",item)) ++verbosity_level;
        else if (!std::strcmp("-verbose-",item) || !std::strcmp("-v-",item)) --verbosity_level;
        else if (!std::strcmp("-verbose",item) || !std::strcmp("-v",item)) {
          int level = 0;
          if (std::sscanf(argument,"%d%c",&level,&end)==1) verbosity_level = level;
          else arg_error("verbose");
        }
      }

      if (is_begin) { print(images,"Start G'MIC instance."); is_begin = false; }
      debug(images,"Current item : '%s', Selection%s, Argument : '%s'.",item,gmic_selection,argument);

      // Begin command interpretation.
      if (*item=='-') {

        //----------------
        // Global options
        //----------------

        // Verbosity (actually, just continue to next command since verbosity has been already processed above).
        if (!std::strcmp("-verbose+",item) || !std::strcmp("-v+",item) ||
            !std::strcmp("-verbose-",item) || !std::strcmp("-v-",item)) continue;
        if (!std::strcmp("-verbose",item) || !std::strcmp("-v",item)) { ++position; continue; }

        // Load command file.
        if (!std::strcmp("-command",item) || !std::strcmp("-m",item)) {
          std::FILE *const file = std::fopen(argument,"r");
          const unsigned int siz = command_names.size();
          if (file) {
            const char *const basename = cimg::basename(argument);
            print(images,"Import G'MIC command(s), from file '%s'",is_fullpath?argument:basename);
            add_commands(file);
            cimg::fclose(file);
          } else {
            print(images,"Import G'MIC command(s), from string '%s'.",argument_text);
            add_commands(argument);
          }
          if (verbosity_level>=0) {
            const unsigned int nb_new = command_names.size()-siz;
            std::fprintf(cimg::output()," (%u command%s added).",nb_new,nb_new!=1?"s":"");
            std::fflush(cimg::output());
          }
          ++position; continue;
        }

        // Switch debug flag.
        if (!std::strcmp("-debug",item)) {
          is_debug = !is_debug;
          print(images,"%s debug mode.",is_debug?"Activate":"Deactivate");
          if (is_debug) verbosity_level = 100;
          continue;
        }

        // Switch fullpath mode.
        if (!std::strcmp("-fullpath",item)) {
          is_fullpath = !is_fullpath;
          print(images,"%s full-path mode.",is_fullpath?"Activate":"Deactivate");
          continue;
        }

        //----------------------
        // Arithmetic operators
        //----------------------

        // Common arithmetic operators.
        gmic_arithmetic_item("-add","-+",
                             operator+=,"Add '%g' to image%s",value,gmic_selection,T,
                             operator+=,"Add image [%d] to image%s",ind[0],gmic_selection,
                             "Add expression %s to image%s",argument_text,gmic_selection,
                             "Add image%s together");

        gmic_arithmetic_item("-sub","--",
                             operator-=,"Substract '%g' to image%s",value,gmic_selection,T,
                             operator-=,"Substract image [%d] to image%s",ind[0],gmic_selection,
                             "Substract expression %s to image%s",argument_text,gmic_selection,
                             "Substract image%s together");

        gmic_arithmetic_item("-mul","-*",
                             operator*=,"Multiply image%s by '%g'",gmic_selection,value,double,
                             mul,"Multiply image%s by image [%d]",gmic_selection,ind[0],
                             "Multiply image%s by expression %s",gmic_selection,argument_text,
                             "Multiply image%s together");

        gmic_arithmetic_item("-div","-/",
                             operator/=,"Divide image%s by '%g'",gmic_selection,value,double,
                             div,"Divide image%s by image [%d]",gmic_selection,ind[0],
                             "Divide image%s by expression %s",gmic_selection,argument_text,
                             "Divide image%s together");

        gmic_arithmetic_item("-pow","-^",
                             pow,"Compute image%s to the power of '%g'",gmic_selection,value,double,
                             pow,"Compute image%s to the power of image [%d]",gmic_selection,ind[0],
                             "Compute image%s to the power of expression %s",gmic_selection,argument_text,
                             "Compute the power of image%s together");

        gmic_arithmetic_item("-min","-min",
                             min,"Compute pointwise minimum between image%s and '%g'",gmic_selection,value,T,
                             min,"Compute pointwise minimum between image%s and image [%d]",gmic_selection,ind[0],
                             "Compute pointwise minimum between image%s and expression %s",gmic_selection,argument_text,
                             "Compute pointwise minimum between image%s together");

        gmic_arithmetic_item("-max","-max",
                             max,"Compute pointwise maximum between image%s and '%g'",gmic_selection,value,T,
                             max,"Compute pointwise maximum between image%s and image [%d]",gmic_selection,ind[0],
                             "Compute pointwise maximum between image%s and expression %s",gmic_selection,argument_text,
                             "Compute pointwise maximum between image%s together");

        gmic_arithmetic_item("-mod","-%",
                             operator%=,"Compute pointwise modulo between image%s and '%g'.",gmic_selection,value,T,
                             operator%=,"Compute pointwise modulo between image%s and image [%d]",gmic_selection,ind[0],
                             "Compute pointwise modulo between image%s and expression %s",gmic_selection,argument_text,
                             "Compute pointwise modulo between image%s together");

        gmic_arithmetic_item("-and","-and",
                             operator&=,"Compute bitwise AND between image%s and '%g'.",gmic_selection,value,T,
                             operator&=,"Compute bitwise AND between image%s and image [%d]",gmic_selection,ind[0],
                             "Compute bitwise AND between image%s and expression %s",gmic_selection,argument_text,
                             "Compute bitwise AND between image%s together");

        gmic_arithmetic_item("-or","-or",
                             operator|=,"Compute bitwise OR between image%s and '%g'.",gmic_selection,value,T,
                             operator|=,"Compute bitwise OR between image%s and image [%d]",gmic_selection,ind[0],
                             "Compute bitwise OR between image%s and expression %s",gmic_selection,argument_text,
                             "Compute bitwise OR between image%s together");

        gmic_arithmetic_item("-xor","-xor",
                             operator^=,"Compute bitwise XOR between image%s and '%g'.",gmic_selection,value,T,
                             operator^=,"Compute bitwise XOR between image%s and image [%d]",gmic_selection,ind[0],
                             "Compute bitwise XOR between image%s and expression %s",gmic_selection,argument_text,
                             "Compute bitwise XOR between image%s together");

        // Other arithmetic operators.
        gmic_simple_item("-cos",cos,"Compute cosine of image%s.");
        gmic_simple_item("-sin",sin,"Compute sine of image%s.");
        gmic_simple_item("-tan",tan,"Compute tangent of image%s.");
        gmic_simple_item("-acos",acos,"Compute arc-cosine of image%s.");
        gmic_simple_item("-asin",asin,"Compute arc-sine of image%s.");
        gmic_simple_item("-atan",atan,"Compute arc-tangent of image%s.");
        gmic_simple_item("-abs",abs,"Compute absolute value of image%s.");
        gmic_simple_item("-sqr",sqr,"Compute square function of image%s.");
        gmic_simple_item("-sqrt",sqrt,"Compute square root of image%s.");
        gmic_simple_item("-exp",exp,"Compute exponential of image%s.");
        gmic_simple_item("-log",log,"Compute logarithm of image%s.");
        gmic_simple_item("-log10",log10,"Compute logarithm-10 of image%s.");

        // Arc-tangent2.
        if (!std::strcmp("-atan2",command_name)) {
          char sep = 0; int ind = no_ind;
          if (std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') {
            gmic_check_indice(ind);
            print(images,"Compute oriented arc-tangent of image%s, using x-argument [%d].",gmic_selection,ind);
            const CImg<T> img0 = images[ind];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],atan2(img0));
          } else arg_error("atan2");
          ++position; continue;
        }

        //---------------------------------------
        // Pointwise operations on pixel values
        //---------------------------------------

        // Type cast.
        if (!std::strcmp("-type",item)) {
          typedef unsigned char uchar;
          typedef unsigned short ushort;
          typedef unsigned int uint;
#ifndef gmic_minimal
          gmic_cast(bool,"bool");
          gmic_cast(uchar,"unsigned char");
          gmic_cast(char,"char");
          gmic_cast(ushort,"unsigned short");
          gmic_cast(short,"short");
          gmic_cast(uint,"unsigned int");
          gmic_cast(int,"int");
          gmic_cast(double,"double");
#endif
          gmic_cast(float,"float");
          arg_error("type");
        }

        // Set pixel (scalar) value.
        if (!std::strcmp("-set",command_name) || !std::strcmp("-=",command_name)) {
          double value = 0; float x = 0, y = 0, z = 0, v = 0; char sepx = 0, sepy = 0, sepz = 0, sepv = 0;
          char argx[4096] = { 0 }, argy[4096] = { 0 }, argz[4096] = { 0 }, argv[4096] = { 0 };
          if ((std::sscanf(argument,"%lf%c",&value,&end)==1 ||
               std::sscanf(argument,"%lf,%4095[0-9.eE%+-]%c",&value,argx,&end)==2 ||
               std::sscanf(argument,"%lf,%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",&value,argx,argy,&end)==3 ||
               std::sscanf(argument,"%lf,%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",&value,argx,argy,argz,&end)==4 ||
               std::sscanf(argument,"%lf,%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",
                           &value,argx,argy,argz,argv,&end)==5) &&
              (!*argx || (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%') || std::sscanf(argx,"%f%c",&x,&end)==1) &&
              (!*argy || (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%') || std::sscanf(argy,"%f%c",&y,&end)==1) &&
              (!*argz || (std::sscanf(argz,"%f%c%c",&z,&sepz,&end)==2 && sepz=='%') || std::sscanf(argz,"%f%c",&z,&end)==1) &&
              (!*argv || (std::sscanf(argv,"%f%c%c",&v,&sepv,&end)==2 && sepv=='%') || std::sscanf(argv,"%f%c",&v,&end)==1)) {
            print(images,"Set scalar value %g at position (%g%s,%g%s,%g%s,%g%s) in image%s.",
                  value,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,sepz=='%'?"%":"",v,sepv=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width()-1)/100:x,1),
                ny = (int)cimg::round(sepy=='%'?y*(img.height()-1)/100:y,1),
                nz = (int)cimg::round(sepz=='%'?z*(img.depth()-1)/100:z,1),
                nv = (int)cimg::round(sepv=='%'?v*(img.spectrum()-1)/100:v,1);
              gmic_apply(images[selection[l]],gmic_set(value,nx,ny,nz,nv));
            }
          } else arg_error("set");
          ++position; continue;
        }

        // Invert endianness.
        gmic_simple_item("-endian",invert_endianness,"Invert data endianness of image%s.");

        // Fill.
        if (!std::strcmp("-fill",command_name) || !std::strcmp("-f",command_name)) {
          char sep = 0; double value = 0; int ind = no_ind;
          if (std::sscanf(argument,"%lf%c",&value,&end)==1) {
            print(images,"Fill image%s with value %g.",gmic_selection,value);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],fill((T)value));
          } else if (std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') {
            gmic_check_indice(ind);
            print(images,"Fill image%s with values from image [%d].",gmic_selection,ind);
            const CImg<T> values = images[ind];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],fill(values));
          } else {
            print(images,"Fill image%s with expression '%s'.",gmic_selection,argument_text);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],fill(argument,true));
          }
          ++position; continue;
        }

        // Threshold.
        if (!std::strcmp("-threshold",command_name) || !std::strcmp("-t",command_name)) {
          char sep = 0; unsigned int soft = 0; double value = 0;
          if (std::sscanf(argument,"%lf%c",&value,&end)==1 ||
              (std::sscanf(argument,"%lf%c%c",&value,&sep,&end)==2 && sep=='%') ||
              std::sscanf(argument,"%lf,%u%c",&value,&soft,&end)==2 ||
              std::sscanf(argument,"%lf%c,%u%c",&value,&sep,&soft,&end)==3) {
            print(images,"%s-threshold image%s with value %g%s.",soft?"Soft":"Hard",gmic_selection,value,sep=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              double vmin = 0, vmax = 0, nvalue = value;
              if (sep=='%') { vmin = img.min_max(vmax); nvalue = vmin + (vmax - vmin)*value/100; }
              gmic_apply(img,threshold((T)nvalue,soft?true:false));
            }
            ++position;
          } else {
#if cimg_display==0
            print(images,"Threshold image%s : run interactive mode (skipped, no display available).",gmic_selection);
#else
            print(images,"Threshold image%s : run interactive mode.",gmic_selection);
            CImgDisplay _disp, &disp = instant_window[0]?instant_window[0]:_disp;
            cimg_forY(selection,l) {
              CImg<T>
                &img = images[selection[l]],
                visu = img.depth()>1?img.get_projections2d(img.width()/2,img.height()/2,img.depth()/2).
                channels(0,cimg::min(3,img.spectrum())-1):img.get_channels(0,cimg::min(3,img.spectrum()-1));
              if (disp) disp.resize(cimg_fitscreen(visu.width(),visu.height(),1),false);
              else disp.assign(cimg_fitscreen(visu.width(),visu.height(),1),0,1);
              double
                vmin = 0, vmax = (double)img.max_min(vmin),
                distmax = std::sqrt(cimg::sqr(disp.width()-1.0) + cimg::sqr(disp.height()-1.0)),
                amount = 50;
              bool stopflag = false, obutt = false;
              int omx = -1, omy = -1;
              CImg<T> res;
              for (disp.show().flush(); !stopflag; ) {
                const unsigned int key = disp.key();
                if (!res) {
                  std::sprintf(tmpstr,"%s : threshold %.3g%%",filenames[selection[l]].data(),amount);
                  disp.display(res=visu.get_threshold((T)(vmin + amount*(vmax-vmin)/100))).
                    set_title("%s",tmpstr).wait();
                }
                const int mx = disp.mouse_x(), my = disp.mouse_y();
                if (disp.button() && mx>=0 && my>=0) {
                  if (omx==mx && omy==my && !obutt) break;
                  omx = mx; omy = my; obutt = true;
                  const double dist = std::sqrt((double)cimg::sqr(mx) + cimg::sqr(my));
                  amount = dist*100/distmax;
                  res.assign();
                } else if (!disp.button()) obutt = false;
                if (disp.is_closed() || (key && key!=cimg::keyCTRLLEFT)) stopflag = true;
                if (key==cimg::keyD && disp.is_keyCTRLLEFT() &&
                    (disp.resize(cimg_fitscreen(3*disp.width()/2,3*disp.height()/2,1),stopflag=false).set_key())==0)
                  disp._is_resized = true;
                if (key==cimg::keyC && disp.is_keyCTRLLEFT() &&
                    (disp.resize(cimg_fitscreen(2*disp.width()/3,2*disp.height()/3,1),stopflag=false).set_key())==0)
                  disp._is_resized = true;
                if (disp.is_resized()) {
                  disp.resize(false).display(res);
                  distmax = std::sqrt(cimg::sqr(disp.width()-1.0) + cimg::sqr(disp.height()-1.0));
                }
              }
              gmic_apply(img,threshold((T)(vmin + amount*(vmax-vmin)/100)));
            }
#endif
          }
          continue;
        }

        // Cut.
        if (!std::strcmp("-cut",command_name) || !std::strcmp("-c",command_name)) {
          char sep0 = 0, sep1 = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          double value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
              ((std::sscanf(arg0,"[%d%c%c",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               std::sscanf(arg0,"%lf%c",&value0,&end)==1) &&
              ((std::sscanf(arg1,"[%d%c%c",&ind1,&sep1,&end)==2 && sep1==']') ||
               (std::sscanf(arg1,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               std::sscanf(arg1,"%lf%c",&value1,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].min(); sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1); value1 = images[ind1].max(); sep1 = 0; }
            print(images,"Cut image%s in value range [%g%s,%g%s].",gmic_selection,value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              double vmin = 0, vmax = 0, nvalue0 = value0, nvalue1 = value1;
              if (sep0=='%') { vmin = img.min_max(vmax); nvalue0 = vmin + (vmax - vmin)*value0/100; }
              if (sep1=='%') { vmin = img.min_max(vmax); nvalue1 = vmin + (vmax - vmin)*value1/100; }
              gmic_apply(img,cut((T)nvalue0,(T)nvalue1));
            }
            ++position;
          } else if (std::sscanf(argument,"[%d%c%c",&(ind0=no_ind),&sep0,&end)==2) {
            if (ind0!=no_ind) gmic_check_indice(ind0);
            value0 = images[ind0].min_max(value1);
            print(images,"Cut image%s in value range [%g,%g].",gmic_selection,value0,value1);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],cut((T)value0,(T)value1));
            ++position;
          } else {
#if cimg_display==0
            print(images,"Cut image%s : run interactive mode (skipped, no display available).",gmic_selection);
#else
            print(images,"Cut image%s : run interactive mode.",gmic_selection);
            CImgDisplay _disp, &disp = instant_window[0]?instant_window[0]:_disp;
            cimg_forY(selection,l) {
              CImg<T>
                &img = images[selection[l]],
                visu = img.depth()>1?img.get_projections2d(img.width()/2,img.height()/2,img.depth()/2).
                channels(0,cimg::min(3,img.spectrum())-1):img.get_channels(0,cimg::min(3,img.spectrum()-1));
              if (disp) disp.resize(cimg_fitscreen(visu.width(),visu.height(),1),false);
              else disp.assign(cimg_fitscreen(visu.width(),visu.height(),1),0,1);
              double vmin = 0, vmax = (double)img.max_min(vmin), amount0 = 0, amount1 = 100;
              bool stopflag = false, obutt = false;
              int omx = -1, omy = -1;
              CImg<T> res;
              for (disp.show().flush(); !stopflag; ) {
                const unsigned int key = disp.key();
                if (!res) {
                  std::sprintf(tmpstr,"%s : cut [%.3g%%,%.3g%%]",
                               filenames[selection[l]].data(),amount0,amount1);
                  disp.display(res = visu.get_cut((T)(vmin + amount0*(vmax-vmin)/100),
                                                  (T)(vmin + amount1*(vmax-vmin)/100))).
                    set_title("%s",tmpstr).wait();
                }
                const int mx = disp.mouse_x(), my = disp.mouse_y();
                if (disp.button() && mx>=0 && my>=0) {
                  if (omx==mx && omy==my && !obutt) break;
                  omx = mx; omy = my; obutt = true;
                  amount0 = mx*100/disp.width(); amount1 = my*100/disp.height();
                  res.assign();
                } else if (!disp.button()) obutt = false;
                if (disp.is_closed() || (key && key!=cimg::keyCTRLLEFT)) stopflag = true;
                if (key==cimg::keyD && disp.is_keyCTRLLEFT() &&
                    (disp.resize(cimg_fitscreen(3*disp.width()/2,3*disp.height()/2,1),stopflag=false).set_key())==0)
                  disp._is_resized = true;
                if (key==cimg::keyC && disp.is_keyCTRLLEFT() &&
                    (disp.resize(cimg_fitscreen(2*disp.width()/3,2*disp.height()/3,1),stopflag=false).set_key())==0)
                  disp._is_resized = true;
                if (disp.is_resized()) disp.resize(false).display(res);
              }
              gmic_apply(img,cut((T)(vmin + amount0*(vmax-vmin)/100),(T)(vmin + amount1*(vmax-vmin)/100)));
            }
#endif
          }
          continue;
        }

        // Normalize.
        if (!std::strcmp("-normalize",command_name) || !std::strcmp("-n",command_name)) {
          char sep0 = 0, sep1 = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          double value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
              ((std::sscanf(arg0,"[%d%c%c",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               std::sscanf(arg0,"%lf%c",&value0,&end)==1) &&
              ((std::sscanf(arg1,"[%d%c%c",&ind1,&sep1,&end)==2 && sep1==']') ||
               (std::sscanf(arg1,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               std::sscanf(arg1,"%lf%c",&value1,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].min(); sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1); value1 = images[ind1].max(); sep1 = 0; }
            print(images,"Normalize image%s in value range [%g%s,%g%s].",gmic_selection,value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              double vmin = 0, vmax = 0, nvalue0 = value0, nvalue1 = value1;
              if (sep0=='%') { vmin = img.min_max(vmax); nvalue0 = vmin + (vmax - vmin)*value0/100; }
              if (sep1=='%') { vmin = img.min_max(vmax); nvalue1 = vmin + (vmax - vmin)*value1/100; }
              gmic_apply(img,normalize((T)nvalue0,(T)nvalue1));
            }
          } else if (std::sscanf(argument,"[%d%c%c",&(ind0=no_ind),&sep0,&end)==2) {
            if (ind0!=no_ind) gmic_check_indice(ind0);
            value0 = images[ind0].min_max(value1);
            print(images,"Normalize image%s in value range [%g,%g].",gmic_selection,value0,value1);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],normalize((T)value0,(T)value1));
          } else arg_error("normalize");
          ++position; continue;
        }

        // Round.
        if (!std::strcmp("-round",command_name)) {
          double rounding_value = 0; int rounding_type = 0;
          if ((std::sscanf(argument,"%lf%c",&rounding_value,&end)==1 ||
               std::sscanf(argument,"%lf,%d%c",&rounding_value,&rounding_type,&end)==2) &&
              rounding_value>=0 && rounding_type>=-1 && rounding_type<=1) {
            print(images,"Round image%s with value %g and %s rounding.",
                  gmic_selection,rounding_value,rounding_type<0?"backward":rounding_type>0?"forward":"nearest");
            cimg_forY(selection,l) gmic_apply(images[selection[l]],round((float)rounding_value,rounding_type));
          } else arg_error("round");
          ++position; continue;
        }

        // Equalize.
        if (!std::strcmp("-equalize",command_name)) {
          int nb_levels = 256; double vmin = 0, vmax = 0; char sep = 0, sepm = 0, sepM = 0;
          if ((std::sscanf(argument,"%d%c",&nb_levels,&end)==1 ||
               (std::sscanf(argument,"%d%c%c",&nb_levels,&sep,&end)==2 && sep=='%') ||
               std::sscanf(argument,"%d,%lf,%lf%c",&nb_levels,&vmin,&vmax,&end)==3 ||
               (std::sscanf(argument,"%d%c,%lf,%lf%c",&nb_levels,&sep,&vmin,&vmax,&end)==4 && sep=='%') ||
               (std::sscanf(argument,"%d,%lf%c,%lf%c",&nb_levels,&vmin,&sepm,&vmax,&end)==4 && sepm=='%') ||
               (std::sscanf(argument,"%d%c,%lf%c,%lf%c",&nb_levels,&sep,&vmin,&sepm,&vmax,&end)==5 && sep=='%' && sepm=='%') ||
               (std::sscanf(argument,"%d,%lf,%lf%c%c",&nb_levels,&vmin,&vmax,&sepM,&end)==4 && sepm=='%') ||
               (std::sscanf(argument,"%d%c,%lf,%lf%c%c",&nb_levels,&sep,&vmin,&vmax,&sepM,&end)==5 && sep=='%' && sepm=='%') ||
               (std::sscanf(argument,"%d,%lf%c,%lf%c%c",&nb_levels,&vmin,&sepm,&vmax,&sepM,&end)==5 && sepm=='%' && sepM=='%') ||
               (std::sscanf(argument,"%d%c,%lf%c,%lf%c%c",&nb_levels,&sep,&vmin,&sepm,&vmax,&sepM,&end)==6 && sep=='%' && sepm=='%' &&
                sepM=='%')) &&
              nb_levels>0) {
            if (vmin==vmax && vmin==0) { vmax = 100; sepM = '%'; }
            print(images,"Equalize histogram of image%s, using %d%s levels in range [%g%s,%g%s].",
                  gmic_selection,nb_levels,sep=='%'?"%":"",vmin,sepm=='%'?"%":"",vmax,sepM=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              double m = vmin, M = vmax;
              int nnb_levels = nb_levels;
              if (sepm=='%') m*=img.min()/100;
              if (sepM=='%') M*=img.max()/100;
              if (sep=='%') nnb_levels = (int)cimg::round(nb_levels*(1+M-m)/100,1);
              gmic_apply(images[selection[l]],equalize(nnb_levels,(T)m,(T)M));
            }
          } else arg_error("equalize");
          ++position; continue;
        }

        // Quantize.
        if (!std::strcmp("-quantize",command_name)) {
          int nb_levels = 0, keep_range = 1;
          if ((std::sscanf(argument,"%d%c",&nb_levels,&end)==1 ||
               std::sscanf(argument,"%d,%d%c",&nb_levels,&keep_range,&end)==2) &&
              nb_levels>0) {
            print(images,"Quantize image%s in %d levels.",gmic_selection,nb_levels);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],quantize(nb_levels,(bool)keep_range));
          } else arg_error("quantize");
          ++position; continue;
        }

        // Add noise.
        if (!std::strcmp("-noise",command_name)) {
          float sigma = 0; char sep = 0; int noise_type = 0;
          if ((std::sscanf(argument,"%f%c",&sigma,&end)==1 ||
               (std::sscanf(argument,"%f%c%c",&sigma,&sep,&end)==2 && sep=='%') ||
               std::sscanf(argument,"%f,%d%c",&sigma,&noise_type,&end)==2 ||
               (std::sscanf(argument,"%f%c,%d%c",&sigma,&sep,&noise_type,&end)==3 && sep=='%')) &&
              sigma>=0 && noise_type>=0 && noise_type<=4) {
            const char *st_type = noise_type==0?"gaussian":noise_type==1?"uniform":noise_type==2?"salt&pepper":
              noise_type==3?"poisson":"rice";
            if (sep=='%') sigma = -sigma;
            print(images,"Add %s noise with standard deviation %g%s to image%s.",
                  st_type,cimg::abs(sigma),sep=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],noise(sigma,noise_type));
          } else arg_error("noise");
          ++position; continue;
        }

        // Rand.
        if (!std::strcmp("-rand",command_name)) {
          double value0 = 0, value1 = 0;
          if (std::sscanf(argument,"%lf,%lf%c",&value0,&value1,&end)==2) {
            print(images,"Fill image%s with random values from range [%g,%g].",gmic_selection,value0,value1);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],rand((T)value0,(T)value1));
          } else arg_error("rand");
          ++position; continue;
        }

        // Compute pointwise norms and orientations.
        gmic_simple_item("-norm",norm,"Compute pointwise L2-norm of pixels in image%s.");
        gmic_simple_item("-orientation",normalize,"Compute pointwise orientation of pixels in image%s.");

        // Map palette.
        if (!std::strcmp("-map",command_name)) {
          unsigned int lut_type = 0; int ind = 0; char sep = 0;
          CImg<T> palette;
          if (std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') {
            gmic_check_indice(ind);
            print(images,"Map palette [%d] on image%s.",ind,gmic_selection);
            palette = images[ind];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],map(palette));
          } else if (std::sscanf(argument,"%u%c",&lut_type,&end)==1 &&
                     lut_type<=2) {
            print(images,"Map %s palette on image%s.",lut_type==0?"default":lut_type==1?"rainbow":"cluster",gmic_selection);
            palette = lut_type==0?CImg<T>::default_LUT256():lut_type==1?CImg<T>::rainbow_LUT256():CImg<T>::contrast_LUT256();
            cimg_forY(selection,l) gmic_apply(images[selection[l]],map(palette));
          } else arg_error("map");
          ++position; continue;
        }

        // Index/cluster image using a palette.
        if (!std::strcmp("-index",command_name)) {
          unsigned int lut_type = 0; int ind = 0, dithering = 0, map_indexes = 0; char sep = 0;
          CImg<T> palette;
          if ((std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') ||
              std::sscanf(argument,"[%d],%d%c",&ind,&dithering,&end)==2 ||
              std::sscanf(argument,"[%d],%d,%d%c",&ind,&dithering,&map_indexes,&end)==3) {
            gmic_check_indice(ind);
            print(images,"Index vector values in image%s by palette [%d], %s dithering%s.",
                  gmic_selection,ind,dithering?"with":"without",map_indexes?" and palette mapping":"");
            palette = images[ind];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],index(palette,dithering?true:false,map_indexes?true:false));
          } else if ((std::sscanf(argument,"%u%c",&lut_type,&end)==1 ||
                      std::sscanf(argument,"%u,%d%c",&lut_type,&dithering,&end)==2 ||
                      std::sscanf(argument,"%u,%d,%d%c",&lut_type,&dithering,&map_indexes,&end)==3) &&
                     lut_type<=2) {
            print(images,"Index vector values in image%s by %s palette, %s dithering%s.",
                  gmic_selection,
                  lut_type==0?"default":lut_type==1?"rainbow":"cluster",dithering?"with":"without",
                  map_indexes?" and index mapping":"");
            palette = lut_type==0?CImg<T>::default_LUT256():lut_type==1?CImg<T>::rainbow_LUT256():CImg<T>::contrast_LUT256();
            cimg_forY(selection,l) gmic_apply(images[selection[l]],index(palette,dithering?true:false,map_indexes?true:false));
          } else arg_error("index");
          ++position; continue;
        }

        //------------------------
        // Color base conversions
        //------------------------
        gmic_simple_item("-rgb2hsv",RGBtoHSV,"Convert image%s from RGB to HSV colorbases.");
        gmic_simple_item("-rgb2hsl",RGBtoHSL,"Convert image%s from RGB to HSL colorbases.");
        gmic_simple_item("-rgb2hsi",RGBtoHSI,"Convert image%s from RGB to HSI colorbases.");
        gmic_simple_item("-rgb2yuv",RGBtoYUV,"Convert image%s from RGB to YUV colorbases.");
        gmic_simple_item("-rgb2ycbcr",RGBtoYCbCr,"Convert image%s from RGB to YCbCr colorbases.");
        gmic_simple_item("-rgb2xyz",RGBtoXYZ,"Convert image%s from RGB to XYZ colorbases.");
        gmic_simple_item("-rgb2lab",RGBtoLab,"Convert image%s from RGB to Lab colorbases.");
        gmic_simple_item("-rgb2cmy",RGBtoCMY,"Convert image%s from RGB to CMY colorbases.");
        gmic_simple_item("-rgb2cmyk",RGBtoCMYK,"Convert image%s from RGB to CMYK colorbases.");
        gmic_simple_item("-cmyk2rgb",CMYKtoRGB,"Convert image%s from CMYK to RGB colorbases.");
        gmic_simple_item("-cmy2rgb",CMYtoRGB,"Convert image%s from CMY to RGB colorbases.");
        gmic_simple_item("-lab2rgb",LabtoRGB,"Convert image%s from Lab to RGB colorbases.");
        gmic_simple_item("-xyz2rgb",XYZtoRGB,"Convert image%s from XYZ to RGB colorbases.");
        gmic_simple_item("-ycbcr2rgb",YCbCrtoRGB,"Convert image%s from YCbCr to RGB colorbases.");
        gmic_simple_item("-yuv2rgb",YUVtoRGB,"Convert image%s from YUV to RGB colorbases.");
        gmic_simple_item("-hsi2rgb",HSItoRGB,"Convert image%s from HSI to RGB colorbases.");
        gmic_simple_item("-hsl2rgb",HSLtoRGB,"Convert image%s from HSL to RGB colorbases.");
        gmic_simple_item("-hsv2rgb",HSVtoRGB,"Convert image%s from HSV to RGB colorbases.");

        //-----------------------
        // Geometric manipulation
        //-----------------------

        // Resize.
        if (!std::strcmp("-resize",command_name) || !std::strcmp("-r",command_name)) {
          char argx[4096] = { 0 }, argy[4096] = { 0 }, argz[4096] = { 0 }, argv[4096] = { 0 };
          char sep = 0, sepx = '%', sepy = '%', sepz = '%', sepv = '%';
          int interpolation = 1, borders = -1, indx = no_ind, indy = no_ind, indz = no_ind, indv = no_ind;
          float valx = 100, valy = 100, valz = 100, valv = 100;
          unsigned int center = 0;
          if (((std::sscanf(argument,"[%d%c%c",&indx,&sep,&end)==2 && sep==']') ||
               std::sscanf(argument,"[%d],%d%c",&indx,&interpolation,&end)==2 ||
               std::sscanf(argument,"[%d],%d,%d%c",&indx,&interpolation,&borders,&end)==3 ||
               std::sscanf(argument,"[%d],%d,%d,%u%c",&indx,&interpolation,&borders,&center,&end)==4) &&
              interpolation>=-1 && interpolation<=5 && borders>=-1 && borders<=2) {
            gmic_check_indice(indx);
            const int
              nvalx = images[indx].width(),
              nvaly = images[indx].height(),
              nvalz = images[indx].depth(),
              nvalv = images[indx].spectrum();
            print(images,"Resize image%s to %dx%dx%dx%d, with %s interpolation.",
                  gmic_selection,nvalx,nvaly,nvalz,nvalv,
                  interpolation<=0?"no":interpolation==1?"nearest neighbor":
                  interpolation==2?"moving average":interpolation==3?"linear":
                  interpolation==4?"grid":"cubic");
            cimg_forY(selection,l) gmic_apply(images[selection[l]],resize(nvalx,nvaly,nvalz,nvalv,interpolation,borders,center?true:false));
            ++position;
          } else if ((std::sscanf(argument,"%4095[][0-9.eE%+-]%c",argx,&end)==1 ||
                      std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",argx,argy,&end)==2 ||
                      std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",argx,argy,argz,&end)==3 ||
                      std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",
                                  argx,argy,argz,argv,&end)==4 ||
                      std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%d%c",
                                  argx,argy,argz,argv,&interpolation,&end)==5 ||
                      std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%d,%d%c",
                                  argx,argy,argz,argv,&interpolation,&borders,&end)==6 ||
                      std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%d,%d,%u%c",
                                  argx,argy,argz,argv,&interpolation,&borders,&center,&end)==7) &&
                     (std::sscanf(argx,"%f%c",&valx,&(sepx=0))==1 ||
                      (std::sscanf(argx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']') ||
                      (std::sscanf(argx,"%f%c%c",&valx,&sepx,&end)==2 && sepx=='%')) &&
                     (!*argy || std::sscanf(argy,"%f%c",&valy,&(sepy=0))==1 ||
                      (std::sscanf(argy,"[%d%c%c",&indy,&sepy,&end)==2 && sepy==']') ||
                      (std::sscanf(argy,"%f%c%c",&valy,&sepy,&end)==2 && sepy=='%')) &&
                     (!*argz || std::sscanf(argz,"%f%c",&valz,&(sepz=0))==1 ||
                      (std::sscanf(argz,"[%d%c%c",&indz,&sepz,&end)==2 && sepz==']') ||
                      (std::sscanf(argz,"%f%c%c",&valz,&sepz,&end)==2 && sepz=='%')) &&
                     (!*argv || std::sscanf(argv,"%f%c",&valv,&(sepv=0))==1 ||
                      (std::sscanf(argv,"[%d%c%c",&indv,&sepv,&end)==2 && sepv==']') ||
                      (std::sscanf(argv,"%f%c%c",&valv,&sepv,&end)==2 && sepv=='%')) &&
                     valx>0 && valy>0 && valz>0 && valv>0 && interpolation>=-1 && interpolation<=5 && borders>=-1 && borders<=2) {
            if (indx!=no_ind) { gmic_check_indice(indx); valx = (float)images[indx].width(); sepx = 0; }
            if (indy!=no_ind) { gmic_check_indice(indy); valy = (float)images[indy].height(); sepy = 0; }
            if (indz!=no_ind) { gmic_check_indice(indz); valz = (float)images[indz].depth(); sepz = 0; }
            if (indv!=no_ind) { gmic_check_indice(indv); valv = (float)images[indv].spectrum(); sepv = 0; }
            print(images,"Resize image%s to %g%s%g%s%g%s%g%s, with %s interpolation.",
                  gmic_selection,valx,sepx=='%'?"%x":"x",valy,sepy=='%'?"%x":"x",valz,
                  sepz=='%'?"%x":"x",valv,sepv=='%'?"% ":"",
                  interpolation<=0?"no":interpolation==1?"nearest neighbor":
                  interpolation==2?"moving average":interpolation==3?"linear":
                  interpolation==4?"grid":"cubic");
            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              const int
                _nvalx = (int)cimg::round(sepx=='%'?valx*img.width()/100:valx,1),
                _nvaly = (int)cimg::round(sepy=='%'?valy*img.height()/100:valy,1),
                _nvalz = (int)cimg::round(sepz=='%'?valz*img.depth()/100:valz,1),
                _nvalv = (int)cimg::round(sepv=='%'?valv*img.spectrum()/100:valv,1),
                nvalx = _nvalx?_nvalx:1,
                nvaly = _nvaly?_nvaly:1,
                nvalz = _nvalz?_nvalz:1,
                nvalv = _nvalv?_nvalv:1;
              gmic_apply(img,resize(nvalx,nvaly,nvalz,nvalv,interpolation,borders,center?true:false));
            }
            ++position;
          } else {
#if cimg_display==0
            print(images,"Resize image%s : run interactive mode (skipped, no display available).",gmic_selection);
#else
            print(images,"Resize image%s : run interactive mode.",gmic_selection);
            CImgDisplay _disp, &disp = instant_window[0]?instant_window[0]:_disp;
            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              if (disp) disp.resize(cimg_fitscreen(img.width(),img.height(),1),false);
              else disp.assign(cimg_fitscreen(img.width(),img.height(),1),0,1);
              std::sprintf(tmpstr,"%s : resize",filenames[selection[l]].data());
              disp.set_title("%s",tmpstr);
              img.get_select(disp,0);
              print(images,"Resize image [%d] to %dx%d, with nearest neighbor interpolation.",selection[l],disp.width(),disp.height());
              gmic_apply(img,resize(disp));
            }
#endif
          }
          continue;
        }

        // Resize2x. and Resize3x.
        gmic_simple_item("-resize2x",resize_doubleXY,"Double size of image%s, using Scale2x algorithm.");
        gmic_simple_item("-resize3x",resize_doubleXY,"Triple size of image%s, using Scale3x algorithm.");

        // Crop.
        if (!std::strcmp("-crop",command_name)) {
          char st0[4096] = { 0 }, st1[4096] = { 0 }, st2[4096] = { 0 }, st3[4096] = { 0 };
          char st4[4096] = { 0 }, st5[4096] = { 0 }, st6[4096] = { 0 }, st7[4096] = { 0 };
          char sep0 = 0, sep1 = 0, sep2 = 0, sep3 = 0, sep4 = 0, sep5 = 0, sep6 = 0, sep7 = 0;
          float a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0, a7 = 0; unsigned int borders = 0;

          if ((std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",st0,st1,&end)==2 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%u%c",st0,st1,&borders,&end)==3) &&
              (std::sscanf(st0,"%f%c",&a0,&end)==1 || (std::sscanf(st0,"%f%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
              (std::sscanf(st1,"%f%c",&a1,&end)==1 || (std::sscanf(st1,"%f%c%c",&a1,&sep1,&end)==2 && sep1=='%'))) {
            print(images,"Crop image%s with (%g%s x (%g%s.",gmic_selection,
                  a0,sep0=='%'?"%)":")",a1,sep1=='%'?"%)":")");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*(img.width()-1)/100:a0,1),
                x1 = (int)cimg::round(sep1=='%'?a1*(img.width()-1)/100:a1,1);
              gmic_apply(img,crop(x0,x1,borders?true:false));
            }
            ++position;
          } else if ((std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",
                                  st0,st1,st2,st3,&end)==4 ||
                      std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%u%c",
                                  st0,st1,st2,st3,&borders,&end)==5) &&
                     (std::sscanf(st0,"%f%c",&a0,&end)==1 || (std::sscanf(st0,"%f%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(st1,"%f%c",&a1,&end)==1 || (std::sscanf(st1,"%f%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (std::sscanf(st2,"%f%c",&a2,&end)==1 || (std::sscanf(st2,"%f%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
                     (std::sscanf(st3,"%f%c",&a3,&end)==1 || (std::sscanf(st3,"%f%c%c",&a3,&sep3,&end)==2 && sep3=='%'))) {
            print(images,"Crop image%s with (%g%s%g%s x (%g%s%g%s.",gmic_selection,
                  a0,sep0=='%'?"%,":",",a1,sep1=='%'?"%)":")",
                  a2,sep2=='%'?"%,":",",a3,sep3=='%'?"%)":")");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*(img.width()-1)/100:a0,1),
                y0 = (int)cimg::round(sep1=='%'?a1*(img.height()-1)/100:a1,1),
                x1 = (int)cimg::round(sep2=='%'?a2*(img.width()-1)/100:a2,1),
                y1 = (int)cimg::round(sep3=='%'?a3*(img.height()-1)/100:a3,1);
              gmic_apply(img,crop(x0,y0,x1,y1,borders?true:false));
            }
            ++position;
          } else if ((std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],"
                                  "%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",
                                  st0,st1,st2,st3,st4,st5,&end)==6 ||
                      std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],"
                                  "%4095[0-9.eE%+-],%4095[0-9.eE%+-],%u%c",
                                  st0,st1,st2,st3,st4,st5,&borders,&end)==7) &&
                     (std::sscanf(st0,"%f%c",&a0,&end)==1 || (std::sscanf(st0,"%f%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(st1,"%f%c",&a1,&end)==1 || (std::sscanf(st1,"%f%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (std::sscanf(st2,"%f%c",&a2,&end)==1 || (std::sscanf(st2,"%f%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
                     (std::sscanf(st3,"%f%c",&a3,&end)==1 || (std::sscanf(st3,"%f%c%c",&a3,&sep3,&end)==2 && sep3=='%')) &&
                     (std::sscanf(st4,"%f%c",&a4,&end)==1 || (std::sscanf(st4,"%f%c%c",&a4,&sep4,&end)==2 && sep4=='%')) &&
                     (std::sscanf(st5,"%f%c",&a5,&end)==1 || (std::sscanf(st5,"%f%c%c",&a5,&sep5,&end)==2 && sep5=='%'))) {
            print(images,"Crop image%s with (%g%s%g%s%g%s x (%g%s%g%s%g%s.",gmic_selection,
                  a0,sep0=='%'?"%,":",",a1,sep1=='%'?"%,":",",a2,sep2=='%'?"%)":")",
                  a3,sep3=='%'?"%,":",",a4,sep4=='%'?"%,":",",a5,sep5=='%'?"%)":")");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*(img.width()-1)/100:a0,1),
                y0 = (int)cimg::round(sep1=='%'?a1*(img.height()-1)/100:a1,1),
                z0 = (int)cimg::round(sep2=='%'?a2*(img.depth()-1)/100:a2,1),
                x1 = (int)cimg::round(sep3=='%'?a3*(img.width()-1)/100:a3,1),
                y1 = (int)cimg::round(sep4=='%'?a4*(img.height()-1)/100:a4,1),
                z1 = (int)cimg::round(sep5=='%'?a5*(img.depth()-1)/100:a5,1);
              gmic_apply(img,crop(x0,y0,z0,x1,y1,z1,borders?true:false));
            }
            ++position;
          } else if ((std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],"
                                  "%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",
                                  st0,st1,st2,st3,st4,st5,st6,st7,&end)==8 ||
                      std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],"
                                  "%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%u%c",
                                  st0,st1,st2,st3,st4,st5,st6,st7,&borders,&end)==9) &&
                     (std::sscanf(st0,"%f%c",&a0,&end)==1 || (std::sscanf(st0,"%f%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(st1,"%f%c",&a1,&end)==1 || (std::sscanf(st1,"%f%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (std::sscanf(st2,"%f%c",&a2,&end)==1 || (std::sscanf(st2,"%f%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
                     (std::sscanf(st3,"%f%c",&a3,&end)==1 || (std::sscanf(st3,"%f%c%c",&a3,&sep3,&end)==2 && sep3=='%')) &&
                     (std::sscanf(st4,"%f%c",&a4,&end)==1 || (std::sscanf(st4,"%f%c%c",&a4,&sep4,&end)==2 && sep4=='%')) &&
                     (std::sscanf(st5,"%f%c",&a5,&end)==1 || (std::sscanf(st5,"%f%c%c",&a5,&sep5,&end)==2 && sep5=='%')) &&
                     (std::sscanf(st6,"%f%c",&a6,&end)==1 || (std::sscanf(st6,"%f%c%c",&a6,&sep6,&end)==2 && sep6=='%')) &&
                     (std::sscanf(st7,"%f%c",&a7,&end)==1 || (std::sscanf(st7,"%f%c%c",&a7,&sep7,&end)==2 && sep7=='%'))) {
            print(images,"Crop image%s with (%g%s%g%s%g%s%g%s x (%g%s%g%s%g%s%g%s.",gmic_selection,
                  a0,sep0=='%'?"%,":",",a1,sep1=='%'?"%,":",",
                  a2,sep2=='%'?"%,":",",a3,sep3=='%'?"%)":")",
                  a4,sep4=='%'?"%,":",",a5,sep5=='%'?"%,":",",
                  a6,sep6=='%'?"%,":",",a7,sep7=='%'?"%)":")");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*(img.width()-1)/100:a0,1),
                y0 = (int)cimg::round(sep1=='%'?a1*(img.height()-1)/100:a1,1),
                z0 = (int)cimg::round(sep2=='%'?a2*(img.depth()-1)/100:a2,1),
                v0 = (int)cimg::round(sep3=='%'?a3*(img.spectrum()-1)/100:a3,1),
                x1 = (int)cimg::round(sep4=='%'?a4*(img.width()-1)/100:a4,1),
                y1 = (int)cimg::round(sep5=='%'?a5*(img.height()-1)/100:a5,1),
                z1 = (int)cimg::round(sep6=='%'?a6*(img.depth()-1)/100:a6,1),
                v1 = (int)cimg::round(sep7=='%'?a7*(img.spectrum()-1)/100:a7,1);
              gmic_apply(img,crop(x0,y0,z0,v0,x1,y1,z1,v1,borders?true:false));
            }
            ++position;
          } else {
#if cimg_display==0
            print(images,"Crop image%s : run interactive mode (skipped, no display available).",gmic_selection);
#else
            print(images,"Crop image%s : run interactive mode.",gmic_selection);
            CImgDisplay _disp, &disp = instant_window[0]?instant_window[0]:_disp;
            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              if (disp) disp.resize(cimg_fitscreen(img.width(),img.height(),1),false);
              else disp.assign(cimg_fitscreen(img.width(),img.height(),1),0,1);
              std::sprintf(tmpstr,"%s : crop",filenames[selection[l]].data());
              disp.set_title("%s",tmpstr);
              const CImg<int> s = img.get_select(disp,2);
              print(images,"Crop image [%d] with (%d,%d,%d) x (%d,%d,%d).",selection[l],s[0],s[1],s[2],s[3],s[4],s[5]);
              gmic_apply(img,crop(s[0],s[1],s[2],s[3],s[4],s[5]));
            }
#endif
          }
          continue;
        }

        // Autocrop.
        if (!std::strcmp("-autocrop",command_name)) {
          print(images,"Auto-crop image%s by color '%s'.",gmic_selection,argument_text);
          cimg_forY(selection,l) {
            CImg<T>& img = images[selection[l]];
            const CImg<T> col = CImg<T>(img.spectrum()).fill(argument,true);
            gmic_apply(img,autocrop(col));
          }
          ++position; continue;
        }

        // Select channels.
        if (!std::strcmp("-channels",command_name)) {
          char sep0 = 0, sep1 = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          float value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%c",arg0,&end)==1 &&
              (std::sscanf(arg0,"%f%c",&value0,&end)==1 ||
               (std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%'))) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].spectrum()-1.0f; sep0 = 0; }
            print(images,"Keep channel %g%s of image%s.",value0,sep0=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.spectrum()-1)/100:value0,1);
              gmic_apply(img,channel(nvalue0));
            }
          } else if (std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
                     (std::sscanf(arg0,"%f%c",&value0,&end)==1 ||
                      (std::sscanf(arg0,"[%d%c%c",&ind0,&sep0,&end)==2 && sep0==']') ||
                      (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(arg1,"%f%c",&value1,&end)==1 ||
                      (std::sscanf(arg1,"[%d%c%c",&ind1,&sep1,&end)==2 && sep1==']') ||
                      (std::sscanf(arg1,"%f%c%c",&value1,&sep1,&end)==2 && sep1=='%'))) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].spectrum()-1.0f; sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1); value1 = images[ind1].spectrum()-1.0f; sep1 = 0; }
            print(images,"Keep channels %g%s..%g%s of image%s.",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.spectrum()-1)/100:value0,1),
                nvalue1 = (int)cimg::round(sep1=='%'?value1*(img.spectrum()-1)/100:value1,1);
              gmic_apply(img,channels(nvalue0,nvalue1));
            }
          } else arg_error("channels");
          ++position; continue;
        }

        // Select slices.
        if (!std::strcmp("-slices",command_name)) {
          char sep0 = 0, sep1 = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          float value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%c",arg0,&end)==1 &&
              (std::sscanf(arg0,"%f%c",&value0,&end)==1 ||
               (std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%'))) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].depth()-1.0f; sep0 = 0; }
            print(images,"Keep slice %g%s of image%s.",value0,sep0=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.depth()-1)/100:value0,1);
              gmic_apply(img,slice(nvalue0));
            }
          } else if (std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
                     (std::sscanf(arg0,"%f%c",&value0,&end)==1 ||
                      (std::sscanf(arg0,"[%d%c%c",&ind0,&sep0,&end)==2 && sep0==']') ||
                      (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(arg1,"%f%c",&value1,&end)==1 ||
                      (std::sscanf(arg1,"[%d%c%c",&ind1,&sep1,&end)==2 && sep1==']') ||
                      (std::sscanf(arg1,"%f%c%c",&value1,&sep1,&end)==2 && sep1=='%'))) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].depth()-1.0f; sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1); value1 = images[ind1].depth()-1.0f; sep1 = 0; }
            print(images,"Keep slices %g%s..%g%s of image%s.",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.depth()-1)/100:value0,1),
                nvalue1 = (int)cimg::round(sep1=='%'?value1*(img.depth()-1)/100:value1,1);
              gmic_apply(img,slices(nvalue0,nvalue1));
            }
          } else arg_error("slices");
          ++position; continue;
        }

        // Select lines.
        if (!std::strcmp("-lines",command_name)) {
          char sep0 = 0, sep1 = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          float value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%c",arg0,&end)==1 &&
              (std::sscanf(arg0,"%f%c",&value0,&end)==1 ||
               (std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%'))) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].height()-1.0f; sep0 = 0; }
            print(images,"Keep line %g%s of image%s.",value0,sep0=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.height()-1)/100:value0,1);
              gmic_apply(img,line(nvalue0));
            }
          } else if (std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
                     (std::sscanf(arg0,"%f%c",&value0,&end)==1 ||
                      (std::sscanf(arg0,"[%d%c%c",&ind0,&sep0,&end)==2 && sep0==']') ||
                      (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(arg1,"%f%c",&value1,&end)==1 ||
                      (std::sscanf(arg1,"[%d%c%c",&ind1,&sep1,&end)==2 && sep1==']') ||
                      (std::sscanf(arg1,"%f%c%c",&value1,&sep1,&end)==2 && sep1=='%'))) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].height()-1.0f; sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1); value1 = images[ind1].height()-1.0f; sep1 = 0; }
            print(images,"Keep lines %g%s..%g%s of image%s.",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.height()-1)/100:value0,1),
                nvalue1 = (int)cimg::round(sep1=='%'?value1*(img.height()-1)/100:value1,1);
              gmic_apply(img,lines(nvalue0,nvalue1));
            }
          } else arg_error("lines");
          ++position; continue;
        }

        // Select columns.
        if (!std::strcmp("-columns",command_name)) {
          char sep0 = 0, sep1 = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          float value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%c",arg0,&end)==1 &&
              (std::sscanf(arg0,"%f%c",&value0,&end)==1 ||
               (std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%'))) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].width()-1.0f; sep0 = 0; }
            print(images,"Keep column %g%s of image%s.",value0,sep0=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.width()-1)/100:value0,1);
              gmic_apply(img,column(nvalue0));
            }
          } else if (std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
                     (std::sscanf(arg0,"%f%c",&value0,&end)==1 ||
                      (std::sscanf(arg0,"[%d%c%c",&ind0,&sep0,&end)==2 && sep0==']') ||
                      (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(arg1,"%f%c",&value1,&end)==1 ||
                      (std::sscanf(arg1,"[%d%c%c",&ind1,&sep1,&end)==2 && sep1==']') ||
                      (std::sscanf(arg1,"%f%c%c",&value1,&sep1,&end)==2 && sep1=='%'))) {
            if (ind0!=no_ind) { gmic_check_indice(ind0); value0 = images[ind0].width()-1.0f; sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1); value1 = images[ind1].width()-1.0f; sep1 = 0; }
            print(images,"Keep columns %g%s..%g%s of image%s.",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.width()-1)/100:value0,1),
                nvalue1 = (int)cimg::round(sep1=='%'?value1*(img.width()-1)/100:value1,1);
              gmic_apply(img,columns(nvalue0,nvalue1));
            }
          } else arg_error("columns");
          ++position; continue;
        }

        // Rotate.
        if (!std::strcmp("-rotate",command_name)) {
          float angle = 0, zoom = 1, cx = 0, cy = 0; int borders = 0, interpolation = 1;
          char argx[4096] = { 0 }, argy[4096] = { 0 }, sepx = 0, sepy = 0;
          if ((std::sscanf(argument,"%f%c",&angle,&end)==1 ||
               std::sscanf(argument,"%f,%d%c",&angle,&borders,&end)==2 ||
               std::sscanf(argument,"%f,%d,%d%c",&angle,&borders,&interpolation,&end)==3 ||
               std::sscanf(argument,"%f,%d,%d,%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",&angle,&borders,&interpolation,argx,argy,&end)==5 ||
               std::sscanf(argument,"%f,%d,%d,%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f%c",&angle,&borders,&interpolation,argx,argy,&zoom,&end)==6) &&
              (!*argx || std::sscanf(argx,"%f%c",&cx,&end)==1 || (std::sscanf(argx,"%f%c%c",&cx,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy || std::sscanf(argy,"%f%c",&cy,&end)==1 || (std::sscanf(argy,"%f%c%c",&cy,&sepy,&end)==2 && sepy=='%')) &&
              borders>=0 && borders<=2 && interpolation>=0 && interpolation<=2) {
            if (*argx) {
              print(images,"Rotate image%s of %g deg., %s borders, %s interpolation and center (%g%s,%g%s).",
                    gmic_selection,angle,
                    borders==0?"black":borders==1?"nearest":"cyclic",
                    interpolation==0?"nearest-neighbor":interpolation==1?"linear":"cubic",
                    cx,sepx=='%'?"%":"",cy,sepy=='%'?"%":"");
              cimg_forY(selection,l) {
                CImg<T> &img = images[selection[l]];
                const float
                  ncx = sepx=='%'?cx*(img.width()-1)/100:cx,
                  ncy = sepy=='%'?cy*(img.height()-1)/100:cy;
                gmic_apply(img,rotate(angle,ncx,ncy,zoom,borders,interpolation));
              }
            } else {
              print(images,"Rotate image%s of %g deg., %s borders and %s interpolation.",
                    gmic_selection,angle,
                    borders==0?"black":borders==1?"nearest":"cyclic",
                    interpolation==0?"nearest-neighbor":interpolation==1?"linear":"cubic");
              cimg_forY(selection,l) gmic_apply(images[selection[l]],rotate(angle,borders,interpolation));
            }
          } else arg_error("rotate");
          ++position; continue;
        }

        // Mirror.
        if (!std::strcmp("-mirror",command_name)) {
          const char axis = cimg::uncase(*argument);
          if (std::strlen(argument)==1 &&
              (axis=='x' || axis=='y' || axis=='z' || axis=='c')) {
            print(images,"Mirror image%s along the %c-axis.",gmic_selection,axis);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],mirror(axis));
          } else arg_error("mirror");
          ++position; continue;
        }

        // Shift.
        if (!std::strcmp("-shift",command_name)) {
          char argx[4096] = { 0 }, argy[4096] = { 0 }, argz[4096] = { 0 }, argv[4096] = { 0 };
          char sepx = 0, sepy = 0, sepz = 0, sepv = 0;
          float dx = 0, dy = 0, dz = 0, dv = 0; unsigned int borders = 0;
          if ((std::sscanf(argument,"%4095[0-9.eE%+-]%c",argx,&end)==1 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",argx,argy,&end)==2 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",argx,argy,argz,&end)==3 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",
                           argx,argy,argz,argv,&end)==4 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%u%c",
                           argx,argy,argz,argv,&borders,&end)==5) &&
              (!*argx || std::sscanf(argx,"%f%c",&dx,&end)==1 || (std::sscanf(argx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy || std::sscanf(argy,"%f%c",&dy,&end)==1 || (std::sscanf(argy,"%f%c%c",&dy,&sepy,&end)==2 && sepy=='%')) &&
              (!*argz || std::sscanf(argz,"%f%c",&dz,&end)==1 || (std::sscanf(argz,"%f%c%c",&dz,&sepz,&end)==2 && sepz=='%')) &&
              (!*argv || std::sscanf(argv,"%f%c",&dv,&end)==1 || (std::sscanf(argv,"%f%c%c",&dv,&sepv,&end)==2 && sepv=='%')) &&
              borders<=2) {
            print(images,"Shift image%s by translation vector (%g%s,%g%s,%g%s,%g%s).",
                  gmic_selection,dx,sepx=='%'?"%":"",dy,sepy=='%'?"%":"",dz,sepz=='%'?"%":"",dv,sepv=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                ndx = (int)cimg::round(sepx=='%'?dx*img.width()/100:dx,1),
                ndy = (int)cimg::round(sepy=='%'?dy*img.height()/100:dy,1),
                ndz = (int)cimg::round(sepz=='%'?dz*img.depth()/100:dz,1),
                ndv = (int)cimg::round(sepv=='%'?dv*img.spectrum()/100:dv,1);
              gmic_apply(images[selection[l]],shift(ndx,ndy,ndz,ndv,borders));
            }
          } else arg_error("shift");
          ++position; continue;
        }

        // Transpose.
        gmic_simple_item("-transpose",transpose,"Transpose image%s.");

        // Invert.
        gmic_simple_item("-invert",invert,"Compute matrix inverse of image%s.");

        // Permute axes.
        if (!std::strcmp("-permute",command_name)) {
          print(images,"Permute axes of image%s, with permutation '%s'.",gmic_selection,argument_text);
          cimg_forY(selection,l) gmic_apply(images[selection[l]],permute_axes(argument));
          ++position; continue;
        }

        // Unroll.
        if (!std::strcmp("-unroll",command_name)) {
          const char axis = cimg::uncase(*argument);
          if (std::strlen(argument)==1 &&
              (axis=='x' || axis=='y' || axis=='z' || axis=='c')) {
            print(images,"Unroll image%s along the %c-axis.",gmic_selection,axis);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],unroll(axis));
          } else arg_error("unroll");
          ++position; continue;
        }

        // Split images.
        if (!std::strcmp("-split",command_name) || !std::strcmp("-s",command_name)) {
          char axis = cimg::uncase(*argument), foo = 0; double value = 0;
          int nb = 0, keep_value = 0, px = 0, py = 0, pz = 0, pv = 0; unsigned int borders = 0;
          if ((std::sscanf(argument,"%c%c",&foo,&end)==1 ||
               std::sscanf(argument,"%c,%d%c",&foo,&nb,&end)==2) &&
              (axis=='x' || axis=='y' || axis=='z' || axis=='c')) {
            if (nb>0) print(images,"Split image%s along the %c-axis into %d parts.",gmic_selection,axis,nb);
            else if (nb<0) print(images,"Split image%s along the %c-axis into blocs of %d pixels.",gmic_selection,axis,-nb);
            else print(images,"Split image%s along the %c-axis.",gmic_selection,axis);
            unsigned int off = 0;
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l] + off;
              const CImg<T>& img = images[ind];
              const CImg<char> filename = filenames[ind].get_mark_filename();
              const char naxis = cimg::uncase(axis);
              if (!((naxis=='x' && img.width()==1) ||
                    (naxis=='y' && img.height()==1) ||
                    (naxis=='z' && img.depth()==1) ||
                    (naxis=='c' && img.spectrum()==1)) || get_version) {
                CImgList<T> split = img.get_split(axis,nb);
                if (get_version) {
                  filenames.insert(split.size(),filename);
                  split.move_to(images,~0U);
                } else {
                  off+=split.size() - 1;
                  filenames.remove(ind); filenames.insert(split.size(),filename,ind);
                  images.remove(ind); split.move_to(images,ind);
                }
              }
            }
          } else if ((std::sscanf(argument,"%d,%u%c",&px,&borders,&end)==2 && px>0) ||
                     (std::sscanf(argument,"%d,%d,%u%c",&px,&py,&borders,&end)==3 && px>0 && py>0) ||
                     (std::sscanf(argument,"%d,%d,%d,%u%c",&px,&py,&pz,&borders,&end)==4 && px>0 && py>0 && pz>0) ||
                     (std::sscanf(argument,"%d,%d,%d,%d,%u%c",&px,&py,&pz,&pv,&borders,&end)==5 && px>0 && py>0 && pz>0 && pv>0)) {
            unsigned int pdim = 0;
            const char *const s_borders = borders?"Neumann":"Dirichlet";
            if (pv) { print(images,"Split image%s into %dx%dx%dx%d patchs, with %s border conditions.",
                            gmic_selection,px,py,pz,pv,s_borders); pdim = 4; }
            else if (pz) { print(images,"Split image%s into %dx%dx%d patchs, with %s border conditions.",
                                 gmic_selection,px,py,pz,s_borders); pdim = 3; }
            else if (py) { print(images,"Split image%s into %dx%d patchs, with %s border conditions.",
                                 gmic_selection,px,py,s_borders); pdim = 2; }
            else { print(images,"Split image%s into %d patchs, with %s border conditions.",
                         gmic_selection,px,s_borders); pdim = 1; }
            unsigned int off = 0;
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l] + off;
              const CImg<char> filename = filenames[ind].get_mark_filename();
              CImgList<T> split = images[ind].get_split_patch(px,py,pz,pv,borders);
              if (get_version) {
                filenames.insert(split.size(),filename);
                split.move_to(images,~0U);
              } else {
                off+=split.size() - 1;
                filenames.remove(ind); filenames.insert(split.size(),filename,ind);
                images.remove(ind); split.move_to(images,ind);
              }
            }
          } else if ((std::sscanf(argument,"%lf%c",&value,&end)==1 || std::sscanf(argument,"%lf,%d%c",&value,&keep_value,&end)==2) &&
                     (keep_value==0 || keep_value==1)) {
            print(images,"Split image%s according to value %g.",gmic_selection,value);
            unsigned int off = 0;
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l] + off;
              const CImg<char> filename = filenames[ind].get_mark_filename();
              CImgList<T> split = images[ind].get_split((T)value,keep_value,false);
              if (get_version) {
                filenames.insert(split.size(),filename);
                split.move_to(images,~0U);
              } else {
                off+=split.size() - 1;
                filenames.remove(ind); filenames.insert(split.size(),filename,ind);
                images.remove(ind); split.move_to(images,ind);
              }
            }
          } else arg_error("split");
          ++position; continue;
        }

        // Append images.
        if (!std::strcmp("-append",command_name) || !std::strcmp("-a",command_name)) {
          char axis = 0, align='p';
          if ((std::sscanf(argument,"%c%c",&axis,&end)==1 ||
               std::sscanf(argument,"%c,%c%c",&axis,&align,&end)==2) &&
              (axis=='x' || axis=='y' || axis=='z' || axis=='c') &&
              (align=='p' || align=='c' || align=='n')) {
            axis = cimg::uncase(axis);
            print(images,"Append image%s along the %c-axis with %s alignment.",
                  gmic_selection,axis,align=='p'?"left":align=='c'?"center":"right");
            CImgList<T> subimages; cimg_forY(selection,l) subimages.insert(images[selection[l]],~0U,true);
            const CImg<char> filename = filenames[selection[0]].get_mark_filename();
            if (get_version) {
              filenames.insert(filename);
              subimages.get_append(axis,align).move_to(images,~0U);
            } else {
              filenames.insert(filename,selection[0]);
              subimages.get_append(axis,align).move_to(images,selection[0]);
              int off = 1;
              cimg_forY(selection,l) {
                const int ind = selection[l] + off;
                images.remove(ind); filenames.remove(ind);
                --off;
              }
            }
          } else arg_error("append");
          ++position; continue;
        }

        // Warp images.
        if (!std::strcmp("-warp",command_name)) {
          int ind0 = no_ind, nb_frames = 1; unsigned int interpolation = 1, relative = 0, borders = 1; char sep = 0;
          if (((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']')||
               std::sscanf(argument,"[%d],%u%c",&ind0,&relative,&end)==2 ||
               std::sscanf(argument,"[%d],%u,%u%c",&ind0,&relative,&interpolation,&end)==3 ||
               std::sscanf(argument,"[%d],%u,%u,%u%c",&ind0,&relative,&interpolation,&borders,&end)==4 ||
               std::sscanf(argument,"[%d],%u,%u,%u,%d%c",&ind0,&relative,&interpolation,&borders,&nb_frames,&end)==5) &&
              borders<=2 && nb_frames>=1) {
            gmic_check_indice(ind0);
            if (nb_frames>1) print(images,"Warp image%s with %s field [%u] and %d frames.",
                                   gmic_selection,relative?"relative":"absolute",ind0,nb_frames);
            else print(images,"Warp image%s with %s field [%u].",gmic_selection,relative?"relative":"absolute",ind0);
            const CImg<T> warp = images[ind0];
            unsigned int off = 0;
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l] + off;
              CImg<T> &img = images[ind];
              const CImg<char> filename = filenames[ind].get_mark_filename();
              CImgList<T> frames(nb_frames);
              cimglist_for(frames,t) frames[t] = img.get_warp(warp*((t+1.0f)/nb_frames),
                                                              relative?true:false,interpolation?true:false,borders);
              if (get_version) {
                filenames.insert(nb_frames,filename);
                frames.move_to(images,~0U);
              } else {
                off+=nb_frames - 1;
                filenames.insert(nb_frames-1,filename,ind);
                images.remove(ind); frames.move_to(images,ind);
              }
            }
          } else arg_error("warp");
          ++position; continue;
        }

        //-----------------------
        // Image filtering
        //-----------------------

        // Quasi-gaussian recursive blur.
        if (!std::strcmp("-blur",command_name)) {
          float sigma = -1; unsigned int borders = 1; char sep = 0;
          if ((std::sscanf(argument,"%f%c",&sigma,&end)==1 ||
               (std::sscanf(argument,"%f%c%c",&sigma,&sep,&end)==2 && sep=='%') ||
               std::sscanf(argument,"%f,%u%c",&sigma,&borders,&end)==2 ||
               (std::sscanf(argument,"%f%c,%u%c",&sigma,&sep,&borders,&end)==3 && sep=='%')) &&
              sigma>=0) {
            print(images,"Blur image%s with standard deviation %g%s.",gmic_selection,sigma,sep=='%'?"%":"");
            if (sep=='%') sigma = -sigma;
            cimg_forY(selection,l) gmic_apply(images[selection[l]],blur(sigma,borders?true:false));
          } else arg_error("blur");
          ++position; continue;
        }

        // Bilateral filter.
        if (!std::strcmp("-bilateral",command_name)) {
          float sigma_s = 0, sigma_r = 0; char sep =  0;
          if ((std::sscanf(argument,"%f,%f%c",&sigma_s,&sigma_r,&end)==2 ||
               (std::sscanf(argument,"%f%c,%f%c",&sigma_s,&sep,&sigma_r,&end)==3 && sep=='%')) &&
              sigma_s>=0 && sigma_r>=0) {
            print(images,"Apply bilateral filter on image%s with standard deviations %g%s and %g.",
                  gmic_selection,sigma_s,sep=='%'?"%":"",sigma_r);
            if (sep=='%') sigma_s = -sigma_s;
            cimg_forY(selection,l) gmic_apply(images[selection[l]],blur_bilateral(sigma_s,sigma_r));
          } else arg_error("bilateral");
          ++position; continue;
        }

        // Patch averaging.
        if (!std::strcmp("-denoise",command_name)) {
          float sigma_s = 10, sigma_r = 10, smoothness = 1; int psize = 5, rsize = 6;
          unsigned int fast_approximation = 0;
          if ((std::sscanf(argument,"%f%c",&sigma_s,&end)==1 ||
               std::sscanf(argument,"%f,%f%c",&sigma_s,&sigma_r,&end)==2 ||
               std::sscanf(argument,"%f,%f,%d%c",&sigma_s,&sigma_r,&psize,&end)==3 ||
               std::sscanf(argument,"%f,%f,%d,%d%c",&sigma_s,&sigma_r,&psize,&rsize,&end)==4 ||
               std::sscanf(argument,"%f,%f,%d,%d,%f%c",&sigma_s,&sigma_r,&psize,&rsize,&smoothness,&end)==5 ||
               std::sscanf(argument,"%f,%f,%d,%d,%f,%u%c",&sigma_s,&sigma_r,&psize,&rsize,&smoothness,&fast_approximation,&end)==6) &&
              sigma_s>=0 && sigma_r>=0 && psize>0 && rsize>0) {
            print(images,"Denoise image%s with %dx%d patchs, standard deviations %lg,%g, lookup size %d and smoothness %g.",
                  gmic_selection,psize,psize,sigma_s,sigma_r,rsize,smoothness);
            cimg_forY(selection,l)
              gmic_apply(images[selection[l]],blur_patch(sigma_s,sigma_r,psize,rsize,smoothness,fast_approximation?true:false));
          } else arg_error("denoise");
          ++position; continue;
        }

        // Smooth.
        if (!std::strcmp("-smooth",command_name)) {
          float amplitude = 0, sharpness = 0.7f, anisotropy = 0.3f, alpha = 0.6f, sigma = 1.1f, dl =0.8f, da = 30.0f, gauss_prec = 2.0f;
          int ind = no_ind; unsigned int interpolation_type = 0, fast_approx = 1; char sep = 0;
          if ((std::sscanf(argument,"%f%c",&amplitude,&end)==1 ||
               std::sscanf(argument,"%f,%f%c",&amplitude,&sharpness,&end)==2 ||
               std::sscanf(argument,"%f,%f,%f%c",&amplitude,&sharpness,&anisotropy,&end)==3 ||
               std::sscanf(argument,"%f,%f,%f,%f%c",&amplitude,&sharpness,&anisotropy,&alpha,&end)==4 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f%c",&amplitude,&sharpness,&anisotropy,&alpha,&sigma,&end)==5 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f,%f%c",&amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&end)==6 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f,%f,%f%c",&amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&da,&end)==7 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f,%f,%f,%f%c",
                           &amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&da,&gauss_prec,&end)==8 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f,%f,%f,%f,%u%c",
                           &amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&da,&gauss_prec,&interpolation_type,&end)==9 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f,%f,%f,%f,%u,%u%c",
                           &amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&da,&gauss_prec,&interpolation_type,&fast_approx,&end)==10) &&
              amplitude>=0 && sharpness>=0 && anisotropy>=0 && anisotropy<=1 && dl>0 && da>=0 && gauss_prec>0 &&
              interpolation_type<=2) {
            if (da>0)
              print(images,"Smooth image%s anisotropically with amplitude %g, sharpness %g, anisotropy %g, alpha %g and sigma %g.",
                    gmic_selection,amplitude,sharpness,anisotropy,alpha,sigma);
            else
              print(images,"Smooth image%s anisotropically with %d iterations, sharpness %g, anisotropy %g, alpha %g and sigma %g.",
                    gmic_selection,(int)amplitude,sharpness,anisotropy,alpha,sigma);
            cimg_forY(selection,l)
              gmic_apply(images[selection[l]],blur_anisotropic(amplitude,sharpness,anisotropy,alpha,sigma,
                                                             dl,da,gauss_prec,interpolation_type,fast_approx?true:false));
          } else if (((std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') ||
                      std::sscanf(argument,"[%d],%f%c",&ind,&amplitude,&end)==2 ||
                      std::sscanf(argument,"[%d],%f,%f%c",&ind,&amplitude,&dl,&end)==3 ||
                      std::sscanf(argument,"[%d],%f,%f,%f%c",&ind,&amplitude,&dl,&da,&end)==4 ||
                      std::sscanf(argument,"[%d],%f,%f,%f,%f%c",&ind,&amplitude,&dl,&da,&gauss_prec,&end)==5 ||
                      std::sscanf(argument,"[%d],%f,%f,%f,%f,%u%c",&ind,&amplitude,&dl,&da,&gauss_prec,&interpolation_type,&end)==5 ||
                      std::sscanf(argument,"[%d],%f,%f,%f,%f,%u,%u%c",
                                  &ind,&amplitude,&dl,&da,&gauss_prec,&interpolation_type,&fast_approx,&end)==6) &&
                     amplitude>=0 && dl>0 && da>=0 && gauss_prec>0 && interpolation_type<=2) {
            gmic_check_indice(ind);
            const CImg<T> tensors = images[ind];
            if (da>0)
              print(images,"Smooth image%s anisotropically with tensor field [%d] and amplitude %g.",gmic_selection,ind,amplitude);
            else
              print(images,"Smooth image%s anisotropically with tensor field [%d] and %d iterations.",gmic_selection,ind,(int)amplitude);
            cimg_forY(selection,l)
              gmic_apply(images[selection[l]],blur_anisotropic(tensors,amplitude,dl,da,gauss_prec,interpolation_type,fast_approx));
          } else arg_error("smooth");
          ++position; continue;
        }

        // Get tensor geometry from image.
        if (!std::strcmp("-edgetensors",command_name)) {
          float sharpness = 0.7f, anisotropy = 0.3f, alpha = 0.6f, sigma = 1.1f; unsigned int is_sqrt = 0;
          if ((std::sscanf(argument,"%f%c",&sharpness,&end)==1 ||
               std::sscanf(argument,"%f,%f%c",&sharpness,&anisotropy,&end)==2 ||
               std::sscanf(argument,"%f,%f,%f%c",&sharpness,&anisotropy,&alpha,&end)==3 ||
               std::sscanf(argument,"%f,%f,%f,%f%c",&sharpness,&anisotropy,&alpha,&sigma,&end)==4 ||
               std::sscanf(argument,"%f,%f,%f,%f,%u%c",&sharpness,&anisotropy,&alpha,&sigma,&is_sqrt,&end)==5) &&
              sharpness>=0 && anisotropy>=0 && anisotropy<=1) {
            print(images,"Compute %stensors for edge-preserving smoothing of image%s, with sharpness %g, anisotropy %g,"
                  "alpha %g and sigma %g.",
                  is_sqrt?"square root of ":"",gmic_selection,sharpness,anisotropy,alpha,sigma);
            cimg_forY(selection,l)
              gmic_apply(images[selection[l]],edge_tensors(sharpness,anisotropy,alpha,sigma,is_sqrt?true:false));
          } else arg_error("edgetensors");
          ++position; continue;
        }

        // Median filter.
        if (!std::strcmp("-median",command_name)) {
          int siz = 3;
          if (std::sscanf(argument,"%d%c",&siz,&end)==1 &&
              siz>=0) {
            print(images,"Apply median filter of size %d on image%s.",siz,gmic_selection);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],blur_median(siz));
          } else arg_error("median");
          ++position; continue;
        }

        // Sharpen.
        if (!std::strcmp("-sharpen",command_name)) {
          float amplitude = 0, edge = 1, alpha = 0, sigma = 0; unsigned int sharpen_type = 0;
          if ((std::sscanf(argument,"%f%c",&amplitude,&end)==1 ||
               std::sscanf(argument,"%f,%u%c",&amplitude,&sharpen_type,&end)==2 ||
               std::sscanf(argument,"%f,%u,%f%c",&amplitude,&sharpen_type,&edge,&end)==3 ||
               std::sscanf(argument,"%f,%u,%f,%f%c",&amplitude,&sharpen_type,&edge,&alpha,&end)==4 ||
               std::sscanf(argument,"%f,%u,%f,%f,%f%c",&amplitude,&sharpen_type,&edge,&alpha,&sigma,&end)==5) &&
              amplitude>=0 && edge>=0) {
            if (sharpen_type) print(images,"Sharpen image%s with shock filters and amplitude %g, edge %g, alpha %g and sigma %g.",
                                    gmic_selection,amplitude,edge,alpha,sigma);
            else print(images,"Sharpen image%s with inverse diffusion and amplitude %g.",gmic_selection,amplitude);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],sharpen(amplitude,sharpen_type?true:false,edge,alpha,sigma));
          } else arg_error("sharpen");
          ++position; continue;
        }

        // Convolve.
        if (!std::strcmp("-convolve",command_name)) {
          int ind = no_ind; unsigned int borders = 1; char sep = 0;
          if ((std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') ||
               std::sscanf(argument,"[%d],%u%c",&ind,&borders,&end)==2) {
            gmic_check_indice(ind);
            print(images,"Convolve image%s with mask [%d].",gmic_selection,ind);
            const CImg<T> mask = images[ind];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],convolve(mask,borders));
          } else arg_error("convolve");
          ++position; continue;
        }

        // Correlate.
        if (!std::strcmp("-correlate",command_name)) {
          int ind = no_ind; unsigned int borders = 1; char sep = 0;
          if ((std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') ||
              std::sscanf(argument,"[%d],%u%c",&ind,&borders,&end)==2) {
            gmic_check_indice(ind);
            print(images,"Correlate image%s with mask [%d].",gmic_selection,ind);
            const CImg<T> mask = images[ind];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],correlate(mask,borders));
          } else arg_error("correlate");
          ++position; continue;
        }

        // Erode.
        if (!std::strcmp("-erode",command_name)) {
          int sx = 3, sy = 3, sz = 1, ind = no_ind; unsigned int borders = 1; char sep = 0;
          if (((std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') ||
               std::sscanf(argument,"[%d],%u%c",&ind,&borders,&end)==2) &&
              sx>=0 && sy>=0 && sz>=0) {
            gmic_check_indice(ind);
            print(images,"Erode image%s with mask [%d].",gmic_selection,ind);
            const CImg<T> mask = images[ind];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],erode(mask,borders));
          } else if ((std::sscanf(argument,"%d%c",&sx,&end)==1) &&
                     sx>=0) {
            print(images,"Erode image%s with a square mask of size %d.",gmic_selection,sx);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],erode((unsigned int)sx));
          } else if ((std::sscanf(argument,"%d,%d%c",&sx,&sy,&end)==2 ||
                      std::sscanf(argument,"%d,%d,%d%c",&sx,&sy,&sz,&end)==3) &&
                     sx>=0 && sy>=0 && sz>=0) {
            print(images,"Erode image%s with a rectangular mask of size %dx%dx%d.",gmic_selection,sx,sy,sz);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],erode((unsigned int)sx,(unsigned int)sy,(unsigned int)sz));
          } else arg_error("erode");
          ++position; continue;
        }

        // Dilate.
        if (!std::strcmp("-dilate",command_name)) {
          int sx = 3, sy = 3, sz = 1, ind = no_ind; unsigned int borders = 1; char sep = 0;
          if (((std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') ||
               std::sscanf(argument,"[%d],%u%c",&ind,&borders,&end)==2) &&
              sx>=0 && sy>=0 && sz>=0) {
            gmic_check_indice(ind);
            print(images,"Dilate image%s with mask [%d].",gmic_selection,ind);
            const CImg<T> mask = images[ind];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],dilate(mask,borders));
          } else if ((std::sscanf(argument,"%d%c",&sx,&end)==1) &&
                     sx>=0) {
            print(images,"Dilate image%s with a square mask of size %d.",gmic_selection,sx);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],dilate((unsigned int)sx));
          } else if ((std::sscanf(argument,"%d,%d%c",&sx,&sy,&end)==2 ||
                      std::sscanf(argument,"%d,%d,%d%c",&sx,&sy,&sz,&end)==3) &&
                     sx>=0 && sy>=0 && sz>=0) {
            print(images,"Dilate image%s with a rectangular mask of size %dx%dx%d.",gmic_selection,sx,sy,sz);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],dilate((unsigned int)sx,(unsigned int)sy,(unsigned int)sz));
          } else arg_error("dilate");
          ++position; continue;
        }

        // Inpaint
        if (!std::strcmp("-inpaint",command_name)) {
          int ind = no_ind; char sep = 0;
          if (std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') {
            gmic_check_indice(ind);
            print(images,"Inpaint image%s with mask [%d].",gmic_selection,ind);
            CImg<T> mask = images[ind];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],inpaint(mask));
          } else arg_error("inpaint");
          ++position; continue;
        }

        // Compute gradient.
        if (!std::strcmp("-gradient",command_name)) {
          char axes[4096] = { 0 }, *naxes = 0; int scheme = 3;
          print(images,"Compute gradient of image%s.",gmic_selection);
          if (std::sscanf(argument,"%4095[xyz]%c",axes,&end)==1 ||
              std::sscanf(argument,"%4095[xyz],%d%c",axes,&scheme,&end)==2) { naxes = axes; ++position; }
          unsigned int off = 0;
          cimg_forY(selection,l) {
            const unsigned int ind = selection[l] + off;
            CImg<T>& img = images[ind];
            const CImg<char> filename = filenames[ind].get_mark_filename();
            CImgList<T> gradient = img.get_gradient(naxes,scheme);
            if (get_version) {
              filenames.insert(gradient.size(),filename);
              gradient.move_to(images,~0U);
            } else {
              off+=gradient.size() - 1;
              filenames.remove(ind); filenames.insert(gradient.size(),filename,ind);
              images.remove(ind); gradient.move_to(images,ind);
            }
          }
          continue;
        }

        // Compute Hessian.
        if (!std::strcmp("-hessian",command_name)) {
          char axes[4096] = { 0 }, *naxes = 0;
          print(images,"Compute Hessian of image%s.",gmic_selection);
          if (std::sscanf(argument,"%4095[xyz]%c",axes,&end)==1) { naxes = axes; ++position; }
          unsigned int off = 0;
          cimg_forY(selection,l) {
            const unsigned int ind = selection[l] + off;
            CImg<T>& img = images[ind];
            const CImg<char> filename = filenames[ind].get_mark_filename();
            CImgList<T> hessian = img.get_hessian(naxes);
            if (get_version) {
              filenames.insert(hessian.size(),filename);
              hessian.move_to(images,~0U);
            } else {
              off+=hessian.size() - 1;
              filenames.remove(ind); filenames.insert(hessian.size(),filename,ind);
              images.remove(ind); hessian.move_to(images,ind);
            }
          }
          continue;
        }

        // Compute Haar transform.
        const bool inv_haar = !std::strcmp("-ihaar",command_name);
        if (!std::strcmp("-haar",command_name) || inv_haar) {
          int nb_scales = 0;
          if (std::sscanf(argument,"%d%c",&nb_scales,&end)==1 &&
              nb_scales>0) {
            print(images,"Compute %sHaar Transform of image%s",inv_haar?"inverse ":"",gmic_selection);
            cimg_forY(selection,l) images[selection[l]].haar(inv_haar,nb_scales);
          } else arg_error(command_name);
          ++position; continue;
        }

        // Compute direct or inverse FFT.
        const bool inv_fft = !std::strcmp("-ifft",command_name);
        if (!std::strcmp("-fft",command_name) || inv_fft) {
          print(images,"Compute %sFourier Transform of complex data",inv_fft?"inverse ":"");
          cimg_forY(selection,l) {
            const unsigned int ind0 = selection[l], ind1 = l+1<selection.height()?selection[l+1]:~0U;
            if (ind1!=~0U) {
              if (verbosity_level>=0) {
                std::fprintf(cimg::output()," ([%u],[%u])%c",ind0,ind1,l==selection.height()-1?'.':',');
                std::fflush(cimg::output());
              }
              CImgList<T> fft(images[ind0],images[ind1],!get_version);
              fft.FFT(inv_fft);
              if (get_version) {
                filenames.insert(filenames[ind0].get_mark_filename());
                filenames.insert(filenames[ind1].get_mark_filename());
                fft.move_to(images,~0U);
              } else {
                filenames[ind0].mark_filename();
                filenames[ind1].mark_filename();
                fft[0].move_to(images[ind0]);
                fft[1].move_to(images[ind1]);
              }
              ++l;
            } else {
              if (verbosity_level>=0) {
                std::fprintf(cimg::output()," ([%u],0)",ind0);
                std::fflush(cimg::output());
              }
              CImgList<T> fft(images[ind0],!get_version);
              fft.insert(fft[0]);
              fft[1].fill(0);
              fft.FFT(inv_fft);
              if (get_version) {
                filenames.insert(2,filenames[ind0].get_mark_filename());
                fft.move_to(images,~0U);
              } else {
                filenames[ind0].mark_filename(); filenames.insert(filenames[ind0],ind0+1);
                fft[0].move_to(images[ind0]);
                images.insert(fft[1],ind0+1);
              }
            }
          }
          continue;
        }

        //-----------------------------
        // Image creation and drawing
        //-----------------------------

        // Histogram.
        if (!std::strcmp("-histogram",command_name)) {
          int nb_levels = 256; double vmin = 0, vmax = 0; char sep = 0, sepm = 0, sepM = 0;
          if ((std::sscanf(argument,"%d%c",&nb_levels,&end)==1 ||
               (std::sscanf(argument,"%d%c%c",&nb_levels,&sep,&end)==2 && sep=='%') ||
               std::sscanf(argument,"%d,%lf,%lf%c",&nb_levels,&vmin,&vmax,&end)==3 ||
               (std::sscanf(argument,"%d%c,%lf,%lf%c",&nb_levels,&sep,&vmin,&vmax,&end)==4 && sep=='%') ||
               (std::sscanf(argument,"%d,%lf%c,%lf%c",&nb_levels,&vmin,&sepm,&vmax,&end)==4 && sepm=='%') ||
               (std::sscanf(argument,"%d%c,%lf%c,%lf%c",&nb_levels,&sep,&vmin,&sepm,&vmax,&end)==5 && sep=='%' && sepm=='%') ||
               (std::sscanf(argument,"%d,%lf,%lf%c%c",&nb_levels,&vmin,&vmax,&sepM,&end)==4 && sepm=='%') ||
               (std::sscanf(argument,"%d%c,%lf,%lf%c%c",&nb_levels,&sep,&vmin,&vmax,&sepM,&end)==5 && sep=='%' && sepm=='%') ||
               (std::sscanf(argument,"%d,%lf%c,%lf%c%c",&nb_levels,&vmin,&sepm,&vmax,&sepM,&end)==5 && sepm=='%' && sepM=='%') ||
               (std::sscanf(argument,"%d%c,%lf%c,%lf%c%c",&nb_levels,&sep,&vmin,&sepm,&vmax,&sepM,&end)==6 && sep=='%' && sepm=='%' &&
                sepM=='%')) &&
              nb_levels>0) {
            if (vmin==vmax && vmin==0) { vmax = 100; sepM = '%'; }
            print(images,"Compute histogram of image%s, using %d%s levels in range [%g%s,%g%s].",
                  gmic_selection,nb_levels,sep=='%'?"%":"",vmin,sepm=='%'?"%":"",vmax,sepM=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              double m = vmin, M = vmax;
              int nnb_levels = nb_levels;
              if (sepm=='%') m*=img.min()/100;
              if (sepM=='%') M*=img.max()/100;
              if (sep=='%') nnb_levels = (int)cimg::round(nb_levels*(1+M-m)/100,1);
              gmic_apply(images[selection[l]],histogram(nnb_levels,(T)m,(T)M));
            }
          } else arg_error("histogram");
          ++position; continue;
        }

        // Distance function.
        if (!std::strcmp("-distance",command_name)) {
          double value = 0; char sep = 0;
          if (std::sscanf(argument,"%lf%c",&value,&end)==1 ||
              (std::sscanf(argument,"%lf%c%c",&value,&sep,&end)==2 && sep=='%')) {
            print(images,"Compute distance map of image%s to isovalue %g%s.",gmic_selection,value,sep=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              double isovalue = value;
              if (sep=='%') { double m, M = img.max_min(m); isovalue = m + value*(M - m)/100; }
              gmic_apply(img,distance((T)isovalue));
            }
          } else arg_error("distance");
          ++position; continue;
        }

        // Apply Eikonal PDE to compute distance to 0.
        if (!std::strcmp("-eikonal",command_name)) {
          int nb_iter = 0; float band_size = 0;
          if ((std::sscanf(argument,"%d%c",&nb_iter,&end)==1 ||
               std::sscanf(argument,"%d,%f%c",&nb_iter,&band_size,&end)==2) &&
              nb_iter>=0 && band_size>=0) {
            print(images,"Apply %d iterations of Eikonal PDE on image%s.",nb_iter,gmic_selection);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],distance_eikonal((unsigned int)nb_iter,band_size));
          } else arg_error("eikonal");
          ++position; continue;
        }

        // Label regions.
        gmic_simple_item("-label",label_regions,"Label regions on image%s.");

        // Displacement field.
        if (!std::strcmp("-displacement",command_name)) {
          float smooth = 0.1f, precision = 0.1f; int ind0 = no_ind, nbscales = 0, itermax = 1000;
          unsigned int backward = 1; char sep = 0;
          if (((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') ||
               std::sscanf(argument,"[%d],%f%c",&ind0,&smooth,&end)==2 ||
               std::sscanf(argument,"[%d],%f,%f%c",&ind0,&smooth,&precision,&end)==3 ||
               std::sscanf(argument,"[%d],%f,%f,%d%c",&ind0,&smooth,&precision,&nbscales,&end)==4 ||
               std::sscanf(argument,"[%d],%f,%f,%d,%d%c",&ind0,&smooth,&precision,&nbscales,&itermax,&end)==5 ||
               std::sscanf(argument,"[%d],%f,%f,%d,%d,%u%c",&ind0,&smooth,&precision,&nbscales,&itermax,&backward,&end)==6) &&
              smooth>=0 && precision>0 && nbscales>=0 && itermax>=0) {
            gmic_check_indice(ind0);
            print(images,"Compute displacement field of image%s with target [%u] and smoothness %g.",gmic_selection,ind0,smooth);
            const CImg<T> target = images[ind0];
            cimg_forY(selection,l) gmic_apply(images[selection[l]],displacement(target,smooth,precision,nbscales,itermax,backward?true:false));
          } else arg_error("displacement");
          ++position; continue;
        }

        // Sort.
        gmic_simple_item("-sort",sort,"Sort values in image%s.");

        // PSNR.
        if (!std::strcmp("-psnr",command_name)) {
          double valmax = 255;
          if (std::sscanf(argument,"%lf%c",&valmax,&end)==1) ++position;
          CImgList<T> subimages;
          cimg_forY(selection,l) subimages.insert(images[l],~0U,true);
          print(images,"Compute %ux%u matrix of PSNR values from image%s, with maximum pixel value %g.",
                subimages.size(),subimages.size(),gmic_selection,valmax);
          CImg<T> res(subimages.size(),subimages.size(),1,1,(T)-1);
          cimg_forXY(res,x,y) if (x>y) res(x,y) = res(y,x) = (T)subimages[x].PSNR(subimages[y],(float)valmax);
          if (get_version) {
            CImg<char>("(PSNR)",7).move_to(filenames);
            res.move_to(images);
          } else {
            if (selection) {
              cimg_forY(selection,l) { const unsigned int ind = selection[l] - l; images.remove(ind); filenames.remove(ind); }
              images.insert(res,selection[0]);
              CImg<char>("(PSNR)",7).move_to(filenames,selection[0]);
            }
          }
          continue;
        }

        // Draw point.
        if (!std::strcmp("-point",command_name)) {
          char argx[4096] = { 0 }, argy[4096] = { 0 }, argz[4096] = { 0 }, color[4096] = { 0 }, sepx = 0, sepy = 0, sepz = 0;
          float x = 0, y = 0, z = 0, opacity = 1;
          if ((std::sscanf(argument,"%4095[0-9.eE%+-]%c",argx,&end)==1 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",argx,argy,&end)==2 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",argx,argy,argz,&end)==3 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f%c",argx,argy,argz,&opacity,&end)==4 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f,%4095[0-9.eE,+-]%c",
                           argx,argy,argz,&opacity,color,&end)==5) &&
              (!*argx || std::sscanf(argx,"%f%c",&x,&end)==1 || (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy || std::sscanf(argy,"%f%c",&y,&end)==1 || (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (!*argz || std::sscanf(argz,"%f%c",&z,&end)==1 || (std::sscanf(argz,"%f%c%c",&z,&sepz,&end)==2 && sepz=='%'))) {
            print(images,"Draw point (%g%s,%g%s,%g%s) with opacity %g and color '%s' on image%s.",
                  x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,sepz=='%'?"%":"",opacity,color[0]?color:"default",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]], col(img.spectrum(),1,1,1,0);
              col.fill(color,true);
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width()-1)/100:x,1),
                ny = (int)cimg::round(sepy=='%'?y*(img.height()-1)/100:y,1),
                nz = (int)cimg::round(sepz=='%'?z*(img.depth()-1)/100:z,1);
              gmic_apply(img,draw_point(nx,ny,nz,col.data(),opacity));
            }
          } else arg_error("point");
          ++position; continue;
        }

        // Draw line.
        if (!std::strcmp("-line",command_name)) {
          char argx0[4096] = { 0 }, argy0[4096] = { 0 }, argx1[4096] = { 0 }, argy1[4096] = { 0 }, color[4096] = { 0 };
          char sepx0 = 0, sepy0 = 0, sepx1 = 0, sepy1 = 0;
          float x0 = 0, y0 = 0, x1 = 0, y1 = 0, opacity = 1;
          if ((std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",
                           argx0,argy0,argx1,argy1,&end)==4 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f%c",
                           argx0,argy0,argx1,argy1,&opacity,&end)==5 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f,%4095[0-9.eE,+-]%c",
                           argx0,argy0,argx1,argy1,&opacity,color,&end)==6) &&
              (std::sscanf(argx0,"%f%c",&x0,&end)==1 || (std::sscanf(argx0,"%f%c%c",&x0,&sepx0,&end)==2 && sepx0=='%')) &&
              (std::sscanf(argy0,"%f%c",&y0,&end)==1 || (std::sscanf(argy0,"%f%c%c",&y0,&sepy0,&end)==2 && sepy0=='%')) &&
              (std::sscanf(argx1,"%f%c",&x1,&end)==1 || (std::sscanf(argx1,"%f%c%c",&x1,&sepx1,&end)==2 && sepx1=='%')) &&
              (std::sscanf(argy1,"%f%c",&y1,&end)==1 || (std::sscanf(argy1,"%f%c%c",&y1,&sepy1,&end)==2 && sepy1=='%'))) {
            print(images,"Draw line (%g%s,%g%s) - (%g%s,%g%s) with opacity %g and color '%s' on image%s.",
                  x0,sepx0=='%'?"%":"",y0,sepy0=='%'?"%":"",x1,sepx1=='%'?"%":"",y1,sepy1=='%'?"%":"",
                  opacity,color[0]?color:"default",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]], col(img.spectrum(),1,1,1,0);
              col.fill(color,true);
              const int
                nx0 = (int)cimg::round(sepx0=='%'?x0*(img.width()-1)/100:x0,1),
                ny0 = (int)cimg::round(sepy0=='%'?y0*(img.height()-1)/100:y0,1),
                nx1 = (int)cimg::round(sepx1=='%'?x1*(img.width()-1)/100:x1,1),
                ny1 = (int)cimg::round(sepy1=='%'?y1*(img.height()-1)/100:y1,1);
              gmic_apply(img,draw_line(nx0,ny0,nx1,ny1,col.data(),opacity));
            }
          } else arg_error("line");
          ++position; continue;
        }

        // Draw polygon.
        if (!std::strcmp("-polygon",command_name)) {
          char arg0[4096] = { 0 }, arg1[4096] = { 0 }, sepx0 = 0, sepy0 = 0;
          int N = 0; float x0 = 0, y0 = 0, opacity = 1;
          if (std::sscanf(argument,"%d%c",&N,&end)==2 && N>2) {
            const char
              *nargument = argument + std::sprintf(tmpstr,"%d",N) + 1,
              *const eargument = argument + std::strlen(argument);
            CImg<float> coords0(N,2,1,1,0);
            CImg<bool> percents(N,2,1,1,0);
            for (int n = 0; n<N; ++n) if (nargument<eargument) {
              if (std::sscanf(nargument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-]",arg0,arg1)==2 &&
                  (std::sscanf(arg0,"%f%c",&x0,&end)==1 || (std::sscanf(arg0,"%f%c%c",&x0,&(sepx0=0),&end)==2 && sepx0=='%')) &&
                  (std::sscanf(arg1,"%f%c",&y0,&end)==1 || (std::sscanf(arg1,"%f%c%c",&y0,&(sepy0=0),&end)==2 && sepy0=='%'))) {
                coords0(n,0) = x0; percents(n,0) = (sepx0=='%');
                coords0(n,1) = y0; percents(n,1) = (sepy0=='%');
                nargument+=std::strlen(arg0) + std::strlen(arg1) + 2;
              } else arg_error("polygon");
            } else arg_error("polygon");
            if (nargument<eargument && std::sscanf(nargument,"%4095[0-9.eE+-]",arg0)==1 &&
                std::sscanf(arg0,"%f",&opacity)==1) nargument+=std::strlen(arg0)+1;
            const char *const color = nargument<eargument?nargument:&(end=0);
            print(images,"Draw %d-vertices polygon with opacity %g and color '%s' on image%s.",
                  N,opacity,color[0]?color:"default",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              CImg<int> coords(coords0);
              cimg_forX(coords,p) {
                if (percents(p,0)) coords(p,0) = (int)cimg::round(coords0(p,0)*(img.width()-1)/100,1);
                if (percents(p,1)) coords(p,1) = (int)cimg::round(coords0(p,1)*(img.height()-1)/100,1);
              }
              CImg<T> col(img.spectrum(),1,1,1,0);
              col.fill(color,true);
              gmic_apply(img,draw_polygon(coords,col.data(),opacity));
            }
          } else arg_error("polygon");
          ++position; continue;
        }

        // Draw spline.
        if (!std::strcmp("-spline",command_name)) {
          char argx0[256] = { 0 }, argy0[256] = { 0 }, argu0[256] = { 0 }, argv0[256] = { 0 },
               argx1[256] = { 0 }, argy1[256] = { 0 }, argu1[256] = { 0 }, argv1[256] = { 0 },
               color[4096] = { 0 };
          char sepx0 = 0, sepy0 = 0, sepu0 = 0, sepv0 = 0, sepx1 = 0, sepy1 = 0, sepu1 = 0, sepv1 = 0;
          float x0 = 0, y0 = 0, u0 = 0, v0 = 0, x1 = 0, y1 = 0, u1 = 0, v1 = 0, opacity = 1;
          if ((std::sscanf(argument,
                           "%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx0,argy0,argu0,argv0,argx1,argy1,argu1,argv1,&end)==8 ||
               std::sscanf(argument,
                           "%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%f%c",
                           argx0,argy0,argu0,argv0,argx1,argy1,argu1,argv1,&opacity,&end)==9 ||
               std::sscanf(argument,
                           "%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%f,%4095[0-9.eE,+-]%c",
                           argx0,argy0,argu0,argv0,argx1,argy1,argu1,argv1,&opacity,color,&end)==10) &&
              (std::sscanf(argx0,"%f%c",&x0,&end)==1 || (std::sscanf(argx0,"%f%c%c",&x0,&sepx0,&end)==2 && sepx0=='%')) &&
              (std::sscanf(argy0,"%f%c",&y0,&end)==1 || (std::sscanf(argy0,"%f%c%c",&y0,&sepy0,&end)==2 && sepy0=='%')) &&
              (std::sscanf(argu0,"%f%c",&u0,&end)==1 || (std::sscanf(argu0,"%f%c%c",&u0,&sepu0,&end)==2 && sepu0=='%')) &&
              (std::sscanf(argv0,"%f%c",&v0,&end)==1 || (std::sscanf(argv0,"%f%c%c",&v0,&sepv0,&end)==2 && sepv0=='%')) &&
              (std::sscanf(argx1,"%f%c",&x1,&end)==1 || (std::sscanf(argx1,"%f%c%c",&x1,&sepx1,&end)==2 && sepx1=='%')) &&
              (std::sscanf(argy1,"%f%c",&y1,&end)==1 || (std::sscanf(argy1,"%f%c%c",&y1,&sepy1,&end)==2 && sepy1=='%')) &&
              (std::sscanf(argu1,"%f%c",&u1,&end)==1 || (std::sscanf(argu1,"%f%c%c",&u1,&sepu1,&end)==2 && sepu1=='%')) &&
              (std::sscanf(argv1,"%f%c",&v1,&end)==1 || (std::sscanf(argv1,"%f%c%c",&v1,&sepv1,&end)==2 && sepv1=='%'))) {
            print(images,"Draw spline from (%g%s,%g%s) [%g%s,%g%s] to (%g%s,%g%s) [%g%s,%g%s], with opacity %g and "
                  "color '%s' on image%s.",
                  x0,sepx0=='%'?"%":"",y0,sepy0=='%'?"%":"",u0,sepu0=='%'?"%":"",v0,sepv0=='%'?"%":"",
                  x1,sepx1=='%'?"%":"",y1,sepy1=='%'?"%":"",u1,sepu1=='%'?"%":"",v1,sepv1=='%'?"%":"",
                  opacity,*color?color:"default",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]], col(img.spectrum(),1,1,1,0);
              col.fill(color,true);
              const int
                nx0 = (int)cimg::round(sepx0=='%'?x0*(img.width()-1)/100:x0,1),
                ny0 = (int)cimg::round(sepy0=='%'?y0*(img.height()-1)/100:y0,1),
                nx1 = (int)cimg::round(sepx1=='%'?x1*(img.width()-1)/100:x1,1),
                ny1 = (int)cimg::round(sepy1=='%'?y1*(img.height()-1)/100:y1,1);
              const float
                nu0 = sepu0=='%'?u0*(img.width()-1)/100:u0,
                nv0 = sepv0=='%'?v0*(img.height()-1)/100:v0,
                nu1 = sepu1=='%'?u1*(img.width()-1)/100:u1,
                nv1 = sepv1=='%'?v1*(img.height()-1)/100:v1;
              gmic_apply(img,draw_spline(nx0,ny0,nu0,nv0,nx1,ny1,nu1,nv1,col.data(),opacity,4));
            }
          } else arg_error("spline");
          ++position; continue;
        }

        // Draw ellipse.
        if (!std::strcmp("-ellipse",command_name)) {
          char argx[256] = { 0 }, argy[256] = { 0 }, argR[256] = { 0 }, argr[256] = { 0 }, color[4096] = { 0 };
          char sepx = 0, sepy = 0, sepR = 0, sepr = 0;
          float x = 0, y = 0, R = 0, r = 0, angle = 0, opacity = 1;
          if ((std::sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx,argy,argR,&end)==3 ||
               std::sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx,argy,argR,argr,&end)==4 ||
               std::sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%f%c",
                           argx,argy,argR,argr,&angle,&end)==5 ||
               std::sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%f,%f%c",
                           argx,argy,argR,argr,&angle,&opacity,&end)==6 ||
               std::sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%f,%f,%4095[0-9.eE,+-]%c",
                           argx,argy,argR,argr,&angle,&opacity,color,&end)==7) &&
              (std::sscanf(argx,"%f%c",&x,&end)==1 || (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (std::sscanf(argy,"%f%c",&y,&end)==1 || (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (std::sscanf(argR,"%f%c",&R,&end)==1 || (std::sscanf(argR,"%f%c%c",&R,&sepR,&end)==2 && sepR=='%')) &&
              (!*argr || std::sscanf(argr,"%f%c",&r,&end)==1 || (std::sscanf(argr,"%f%c%c",&r,&sepr,&end)==2 && sepr=='%'))) {
            if (!*argr) r = R;
            print(images,"Draw ellipse centered at (%g%s,%g%s) with radii (%g%s,%g%s), orientation %g deg, opacity %g and "
                  "color '%s' on image%s.",
                  x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",R,sepR=='%'?"%":"",r,sepr=='%'?"%":"",angle,
                  opacity,color[0]?color:"default",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]], col(img.spectrum(),1,1,1,0);
              col.fill(color,true);
              const float rmax = std::sqrt((float)cimg::sqr(img.width()) + cimg::sqr(img.height()))/2;
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width()-1)/100:x,1),
                ny = (int)cimg::round(sepy=='%'?y*(img.height()-1)/100:y,1);
              const float
                nR = (float)cimg::round(sepR=='%'?R*rmax/100:R,1),
                nr = (float)cimg::round(sepr=='%'?r*rmax/100:r,1);
                gmic_apply(img,draw_ellipse(nx,ny,nR,nr,angle,col.data(),opacity));
            }
          } else arg_error("ellipse");
          ++position; continue;
        }

        // Draw text.
        if (!std::strcmp("-text",command_name)) {
          char argx[4096] = { 0 }, argy[4096] = { 0 }, color[4096] = { 0 }, text[4096] = { 0 }, sepx = 0, sepy = 0;
          float x = 0, y = 0, opacity = 1; int siz = 11;
          if ((std::sscanf(argument,"%4095[^,]%c",text,&end)==1 ||
               std::sscanf(argument,"%4095[^,],%4095[0-9.eE%+-]%c",text,argx,&end)==2 ||
               std::sscanf(argument,"%4095[^,],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",text,argx,argy,&end)==3 ||
               std::sscanf(argument,"%4095[^,],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%d%c",text,argx,argy,&siz,&end)==4 ||
               std::sscanf(argument,"%4095[^,],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%d,%f%c",text,argx,argy,&siz,&opacity,&end)==5 ||
               std::sscanf(argument,"%4095[^,],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%d,%f,%4095[0-9.eE,+-]%c",
                           text,argx,argy,&siz,&opacity,color,&end)==6) &&
              (!*argx || std::sscanf(argx,"%f%c",&x,&end)==1 || (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy || std::sscanf(argy,"%f%c",&y,&end)==1 || (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              siz>0) {
            for (char *s = text; *s; ++s) if (*s==29) *s=',';
            print(images,"Draw text '%s' at position (%g%s,%g%s) with font size %d, opacity %g and color '%s' on image%s.",
                  text,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",siz,opacity,color[0]?color:"default",gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]], col(img.spectrum(),1,1,1,0);
              col.fill(color,true);
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width()-1)/100:x,1),
                ny = (int)cimg::round(sepy=='%'?y*(img.height()-1)/100:y,1);
              gmic_apply(img,draw_text(nx,ny,text,col.data(),0,opacity,siz));
            }
          } else arg_error("text");
          ++position; continue;
        }

        // Draw image.
        if (!std::strcmp("-image",command_name)) {
          char argx[4096] = { 0 }, argy[4096] = { 0 }, argz[4096] = { 0 }, sep = 0, sepx = 0, sepy = 0, sepz = 0;
          int ind = no_ind, indm = no_ind; float x = 0, y = 0, z = 0, opacity = 1;
          if (((std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==1 && sep==']') ||
               std::sscanf(argument,"[%d],%4095[0-9.eE%+-]%c",&ind,argx,&end)==2 ||
               std::sscanf(argument,"[%d],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",&ind,argx,argy,&end)==3 ||
               std::sscanf(argument,"[%d],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",&ind,argx,argy,argz,&end)==4 ||
               std::sscanf(argument,"[%d],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f%c",&ind,argx,argy,argz,&opacity,&end)==5 ||
               (std::sscanf(argument,"[%d],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f,[%d%c%c",
                            &ind,argx,argy,argz,&opacity,&indm,&sep,&end)==7 && sep==']')) &&
              (!*argx || std::sscanf(argx,"%f%c",&x,&end)==1 || (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy || std::sscanf(argy,"%f%c",&y,&end)==1 || (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (!*argz || std::sscanf(argz,"%f%c",&z,&end)==1 || (std::sscanf(argz,"%f%c%c",&z,&sepz,&end)==2 && sepz=='%'))) {
            gmic_check_indice(ind);
            const CImg<T> sprite = images[ind];
            CImg<T> mask;
            if (indm!=no_ind) {
              gmic_check_indice(indm);
              mask = images[indm];
              print(images,"Draw image [%d] at (%g%s,%g%s,%g%s), with opacity %g and mask [%d] on image%s.",
                    ind,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,sepz=='%'?"%":"",opacity,indm,gmic_selection);
            } else print(images,"Draw image [%d] at (%g%s,%g%s,%g%s) with opacity %g on image%s.",
                         ind,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,sepz=='%'?"%":"",opacity,gmic_selection);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width()-1)/100:x,1),
                ny = (int)cimg::round(sepy=='%'?y*(img.height()-1)/100:y,1),
                nz = (int)cimg::round(sepz=='%'?z*(img.depth()-1)/100:z,1);
              if (indm!=no_ind) { gmic_apply(img,draw_image(nx,ny,nz,sprite,mask,opacity)); }
              else { gmic_apply(img,draw_image(nx,ny,nz,sprite,opacity)); }
            }
          } else arg_error("image");
          ++position; continue;
        }

        // Draw 3D object.
        if (!std::strcmp("-object3d",command_name)) {
          char argx[4096] = { 0 }, argy[4096] = { 0 }, sep = 0, sepx = 0, sepy = 0;
          float x = 0, y = 0, z = 0, opacity = 1; int ind = no_ind, is_zbuffer = 1;
          if (((std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') ||
               std::sscanf(argument,"[%d],%4095[0-9.eE%+-]%c",&ind,argx,&end)==2 ||
               std::sscanf(argument,"[%d],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",&ind,argx,argy,&end)==3 ||
               std::sscanf(argument,"[%d],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f%c",&ind,argx,argy,&z,&end)==4 ||
               std::sscanf(argument,"[%d],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f,%f%c",&ind,argx,argy,&z,&opacity,&end)==5 ||
               std::sscanf(argument,"[%d],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f,%f,%d%c",&ind,argx,argy,&z,&opacity,&is_zbuffer,&end)==6) &&
              (!*argx || std::sscanf(argx,"%f%c",&x,&end)==1 || (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy || std::sscanf(argy,"%f%c",&y,&end)==1 || (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%'))) {
            gmic_check_indice(ind);
            if (!images[ind].is_CImg3d())
              error(images,"Command 'object3d' : Invalid 3D object [%d] specified in argument '%s'.",ind,argument_text);
            print(images,"Draw 3D object [%d] at (%g%s,%g%s,%g) with opacity %g on image%s.",
                  ind,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,opacity,gmic_selection);
            CImgList<unsigned int> primitives;
            CImgList<unsigned char> colors;
            CImg<float> opacities, vertices(images[ind]);
            vertices.CImg3dtoobject3d(primitives,colors,opacities);
            opacities*=opacity;
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const float
                nx = sepx=='%'?x*(img.width()-1)/100:x,
                ny = sepy=='%'?y*(img.height()-1)/100:y;
              CImg<float> zbuffer(is_zbuffer?img.width():0,is_zbuffer?img.height():0,1,1,0);
              gmic_apply(img,draw_object3d(nx,ny,z,vertices,primitives,colors,opacities,
                                           render3d,is_double3d,focale3d,light3d_x,light3d_y,light3d_z,specular_light3d,specular_shine3d,
                                           zbuffer));
            }
          } else arg_error("object3d");
          ++position; continue;
        }

        // Draw plasma fractal.
        if (!std::strcmp("-plasma",command_name)) {
          float alpha = 1, beta = 1, opacity = 1;
          if (std::sscanf(argument,"%f%c",&alpha,&end)==1 ||
              std::sscanf(argument,"%f,%f%c",&alpha,&beta,&end)==2 ||
              std::sscanf(argument,"%f,%f,%f%c",&alpha,&beta,&opacity,&end)==3) {
            print(images,"Draw plasma in image%s with alpha %g, beta %g and opacity %g.",gmic_selection,alpha,beta,opacity);
            cimg_forY(selection,l) gmic_apply(images[selection[l]],draw_plasma(alpha,beta,opacity));
          } else arg_error("plasma");
          ++position; continue;
        }

        // Draw Mandelbrot/Julia fractal.
        if (!std::strcmp("-mandelbrot",command_name)) {
          double z0r = -2, z0i = -2, z1r = 2, z1i = 2, paramr = 0, parami = 0;
          float opacity = 1; int itermax = 100; unsigned int julia = 0;
          if ((std::sscanf(argument,"%lf,%lf,%lf,%lf%c",&z0r,&z0i,&z1r,&z1i,&end)==4 ||
               std::sscanf(argument,"%lf,%lf,%lf,%lf,%d%c",&z0r,&z0i,&z1r,&z1i,&itermax,&end)==5 ||
               std::sscanf(argument,"%lf,%lf,%lf,%lf,%d,%u,%lf,%lf%c",&z0r,&z0i,&z1r,&z1i,&itermax,&julia,&paramr,&parami,&end)==8 ||
               std::sscanf(argument,"%lf,%lf,%lf,%lf,%d,%u,%lf,%lf,%f%c",
                           &z0r,&z0i,&z1r,&z1i,&itermax,&julia,&paramr,&parami,&opacity,&end)==9) &&
              itermax>=0) {
            print(images,"Draw %s fractal in image%s from complex area (%g,%g)-(%g,%g) with c0 = (%g,%g) (%d iterations).",
                  julia?"Julia":"Mandelbrot",gmic_selection,z0r,z0i,z1r,z1i,paramr,parami,itermax);
            cimg_forY(selection,l)
              gmic_apply(images[selection[l]],draw_mandelbrot(CImg<T>(),opacity,z0r,z0i,z1r,z1i,itermax,true,
                                                            julia?true:false,paramr,parami));
          } else arg_error("mandelbrot");
          ++position; continue;
        }

        // Draw quiver.
        if (!std::strcmp("-quiver",command_name)) {
          int ind = no_ind, sampling = 25; unsigned int arrows = 1; float factor = -20, opacity = 1; char color[4096] = { 0 };
          if ((std::sscanf(argument,"[%d]%c",&ind,&end)==1 ||
               std::sscanf(argument,"[%d],%d%c",&ind,&sampling,&end)==2 ||
               std::sscanf(argument,"[%d],%d,%f%c",&ind,&sampling,&factor,&end)==3 ||
               std::sscanf(argument,"[%d],%d,%f,%u%c",&ind,&sampling,&factor,&arrows,&end)==4 ||
               std::sscanf(argument,"[%d],%d,%f,%u,%f%c",&ind,&sampling,&factor,&arrows,&opacity,&end)==5 ||
               std::sscanf(argument,"[%d],%d,%f,%u,%f,%4095[0-9.eE,+-]%c",&ind,&sampling,&factor,&arrows,&opacity,color,&end)==6) &&
              sampling>0) {
            gmic_check_indice(ind);
            print(images,"Draw 2D vector field on image%s with sampling %d, factor %g, opacity %g and color '%s'.",
                  gmic_selection,sampling,factor,opacity,color);
            const CImg<T> flow = images[ind];
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]], col(img.spectrum(),1,1,1,0);
              col.fill(color,true);
              gmic_apply(img,draw_quiver(flow,col.data(),opacity,sampling,factor,arrows?true:false));
            }
          } else arg_error("quiver");
          ++position; continue;
        }

        // Flood fill.
        if (!std::strcmp("-flood",command_name)) {
          char argx[4096] = { 0 }, argy[4096] = { 0 }, argz[4096] = { 0 }, color[4096] = { 0 }, sepx = 0, sepy = 0, sepz = 0;
          float x = 0, y = 0, z = 0, tolerance = 0, opacity = 1;
          if ((std::sscanf(argument,"%4095[0-9.eE%+-]%c",argx,&end)==1 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",argx,argy,&end)==2 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-]%c",argx,argy,argz,&end)==3 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f%c",argx,argy,argz,&tolerance,&end)==4 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f,%f%c",
                           argx,argy,argz,&tolerance,&opacity,&end)==5 ||
               std::sscanf(argument,"%4095[0-9.eE%+-],%4095[0-9.eE%+-],%4095[0-9.eE%+-],%f,%f,%4095[0-9.eE,+-]%c",
                           argx,argy,argz,&tolerance,&opacity,color,&end)==6) &&
              (!*argx || std::sscanf(argx,"%f%c",&x,&end)==1 || (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy || std::sscanf(argy,"%f%c",&y,&end)==1 || (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (!*argz || std::sscanf(argz,"%f%c",&z,&end)==1 || (std::sscanf(argz,"%f%c%c",&z,&sepz,&end)==2 && sepz=='%')) &&
              tolerance>=0) {
            print(images,"Flood fill image%s from (%g%s,%g%s,%g%s) with tolerance %g, opacity %g and color '%s'.",
                  gmic_selection,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,sepz=='%'?"%":"",tolerance,opacity,color);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]], col(img.spectrum(),1,1,1,0);
              col.fill(color,true);
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width()-1)/100:x,1),
                ny = (int)cimg::round(sepy=='%'?y*(img.height()-1)/100:y,1),
                nz = (int)cimg::round(sepz=='%'?z*(img.depth()-1)/100:z,1);
              gmic_apply(img,draw_fill(nx,ny,nz,col.data(),opacity,tolerance));
            }
          } else arg_error("flood");
          ++position; continue;
        }

        //-------------------------
        // Image list manipulation
        //-------------------------

        // Remove selected images.
        if (!std::strcmp("-remove",command_name) || !std::strcmp("-rm",command_name)) {
          print(images,"Remove image%s",gmic_selection);
          unsigned int off = 0;
          cimg_forY(selection,l) {
            const unsigned int ind = selection[l] - off;
            images.remove(ind); filenames.remove(ind);
            ++off;
          }
          if (verbosity_level>=0) {
            std::fprintf(cimg::output()," (%u image%s left).",images.size(),images.size()==1?"":"s");
            std::fflush(cimg::output());
          }
          continue;
        }

        // Keep selected images.
        if (!std::strcmp("-keep",command_name) || !std::strcmp("-k",command_name)) {
          print(images,"Keep image%s",gmic_selection);
          CImgList<T> nimages(selection.size());
          CImgList<char> nfilenames(selection.size());
          cimg_forY(selection,l) {
            nimages[l].swap(images[selection[l]]);
            nfilenames[l].swap(filenames[selection[l]]);
          }
          nimages.move_to(images);
          nfilenames.move_to(filenames);
          if (verbosity_level>=0) {
            std::fprintf(cimg::output()," (%u image%s left).",images.size(),images.size()==1?"":"s");
            std::fflush(cimg::output());
          }
          continue;
        }

        // Move images to specified position.
        if (!std::strcmp("-move",command_name) || !std::strcmp("-mv",command_name)) {
          float pos = 0; char sep = 0;
          if (std::sscanf(argument,"%f%c",&pos,&end)==1 ||
              (std::sscanf(argument,"%f%c%c",&pos,&sep,&end)==2 && sep=='%')) {
            int ind0 = (int)(sep=='%'?pos*images.size()/100:pos);
            if (ind0<0) ind0+=images.size();
            if (ind0<0) ind0 = 0;
            if (ind0>(int)images.size()) ind0 = images.size();
            print(images,"Move image%s to position %d.",gmic_selection,ind0);
            CImgList<T> nimages;
            CImgList<char> nfilenames;
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              images[ind].move_to(nimages);
              filenames[ind].move_to(nfilenames);
            }
            nimages.move_to(images,ind0);
            nfilenames.move_to(filenames,ind0);
            cimglist_for(images,l) if (!images[l]) { images.remove(l); filenames.remove(l--); }
          } else arg_error("move");
          ++position; continue;
        }

        // Reverse images order.
        if (!std::strcmp("-reverse",command_name)) {
          print(images,"Reverse positions of image%s.",gmic_selection);
          CImgList<T> nimages(selection.size());
          CImgList<char> nfilenames(selection.size());
          cimg_forY(selection,l) { nimages[l].swap(images[selection[l]]); nfilenames[l].swap(filenames[selection[l]]); }
          nimages.reverse(); nfilenames.reverse();
          cimg_forY(selection,l) { nimages[l].swap(images[selection[l]]); nfilenames[l].swap(filenames[selection[l]]); }
          continue;
        }

        // Set image name.
        if (!std::strcmp("-name",command_name)) {
          int ind = no_ind; char sep = 0;
          if (std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') {
            gmic_check_indice(ind);
            print(images,"Set name of image%s to '%s'.",gmic_selection,filenames[ind].data());
            cimg_forY(selection,l) filenames[selection[l]].assign(filenames[ind]);
          } else {
            print(images,"Set name of image%s to '%s'.",gmic_selection,argument_text);
            cimg_forY(selection,l) filenames[selection[l]].assign(argument,std::strlen(argument)+1);
          }
          ++position; continue;
        }

        //-------------------------
        // 3D objects manipulation
        //-------------------------

        // Generate 3D line.
        if (!std::strcmp("-line3d",item)) {
          float x0 = 0, y0 = 0, z0 = 0, x1 = 0, y1 = 0, z1 = 0;
          if (std::sscanf(argument,"%f,%f,%f,%f,%f,%f%c",&x0,&y0,&z0,&x1,&y1,&z1,&end)==6) {
            print(images,"Generates 3D line (%g,%g,%g)-(%g,%g,%g).",x0,y0,z0,x1,y1,z1);
            CImg<T>(1,21,1,1,'C'+0.5,'I'+0.5,'m'+0.5,'g'+0.5,'3'+0.5,'d'+0.5,
                    2.,1.,x0,y0,z0,x1,y1,z1,2.,0.,1.,
                    200.,200.,200.,1.).move_to(images);
            CImg<char>("(3D line)",10).move_to(filenames);
          } else arg_error("line3d");
          ++position; continue;
        }

        // Generate 3D triangle.
        if (!std::strcmp("-triangle3d",item)) {
          float x0 = 0, y0 = 0, z0 = 0, x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0;
          if (std::sscanf(argument,"%f,%f,%f,%f,%f,%f,%f,%f,%f%c",&x0,&y0,&z0,&x1,&y1,&z1,&x2,&y2,&z2,&end)==9) {
            print(images,"Generates 3D triangle (%g,%g,%g)-(%g,%g,%g)-(%g,%g,%g).",x0,y0,z0,x1,y1,z1,x2,y2,z2);
            CImg<T>(1,25,1,1,'C'+0.5,'I'+0.5,'m'+0.5,'g'+0.5,'3'+0.5,'d'+0.5,
                    3.,1.,x0,y0,z0,x1,y1,z1,x2,y2,z2,3.,0.,1.,2.,
                    200.,200.,200.,1.).move_to(images);
            CImg<char>("(3D triangle)",14).move_to(filenames);
          } else arg_error("triangle3d");
          ++position; continue;
        }

        // Generate 3D quadrangle.
        if (!std::strcmp("-quadrangle3d",item)) {
          float x0 = 0, y0 = 0, z0 = 0, x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0, x3 = 0, y3 = 0, z3 = 0;
          if (std::sscanf(argument,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f%c",&x0,&y0,&z0,&x1,&y1,&z1,&x2,&y2,&z2,&x3,&y3,&z3,&end)==12) {
            print(images,"Generates 3D quadrangle (%g,%g,%g)-(%g,%g,%g)-(%g,%g,%g)-(%g,%g,%g).",x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);
            CImg<T>(1,29,1,1,'C'+0.5,'I'+0.5,'m'+0.5,'g'+0.5,'3'+0.5,'d'+0.5,
                    4.,1.,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,4.,0.,1.,2.,3.,
                    200.,200.,200.,1.).move_to(images);
            CImg<char>("(3D quadrangle)",16).move_to(filenames);
          } else arg_error("quadrangle3d");
          ++position; continue;
        }

        // Generate 3D box.
        if (!std::strcmp("-box3d",item)) {
          float sx = 100, sy = 100, sz = 100;
          if ((std::sscanf(argument,"%f%c",&sx,&end)==1 && ((sz=sy=sx),1)) ||
              std::sscanf(argument,"%f,%f,%f%c",&sx,&sy,&sz,&end)==3) {
            print(images,"Generate 3D box with size (%g,%g,%g).",sx,sy,sz);
            CImgList<unsigned int> primitives;
            CImg<float> vertices = CImg<T>::box3d(primitives,sx,sy,sz);
            CImgList<unsigned char> colors(primitives.size(),1,3,1,1,200);
            CImg<float> opacities(1,primitives.size(),1,1,1);
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            CImg<char>("(3D cube)",10).move_to(filenames);
          } else arg_error("box3d");
          ++position; continue;
        }

        // Generate 3D cone.
        if (!std::strcmp("-cone3d",item)) {
          float radius = 100, height = 200; int subdivisions = 24;
          if ((std::sscanf(argument,"%f%c",&radius,&end)==1 ||
               std::sscanf(argument,"%f,%f%c",&radius,&height,&end)==2 ||
               std::sscanf(argument,"%f,%f,%d%c",&radius,&height,&subdivisions,&end)==3) &&
              subdivisions>0) {
            print(images,"Generate 3D cone with radius %g, height %g and %d subdivisions.",radius,height,subdivisions);
            CImgList<unsigned int> primitives;
            CImg<float> vertices = CImg<T>::cone3d(primitives,radius,height,subdivisions);
            CImgList<unsigned char> colors(primitives.size(),1,3,1,1,200);
            CImg<float> opacities(1,primitives.size(),1,1,1);
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            CImg<char>("(3D cone)",10).move_to(filenames);
          } else arg_error("cone3d");
          ++position; continue;
        }

        // Generate 3D cylinder.
        if (!std::strcmp("-cylinder3d",item)) {
          float radius = 100, height = 200; int subdivisions = 24;
          if ((std::sscanf(argument,"%f%c",&radius,&end)==1 ||
               std::sscanf(argument,"%f,%f%c",&radius,&height,&end)==2 ||
               std::sscanf(argument,"%f,%f,%d%c",&radius,&height,&subdivisions,&end)==3) &&
              subdivisions>0) {
            print(images,"Generate 3D cylinder with radius %g, height %g and %d subdivisions.",radius,height,subdivisions);
            CImgList<unsigned int> primitives;
            CImg<float> vertices = CImg<T>::cylinder3d(primitives,radius,height,subdivisions);
            CImgList<unsigned char> colors(primitives.size(),1,3,1,1,200);
            CImg<float> opacities(1,primitives.size(),1,1,1);
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            CImg<char>("(3D cylinder)",14).move_to(filenames);
          } else arg_error("cylinder3d");
          ++position; continue;
        }

        // Generate 3D torus.
        if (!std::strcmp("-torus3d",item)) {
          float radius1 = 100, radius2 = 30; int subdivisions1 = 24, subdivisions2 = 12;
          if ((std::sscanf(argument,"%f%c",&radius1,&end)==1 ||
               std::sscanf(argument,"%f,%f%c",&radius1,&radius2,&end)==2 ||
               std::sscanf(argument,"%f,%f,%d,%c",&radius1,&radius2,&subdivisions1,&end)==3 ||
               std::sscanf(argument,"%f,%f,%d,%d%c",&radius1,&radius2,&subdivisions1,&subdivisions2,&end)==4) &&
              subdivisions1>0 && subdivisions2>0) {
            print(images,"Generate 3D torus with radii %g and %g, and subdivisions %d and %d.",radius1,radius2,subdivisions1,subdivisions2);
            CImgList<unsigned int> primitives;
            CImg<float> vertices = CImg<T>::torus3d(primitives,radius1,radius2,subdivisions1,subdivisions2);
            CImgList<unsigned char> colors(primitives.size(),1,3,1,1,200);
            CImg<float> opacities(1,primitives.size(),1,1,1);
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            CImg<char>("(3D torus)",11).move_to(filenames);
          } else arg_error("torus3d");
          ++position; continue;
        }

        // Generate 3D plane.
        if (!std::strcmp("-plane3d",item)) {
          float size_x = 100, size_y = 30; int subdivisions_x = 24, subdivisions_y = 12;
          if ((std::sscanf(argument,"%f%c",&size_x,&end)==1 ||
               std::sscanf(argument,"%f,%f%c",&size_x,&size_y,&end)==2 ||
               std::sscanf(argument,"%f,%f,%d%c",&size_x,&size_y,&subdivisions_x,&end)==3 ||
               std::sscanf(argument,"%f,%f,%d,%d%c",&size_x,&size_y,&subdivisions_x,&subdivisions_y,&end)==4) &&
              subdivisions_x>0 && subdivisions_y>0) {
            print(images,"Generate 3D plane with dimensions (%g,%g), and subdivisions (%d,%d).",size_x,size_y,subdivisions_x,subdivisions_y);
            CImgList<unsigned int> primitives;
            CImg<float> vertices = CImg<T>::plane3d(primitives,size_x,size_y,subdivisions_x,subdivisions_y);
            CImgList<unsigned char> colors(primitives.size(),1,3,1,1,200);
            CImg<float> opacities(1,primitives.size(),1,1,1);
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            CImg<char>("(3D plane)",11).move_to(filenames);
          } else arg_error("plane3d");
          ++position; continue;
        }

        // Generate 3D sphere.
        if (!std::strcmp("-sphere3d",item)) {
          float radius = 100; int recursions = 3;
          if ((std::sscanf(argument,"%f%c",&radius,&end)==1 ||
               std::sscanf(argument,"%f,%d%c",&radius,&recursions,&end)==2) &&
              recursions>=0) {
            print(images,"Generate 3D sphere with radius %g and %d recursions.",radius,recursions);
            CImgList<unsigned int> primitives;
            CImg<float> vertices = CImg<T>::sphere3d(primitives,radius,recursions);
            CImgList<unsigned char> colors(primitives.size(),1,3,1,1,200);
            CImg<float> opacities(1,primitives.size(),1,1,1);
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            CImg<char>("(3D sphere)",12).move_to(filenames);
          } else arg_error("sphere3d");
          ++position; continue;
        }

        // Build 3D elevation.
        if (!std::strcmp("-elevation3d",command_name)) {
          float x0 = -3, y0 = -3, x1 = 3, y1 = 3; int dx = 256, dy = 256;
          char sep = 0, sepx = 0, sepy = 0, formula[4096] = { 0 };
          if ((std::sscanf(argument,"'%4095[^']'%c",formula,&end)==1 ||
               std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f%c",formula,&x0,&y0,&x1,&y1,&end)==5 ||
               std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%d,%d%c",formula,&x0,&y0,&x1,&y1,&dx,&dy,&end)==7 ||
               (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%d%c,%d%c",formula,&x0,&y0,&x1,&y1,&dx,&sepx,&dy,&end)==8 &&
                sepx=='%') ||
               (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%d,%d%c%c",formula,&x0,&y0,&x1,&y1,&dx,&dy,&sepy,&end)==8 &&
                sepy=='%')||
               (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%d%c,%d%c%c",formula,&x0,&y0,&x1,&y1,&dx,&sepx,&dy,&sepy,&end)==9 &&
                sepx=='%' && sepy=='%')) &&
              dx>0 && dy>0) {
            print(images,"Build 3D elevation of mathematical expression '%s', range (%g,%g)-(%g,%g) and size %d%sx%d%s.",
                  formula,x0,y0,x1,y1,dx,sepx=='%'?"%":"",dy,sepy=='%'?"%":"");
            if (sepx=='%') dx = -dx;
            if (sepy=='%') dy = -dy;
            CImgList<unsigned int> primitives;
            CImg<T> vertices = CImg<T>::elevation3d(primitives,(const char*)formula,x0,y0,x1,y1,dx,dy);
            CImgList<unsigned char> colors(primitives.size(),1,3,1,1,200);
            CImg<float> opacities(1,primitives.size(),1,1,1);
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            std::sprintf(tmpstr,"(3D elevation of '%s')",formula); CImg<char>(tmpstr,std::strlen(tmpstr)+1).move_to(filenames);
          } else {
            int ind = 0; float fact = 1;
            CImg<typename CImg<T>::Tfloat> elev;
            if (std::sscanf(argument,"[%d%c%c",&ind,&sep,&end)==2 && sep==']') {
              gmic_check_indice(ind);
              print(images,"Build 3D elevation of image%s with elevation map [%d].",gmic_selection,ind);
              images[ind].get_norm().move_to(elev);
              cimg_forY(selection,l) {
                CImg<T>& img = images[selection[l]];
                CImgList<unsigned int> primitives;
                CImgList<unsigned char> colors;
                CImg<float>
                  vertices = img.get_elevation3d(primitives,colors,elev),
                  opacities(1,primitives.size(),1,1,1);
                vertices.object3dtoCImg3d(primitives,colors,opacities);
                gmic_apply(img,replace(vertices));
              }
            } else {
              if (std::sscanf(argument,"%f%c",&fact,&end)==1)
                print(images,"Build 3D elevation of image%s with elevation factor %g.",gmic_selection,fact);
              else print(images,"Build 3D elevation of image%s.",gmic_selection);
              cimg_forY(selection,l) {
                CImg<T>& img = images[selection[l]];
                CImgList<unsigned int> primitives;
                CImgList<unsigned char> colors;
                if (fact==1 && img.spectrum()==1) elev = img.get_shared();
                else (img.get_norm().move_to(elev))*=fact;
                CImg<float>
                  vertices = img.get_elevation3d(primitives,colors,elev),
                  opacities(1,primitives.size(),1,1,1);
                vertices.object3dtoCImg3d(primitives,colors,opacities);
                gmic_apply(img,replace(vertices));
              }
            }
          }
          ++position; continue;
        }

        // Build 3D isoline.
        if (!std::strcmp("-isoline3d",command_name)) {
          float x0 = -3, y0 = -3, x1 = 3, y1 = 3, value = 0; int dx = 256, dy = 256;
          char sep = 0, sepx = 0, sepy = 0, formula[4096] = { 0 };
          if (std::sscanf(argument,"%f%c",&value,&end)==1 ||
              std::sscanf(argument,"%f%c%c",&value,&sep,&end)==2) {
            print(images,"Build 3D isoline of image%s using isovalue %g%s.",gmic_selection,value,sep=='%'?"%":"");
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              CImg<T>& img = images[ind];
              CImg<float> vertices, opacities;
              CImgList<unsigned int> primitives;
              CImgList<unsigned char> colors;
              CImg<unsigned char> palette;
              palette.assign(3,img.spectrum(),1,1,220).noise(35,1);
              if (img.spectrum()==1) palette(0) = palette(1) = palette(2) = 255;
              else {
                palette(0,0) = 255; palette(1,0) = 30; palette(2,0) = 30;
                palette(0,1) = 30; palette(1,1) = 255; palette(2,1) = 30;
                if (img.spectrum()>=3) palette(0,2) = 30; palette(1,2) = 30; palette(2,2) = 255;
              }
              cimg_forC(img,k) {
                const CImg<T> channel = img.get_shared_channel(k);
                float nvalue = value;
                if (sep=='%') { float M = 0, m = (float)channel.min_max(M); nvalue = m + (M-m)*value/100; }
                CImgList<unsigned int> prims;
                const CImg<float> pts = img.get_shared_channel(k).get_isoline3d(prims,nvalue);
                vertices.append_object3d(primitives,pts,prims);
                colors.insert(prims.size(),CImg<unsigned char>::vector(palette(0,k),palette(1,k),palette(2,k)));
              }
              opacities.assign(1,primitives.size(),1,1,1);
              if (!vertices) warning(images,"Command 'isoline3d' : Isovalue %g not found in image [%u].",value,ind);
              vertices.object3dtoCImg3d(primitives,colors,opacities);
              gmic_apply(img,replace(vertices));
            }
          } else if ((std::sscanf(argument,"'%4095[^']',%f%c",formula,&value,&end)==2 ||
                      std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f%c",formula,&value,&x0,&y0,&x1,&y1,&end)==6 ||
                      std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%d,%d%c",formula,&value,&x0,&y0,&x1,&y1,&dx,&dy,&end)==8 ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%d%c,%d%c",
                                   formula,&value,&x0,&y0,&x1,&y1,&dx,&sepx,&dy,&end)==9 &&
                       sepx=='%') ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%d,%d%c%c",
                                   formula,&value,&x0,&y0,&x1,&y1,&dx,&dy,&sepy,&end)==9 &&
                       sepy=='%') ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%d%c,%d%c%c",
                                   formula,&value,&x0,&y0,&x1,&y1,&dx,&sepx,&dy,&sepy,&end)==10 &&
                       sepx=='%' && sepy=='%')) &&
                     dx>0 && dy>0) {
            print(images,"Build 3D isoline %g of mathematical expression '%s' in range (%g,%g)-(%g,%g) with size %d%sx%d%s.",
                  value,formula,x0,y0,x1,y1,dx,sepx=='%'?"%":"",dy,sepy=='%'?"%":"");
            if (sepx=='%') dx = -dx;
            if (sepy=='%') dy = -dy;
            CImgList<unsigned int> primitives;
            CImg<T> vertices = CImg<T>::isoline3d(primitives,(const char*)formula,value,x0,y0,x1,y1,dx,dy);
            CImgList<unsigned char> colors(primitives.size(),1,3,1,1,200);
            CImg<float> opacities(1,primitives.size(),1,1,1);
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            std::sprintf(tmpstr,"(3D isoline %g of '%s')",value,formula); CImg<char>(tmpstr,std::strlen(tmpstr)+1).move_to(filenames);
          } else arg_error("isoline3d");
          ++position; continue;
        }

        // Build 3D isosurface.
        if (!std::strcmp("-isosurface3d",command_name)) {
          float x0 = -3, y0 = -3, z0 = -3, x1 = 3, y1 = 3, z1 = 3, value = 0;
          int dx = 32, dy = 32, dz = 32; char sep = 0, sepx = 0, sepy = 0, sepz = 0, formula[4096] = { 0 };
          if (std::sscanf(argument,"%f%c",&value,&end)==1 ||
              std::sscanf(argument,"%f%c%c",&value,&sep,&end)==2) {
            print(images,"Build 3D isosurface of image%s using isovalue %g%s.",gmic_selection,value,sep=='%'?"%":"");
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              CImg<T>& img = images[ind];
              CImg<float> vertices, opacities;
              CImgList<unsigned int> primitives;
              CImgList<unsigned char> colors;
              CImg<unsigned char> palette;
              palette.assign(3,img.spectrum(),1,1,220).noise(35,1);
              if (img.spectrum()==1) palette(0) = palette(1) = palette(2) = 255;
              else {
                palette(0,0) = 255; palette(1,0) = 30; palette(2,0) = 30;
                palette(0,1) = 30; palette(1,1) = 255; palette(2,1) = 30;
                if (img.spectrum()>=3) palette(0,2) = 30; palette(1,2) = 30; palette(2,2) = 255;
              }
              cimg_forC(img,k) {
                const CImg<T> channel = img.get_shared_channel(k);
                float nvalue = value;
                if (sep=='%') { float M = 0, m = (float)channel.min_max(M); nvalue = m + (M-m)*value/100; }
                CImgList<unsigned int> prims;
                const CImg<float> pts = channel.get_isosurface3d(prims,nvalue);
                vertices.append_object3d(primitives,pts,prims);
                colors.insert(prims.size(),CImg<unsigned char>::vector(palette(0,k),palette(1,k),palette(2,k)));
              }
              opacities.assign(1,primitives.size(),1,1,1);
              if (!vertices) warning(images,"Command 'isosurface3d' : Isovalue %g not found in image [%u].",value,ind);
              vertices.object3dtoCImg3d(primitives,colors,opacities);
              gmic_apply(img,replace(vertices));
            }
          } else if ((std::sscanf(argument,"'%4095[^']',%f%c",formula,&value,&end)==2 ||
                      std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%f,%f%c",formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&end)==8 ||
                      std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%f,%f,%d,%d,%d%c",
                                  formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&dx,&dy,&dz,&end)==11 ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%f,%f,%d%c,%d,%d%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&dx,&sepx,&dy,&dz,&end)==12 &&
                       sepx=='%') ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%f,%f,%d,%d%c,%d%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&dx,&dy,&sepy,&dz,&end)==12 &&
                       sepy=='%') ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%f,%f,%d,%d,%d%c%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&dx,&dy,&dz,&sepz,&end)==12 &&
                       sepz=='%') ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%f,%f,%d%c,%d%c,%d%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&dx,&sepx,&dy,&sepy,&dz,&end)==13 &&
                       sepx=='%' && sepy=='%') ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%f,%f,%d%c,%d,%d%c%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&dx,&sepx,&dy,&dz,&sepz,&end)==13 &&
                       sepx=='%' && sepz=='%') ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%f,%f,%d,%d%c,%d%c%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&dx,&dy,&sepy,&dz,&sepz,&end)==13 &&
                       sepy=='%' && sepz=='%') ||
                      (std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%f,%f,%d%c,%d%c,%d%c%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&dx,&sepx,&dy,&sepy,&dz,&sepz,&end)==14 &&
                       sepx=='%' && sepy=='%' && sepz=='%')) &&
                     dx>0 && dy>0 && dz>0) {
            print(images,"Build 3D isosurface %g of mathematical expression '%s' in range (%g,%g,%g)-(%g,%g,%g) with size %d%sx%d%sx%d%s.",
                  value,formula,x0,y0,z0,x1,y1,z1,dx,sepx=='%'?"%":"",dy,sepy=='%'?"%":"",dz,sepz=='%'?"%":"");
            if (sepx=='%') dx = -dx;
            if (sepy=='%') dy = -dy;
            if (sepz=='%') dz = -dz;
            CImgList<unsigned int> primitives;
            CImg<T> vertices = CImg<T>::isosurface3d(primitives,(const char*)formula,value,x0,y0,z0,x1,y1,z1,dx,dy,dz);
            CImgList<unsigned char> colors(primitives.size(),1,3,1,1,200);
            CImg<float> opacities(1,primitives.size(),1,1,1);
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            std::sprintf(tmpstr,"(3D isosurface %g of '%s')",value,formula); CImg<char>(tmpstr,std::strlen(tmpstr)+1).move_to(filenames);
          } else arg_error("isosurface3d");
          ++position; continue;
        }

        // Build 3D streamline.
        if (!std::strcmp("-streamline3d",command_name)) {
          float x0 = 0, y0 = 0, z0 = 0, L = 100, dl = 0.1f;
          unsigned int interp = 2, is_backward = 0, is_oriented = 0; char formula[4096] = { 0 };
          if ((std::sscanf(argument,"%f,%f,%f%c",&x0,&y0,&z0,&end)==3 ||
               std::sscanf(argument,"%f,%f,%f,%f%c",&x0,&y0,&z0,&L,&end)==4 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f%c",&x0,&y0,&z0,&L,&dl,&end)==5 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f,%u%c",&x0,&y0,&z0,&L,&dl,&interp,&end)==6 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f,%u,%u%c",&x0,&y0,&z0,&L,&dl,&interp,&is_backward,&end)==7 ||
               std::sscanf(argument,"%f,%f,%f,%f,%f,%u,%u,%u%c",&x0,&y0,&z0,&L,&dl,&interp,&is_backward,&is_oriented,&end)==8) &&
              L>=0 && dl>0 && interp<4) {
            print(images,"Build 3D streamline of image%s at (%g,%g,%g).",gmic_selection,x0,y0,z0);
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              CImg<T>& img = images[ind];
              CImg<T> vertices = img.get_streamline(x0,y0,z0,L,dl,interp,is_backward?true:false,is_oriented?true:false), opacities;
              CImgList<unsigned int> primitives;
              CImgList<unsigned char> colors;
              if (vertices.width()>1) {
                primitives.assign(vertices.width()-1,1,2);
                cimglist_for(primitives,l) { primitives(l,0) = l; primitives(l,1) = l+1; }
                colors.assign(primitives.size(),1,3,1,1,200);
                opacities.assign(1,primitives.size(),1,1,1);
              } else {
                vertices.assign();
                warning(images,"Command 'streamline3d' : Empty streamline starting from (%g,%g,%g) in image [%u].",x0,y0,z0,ind);
              }
              vertices.object3dtoCImg3d(primitives,colors,opacities);
              gmic_apply(img,replace(vertices));
            }
          } else if ((std::sscanf(argument,"'%4095[^']',%f,%f,%f%c",formula,&x0,&y0,&z0,&end)==4 ||
                      std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f%c",formula,&x0,&y0,&z0,&L,&end)==5 ||
                      std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f%c",formula,&x0,&y0,&z0,&L,&dl,&end)==6 ||
                      std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%u%c",formula,&x0,&y0,&z0,&L,&dl,&interp,&end)==7 ||
                      std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%u,%u%c",formula,&x0,&y0,&z0,&L,&dl,&interp,&is_backward,&end)==8 ||
                      std::sscanf(argument,"'%4095[^']',%f,%f,%f,%f,%f,%u,%u,%u%c",
                                  formula,&x0,&y0,&z0,&L,&dl,&interp,&is_backward,&is_oriented,&end)==9) &&
                     dl>0 && interp<4) {
            print(images,"Build 3D streamline of mathematical expression '%s' at (%g,%g,%g).",formula,x0,y0,z0);
            CImg<T> vertices = CImg<T>::streamline((const char *)formula,x0,y0,z0,L,dl,interp,is_backward?true:false,
                                                   is_oriented?true:false), opacities;
            CImgList<unsigned int> primitives;
            CImgList<unsigned char> colors;
            if (vertices.width()>1) {
              primitives.assign(vertices.width()-1,1,2);
              cimglist_for(primitives,l) { primitives(l,0) = l; primitives(l,1) = l+1; }
              colors.assign(primitives.size(),1,3,1,1,200);
              opacities.assign(1,primitives.size(),1,1,1);
            } else { vertices.assign(); warning(images,"Command 'streamline3d' : Empty streamline starting from (%g,%g,%g) in expression '%s'.",x0,y0,z0,formula); }
            vertices.object3dtoCImg3d(primitives,colors,opacities).move_to(images);
            std::sprintf(tmpstr,"(3D streamline of '%s' at (%g,%g,%g))",formula,x0,y0,z0);
            CImg<char>(tmpstr,std::strlen(tmpstr)+1).move_to(filenames);
          } else arg_error("streamline3d");
          ++position; continue;
        }

        // Add 3D objects together or shift a 3D object.
        if (!std::strcmp("-add3d",command_name) || !std::strcmp("-+3d",command_name)) {
          float tx = 0, ty = 0, tz = 0; int ind0 = no_ind; char sep = 0;
          if (std::sscanf(argument,"%f%c",&tx,&end)==1 ||
              std::sscanf(argument,"%f,%f%c",&tx,&ty,&end)==2 ||
              std::sscanf(argument,"%f,%f,%f%c",&tx,&ty,&tz,&end)==3) {
            print(images,"Shift 3D object%s with vector (%g,%g,%g).",gmic_selection,tx,ty,tz);
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              if (!images[ind].is_CImg3d())
                error(images,"Command 'add3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
              gmic_apply(images[ind],shiftCImg3d(tx,ty,tz));
            }
            ++position;
          } else if (std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') {
            gmic_check_indice(ind0);
            const CImg<T> img0 = images[ind0];
            if (!img0.is_CImg3d())
              error(images,"Command 'add3d' : Invalid 3D object [%d] in specified argument '%s'.",ind0,argument_text);
            print(images,"Merge object [%d] with 3D object%s.",ind0,gmic_selection);
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              const CImg<T> &img = images[ind];
              if (!img.is_CImg3d())
                error(images,"Command 'add3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
              gmic_apply(images[ind],appendCImg3d(img0));
            }
            ++position;
          } else {
            print(images,"Merge 3D object%s together.",gmic_selection);
            if (selection) {
              const unsigned int ind0 = selection[0];
              if (!images[ind0].is_CImg3d())
                error(images,"Command 'add3d' : Invalid 3D object [%d] in selected image%s.",ind0,gmic_selection);
              if (get_version) {
                CImg<T> img0 = images[ind0];
                for (unsigned int siz = selection.size(), l = 1; l<siz; ++l) img0.appendCImg3d(images[selection[l]]);
                filenames.insert(filenames[ind0].get_mark_filename());
                img0.move_to(images);
              } else for (unsigned int siz = selection.size(), off = 0, l = 1; l<siz; ++l) {
                const unsigned int ind = selection[l] - off;
                if (!images[ind].is_CImg3d())
                  error(images,"Command 'add3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
                images[ind0].appendCImg3d(images[ind]);
                filenames[ind0].mark_filename();
                images.remove(ind); filenames.remove(ind);
                ++off;
              }
            }
          }
          continue;
        }

        // Shift 3D object with opposite vector.
        if (!std::strcmp("-sub3d",command_name) || !std::strcmp("--3d",command_name)) {
          float tx = 0, ty = 0, tz = 0;
          if (std::sscanf(argument,"%f%c",&tx,&end)==1 ||
              std::sscanf(argument,"%f,%f%c",&tx,&ty,&end)==2 ||
              std::sscanf(argument,"%f,%f,%f%c",&tx,&ty,&tz,&end)==3) {
            print(images,"Shift 3D object%s with vector -(%g,%g,%g).",gmic_selection,tx,ty,tz);
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              if (!images[ind].is_CImg3d())
                error(images,"Command 'sub3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
              gmic_apply(images[ind],shiftCImg3d(-tx,-ty,-tz));
            }
          } else arg_error("sub3d");
          ++position; continue;
        }

        // Scale a 3D object.
        const bool divide3d = !std::strcmp("-div3d",command_name) || !std::strcmp("-/3d",command_name);
        if (!std::strcmp("-mul3d",command_name) || !std::strcmp("-*3d",command_name) || divide3d) {
          float sx = 0, sy = 1, sz = 1;
          if ((std::sscanf(argument,"%f%c",&sx,&end)==1 && ((sz=sy=sx),1)) ||
              std::sscanf(argument,"%f,%f%c",&sx,&sy,&end)==2 ||
              std::sscanf(argument,"%f,%f,%f%c",&sx,&sy,&sz,&end)==3) {
            if (divide3d) print(images,"Scale 3D object%s with factors (1/%g,1/%g,1/%g).",gmic_selection,sx,sy,sz);
            else print(images,"Scale 3D object%s with factors (%g,%g,%g).",gmic_selection,sx,sy,sz);
            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              CImgList<unsigned int> primitives;
              CImgList<unsigned char> colors;
              CImg<float> opacities;
              CImg<T> vertices;
              if (get_version) vertices.assign(img); else img.move_to(vertices);
              vertices.CImg3dtoobject3d(primitives,colors,opacities);
              if (divide3d) {
                vertices.get_shared_line(0)/=sx;
                vertices.get_shared_line(1)/=sy;
                vertices.get_shared_line(2)/=sz;
              } else {
                vertices.get_shared_line(0)*=sx;
                vertices.get_shared_line(1)*=sy;
                vertices.get_shared_line(2)*=sz;
              }
              vertices.object3dtoCImg3d(primitives,colors,opacities);
              if (get_version) {
                filenames.insert(filenames[selection[l]].get_mark_filename());
                vertices.move_to(images);
              } else {
                filenames[selection[l]].mark_filename();
                vertices.move_to(images[selection[l]]);
              }
            }
          } else { if (divide3d) arg_error("div3d"); else arg_error("mul3d"); }
          ++position; continue;
        }

        // Center a 3D object.
        if (!std::strcmp("-center3d",command_name) || !std::strcmp("-c3d",command_name)) {
          print(images,"Center 3D object%s.",gmic_selection);
          cimg_forY(selection,l) {
            const unsigned int ind = selection[l];
            if (!images[ind].is_CImg3d())
              error(images,"Command 'center3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
            gmic_apply(images[ind],centerCImg3d());
          }
          continue;
        }

        // Normalize a 3D object.
        if (!std::strcmp("-normalize3d",command_name) || !std::strcmp("-n3d",command_name)) {
          print(images,"Normalize 3D object%s.",gmic_selection);
          cimg_forY(selection,l) {
            const unsigned int ind = selection[l];
            if (!images[ind].is_CImg3d())
              error(images,"Command 'normalize3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
            gmic_apply(images[ind],normalizeCImg3d());
          }
          continue;
        }

        // Rotate a 3D object.
        if (!std::strcmp("-rotate3d",command_name) || !std::strcmp("-rot3d",command_name)) {
          float u = 0, v = 0, w = 1, angle = 0;
          if (std::sscanf(argument,"%f,%f,%f,%f%c",&u,&v,&w,&angle,&end)==4) {
            print(images,"Rotate 3D object%s around axis (%g,%g,%g) with angle %g.",gmic_selection,u,v,w,angle);
            const CImg<float> rot = CImg<float>::rotation_matrix(u,v,w,(float)(angle*cimg::PI/180));
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              if (!images[ind].is_CImg3d())
                error(images,"Command 'rotate3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
              gmic_apply(images[ind],rotateCImg3d(rot));
            }
          } else arg_error("rotate3d");
          ++position; continue;
        }

        // Set color of 3D objects.
        if (!std::strcmp("-color3d",command_name) || !std::strcmp("-col3d",command_name)) {
          float R = 200, G = 200, B = 200, opacity = -1;
          if (std::sscanf(argument,"%f,%f,%f%c",&R,&G,&B,&end)==3 ||
              std::sscanf(argument,"%f,%f,%f,%f%c",&R,&G,&B,&opacity,&end)==4) {
            const bool set_opacity = (opacity>=0);
            R = (float)cimg::round(R,1); G = (float)cimg::round(G,1); B = (float)cimg::round(B,1);
            if (R<0) R = 0; if (R>255) R = 255;
            if (G<0) G = 0; if (G>255) G = 255;
            if (B<0) B = 0; if (B>255) B = 255;
            if (set_opacity) print(images,"Set colors of 3D object%s to (%g,%g,%g) and opacity to %g.",gmic_selection,R,G,B,opacity);
            else print(images,"Set color of 3D object%s to (%g,%g,%g).",gmic_selection,R,G,B);
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              if (!images[ind].is_CImg3d())
                error(images,"Command 'color3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
              gmic_apply(images[ind],coloropacityCImg3d(R,G,B,opacity,true,set_opacity));
            }
          } else arg_error("color3d");
          ++position; continue;
        }

        // Set opacity of 3D objects.
        if (!std::strcmp("-opacity3d",command_name) || !std::strcmp("-o3d",command_name)) {
          float opacity = 1;
          if (std::sscanf(argument,"%f%c",&opacity,&end)==1) {
            print(images,"Set opacity of 3D object%s to %g.",gmic_selection,opacity);
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              if (!images[ind].is_CImg3d())
                error(images,"Command 'opacity3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
              gmic_apply(images[ind],coloropacityCImg3d(0,0,0,opacity,false,true));
            }
          } else arg_error("opacity3d");
          ++position; continue;
        }

        // Reverse 3D orientation.
        if (!std::strcmp("-reverse3d",command_name) || !std::strcmp("-r3d",command_name)) {
          print(images,"Reverse orientations of 3D object%s.",gmic_selection);
          cimg_forY(selection,l) {
            const unsigned int ind = selection[l];
            CImg<T> &img = images[ind];
            if (!img.is_CImg3d())
              error(images,"Command 'reverse3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
            CImgList<unsigned int> primitives;
            CImgList<unsigned char> colors;
            CImg<float> opacities;
            CImg<T> vertices;
            if (get_version) vertices.assign(img); else img.move_to(vertices);
            vertices.CImg3dtoobject3d(primitives,colors,opacities);
            if (primitives) primitives.reverse_object3d();
            vertices.object3dtoCImg3d(primitives,colors,opacities);
            if (get_version) {
              filenames.insert(filenames[selection[l]].get_mark_filename());
              vertices.move_to(images);
            } else {
              filenames[selection[l]].mark_filename();
              vertices.move_to(images[selection[l]]);
            }
          }
          continue;
        }

        // Set primitive mode for a selected 3D objects.
        if (!std::strcmp("-primitives3d",command_name) || !std::strcmp("-p3d",command_name)) {
          int mode = 0;
          if (std::sscanf(argument,"%d%c",&mode,&end)==1 &&
              mode>=0 && mode<=1) {
            print(images,"Set primitive mode '%s' for 3D object%s.",mode==0?"points":"segments",gmic_selection);
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              CImg<T> &img = images[ind];
              if (!img.is_CImg3d())
                error(images,"Command 'primitives3d' : Invalid 3D object [%d] in selected image%s.",ind,gmic_selection);
              CImgList<unsigned int> primitives;
              CImgList<unsigned char> colors;
              CImg<float> opacities;
              CImg<T> vertices;
              if (get_version) vertices.assign(img); else img.move_to(vertices);
              vertices.CImg3dtoobject3d(primitives,colors,opacities);
              CImgList<float> nopacities;
              cimglist_for(primitives,k) {
                CImg<unsigned int>& p = primitives[0];
                CImg<unsigned char>& c = colors[0];
                float opacity = opacities[k];
                switch (p.size()) {
                case 2 : // Segment.
                  if (mode==0) {
                    primitives.insert(CImg<unsigned int>::vector(p[0])).insert(CImg<unsigned int>::vector(p[1]));
                    colors.insert(2,CImg<unsigned char>::vector(c[0],c[1],c[2]));
                    nopacities.insert(2,CImg<float>::vector(opacity));
                  } else {
                    primitives.insert(CImg<unsigned int>::vector(p[0],p[1]));
                    colors.insert(CImg<unsigned char>::vector(c[0],c[1],c[2]));
                    nopacities.insert(CImg<float>::vector(opacity));
                  }
                  break;
                case 3 : // Triangle.
                  if (mode==0) {
                    primitives.insert(CImg<unsigned int>::vector(p[0])).insert(CImg<unsigned int>::vector(p[1])).
                      insert(CImg<unsigned int>::vector(p[2]));
                    colors.insert(3,CImg<unsigned char>::vector(c[0],c[1],c[2]));
                    nopacities.insert(3,CImg<float>::vector(opacity));
                  } else {
                    primitives.insert(CImg<unsigned int>::vector(p[0],p[1])).insert(CImg<unsigned int>::vector(p[1],p[2])).
                      insert(CImg<unsigned int>::vector(p[2],p[0]));
                    colors.insert(3,CImg<unsigned char>::vector(c[0],c[1],c[2]));
                    nopacities.insert(3,CImg<float>::vector(opacity));
                  }
                  break;
                case 4 : // Quadrangle.
                  if (mode==0) {
                    primitives.insert(CImg<unsigned int>::vector(p[0])).insert(CImg<unsigned int>::vector(p[1])).
                      insert(CImg<unsigned int>::vector(p[2])).insert(CImg<unsigned int>::vector(p[3]));
                    colors.insert(4,CImg<unsigned char>::vector(c[0],c[1],c[2]));
                    nopacities.insert(4,CImg<float>::vector(opacity));
                  } else {
                    primitives.insert(CImg<unsigned int>::vector(p[0],p[1])).insert(CImg<unsigned int>::vector(p[1],p[2])).
                      insert(CImg<unsigned int>::vector(p[2],p[3])).insert(CImg<unsigned int>::vector(p[3],p[0]));
                    colors.insert(4,CImg<unsigned char>::vector(c[0],c[1],c[2]));
                    nopacities.insert(4,CImg<float>::vector(opacity));
                  }
                  break;
                default : // Other primitives.
                  primitives.insert(p);
                  colors.insert(c);
                  nopacities.insert(CImg<float>::vector(opacity));
                }
                primitives.remove(0); colors.remove(0);
              }
              opacities = nopacities>'y';
              vertices.object3dtoCImg3d(primitives,colors,opacities);
              if (get_version) {
                filenames.insert(filenames[selection[l]].get_mark_filename());
                vertices.move_to(images);
              } else {
                filenames[selection[l]].mark_filename();
                vertices.move_to(images[selection[l]]);
              }
            }
          } else arg_error("primitives3d");
          ++position; continue;
        }

        // Split 3D objects into 6 vector images {header,N,vertices,primitives,colors,opacities}
        if (!std::strcmp("-split3d",command_name) || !std::strcmp("-s3d",command_name)) {
          print(images,"Split 3D object%s into its different characteristics.",gmic_selection);
          unsigned int off = 0;
          cimg_forY(selection,l) {
            const unsigned int ind = selection[l] + off;
            CImg<T> &img = images[ind];
            if (!img.is_CImg3d())
              error(images,"Command 'split3d' : Invalid 3D object [%d] in selected image%s.",ind-off,gmic_selection);
            const CImg<char> filename = filenames[ind].get_mark_filename();
            CImgList<unsigned int> primitives;
            CImgList<unsigned char> colors;
            CImg<float> opacities;
            CImg<T> vertices;
            if (get_version) vertices.assign(img); else img.move_to(vertices);
            vertices.CImg3dtoobject3d(primitives,colors,opacities);
            CImgList<T> split;
            (CImg<T>("CImg3d",1,6)+=0.5f).move_to(split);
            CImg<T>::vector((T)vertices.width(),(T)primitives.size()).move_to(split);
            vertices.resize(-100,3,1,1,0).transpose().unroll('y').move_to(split);
            CImgList<T> _prims;
            cimglist_for(primitives,p) {
              CImg<T>::vector((T)primitives[p].size()).move_to(_prims);
              primitives[p].unroll('y').move_to(_prims);
            }
            primitives.assign();
            (_prims>'y').move_to(split);
            (colors>'x').transpose().unroll('y').move_to(split); colors.assign();
            opacities.move_to(split);
            if (get_version) {
              filenames.insert(split.size(),filename);
              split.move_to(images,~0U);
            } else {
              off+=split.size() - 1;
              filenames.remove(ind); filenames.insert(split.size(),filename,ind);
              images.remove(ind); split.move_to(images,ind);
            }
          }
          continue;
        }

        // Set 3D light position.
        if (!std::strcmp("-light3d",item) || !std::strcmp("-l3d",item)) {
          float lx = 0, ly = 0, lz = -5000;
          if (std::sscanf(argument,"%f,%f,%f%c",&lx,&ly,&lz,&end)==3) {
            print(images,"Set 3D light position at (%g,%g,%g).",lx,ly,lz);
            light3d_x = lx;
            light3d_y = ly;
            light3d_z = lz;
          } else arg_error("light3d");
          ++position; continue;
        }

        // Set 3D focale.
        if (!std::strcmp("-focale3d",item) || !std::strcmp("-f3d",item)) {
          float focale = 500;
          if (std::sscanf(argument,"%f%c",&focale,&end)==1 && focale>0) {
            focale3d = focale;
            print(images,"Set 3D focale to %g.",focale);
          } else arg_error("focale3d");
          ++position; continue;
        }

        // Set 3D specular light parameters.
        if (!std::strcmp("-specl3d",item) || !std::strcmp("-sl3d",item)) {
          float value = 0;
          if (std::sscanf(argument,"%f%c",&value,&end)==1 && value>=0) {
            specular_light3d = value;
            print(images,"Set amount of 3D specular light to %g.",specular_light3d);
          } else arg_error("specl3d");
          ++position; continue;
        }

        if (!std::strcmp("-specs3d",item) || !std::strcmp("-ss3d",item)) {
          float value = 0;
          if (std::sscanf(argument,"%f%c",&value,&end)==1 && value>=0) {
            specular_shine3d = value;
            print(images,"Set shininess of 3D specular light to %g.",specular_shine3d);
          }
          else arg_error("specs3d");
          ++position; continue;
        }

        // Switch double-sided mode for 3D rendering.
        if (!std::strcmp("-double3d",item) || !std::strcmp("-db3d",item)) {
          int value = 0;
          if (std::sscanf(argument,"%d%c",&value,&end)==1) {
            is_double3d = value?true:false;
            print(images,"%s double-sided mode for 3D rendering.",is_double3d?"Enable":"Disable");
          } else arg_error("double3d");
          ++position; continue;
        }

        // Set 3D rendering mode.
        if (!std::strcmp("-mode3d",item) || !std::strcmp("-m3d",item)) {
          int value = 0;
          if (std::sscanf(argument,"%d%c",&value,&end)==1 && value>=-1 && value<=5) {
            render3d = value;
            print(images,"Set static 3D render mode to %s.",
                  render3d==-1?"bounding-box":
                  render3d==0?"pointwise":render3d==1?"linear":render3d==2?"flat":
                  render3d==3?"flat-shaded":render3d==4?"Gouraud-shaded":
                  render3d==5?"Phong-shaded":"none");
          } else arg_error("mode3d");
          ++position; continue;
        }

        if (!std::strcmp("-moded3d",item) || !std::strcmp("-md3d",item)) {
          int value = 0;
          if (std::sscanf(argument,"%d%c",&value,&end)==1 && value>=-1 && value<=5) {
            renderd3d = value;
            print(images,"Set dynamic 3D render mode to %s.",
                  renderd3d==-1?"bounding-box":
                  renderd3d==0?"pointwise":renderd3d==1?"linear":renderd3d==2?"flat":
                  renderd3d==3?"flat-shaded":renderd3d==4?"Gouraud-shaded":
                  renderd3d==5?"Phong-shaded":"none");
          } else arg_error("moded3d");
          ++position; continue;
        }

        // Set 3D background color.
        if (!std::strcmp("-background3d",item) || !std::strcmp("-b3d",item)) {
          int R = 0, G = 0, B = 0;
          const int nb = std::sscanf(argument,"%d,%d,%d%c",&R,&G,&B,&end);
          if (nb>=1 && nb<=3) {
            switch (nb) {
            case 1 : background3d[0] = background3d[1] = background3d[2] = R; break;
            case 2 : background3d[0] = R; background3d[1] = background3d[2] = G; break;
            case 3 : background3d[0] = R; background3d[1] = G; background3d[2] = B; break;
            }
            print(images,"Set 3D background color to (%d,%d,%d).",
                  (int)background3d[0],(int)background3d[1],(int)background3d[2]);
          } else arg_error("background3d");
          ++position; continue;
        }

        //----------------------
        // Procedural commands.
        //----------------------

        // No operations : do nothing
        if (!std::strcmp("-nop",item)) {
          if (verbosity_level>0) print(images,"Do nothing.");
          continue;
        }

        // Skip next argument;
        if (!std::strcmp("-skip",item)) {
          if (verbosity_level>0) print(images,"Skip argument '%s'.",argument_text);
          ++position;
          continue;
        }

        // Echo.
        if (!std::strcmp("-echo",item) || !std::strcmp("-e",item)) {
          print(images,"%s",argument);
          ++position; continue;
        }

        // Print.
        if (!std::strcmp("-print",command_name)) {
          if (images.size()) {
            print(images,"Print image%s.\n\n",gmic_selection);
            if (verbosity_level>=0) cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              std::sprintf(tmpstr,"image [%u] = '%s'",ind,filenames[ind].data());
              images[ind].print(tmpstr);
            }
            is_released = true;
          } else print(images,"Print image[].");
          continue;
        }

        // Return.
        if (!std::strcmp("-return",item)) {
          if (verbosity_level>0) print(images,"Return.");
          dowhiles.assign();
          repeatdones.assign();
          locals.assign();
          position = command_line.size();
          continue;
        }

        // Quit.
        if (!std::strcmp("-quit",item) || !std::strcmp("-q",item) || *cancel) {
          print(images,"Quit.");
          is_released = is_quit = true;
          dowhiles.assign();
          repeatdones.assign();
          locals.assign();
          position = command_line.size();
          continue;
        }

        // Exec.
        if (!std::strcmp("-exec",item)) {
          print(images,"Execute external command '%s'.\n",argument);
          int err = cimg::system(argument);
          err = ++position; continue;
        }

        // Do..while.
        if (!std::strcmp("-do",item)) {
          if (verbosity_level>0) print(images,"Do : Start do..while loop.");
          CImg<unsigned int>::vector(position).move_to(dowhiles);
          continue;
        }

        if (!std::strcmp("-while",item)) {
          if (!dowhiles) error(images,"Command 'while' : Missing associated 'do' command.");
          float number = 0; bool cond = false;
          if (std::sscanf(argument,"%f%c",&number,&end)==1) {
            cond = (bool)number;
            if (verbosity_level>0) print(images,"While : Check '%s' -> %s.",argument_text,cond?"true":"false");
            if (cond) { position = dowhiles.back()(0); continue; }
            else dowhiles.remove();
          } else arg_error("while");
          ++position; continue;
        }

        // If..[elif]..[else]..endif
        if (!std::strcmp("-if",item) || (!std::strcmp("-elif",item) && check_elif)) {
          check_elif = false;
          float number = 0; bool cond = false;
          if (std::sscanf(argument,"%f%c",&number,&end)==1) {
            cond = (bool)number;
            if (verbosity_level>0) print(images,"If : Check '%s' -> %s.",argument_text,cond?"true":"false");
            if (!cond) {
              for (int nbifs = 1; nbifs && position<command_line.size(); ++position) {
                const char *const it = command_line[position].data();
                if (!std::strcmp("-if",it)) ++nbifs;
                if (!std::strcmp("-endif",it)) --nbifs;
                if (nbifs==1) {
                  if (!std::strcmp("-else",it)) --nbifs;
                  if (!std::strcmp("-elif",it)) { --nbifs; check_elif = true; --position;}
                }
              }
              continue;
            }
          } else arg_error("if");
          ++position; continue;
        }

        if (!std::strcmp("-else",item) || !std::strcmp("-elif",item)) {
          check_elif = false;
          for (int nbifs = 1; nbifs && position<command_line.size(); ++position) {
            if (!std::strcmp("-if",command_line[position].data())) ++nbifs;
            if (!std::strcmp("-endif",command_line[position].data())) --nbifs;
          }
          continue;
        }

        if (!std::strcmp("-endif",item)) {
          check_elif = false;
          if (verbosity_level>0) print(images,"End if.");
          continue;
        }

        // Repeat..[break]..done
        if (!std::strcmp("-repeat",item)) {
          float number = 0;
          if (std::sscanf(argument,"%f%c",&number,&end)==1) {
            const int nb = (int)number;
            if (verbosity_level>0) print(images,"Start %d iteration%s of a repeat..done loop.",nb,nb>1?"s":"");
            if (nb>0) CImg<unsigned int>::vector(position+1,nb,0).move_to(repeatdones);
            else {
              int nbrepeats = 0;
              for (nbrepeats = 1; nbrepeats && position<command_line.size(); ++position) {
                const char *it = command_line[position].data();
                if (!std::strcmp("-repeat",it)) ++nbrepeats;
                if (!std::strcmp("-done",it)) --nbrepeats;
              }
              if (nbrepeats && position>=command_line.size())
                error(images,"Command 'repeat' : Missing associated 'done' command.");
              continue;
            }
          } else arg_error("repeat");
          ++position; continue;
        }

        if (!std::strcmp("-done",item)) {
          if (!repeatdones) error(images,"Command 'done' : Missing associated 'repeat' command.");
          if (--repeatdones.back()(1)) {
            ++repeatdones.back()(2);
            position = repeatdones.back()(0);
          }
          else repeatdones.remove();
          continue;
        }

        if (!std::strcmp("-error",item)) {
          error(images,argument);
        }

        if (!std::strcmp("-warning",item)) {
          warning(images,argument);
          ++position; continue;
        }

        if (!std::strcmp("-check",item)) {
          bool cond = false;
          try { cond = cimg::eval(argument); } catch (CImgException&) { cond = false; }
          if (!cond) error(images,"Command 'check' : Expression '%s' is not verified.",argument_text);
          ++position; continue;
        }

        // Set progress bar.
        if (!std::strcmp("-progress",item)) {
          float value = -1;
          if (std::sscanf(argument,"%f%c",&value,&end)==1) {
            if (value<0) value = -1; else if (value>100) value = 100;
            *progress = value;
          } else arg_error("progress");
          ++position; continue;
        }

        // Push/pop variable to/from global stack.
        if (!std::strcmp("-push",command_name) || !std::strcmp("-p",command_name)) {
          if (!is_restriction) selection.assign(1,1,1,1,stack.size());
          print(images,"Push item '%s' on the global stack at position%s.",argument_text,gmic_selection);
          cimg_forY(selection,l) CImg<char>(argument,std::strlen(argument)+1).move_to(stack,selection[l]);
          ++position; continue;
        }

        if (!std::strcmp("-pop",command_name) || !std::strcmp("-pp",command_name)) {
          if (!is_restriction) CImg<unsigned int>::sequence(stack.size(),0,stack.size()-1).move_to(selection);
          print(images,"Pop item%s from the global stack",gmic_selection);
          unsigned int off = 0;
          cimg_forY(selection,l) { const unsigned int ind = selection[l] - off; stack.remove(ind); ++off; }
          if (verbosity_level>=0) {
            std::fprintf(cimg::output()," (%u item%s left).",stack.size(),stack.size()==1?"":"s");
            std::fflush(cimg::output());
          }
          continue;
        }

        // Handle local environnements.
        if (!std::strcmp("-local",command_name) || !std::strcmp("-l",command_name)) {
          print(images,"Start local environment with image%s.",gmic_selection);
          locals.insert(selection);
          const unsigned int siz = selection.size();
          CImgList<T> nimages(siz);
          CImgList<char> nfilenames(siz);
          if (get_version) {
            cimg_forY(selection,l) { nimages[l].assign(images[selection[l]]); nfilenames[l].assign(filenames[selection[l]]); }
            parse(command_line,position,nimages,nfilenames,dowhiles,repeatdones,locals,false);
            nimages.move_to(images,~0U); nfilenames.move_to(filenames,~0U);
          } else {
            cimg_forY(selection,l) { nimages[l].swap(images[selection[l]]); nfilenames[l].swap(filenames[selection[l]]); }
            parse(command_line,position,nimages,nfilenames,dowhiles,repeatdones,locals,false);
            const unsigned int nb = cimg::min(siz,nimages.size());
            for (unsigned int i = 0; i<nb; ++i) {
              images[selection[i]].swap(nimages[0]); filenames[selection[i]].swap(nfilenames[0]);
              nimages.remove(0); nfilenames.remove(0);
            }
            if (nb<siz) for (unsigned int off = 0, l = nb; l<siz; ++l, ++off) {
                const unsigned int ind = selection[l] - off;
                images.remove(ind); filenames.remove(ind);
              } else if (nimages) {
              const unsigned int ind0 = siz?selection[siz-1]+1:images.size();
              images.insert(nimages,ind0); filenames.insert(nimages.size(),CImg<char>("(unnamed)",10),ind0);
            }
          }
          continue;
        }

        if (!std::strcmp("-endlocal",item) || !std::strcmp("-endl",item)) {
          if (locals) {
            print(images,"End local environment with image%s.",selection2string(locals.back(),filenames,true));
            locals.remove();
            return *this;
          } else error(images,"Command 'endlocal' : No local environment has been previously started.");
          continue;
        }

        // Handle patch processing mode.
        if (!std::strcmp("-patch",item)) {
          int sx = 0, sy = 0, sz = 0, sv = 0; unsigned int borders = 0;
          if (std::sscanf(argument,"%d,%d,%d,%d,%u%c",&sx,&sy,&sz,&sv,&borders,&end)==5 &&
              sx>0 && sy>0 && sz>0 && sv>0) {
            print(images,"Enable %dx%dx%dx%d patch processing, with %s border conditions.",sx,sy,sz,sv,borders?"Neumann":"Dirichlet");
            patch_borders = borders;
            patch_w = (unsigned int)sx; patch_h = (unsigned int)sy;
            patch_d = (unsigned int)sz; patch_c = (unsigned int)sv;
          } else if (std::sscanf(argument,"%d,%d,%d,%u%c",&sx,&sy,&sz,&borders,&end)==4 &&
                     sx>0 && sy>0 && sz>0) {
            print(images,"Enable %dx%dx%d patch processing, with %s border conditions.",sx,sy,sz,borders?"Neumann":"Dirichlet");
            patch_borders = borders;
            patch_w = (unsigned int)sx; patch_h = (unsigned int)sy;
            patch_d = (unsigned int)sz; patch_c = 0;
          } else if (std::sscanf(argument,"%d,%d,%u%c",&sx,&sy,&borders,&end)==3 &&
                     sx>0 && sy>0) {
            print(images,"Enable %dx%d patch processing, with %s border conditions.",sx,sy,borders?"Neumann":"Dirichlet");
            patch_borders = borders;
            patch_w = (unsigned int)sx; patch_h = (unsigned int)sy;
            patch_d = patch_c = 0;
          } else if (std::sscanf(argument,"%d,%u%c",&sx,&borders,&end)==2 &&
                     sx>0) {
            print(images,"Enable %d patch processing, with %s border conditions.",sx,borders?"Neumann":"Dirichlet");
            patch_borders = borders;
            patch_w = (unsigned int)sx;
            patch_h = patch_d = patch_c = 0;
          } else arg_error("patch");
          ++position; continue;
        }

        if (!std::strcmp("-endpatch",item)) {
          print(images,"Disable patch processing.");
          patch_w = patch_h = patch_d = patch_c = 0; patch_borders = false;
          continue;
        }

        //--------------------------
        // Input/Output and Display
        //--------------------------

        // Display.
        if (!std::strcmp("-display",command_name) || !std::strcmp("-d",command_name)) {
          is_released |= display_images(images,filenames,selection,true);
          continue;
        }

        // Display images in the instant display window.
        unsigned int wind = 0;
        if ((!std::strcmp("-window",command_name) || !std::strcmp("-w",command_name) ||
             std::sscanf(command_name,"-window%u%c",&wind,&end)==1 ||
             std::sscanf(command_name,"-w%u%c",&wind,&end)==1) &&
            wind<10) {
          int dimw = -1, dimh = -1, norm = -1, fs = -1;
          if (std::sscanf(argument,"%d%c",&dimw,&end)==1) { ++position; dimh = dimw; }
          else if ((std::sscanf(argument,"%d,%d%c",&dimw,&dimh,&end)==2 ||
                    std::sscanf(argument,"%d,%d,%d%c",&dimw,&dimh,&norm,&end)==3 ||
                    std::sscanf(argument,"%d,%d,%d,%d%c",&dimw,&dimh,&norm,&fs,&end)==4) &&
                   dimw>=-1 && dimh>=-1 && norm>=-1 && norm<=3) ++position;
          else dimw = dimh = norm = -1;
          if (dimh==0) dimw = 0; else if (dimh==-1) dimw = -1;
#if cimg_display==0
          print(images,"Display image%s in instant window [%d] (skipped, no display available).",gmic_selection,wind);
#else
          if (!dimw) { // Close.
            print(images,"Close instant window.");
            instant_window[wind].assign();
          } else {
            CImgList<T> subimages;
            cimg_forY(selection,l) subimages.insert(images[selection[l]],~0U,true);
            const char *const title = selection2string(selection,filenames,false);
            if (instant_window[wind]) { // Resize and refresh.
              if (dimw>0) instant_window[wind].resize(dimw,dimh,false); else instant_window[wind].resize(false);
              if (norm>=0) instant_window[wind]._normalization = norm;
              if (std::strcmp(instant_window[wind].title(),title)) instant_window[wind].set_title(title);
              if (fs>=0 && (bool)fs!=instant_window[wind].is_fullscreen()) instant_window[wind].toggle_fullscreen(false);
            } else { // Create.
              if (dimw<0) {
                dimw = dimh = 0;
                if (selection) cimg_forY(selection,l) {
                    const CImg<T>& img = images[selection[l]];
                    dimw+=img.width(); if (img.height()>dimh) dimh = img.height();
                  } else dimw = dimh = 256;
              }
              instant_window[wind].assign(dimw,dimh,title,norm<0?3:norm,fs<0?false:(bool)fs);
              if (norm==2) {
                if (subimages) instant_window[wind]._min = (float)subimages.min_max(instant_window[wind]._max);
                else { instant_window[wind]._min = 0; instant_window[wind]._max = 255; }
              }
            }
            print(images,"Display image%s in %dx%d %sinstant window [%d], with%snormalization.",
                  gmic_selection,instant_window[wind].width(),instant_window[wind].height(),
                  instant_window[wind].is_fullscreen()?"fullscreen ":"",wind,
                  instant_window[wind].normalization()==0?"out ":
                  instant_window[wind].normalization()==1?" ":
                  instant_window[wind].normalization()==2?" 1st-time ":" auto-");
            if (subimages) subimages.display(instant_window[wind]);
          }
#endif
          is_released = true;
          continue;
        }

        // Wait a given delay of for a user event from the instant display window.
        if (!std::strcmp("-wait",command_name)) {
          if (!is_restriction) CImg<unsigned int>::vector(0,1,2,3,4,5,6,7,8,9).move_to(selection);
          int delay = 0;
          if (std::sscanf(argument,"%d%c",&delay,&end)==1) ++position; else delay = 0;
#if cimg_display==0
          if (!delay) print(images,"Wait for a user event from the instant window%s (skipped, no display available).",gmic_selection);
          else { print(images,"Wait for %d milliseconds from instant window%s.",delay,gmic_selection); cimg::wait(delay<0?-delay:delay); }
#else
          if (!delay) {
            print(images,"Wait for any user event from the instant window%s.",gmic_selection);
            CImgDisplay *const iw = instant_window;
            switch (selection.size()) {
            case 1 : CImgDisplay::wait(iw[selection[0]]); break;
            case 2 : CImgDisplay::wait(iw[selection[0]],iw[selection[1]]); break;
            case 3 : CImgDisplay::wait(iw[selection[0]],iw[selection[1]],iw[selection[2]]); break;
            case 4 : CImgDisplay::wait(iw[selection[0]],iw[selection[1]],iw[selection[2]],iw[selection[3]]); break;
            case 5 : CImgDisplay::wait(iw[selection[0]],iw[selection[1]],iw[selection[2]],iw[selection[3]],iw[selection[4]]); break;
            case 6 : CImgDisplay::wait(iw[selection[0]],iw[selection[1]],iw[selection[2]],iw[selection[3]],iw[selection[4]],
                                       iw[selection[5]]); break;
            case 7 : CImgDisplay::wait(iw[selection[0]],iw[selection[1]],iw[selection[2]],iw[selection[3]],iw[selection[4]],
                                       iw[selection[5]],iw[selection[6]]); break;
            case 8 : CImgDisplay::wait(iw[selection[0]],iw[selection[1]],iw[selection[2]],iw[selection[3]],iw[selection[4]],
                                       iw[selection[5]],iw[selection[6]],iw[selection[7]]); break;
            case 9 : CImgDisplay::wait(iw[selection[0]],iw[selection[1]],iw[selection[2]],iw[selection[3]],iw[selection[4]],
                                       iw[selection[5]],iw[selection[6]],iw[selection[7]],iw[selection[8]]); break;
            case 10 : CImgDisplay::wait(iw[selection[0]],iw[selection[1]],iw[selection[2]],iw[selection[3]],iw[selection[4]],
                                        iw[selection[5]],iw[selection[6]],iw[selection[7]],iw[selection[8]],iw[selection[9]]); break;
            }
          } else if (delay<0) {
            print(images,"Wait for %d milliseconds and flush display events of instant window%s.",-delay,gmic_selection);
            if (selection && instant_window[selection[0]]) instant_window[selection[0]].wait(-delay); else cimg::wait(-delay);
            cimg_forY(selection,l) instant_window[selection[l]].flush();
          } else {
            print(images,"Wait for %d milliseconds from instant window%s",delay,gmic_selection);
            if (selection && instant_window[selection[0]]) instant_window[selection[0]].wait(delay); else cimg::wait(delay);
          }
#endif
          continue;
        }

        // Display 3D object.
        if (!std::strcmp("-display3d",command_name) || !std::strcmp("-d3d",command_name)) {
          is_released |= display_objects3d(images,filenames,selection,true);
          continue;
        }

        // Display as a graph plot.
        if (!std::strcmp("-plot",command_name)) {
          unsigned int plot_type = 1, vertex_type = 1; int resolution = 65536;
          double ymin = 0, ymax = 0, xmin = 0, xmax = 0; char sep = 0, formula[4096] = { 0 };
          if (((std::sscanf(argument,"'%1023[^']%c%c",formula,&sep,&end)==2 && sep=='\'') ||
               std::sscanf(argument,"'%1023[^']',%lf,%lf%c",formula,&xmin,&xmax,&end)==3 ||
               std::sscanf(argument,"'%1023[^']',%lf,%lf,%lf,%lf%c",formula,&xmin,&xmax,&ymin,&ymax,&end)==5 ||
               std::sscanf(argument,"'%1023[^']',%lf,%lf,%lf,%lf,%d%c",
                           formula,&xmin,&xmax,&ymin,&ymax,&resolution,&end)==6 ||
               std::sscanf(argument,"'%1023[^']',%lf,%lf,%lf,%lf,%d,%u%c",
                           formula,&xmin,&xmax,&ymin,&ymax,&resolution,&plot_type,&end)==7 ||
               std::sscanf(argument,"'%1023[^']',%lf,%lf,%lf,%lf,%d,%u,%u%c",
                           formula,&xmin,&xmax,&ymin,&ymax,&resolution,&plot_type,&vertex_type,&end)==8) &&
              resolution>0 && plot_type<=3 && vertex_type<=7) {
            if (xmin==0 && xmax==0) { xmin = -4; xmax = 4; }
            if (!plot_type && !vertex_type) plot_type = 1;
            if (resolution<1) resolution = 65536;
            CImgList<double> tmp_img(1);
            CImg<double> &img = tmp_img[0];
            img.assign(resolution--).eval(formula);
            const double dx = xmax - xmin;
            cimg_forX(img,X) img(X) = img.eval(0,xmin+X*dx/resolution);
            CImgList<char> tmp_filename;
            CImg<char>(formula,std::strlen(formula)+1).move_to(tmp_filename);
            is_released |= display_plots(tmp_img,tmp_filename,
                                         CImg<unsigned int>::vector(0),plot_type,vertex_type,xmin,xmax,ymin,ymax,true);
            ++position;
          } else {
            plot_type = 1; vertex_type = 0; ymin = ymax = xmin = xmax = 0;
            if ((std::sscanf(argument,"%u%c",&plot_type,&end)==1 ||
                 std::sscanf(argument,"%u,%u%c",&plot_type,&vertex_type,&end)==2 ||
                 std::sscanf(argument,"%u,%u,%lf,%lf%c",&plot_type,&vertex_type,&xmin,&xmax,&end)==4 ||
                 std::sscanf(argument,"%u,%u,%lf,%lf,%lf,%lf%c",&plot_type,&vertex_type,&xmin,&xmax,&ymin,&ymax,&end)==6) &&
                plot_type<=3 && vertex_type<=7) ++position;
            if (!plot_type && !vertex_type) plot_type = 1;
            is_released |= display_plots(images,filenames,selection,plot_type,vertex_type,xmin,xmax,ymin,ymax,true);
          }
          continue;
        }

        // Select image feature.
        if (!std::strcmp("-select",command_name)) {
          unsigned int select_type = 0;
          if (std::sscanf(argument,"%u%c",&select_type,&end)==1 &&
              select_type<=3) {
            print(images,"Select %s in image%s.",select_type==0?"point":select_type==1?"segment":select_type==2?"rectangle":"ellipse",gmic_selection);
#if cimg_display==0
            cimg_forY(selection,l) gmic_apply(images[selection[l]],select(filenames[selection[l]].data(),select_type));
#else
            if (instant_window[0]) { cimg_forY(selection,l) gmic_apply(images[selection[l]],select(instant_window[0],select_type)); }
            else { cimg_forY(selection,l) gmic_apply(images[selection[l]],select(filenames[selection[l]].data(),select_type)); }
#endif
          } else arg_error("select");
          ++position; continue;
        }

        // Output.
        if (!std::strcmp("-output",command_name) || !std::strcmp("-o",command_name)) {
          char filename[4096] = { 0 }; char options[4096] = { 0 };
          if (std::sscanf(argument,"%4095[^,],%s",filename,options)!=2) std::strcpy(filename,argument);
          for (char *s = filename; *s; ++s) if (*s==29) *s=',';
          const char *const ext = cimg::split_filename(filename);
          if (!cimg::strcasecmp("off",ext)) {
            char nfilename[4096] = { 0 };
            std::strcpy(nfilename,filename);
            const unsigned int siz = selection.size();
            cimg_forY(selection,l) {
              const unsigned int ind = selection[l];
              if (siz!=1) cimg::number_filename(filename,l,6,nfilename);
              if (!images[ind].is_CImg3d())
                error(images,"Command 'output' : 3D object file '%s', invalid 3D object [%u] in selected image%s.",nfilename,ind,gmic_selection);
              print(images,"Output 3D object [%u] as file '%s'.",ind,nfilename);
              CImgList<unsigned int> primitives;
              CImgList<unsigned char> colors;
              CImg<float> opacities, vertices(images[ind]);
              vertices.CImg3dtoobject3d(primitives,colors,opacities).save_off(nfilename,primitives,colors);
            }
          } else if (!cimg::strcasecmp("jpeg",ext) || !cimg::strcasecmp("jpg",ext)) {
            int quality = 100;
            if (std::sscanf(options,"%d%c",&quality,&end)!=1) quality = 100;
            if (quality<0) quality = 0; else if (quality>100) quality = 100;
            CImgList<T> output_images;
            cimg_forY(selection,l) output_images.insert(images[selection[l]],~0U,true);
            print(images,"Output image%s as file '%s', with quality %u%%",gmic_selection,filename,quality);
            if (!output_images) throw CImgInstanceException("CImgList<%s>::save() : File '%s, instance list (%u,%p) is empty.",
                                                            output_images.pixel_type(),filename,
                                                            output_images.size(),output_images.data());
            if (output_images.size()==1) output_images[0].save_jpeg(filename,quality);
            else {
              char nfilename[1024];
              cimglist_for(output_images,l) {
                cimg::number_filename(filename,l,6,nfilename);
                output_images[l].save_jpeg(nfilename,quality);
              }
            }
          } else {
            CImgList<T> output_images;
            cimg_forY(selection,l) output_images.insert(images[selection[l]],~0U,true);
            print(images,"Output image%s as file '%s'.",gmic_selection,filename);
            output_images.save(filename);
          }
          is_released = true; ++position; continue;
        }

        // Run custom command, if found.
        if (std::strcmp("-i",command_name) && std::strcmp("-input",command_name)) {
          const char *custom_command_name = 0;
          bool custom_command_found = false, has_arguments = false;
          CImg<char> substituted_command;
          cimglist_for(command_names,l) {
            custom_command_name = command_names[l].data();
            const char *const command = command_definitions[l].data();

            if (!std::strcmp(command_name+1,custom_command_name) && *command) {
              CImgList<char> arguments(256);
              unsigned int nb_arguments = 0;
              char s_argument[4096] = { 0 };
              custom_command_found = true;
              debug(images,"Substitute command '%s' by '%s'.",custom_command_name,command);

              // Get command-line values of custom command arguments.
              if (argument)
                for (const char *nargument = argument; nb_arguments<255 && *nargument; ) {
                  char *ns_argument = s_argument;
                  for (bool is_dquoted = false; *nargument && (*nargument!=',' || is_dquoted); ++nargument) {
                    if (*nargument=='\"') is_dquoted = !is_dquoted;
                    *(ns_argument++) = *nargument;
                  }
                  if (ns_argument!=s_argument) {
                    *ns_argument = 0;
                    CImg<char>(s_argument,std::strlen(s_argument)+1).move_to(arguments[++nb_arguments]);
                  }
                  if (*nargument) ++nargument;
                }

              // Substitute arguments in custom command expression.
              CImgList<char> lreplacement;
              for (const char *ncommand = command; *ncommand;) if (*ncommand=='$') {
                char *replace_text = &(s_argument[0] = 0), sep = 0; int ind = 0, ind1 = 0;

                // Substitute $?.
                if (ncommand[1]=='?') {
                  std::sprintf(s_argument,"%s",gmic_selection);
                  ncommand+=2;

                // Substitute $#.
                } else if (ncommand[1]=='#') {
                  std::sprintf(s_argument,"%u",nb_arguments);
                  ncommand+=2;
                  has_arguments = true;

                  // Substitute $*.
                } else if (ncommand[1]=='*') {
                  for (unsigned int j = 1; j<=nb_arguments; ++j) {
                    replace_text+=std::sprintf(replace_text,"%s",arguments[j].data());
                    if (j<nb_arguments) *(replace_text++) = ',';
                  }
                  replace_text = s_argument;
                  ncommand+=2;
                  has_arguments = true;

                  // Substitute ${i*}.
                } else if (std::sscanf(ncommand,"${%d*%c",&ind,&sep)==2 &&
                           ind>0 && ind<256 && sep=='}') {
                  for (unsigned int j = (unsigned int)ind; j<=nb_arguments; ++j) {
                    replace_text+=std::sprintf(replace_text,"%s",arguments[j].data());
                    if (j<nb_arguments) *(replace_text++) = ',';
                  }
                  replace_text = s_argument;
                  ncommand+=std::sprintf(tmpstr,"${%d*}",ind);
                  has_arguments = true;

                  // Substitute $i and ${i}.
                } else if ((std::sscanf(ncommand,"$%d",&ind)==1 ||
                            (std::sscanf(ncommand,"${%d%c",&ind,&sep)==2 && sep=='}')) &&
                           ind>0 && ind<256) {
                  if (!arguments[ind]) {
                    if (sep=='}') error(images,"Command '%s' : Argument '$%d' is undefined (in expression '${%d}').",
                                        custom_command_name,ind,ind);
                    else error(images,"Command '%s' : Argument '$%d' is undefined (in expression '$%d').",
                               custom_command_name,ind,ind);
                  }
                  replace_text = arguments[ind].data();
                  ncommand+=std::sprintf(tmpstr,"$%d",ind) + (sep=='}'?2:0);
                  if (ind>0) has_arguments = true;

                  // Substitute ${i=$j}.
                } else if (std::sscanf(ncommand,"${%d=$%d%c",&ind,&ind1,&sep)==3 && sep=='}' &&
                           ind>0 && ind<256 && ind1>0 && ind1<256) {
                  if (!arguments[ind1])
                    error(images,"Command '%s' : Argument '$%d' is undefined (in expression '${%d=$%d}').",
                          custom_command_name,ind1,ind,ind1);
                  if (!arguments[ind]) arguments[ind] = arguments[ind1];
                  replace_text = arguments[ind].data();
                  ncommand+=std::sprintf(tmpstr,"${%d=$%d}",ind,ind1);
                  has_arguments = true;

                  // Substitute ${i=$#}.
                } else if (std::sscanf(ncommand,"${%d=$#%c",&ind,&sep)==2 && sep=='}' &&
                           ind>0 && ind<256) {
                  if (!arguments[ind]) {
                    std::sprintf(s_argument,"%u",nb_arguments);
                    CImg<char>(s_argument,std::strlen(s_argument)+1).move_to(arguments[ind]);
                  }
                  replace_text = arguments[ind].data();
                  ncommand+=std::sprintf(tmpstr,"${%d=$#}",ind);
                  has_arguments = true;

                  // Substitute ${i=default}.
                } else if (std::sscanf(ncommand,"${%d=%4095[^}]%c",&ind,tmpstr,&sep)==3 && sep=='}' &&
                           ind>0 && ind<256) {
                  ncommand+=std::strlen(tmpstr) + 4;
                  if (!arguments[ind]) CImg<char>(tmpstr,std::strlen(tmpstr)+1).move_to(arguments[ind]);
                  ncommand+=std::sprintf(tmpstr,"%d",ind);
                  replace_text = arguments[ind].data();
                  has_arguments = true;

                  // Substitute any other expression starting by '$'.
                } else {
                  s_argument[0] = '$';
                  if (std::sscanf(ncommand,"%4095[^$]",s_argument+1)!=1) { s_argument[1] = 0; ++ncommand; }
                  else ncommand+=std::strlen(s_argument);
                }

                const int replace_length = std::strlen(replace_text);
                if (replace_length)
                  CImg<char>(replace_text,replace_length).move_to(lreplacement);

              } else {
                std::sscanf(ncommand,"%4095[^$]",s_argument);
                const int replace_length = std::strlen(s_argument);
                if (replace_length) {
                  CImg<char>(s_argument,replace_length).move_to(lreplacement);
                  ncommand+=std::strlen(s_argument);
                }
              }
              const CImg<char> zero(1,1,1,1,0);
              (lreplacement.insert(zero)>'x').move_to(substituted_command);
              debug(images,"Expand command '%s' to '%s'.",custom_command_name,substituted_command.data());
              break;
            }
          }

          if (custom_command_found) {
            const CImgList<char> ncommand_line = command_line_to_CImgList(substituted_command.data());
            CImgList<unsigned int> ndowhiles, nrepeatdones, nlocals;
            const unsigned int siz = selection.size();
            CImgList<char> nfilenames(siz);
            CImgList<T> nimages(siz);
            unsigned int nposition = 0;

            if (get_version) {
              cimg_forY(selection,l) { nimages[l].assign(images[selection[l]]); nfilenames[l].assign(filenames[selection[l]]); }
              CImg<char>(custom_command_name,std::strlen(custom_command_name)+1).move_to(scope);
              parse(ncommand_line,nposition,nimages,nfilenames,ndowhiles,nrepeatdones,nlocals,false);
              scope.remove();
              nimages.move_to(images,~0U); nfilenames.move_to(filenames,~0U);
            } else {
              cimg_forY(selection,l) { nimages[l].swap(images[selection[l]]); nfilenames[l].swap(filenames[selection[l]]); }
              CImg<char>(custom_command_name,std::strlen(custom_command_name)+1).move_to(scope);
              parse(ncommand_line,nposition,nimages,nfilenames,ndowhiles,nrepeatdones,nlocals,false);
              scope.remove();
              const unsigned int nb = cimg::min(siz,nimages.size());
              for (unsigned int i = 0; i<nb; ++i) {
                images[selection[i]].swap(nimages[0]); filenames[selection[i]].swap(nfilenames[0]);
                nimages.remove(0); nfilenames.remove(0);
              }
              if (nb<siz) for (unsigned int off = 0, l = nb; l<siz; ++l, ++off) {
                  const unsigned int ind = selection[l] - off;
                  images.remove(ind); filenames.remove(ind);
                } else if (nimages) {
                const unsigned int ind0 = siz?selection[siz-1]+1:images.size();
                filenames.insert(nimages.size(),CImg<char>("(unnamed)",10),ind0);
                nimages.move_to(images,ind0);
              }
            }
            if (has_arguments) ++position;
            continue;
          }
        }
      }

      // Input.
      if (!std::strcmp("-input",command_name) || !std::strcmp("-i",command_name)) ++position;
      else {
        if (get_version) --item;
        argument = item;
        if (std::strlen(argument)>=64) {
          std::memcpy(argument_text,argument,60*sizeof(char));
          argument_text[60] = argument_text[61] = argument_text[62] = '.'; argument_text[63] = 0;
        } else std::strcpy(argument_text,argument);
        command_restriction[0] = 0;
      }
      if (!is_restriction || !selection) selection.assign(1,1,1,1,images.size());
      CImgList<T> input_images;
      CImgList<char> input_filenames;
      bool obj3d = false;
      char st_inds[4096] = { 0 }, stx[4096] = { 0 }, sty[4096] = { 0 }, stz[4096] = { 0 }, stv[4096] = { 0 }, st_values[4096] = { 0 };
      char sep = 0, sepx = 0, sepy = 0, sepz = 0, sepv = 0;
      int nb = 1, indx = no_ind, indy = no_ind, indz = no_ind, indv = no_ind;
      float dx = 0, dy = 1, dz = 1, dv = 1;

      if ((std::sscanf(argument,"[%4095[0-9%,:-]%c%c",st_inds,&sep,&end)==2 && sep==']') ||
          std::sscanf(argument,"[%4095[0-9%,:-]]x%d%c",st_inds,&nb,&end)==2) { // Nb copies of existing images.

        const CImg<unsigned int> indices = selection2cimg(st_inds,images.size(),"-input",false);
        char st_tmp[4096] = { 0 }; std::strcpy(st_tmp,selection2string(indices,filenames,true));
        if (nb<=0) arg_error("input");
        if (nb!=1) print(images,"Input %d copies of image%s at position%s",nb,st_tmp,gmic_selection);
        else print(images,"Input copy of image%s at position%s",st_tmp,gmic_selection);
        for (int i = 0; i<nb; ++i) cimg_foroff(indices,l) {
          input_images.insert(images[indices[l]]);
          input_filenames.insert(filenames[indices[l]]);
        }
      } else if ((std::sscanf(argument,"%4095[][0-9.eE%+-]%c",stx,&end)==1 ||
                  std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",stx,sty,&end)==2 ||
                  std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",stx,sty,stz,&end)==3 ||
                  std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-]%c",
                              stx,sty,stz,stv,&end)==4 ||
                  std::sscanf(argument,"%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[][0-9.eE%+-],%4095[^\n]",
                              stx,sty,stz,stv,&(st_values[0]=0))==5) &&
                 (!*stx || std::sscanf(stx,"%f%c",&dx,&end)==1 ||
                  (std::sscanf(stx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%') ||
                  (std::sscanf(stx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']')) &&
                 (!*sty || std::sscanf(sty,"%f%c",&dy,&end)==1 ||
                  (std::sscanf(sty,"%f%c%c",&dy,&sepy,&end)==2 && sepy=='%') ||
                  (std::sscanf(sty,"[%d%c%c",&indy,&sepy,&end)==2 && sepy==']')) &&
                 (!*stz || std::sscanf(stz,"%f%c",&dz,&end)==1 ||
                  (std::sscanf(stz,"%f%c%c",&dz,&sepz,&end)==2 && sepz=='%') ||
                  (std::sscanf(stz,"[%d%c%c",&indz,&sepz,&end)==2 && sepz==']')) &&
                 (!*stv || std::sscanf(stv,"%f%c",&dv,&end)==1 ||
                  (std::sscanf(stv,"%f%c%c",&dv,&sepv,&end)==2 && sepv=='%') ||
                  (std::sscanf(stv,"[%d%c%c",&indv,&sepv,&end)==2 && sepv==']'))) { // New image with specified dimensions and values.
        if (indx!=no_ind) { gmic_check_indice(indx); dx = (float)images[indx].width(); sepx = 0; }
        if (indy!=no_ind) { gmic_check_indice(indy); dy = (float)images[indy].height(); sepy = 0; }
        if (indz!=no_ind) { gmic_check_indice(indz); dz = (float)images[indz].depth(); sepz = 0; }
        if (indv!=no_ind) { gmic_check_indice(indv); dv = (float)images[indv].spectrum(); sepv = 0; }
        int idx = (int)dx, idy = (int)dy, idz = (int)dz, idv = (int)dv;
        if (sepx=='%') { idx = (int)cimg::round(dx*last_image.width()/100,1); if (!idx) ++idx; }
        if (sepy=='%') { idy = (int)cimg::round(dy*last_image.height()/100,1); if (!idy) ++idy; }
        if (sepz=='%') { idz = (int)cimg::round(dz*last_image.depth()/100,1); if (!idz) ++idz; }
        if (sepv=='%') { idv = (int)cimg::round(dv*last_image.spectrum()/100,1); if (!idv) ++idv; }
        if (idx<=0 || idy<=0 || idz<=0 || idv<=0) arg_error("input");
        if (*st_values) {
          print(images,"Input image at position%s, with values '%s'",gmic_selection,st_values);
        } else print(images,"Input black image at position%s",gmic_selection);
        CImg<T> new_image(idx,idy,idz,idv,0);
        if (*st_values) new_image.fill(st_values,true);
        new_image.move_to(input_images);
        filenames.insert(input_images.size(),CImg<char>("(unnamed)",10));
      } else if (*argument=='(' && argument[std::strlen(argument)-1]==')' &&
                 std::sscanf(argument,"(%4095[^\n]",stx)==1) { // New IxJxKxL image specified as array.
        stx[std::strlen(stx)-1] = 0;
        unsigned int cx = 0, cy = 0, cz = 0, cv = 0, maxcx = 0, maxcy = 0, maxcz = 0;
        const char *nargument = 0;
        for (nargument = stx; *nargument; ) {
          char s_value[256] = { 0 }, separator = 0; double value = 0;
          if (std::sscanf(nargument,"%255[0-9.eE+-]%c",s_value,&separator)>0 &&
              std::sscanf(s_value,"%lf%c",&value,&end)==1) {
            if (cx>maxcx) maxcx = cx;
            if (cy>maxcy) maxcy = cy;
            if (cz>maxcz) maxcz = cz;
            switch (separator) {
            case '^' : cx = cy = cz = 0; ++cv; break;
            case '/' : cx = cy = 0; ++cz; break;
            case ';' : cx = 0; ++cy; break;
            default : ++cx;
            }
            nargument+=std::strlen(s_value) + (separator?1:0);
          } else break;
        }
        if (*nargument) arg_error("input");

        CImg<T> img(maxcx+1,maxcy+1,maxcz+1,cv+1,0);
        cx = cy = cz = cv = 0;
        for (nargument = stx; *nargument; ) {
          char s_value[256] = { 0 }, separator = 0; double value = 0;
          if (std::sscanf(nargument,"%255[0-9.eE+-]%c",s_value,&separator)>0 &&
              std::sscanf(s_value,"%lf%c",&value,&end)==1) {
            img(cx,cy,cz,cv) = (T)value;
            switch (separator) {
            case '^' : cx = cy = cz = 0; ++cv; break;
            case '/' : cx = cy = 0; ++cz; break;
            case ';' : cx = 0; ++cy; break;
            default : ++cx;
            }
            nargument+=std::strlen(s_value) + (separator?1:0);
          } else break;
        }
        print(images,"Input image at position%s, with values '%s'",gmic_selection,argument_text);
        img.move_to(input_images); filenames.insert(CImg<char>("(unnamed)",10));
      } else { // Filename.
        char filename[4096] = { 0 }, options[4096] = { 0 };
        if (argument[0]!='-' || (argument[1] && argument[1]!='.')) {
          std::strcpy(filename,argument);
          std::FILE *file = std::fopen(filename,"r");
          if (file) std::fclose(file);
          else {
            std::sscanf(argument,"%4095[^,],%s",filename,options);
            if (!(file=std::fopen(filename,"r"))) {
              if (filename[0]=='-') error(images,"Command '%s' : Command not found.",filename+1);
              else error(images,"Command 'input' : File '%s' not found.",argument_text);
            }
            std::fclose(file);
          }
        } else std::strcpy(filename,argument);

        const char *const basename = cimg::basename(filename), *const ext = cimg::split_filename(filename);
        if (!cimg::strcasecmp("off",ext)) {

          // 3D object file.
          print(images,"Input 3D object '%s' at position%s",is_fullpath?filename:basename,gmic_selection);
          CImgList<unsigned int> primitives;
          CImgList<unsigned char> colors;
          CImg<float> opacities, vertices = CImg<float>::get_load_off(filename,primitives,colors);
          opacities.assign(1,primitives.size(),1,1,1);
          vertices.object3dtoCImg3d(primitives,colors,opacities);
          vertices.move_to(input_images);
          CImg<char>(is_fullpath?filename:basename,std::strlen(is_fullpath?filename:basename)+1).move_to(input_filenames);
          obj3d = true;
        } else if (!cimg::strcasecmp(ext,"avi") ||
                   !cimg::strcasecmp(ext,"mov") ||
                   !cimg::strcasecmp(ext,"asf") ||
                   !cimg::strcasecmp(ext,"divx") ||
                   !cimg::strcasecmp(ext,"flv") ||
                   !cimg::strcasecmp(ext,"mpg") ||
                   !cimg::strcasecmp(ext,"m1v") ||
                   !cimg::strcasecmp(ext,"m2v") ||
                   !cimg::strcasecmp(ext,"m4v") ||
                   !cimg::strcasecmp(ext,"mjp") ||
                   !cimg::strcasecmp(ext,"mkv") ||
                   !cimg::strcasecmp(ext,"mpe") ||
                   !cimg::strcasecmp(ext,"movie") ||
                   !cimg::strcasecmp(ext,"ogm") ||
                   !cimg::strcasecmp(ext,"qt") ||
                   !cimg::strcasecmp(ext,"rm") ||
                   !cimg::strcasecmp(ext,"vob") ||
                   !cimg::strcasecmp(ext,"wmv") ||
                   !cimg::strcasecmp(ext,"xvid") ||
                   !cimg::strcasecmp(ext,"mpeg")) {

          // Image sequence file.
          unsigned int value0 = 0, value1 = 0, step = 1; char sep0 = 0, sep1 = 0;
          if ((std::sscanf(options,"%u%c,%u%c,%u%c",&value0,&sep0,&value1,&sep1,&step,&end)==5 && sep0=='%' && sep1=='%') ||
              (std::sscanf(options,"%u%c,%u,%u%c",&value0,&sep0,&value1,&step,&end)==4 && sep0=='%') ||
              (std::sscanf(options,"%u,%u%c,%u%c",&value0,&value1,&sep1,&step,&end)==4 && sep1=='%') ||
              (std::sscanf(options,"%u,%u,%u%c",&value0,&value1,&step,&end)==3) ||
              (std::sscanf(options,"%u%c,%u%c%c",&value0,&sep0,&value1,&sep1,&end)==4 && sep0=='%' && sep1=='%') ||
              (std::sscanf(options,"%u%c,%u%c",&value0,&sep0,&value1,&end)==3 && sep0=='%') ||
              (std::sscanf(options,"%u,%u%c%c",&value0,&value1,&sep1,&end)==3 && sep1=='%') ||
              (std::sscanf(options,"%u,%u%c",&value0,&value1,&end)==2)) { // Read several frames
            print(images,"Input frames %u%s..%u%s with step %u of file '%s' at position%s",
                  value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",step,is_fullpath?filename:basename,gmic_selection);
            if (sep0=='%' || sep1=='%') {
              const unsigned int nb_frames = CImg<unsigned int>::get_load_ffmpeg(filename,0,0,0)[0];
              if (sep0=='%') value0 = (unsigned int)cimg::round(value0*nb_frames/100,1);
              if (sep1=='%') value1 = (unsigned int)cimg::round(value1*nb_frames/100,1);
            }
          } else if ((std::sscanf(options,"%u%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
                     (std::sscanf(options,"%u%c",&value0,&end)==1)) { // Read one frame
            print(images,"Input frame %u%s of file '%s' at position%s",value0,sep0=='%'?"%":"",is_fullpath?filename:basename,gmic_selection);
            if (sep0=='%') {
              const unsigned int nb_frames = CImg<unsigned int>::get_load_ffmpeg(filename,0,0,0)[0];
              value0 = (unsigned int)cimg::round(value0*nb_frames/100,1);
            }
            value1 = value0; step = 1;
          } else { // Read all frames
            print(images,"Input all frames of file '%s' at position%s",is_fullpath?filename:basename,gmic_selection);
            value0 = 0; value1 = ~0U; sep0 = sep1 = 0; step = 1;
          }
          input_images.load_ffmpeg(filename,value0,value1,step);
          if (input_images)
            input_filenames.insert(input_images.size(),
                                   CImg<char>(is_fullpath?filename:basename,std::strlen(is_fullpath?filename:basename)+1));
        } else if (!cimg::strcasecmp("raw",ext)) {

          // Raw file.
          int dx = 0, dy = 1, dz = 1, dv = 1;
          if (std::sscanf(options,"%d,%d,%d,%d",&dx,&dy,&dz,&dv)>0) {
            if (dx<=0 || dy<=0 || dz<=0 || dv<=0)
              error(images,"Command 'input' : RAW file '%s', invalid specified dimensions %dx%dx%dx%d.",filename,dx,dy,dz,dv);
            print(images,"Input RAW file '%s' at position%s",is_fullpath?filename:basename,gmic_selection);
            CImg<T>::get_load_raw(filename,dx,dy,dz,dv).move_to(input_images);
            input_filenames.insert(CImg<char>(is_fullpath?filename:basename,std::strlen(is_fullpath?filename:basename)+1));
          } else error(images,"Command 'input' : RAW file '%s', image dimensions must be specified as a file option.",filename);
        } else if (!cimg::strcasecmp("yuv",ext)) {

          // YUV file.
          int dx = 0, dy = 0; unsigned int first = 0, last = ~0U, step = 1;
          if (std::sscanf(options,"%d,%d,%u,%u,%u",&dx,&dy,&first,&last,&step)>0) {
            if (dx<=0 || dy<=0)
              error(images,"Command 'input' : YUV file '%s', invalid specified dimensions %dx%d.",filename,dx,dy);
            print(images,"Input YUV file '%s' at position%s",is_fullpath?filename:basename,gmic_selection);
            input_images.load_yuv(filename,dx,dy,first,last,step);
            input_filenames.insert(input_images.size(),
                                   CImg<char>(is_fullpath?filename:basename,std::strlen(is_fullpath?filename:basename)+1));
          } else error(images,"Command 'input' : YUV file '%s', image dimensions must be specified as a file option.",filename);
        } else if (!cimg::strcasecmp("gmic",ext)) {

          // G'MIC custom command file
          print(images,"Load command file '%s'",is_fullpath?filename:basename);
          const unsigned int siz = command_names.size();
          std::FILE *const file = cimg::fopen(argument,"r");
          add_commands(file);
          cimg::fclose(file);
          if (verbosity_level>=0) {
            std::fprintf(cimg::output()," (%u commands added).",command_names.size()-siz);
            std::fflush(cimg::output());
          }
          continue;
        } else {

          // Other file type.
          print(images,"Input file '%s' at position%s",is_fullpath?filename:basename,gmic_selection);
          input_images.load(filename);
          input_filenames.insert(input_images.size(),CImg<char>(is_fullpath?filename:basename,std::strlen(is_fullpath?filename:basename)+1));
        }
      }

      if (verbosity_level>=0) {
        if (input_images) {
          const unsigned int last = input_images.size()-1;
          if (obj3d) {
            std::fprintf(cimg::output()," (%d vertices, %u primitives).",
                         (unsigned int)input_images(0,6),
                         (unsigned int)input_images(0,7));
            std::fflush(cimg::output());
          } else if (input_images.size()==1) {
            std::fprintf(cimg::output()," (1 image %dx%dx%dx%d).",
                         input_images[0].width(),input_images[0].height(),input_images[0].depth(),input_images[0].spectrum());
            std::fflush(cimg::output());
          } else {
            std::fprintf(cimg::output()," (%u images [0] = %dx%dx%dx%d, %s[%u] = %dx%dx%dx%d).",
                         input_images.size(),
                         input_images[0].width(),input_images[0].height(),input_images[0].depth(),input_images[0].spectrum(),
                         last==1?"":"..,",last,
                         input_images[last].width(),input_images[last].height(),input_images[last].depth(),input_images[last].spectrum());
            std::fflush(cimg::output());
          }
        } else {
          std::fprintf(cimg::output()," (no available data).");
          std::fflush(cimg::output());
        }
      }

      for (unsigned int l = 0, siz = selection.size()-1, off = 0; l<=siz; ++l) {
        const unsigned int ind = selection[l] + off;
        off+=input_images.size();
        filenames.insert(input_filenames,ind);
        if (l!=siz) images.insert(input_images,ind);
        else input_images.move_to(images,ind);
      }

    }

    // Post-checking.
    if (filenames.size()!=images.size())
      error("Internal error : Images (%u) and filenames (%u) are in an inconsistent state, when returning.",filenames.size(),images.size());
    if (!scope.size())
      error("Internal error : Scope is empty, when returning.");
    if (dowhiles)
      warning(images,"A 'while' command is probably missing.");
    if (repeatdones)
      warning(images,"A 'done' command is probably missing.");
    if (locals)
      error(images,"A 'endlocal' directive is probably missing.");
    if (initial_call && !is_end && stack) {
      warning(images,"A 'pop' directive is probably missing (global stack contains %u element%s at exit).",
              stack.size(),stack.size()>1?"s":"");
      stack.assign();
    }

    // Display final result if necessary (not 'released' before).
    if (initial_call && !is_end) {
      if (images.size() && !is_released) {
        if (!display_objects3d(images,filenames,CImg<unsigned int>::sequence(images.size(),0,images.size()-1),false))
          display_images(images,filenames,CImg<unsigned int>::sequence(images.size(),0,images.size()-1),true);
      }
      print(images,"End G'MIC instance.\n");
      is_end = true;
    }

  } catch (CImgException &e) {
    const char *error_message = e.message();
    char tmp[4096] = { 0 }, sep = 0;
    if (std::sscanf(error_message,"%4095[^>]>:%c",tmp,&sep)==2 && sep==':') error_message+=std::strlen(tmp)+3;
    error(images,error_message);
  }
  return *this;
}

// Small hack to separate the compilation of G'MIC in different pixel types.
// (only intended to save computer memory when compiling !)
//--------------------------------------------------------------------------
#ifdef gmic_minimal
gmic& gmic::parse_float(const CImgList<char>& command_line, unsigned int& position,CImgList<float>& images,
                        CImgList<char>& filenames,
                        CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                        CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<float>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<float>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#else
#if defined(gmic_bool) || !defined(gmic_separate_compilation)
gmic& gmic::parse_bool(const CImgList<char>& command_line, unsigned int& position, CImgList<bool>& images,
                       CImgList<char>& filenames,
                       CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                       CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<bool>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<bool>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#endif
#if defined(gmic_uchar) || !defined(gmic_separate_compilation)
gmic& gmic::parse_uchar(const CImgList<char>& command_line, unsigned int& position, CImgList<unsigned char>& images,
                        CImgList<char>& filenames,
                        CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                        CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<unsigned char>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<unsigned char>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#endif
#if defined(gmic_char) || !defined(gmic_separate_compilation)
gmic& gmic::parse_char(const CImgList<char>& command_line, unsigned int& position, CImgList<char>& images,
                       CImgList<char>& filenames,
                       CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                       CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<char>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<char>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#endif
#if defined(gmic_ushort) || !defined(gmic_separate_compilation)
gmic& gmic::parse_ushort(const CImgList<char>& command_line, unsigned int& position, CImgList<unsigned short>& images,
                         CImgList<char>& filenames,
                         CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                         CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<unsigned short>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<unsigned short>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#endif
#if defined(gmic_short) || !defined(gmic_separate_compilation)
gmic& gmic::parse_short(const CImgList<char>& command_line, unsigned int& position, CImgList<short>& images,
                        CImgList<char>& filenames,
                        CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                        CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<short>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<short>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#endif
#if defined(gmic_uint) || !defined(gmic_separate_compilation)
gmic& gmic::parse_uint(const CImgList<char>& command_line, unsigned int& position, CImgList<unsigned int>& images,
                       CImgList<char>& filenames,
                       CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                       CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<unsigned int>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<unsigned int>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#endif
#if defined(gmic_int) || !defined(gmic_separate_compilation)
gmic& gmic::parse_int(const CImgList<char>& command_line, unsigned int& position, CImgList<int>& images,
                      CImgList<char>& filenames,
                      CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                      CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<int>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<int>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#endif
#if defined(gmic_float) || !defined(gmic_separate_compilation)
gmic& gmic::parse_float(const CImgList<char>& command_line, unsigned int& position, CImgList<float>& images,
                        CImgList<char>& filenames,
                        CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                        CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<float>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<float>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#endif
#if defined(gmic_double) || !defined(gmic_separate_compilation)
gmic& gmic::parse_double(const CImgList<char>& command_line, unsigned int& position, CImgList<double>& images,
                         CImgList<char>& filenames,
                         CImgList<unsigned int>& dowhiles, CImgList<unsigned int>& repeatdones,
                         CImgList<unsigned int>& locals, const bool initial_call) {
  return parse(command_line,position,images,filenames,dowhiles,repeatdones,locals,initial_call);
}
template gmic::gmic(const int, const char *const *const, CImgList<double>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
template gmic::gmic(const char *const, CImgList<double>&,
                    const char *const custom_commands, const bool default_commands, float *const p_progress, int *const p_cancel);
#endif
#endif
#endif

//------------------------
// Documentation and help
//------------------------
#if defined(gmic_main) || (!defined(gmic_separate_compilation) && !defined(gmic_minimal))
extern char data_gmic_def[];

#define gmic_usage(usage) \
  if (!command_name) cimg_usage(usage); else { \
    cimg::option(0,0,0,0,0,true); if (display_usage) cimg_usage(usage); \
  }
#define gmic_section(str) if (!command_name) cimg_help(str); else is_command_name = false;
#define gmic_help(str) if (!command_name || is_command_name) cimg_help(str);
#define gmic_option(name,defaut,usage) \
  if (!command_name) cimg_option(name,defaut,usage); \
  else { \
    if (is_help_displayed) return true; \
    is_command_name = !std::strcmp(command_name,name); \
    if (is_command_name) { cimg_option(name,defaut,usage); is_help_displayed = true; } \
    }
#define ___ "                     "
#define _   "        "

bool help(const int argc, const char *const *const argv, const char *const command_name=0, const bool display_usage=true) {
  bool is_command_name = true, is_help_displayed;
  gmic_usage("GREYC's Magic Image Converter");
  if (display_usage) {
    char version[1024] = { 0 };
    std::sprintf(version,"        Version %d.%d.%d.%d, Copyright (C) 2008-2009, David Tschumperle (http://gmic.sourceforge.net)\n",
                 gmic_version/1000,(gmic_version/100)%10,(gmic_version/10)%10,gmic_version%10);
    gmic_help(version);
  }
  is_help_displayed = is_command_name = false;

  gmic_section(" Usage\n"
               " -----\n");

  gmic_help(" gmic [file1 | instruction1 [arg1_1,arg1_2,..]] .. [fileN | instructionN [argN_1,argN_2,..]]\n");
  gmic_help(" 'gmic' is an open-source interpreter of the G'MIC language, a script-based");
  gmic_help(" programming language dedicated to the design of complex image processing pipelines.");
  gmic_help(" It can be typically used to convert, manipulate, and visualize datasets composed of");
  gmic_help(" one or several 1D/2D/3D multi-spectral images.\n");
  gmic_help(" The G'MIC language is small and quite easy to learn. It follows these simple rules :\n");
  gmic_help("   - At any time, one considers a numbered list of images which are all stored in");
  gmic_help("     computer memory.");
  gmic_help("   - Pixels of all these images have the same datatype which is one of the following types:");
  gmic_help("     { bool | uchar | char | ushort | short | uint | int | float | double }.");
  gmic_help("   - The first image of the list has indice '0' and is denoted by [0].");
  gmic_help("   - Negative indices are treated in a cyclic way (i.e. image [-1] stands for the");
  gmic_help("     last image, [-2] the penultimate one, and so on..).");
  gmic_help("   - An image processing pipeline is described as a sequence of G'MIC items,");
  gmic_help("     separated by spaces ' ', read and interpreted from the left to the right.");
  gmic_help("   - Items can be either commands, command arguments, filenames or input strings.");
  gmic_help("   - On the command line, any string following 'gmic' is considered as a G'MIC item.");
  gmic_help("   - An item starting by '-' usually designates a G'MIC command.");
  gmic_help("   - Some command may have two equivalent names (regular and short, for instance");
  gmic_help("     command items '-resize' and '-r' are strictly equivalent).");
  gmic_help("   - A command may have mandatory or optional arguments.");
  gmic_help("   - The command argument is the item directly following the command name.");
  gmic_help("   - When multiple arguments are required, they are separated by commas ','.");
  gmic_help("   - When an input filename or an input string is encountered, the corresponding");
  gmic_help("     image data are loaded/created and inserted at the end of the image list");
  gmic_help("     (this is actually similar to giving an argument to the command '-input' or '-i').");
  gmic_help("   - Filenames '-' and '-.ext' stand for the standard input/output streams");
  gmic_help("     (optionally, with data forced to be in the specified 'ext' file format).");
  gmic_help("   - Input strings can be used to insert new images \"from scratch\" at the end of the list :");
  gmic_help("       _ 'width[%],_height[%],_depth[%],_spectrum[%],_values' : ");
  gmic_help("         Insert a new image with specified dimensions and values (adding '%' to a");
  gmic_help("         dimension means 'percentage of the same dimension, get from the last");
  gmic_help("         available image'). A dimension 'width','height','depth' or 'spectrum' can be");
  gmic_help("         also written as '[indice]' in which case its value is taken from the same");
  gmic_help("         dimension of the specified existing image [indice].");
  gmic_help("       _ '[indice]' or '[indice]xN' : Insert 1 or N copies of the image [indice].");
  gmic_help("       _ '(v1,v2,..)' : Create a new image containing the specified values.");
  gmic_help("         Value separators inside parentheses can be ',' (column), ';' (line), ");
  gmic_help("         '/' (slice) or '^' (channel).");
  gmic_help("   - The execution of a command may be restricted to a sub-selection of the image");
  gmic_help("     list, by appending '[selection]' to the command name.");
  gmic_help("     Selections can be defined in many different ways. For instance : ");
  gmic_help("       _ '-command[0,1,3]' : Apply command only on images [0],[1] and [3].");
  gmic_help("       _ '-command[3-5]' : Apply command only on images [3] to [5]");
  gmic_help("                           (i.e. on images [3], [4] and [5]).");
  gmic_help("       _ '-command[50%-100%]' : Apply command on the second half of the image list.");
  gmic_help("       _ '-command[0,-4--1]' : Apply command on the first and the four latest images.");
  gmic_help("       _ '-command[0-9:3]' : Apply command only on images 0 to 9, with a step of 3");
  gmic_help("                             (i.e. on images [0], [3], [6] and [9]).");
  gmic_help("       _ '-command[0,2-4,50%--1]' : Apply command on images [0],[2],[3],[4] and on");
  gmic_help("                                    the second half of the image list.");
  gmic_help("   - When no selection is specified, a command is applied by default on all images.");
  gmic_help("   - A command starting with '--' instead of '-' does not act 'in-place' but inserts");
  gmic_help("     its result as one or several new image(s), at the end of the list.");
  gmic_help("   - Any expression starting with '@' in an item is substituted before item interpretation :");
  gmic_help("       _ '@#' is substituted by the current number of images in the list.");
  gmic_help("       _ '@*' is substituted by the current number of items in the global stack.");
  gmic_help("       _ '@{*}' or '@{*,subset}' are substituted by the stack content, or a subset of it.");
  gmic_help("       _ '@>' and '@<' are substituted by the current number of running 'repeat-done' loops.");
  gmic_help("       _ '@{>}' or '@{>,subset}' are substituted by the indices (or a subset of them) of");
  gmic_help("         the running 'repeat-done' loops, given in ascending order, from 0 to N-1.");
  gmic_help("       _ '@{<}' or '@{<,subset}' do the same in descending order, from N-1 to 0.");
  gmic_help("       _ '@!' is substituted by the visibility state of the instant display window [0]");
  gmic_help("         (can be { 0=closed | 1=visible }).");
  gmic_help("       _ '@{!,feature}' or '@{!indice,feature}' is substituted by a specific feature of the");
  gmic_help("         instant display window [0] (or [indice], if specified). The retrieved 'feature' can be");
  gmic_help("         one of the followings :");
  gmic_help("            . 'w' : get display width.");
  gmic_help("            . 'h' : get display height.");
  gmic_help("            . 'u' : get screen width (do not depend on the instant window).");
  gmic_help("            .' v' : get screen height (do not depend on the instant window).");
  gmic_help("            . 'x' : get X-coordinate of the mouse position.");
  gmic_help("            . 'y' : get Y-coordinate of the mouse position.");
  gmic_help("            . 'n' : get type of normalization.");
  gmic_help("            . 'b' : get state of the mouse buttons.");
  gmic_help("            . 'o' : get state of the mouse wheel.");
  gmic_help("            . 'c' : get boolean telling if the instant display has been closed.");
  gmic_help("            . 'r' : get boolean telling if the instant display has been resized.");
  gmic_help("            . 'm' : get boolean telling if the instant display window has been moved.");
  gmic_help("            . Any other feature describe a key name whose state { 0=pressed | 1=released } is returned.");
  gmic_help("       _ '@indice' or '@{indice,feature}' is substituted by the pixel values of the image [indice],");
  gmic_help("         or by a specific feature (or subset) of this image.");
  gmic_help("         The retrieved 'feature' can be one of the followings :");
  gmic_help("            . 'w' : get image width (number of columns).");
  gmic_help("            . 'h' : get image height (number of lines).");
  gmic_help("            . 'd' : get image depth (number of slices).");
  gmic_help("            . 's' : get image spectrum (number of channels).");
  gmic_help("            . '#' : get number of image values (width x height x depth x spectrum).");
  gmic_help("            . '+' : get sum of all pixel values.");
  gmic_help("            . '-' : get difference of all pixel values.");
  gmic_help("            . '*' : get product of all pixel values.");
  gmic_help("            . '/' : get quotient of all pixel values.");
  gmic_help("            . 'm' : get minimum pixel value.");
  gmic_help("            . 'M' : get maximum pixel value.");
  gmic_help("            . 'a' : get average pixel value.");
  gmic_help("            . 'v' : get variance of pixel values.");
  gmic_help("            . 't' : get text string from the image values, regarded as ASCII codes.");
  gmic_help("            . 'c' : get (x,y,z,c) coordinates of the minimum value.");
  gmic_help("            . 'C' : get (x,y,z,c) coordinates of the maximum value.");
  gmic_help("            . '(x,_y,_z,_c,_borders)' : get pixel value at coordinates (x,y,z,c).");
  gmic_help("         Any other 'feature' is considered as a desired subset of image values, for")
  gmic_help("         instance, expression '@{-1,0-50%}' is substituted by all values (separated");
  gmic_help("         by commas ',') coming from the first half of the last image.");
  gmic_help("   - Any other expression involving braces (as '{expression}') is considered as a mathematical");
  gmic_help("     expression and is evaluated. If the string is not evaluable (invalid expression), it is replaced by");
  gmic_help("     the sequence of ASCII codes that compose the string, separated by commas ','.");
  gmic_help("   - \"Items\" delimited by double quotes '\"' may contain spaces, commas or ESC sequences.");
  gmic_help("   - Some G'MIC commands may result to the generation of 3D objects.");
  gmic_help("   - In G'MIC, a 3D object is stored as a one-column image, containing all object data, in the");
  gmic_help("     following order : [ header, vertices, faces, colors, opacities ].");
  gmic_help("   - Custom G'MIC commands can be defined by the user, through a command file.");
  gmic_help("   - A command file is a simple ASCII text file, where each line starts either by");
  gmic_help("     'instruction_name : substitution' or 'substitution (continuation)' or '# comment'.");
  gmic_help("   - A default command file 'gmic_def.raw' is distributed within the G'MIC package.");
  gmic_help("     Looking at it is a good start to learn how to create your own custom commands.");
  gmic_help("   - Commands defined in file 'gmic_def.raw' are already included by default in the");
  gmic_help("     interpreter. Their explicit inclusion (using command '-m') is useless.");
  gmic_help("   - In custom commands, expressions starting with '$' are substituted this way :");
  gmic_help("       _ '$#' is substituted by the number of specified arguments.");
  gmic_help("       _ '$*' is substituted by all specified arguments, separated by commas ','.");
  gmic_help("       _ '$i' and '${i}' are substituted by the i-th specified argument");
  gmic_help("          ('i' starts from '1' to '$#').");
  gmic_help("       _ '${i*}' is substituted by all arguments whose indice is higher or equal to i.");
  gmic_help("       _ '${i=default}' is substituted by the value of $i (if defined) or by its new");
  gmic_help("         default value 'default' else ('default' can be a $-expression as well).");
  gmic_help("       _ '$?' is substituted by a string telling about the image selection");
  gmic_help("         (should be used in command descriptions only).\n");
  gmic_help(" All currently recognized G'MIC commands are listed below.");
  gmic_help(" Possible formats for required command arguments (if any) are separated by '|'.");
  gmic_help(" An argument specified in '[]' or starting by '_' is optional except when standing for");
  gmic_help(" an existing image [indice] of the current image list. In this case, the characters");
  gmic_help(" '[' and ']' are mandatory when writting the item.\n");

  gmic_section(" Global options\n"
               " --------------\n");

  gmic_option("-help","_command","");
  gmic_help(_"Display help (optionally for specified command only) and quit.");
  gmic_help(_"(eq. to '-h').\n");

  gmic_option("-verbose","level","");
  gmic_help(_"Set verbosity level to 'level'.");
  gmic_help(_"(eq. to '-v').\n");
  gmic_help(_"When 'level>=0', execution messages are send to the standard output stream.");
  gmic_help(_"Default value for 'level' is '0'.\n");

  gmic_option("-verbose+","","");
  gmic_help(_"Increment verbosity level.");
  gmic_help(_"(eq. to '-v+').\n");

  gmic_option("-verbose-","","");
  gmic_help(_"Decrement verbosity level.");
  gmic_help(_"(eq. to '-v-').\n");

  gmic_option("-command","filename |","");
  gmic_help(___"\"string\"");
  gmic_help(_"Import G'MIC command(s) from specified file or string.");
  gmic_help(_"(eq. to '-m').\n");
  gmic_help(_"Custom commands can be used directly after the command execution.\n");

  gmic_option("-debug","","");
  gmic_help(_"Switch debug flag.\n");
  gmic_help(_"When activated, the debug mode outputs additionnal log message describing the");
  gmic_help(_"internal state of the interpreter. Should be used by developers only.\n");

  gmic_option("-fullpath","","");
  gmic_help(_"Switch full path flag.\n");
  gmic_help(_"When activated, the displayed image names contains the full path filename");
  gmic_help(_"including the location folders.\n");

  gmic_section(" Arithmetic operators\n"
               " --------------------\n");

  gmic_option("-add","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Add value 'value', image [indice] or 'filename', mathematical expression 'formula'");
  gmic_help(_"to selected images, or add all selected images together.");
  gmic_help(_"(eq. to '-+').\n");

  gmic_option("-sub","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Substract value 'value', image [indice] or 'filename', mathematical expression 'formula'");
  gmic_help(_"to selected images, or substract all selected images together.");
  gmic_help(_"(eq. to '--').\n");

  gmic_option("-mul","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Multiply selected images by value 'value', image [indice] or 'filename', mathematical");
  gmic_help(_"expression 'formula', or multiply all selected images together.");
  gmic_help(_"(eq. to '-*').\n");

  gmic_option("-div","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Divide selected image by value 'value', image [indice] or 'filename', mathematical");
  gmic_help(_"expression 'formula', or divide all selected images together.");
  gmic_help(_"(eq. to '-/').\n");

  gmic_option("-pow","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Compute selected image to the power of value 'value', image [indice] or 'filename',");
  gmic_help(_"mathematical expression 'formula', or compute power of all selected images together.");
  gmic_help(_"(eq. to '-^').\n");

  gmic_option("-min","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Compute minimum between selected images and value 'value', image [indice] or 'filename',");
  gmic_help(_"mathematical expression 'formula', or compute minimum of all selected images together.\n");

  gmic_option("-max","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Compute maximum between selected images and value 'value', image [indice] or 'filename',");
  gmic_help(_"mathematical expression 'formula', or compute maximum of all selected images together.\n");

  gmic_option("-mod","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Compute modulo of selected images with value 'value', image [indice] or 'filename',");
  gmic_help(_"mathematical expression 'formula', or compute modulo of all selected images together.\n");

  gmic_option("-and","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Compute bitwise AND of selected images with value 'value', image [indice] or 'filename',");
  gmic_help(_"mathematical expression 'formula', or compute bitwise AND of all selected images together.\n");

  gmic_option("-or","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Compute bitwise OR of selected images with value 'value', image [indice] or 'filename',");
  gmic_help(_"mathematical expression 'formula', or compute bitwise OR of all selected images together.\n");

  gmic_option("-xor","value |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'filename' |");
  gmic_help(___"'formula' |");
  gmic_help(___"(no args)");
  gmic_help(_"Compute bitwise XOR of selected images with value 'value', image [indice] or 'filename',");
  gmic_help(_"mathematical expression 'formula', or compute bitwise XOR of all selected images together.\n");

  gmic_option("-cos","","");
  gmic_help(_"Compute pointwise cosine values of selected images.\n");

  gmic_option("-sin","","");
  gmic_help(_"Compute pointwise sine values of selected images.\n");

  gmic_option("-tan","","");
  gmic_help(_"Compute pointwise tangent values of selected images.\n");

  gmic_option("-acos","","");
  gmic_help(_"Compute pointwise arc-cosine values of selected images.\n");

  gmic_option("-asin","","");
  gmic_help(_"Compute pointwise arc-sine values of selected images.\n");

  gmic_option("-atan","","");
  gmic_help(_"Compute pointwise arc-tangent values of selected images.\n");

  gmic_option("-atan2","[indice]","");
  gmic_help(_"Compute pointwise oriented arc-tangent values of selected images.\n");
  gmic_help(_"The selected images are regarded as containing the y-arguments while the specified");
  gmic_help(_"argument gives the corresponding x-arguments of the atan2() function.\n");

  gmic_option("-abs","","");
  gmic_help(_"Compute pointwise absolute values of selected images.\n");

  gmic_option("-sqr","","");
  gmic_help(_"Compute pointwise square values of selected images.\n");

  gmic_option("-sqrt","","");
  gmic_help(_"Compute pointwise square root values of selected images.\n");

  gmic_option("-exp","","");
  gmic_help(_"Compute pointwise exponential values of selected images.\n");

  gmic_option("-log","","");
  gmic_help(_"Compute pointwise logarithm values of selected images.\n");

  gmic_option("-log10","","");
  gmic_help(_"Compute pointwise logarithm_10 values of selected images.\n");

  gmic_section(" Basic pixel manipulation\n"
               " ------------------------\n");

  gmic_option("-type","datatype","");
  gmic_help(_"Cast all images into specified 'datatype'.\n");
#ifndef gmic_minimal
  gmic_help(_"'datatype' can be { bool | uchar | char | ushort | short | uint |");
  gmic_help(_"                     int | float | double }.\n");
#else
  gmic_help(_"'datatype' can be only 'float' in current minimal mode.\n");
#endif

  gmic_option("-set","value,_x,_y,_z,_c","");
  gmic_help(_"Set scalar value at specified location in selected images.");
  gmic_help(_"(eq. to '-=').\n");
  gmic_help(_"If specified location is outside image bounds, nothing happens.");
  gmic_help(_"Default values for '_x','_y','_z','_v' are '0'.\n");

  gmic_option("-endian","","");
  gmic_help(_"Invert data endianness of selected image buffers.\n");

  gmic_option("-fill","value1,value2,.. |","");
  gmic_help(___"[indice] |");
  gmic_help(___"formula");
  gmic_help(_"Fill selected images with values taken from specified value list, existing image");
  gmic_help(_"or mathematical expression.");
  gmic_help(_"(eq. to '-f').\n");

  gmic_option("-threshold","value[%],_soft |","");
  gmic_help(___"(no args)");
  gmic_help(_"Threshold pixel values of selected images.");
  gmic_help(_"(eq. to '-t').\n");
  gmic_help(_"'_soft' can be { 0=hard thresholding | 1=soft thresholding }.");
  gmic_help(_"(noargs) runs interactive mode (uses the instant window [0] if opened).\n");

  gmic_option("-cut","{ value_min[%] | [indice_min] },{ value_max[%] | [indice_max] } |","");
  gmic_help(___"[indice] |");
  gmic_help(___"(no args)");
  gmic_help(_"Cut pixel values of selected images in specified range.");
  gmic_help(_"(eq. to '-c').\n");
  gmic_help(_"(noargs) runs interactive mode (uses the instant window [0] if opened).\n");

  gmic_option("-normalize","{ value_min[%] | [indice] },{ value_max[%] | [indice] }","");
  gmic_help(___"[indice]");
  gmic_help(_"Linearly normalize pixel values of selected images in specified range.");
  gmic_help(_"(eq. to '-n').\n");

  gmic_option("-round","rounding_value>=0,_rounding_type","");
  gmic_help(_"Round pixel values of selected images.\n");
  gmic_help(_"'_rounding_type' can be { -1=backward | 0=nearest | 1=forward }.\n");

  gmic_option("-equalize","nb_levels>0[%],_value_min[%],_value_max[%]","");
  gmic_help(_"Equalize histograms of selected images.\n");
  gmic_help(_"When range [_valmin,_valmax] is specified the equalization is done only\n");
  gmic_help(_"on the specified value range.\n");

  gmic_option("-quantize","nb_levels>0,_preserve_value_range={0|1}","");
  gmic_help(_"Uniformly quantize selected images.\n");

  gmic_option("-noise","std_variation>=0[%],_noise_type","");
  gmic_help(_"Add random noise to selected images.\n");
  gmic_help(_"'noise_type' can be { 0=gaussian | 1=uniform | 2=salt&pepper | 3=poisson | 4=rice }.\n");

  gmic_option("-rand","value_min,value_max","");
  gmic_help(_"Fill selected images with random values in specified range.\n");

  gmic_option("-norm","","");
  gmic_help(_"Compute pointwise L2-norm of pixels in selected images.\n");

  gmic_option("-orientation","","");
  gmic_help(_"Compute pointwise orientation of pixels in selected images.\n");

  gmic_option("-map","[indice] |","");
  gmic_help(___"predefined_palette");
  gmic_help(_"Map vector-valued palette to selected indexed scalar images.\n");
  gmic_help(_"'predefined_palette' can be { 0=default | 1=rainbow | 2=cluster }.\n");

  gmic_option("-index","{ [indice] | predefined_palette },_is_dithered={0|1},_map_palette={0|1}","");
  gmic_help(_"Index selected vector-valued images by specified palette.\n");
  gmic_help(_"'predefined_palette' can be { 0=default | 1=rainbow | 2=cluster }.\n");

  gmic_section(" Color bases conversions\n"
               " -----------------------\n");

  gmic_option("-rgb2hsv","","");
  gmic_help(_"Convert selected images from RGB to HSV colorbases.\n");

  gmic_option("-rgb2hsl","","");
  gmic_help(_"Convert selected images from RGB to HSL colorbases.\n");

  gmic_option("-rgb2hsi","","");
  gmic_help(_"Convert selected images from RGB to HSI colorbases.\n");

  gmic_option("-rgb2yuv","","");
  gmic_help(_"Convert selected images from RGB to YUV colorbases.\n");

  gmic_option("-rgb2ycbcr","","");
  gmic_help(_"Convert selected images from RGB to YCbCr colorbases.\n");

  gmic_option("-rgb2xyz","","");
  gmic_help(_"Convert selected images from RGB to XYZ colorbases.\n");

  gmic_option("-rgb2lab","","");
  gmic_help(_"Convert selected images from RGB to Lab colorbases.\n");

  gmic_option("-rgb2cmy","","");
  gmic_help(_"Convert selected images from RGB to CMY colorbases.\n");

  gmic_option("-rgb2cmyk","","");
  gmic_help(_"Convert selected images from RGB to CMYK colorbases.\n");

  gmic_option("-hsv2rgb","","");
  gmic_help(_"Convert selected images from HSV to RGB colorbases.\n");

  gmic_option("-hsl2rgb","","");
  gmic_help(_"Convert selected images from HSL to RGB colorbases.\n");

  gmic_option("-hsi2rgb","","");
  gmic_help(_"Convert selected images from HSI to RGB colorbases.\n");

  gmic_option("-yuv2rgb","","");
  gmic_help(_"Convert selected images from YUV to RGB colorbases.\n");

  gmic_option("-ycbcr2rgb","","");
  gmic_help(_"Convert selected images from YCbCr to RGB colorbases.\n");

  gmic_option("-xyz2rgb","","");
  gmic_help(_"Convert selected images from XYZ to RGB colorbases.\n");

  gmic_option("-lab2rgb","","");
  gmic_help(_"Convert selected images from Lab to RGB colorbases.\n");

  gmic_option("-cmy2rgb","","");
  gmic_help(_"Convert selected images from CMY to RGB colorbases.\n");

  gmic_option("-cmyk2rgb","","");
  gmic_help(_"Convert selected images from CMYK to RGB colorbases.\n");

  gmic_section(" Geometric manipulation\n"
               " ----------------------\n");

  gmic_option("-resize","[indice],_interpolation,_borders,_center={0|1} |","");
  gmic_help(___"{[indice] | width>0[%]},_{[indice] | height>0[%]},_{[indice] | depth>0[%]},");
  gmic_help(___"  _{[indice] | spectrum>0[%]},_interpolation,_borders,_center |");
  gmic_help(___"(noargs)");
  gmic_help(_"Resize selected images with specified geometry and interpolation.");
  gmic_help(_"(eq. to '-r').\n");
  gmic_help(_"'interpolation' can be { -1=none (memory) | 0=none | 1=nearest | 2=average |");
  gmic_help(_"                          3=linear | 4=grid | 5=cubic }.");
  gmic_help(_"'borders' can be { -1=none, 2=zero, 3=nearest, 2=repeat }.");
  gmic_help(_"(noargs) runs interactive mode (uses the instant window [0] if opened).\n");

  gmic_option("-resize2x","","");
  gmic_help(_"Resize selected images using the Scale2x algorithm.\n");

  gmic_option("-resize3x","","");
  gmic_help(_"Resize selected images using the Scale3x algorithm.\n");

  gmic_option("-crop","x0[%],x1[%],_borders={0|1} |","");
  gmic_help(___"x0[%],y0[%],x1[%],y1[%],_borders |");
  gmic_help(___"x0[%],y0[%],z0[%],x1[%],y1[%],z1[%],_borders |");
  gmic_help(___"x0[%],y0[%],z0[%],v0[%],x1[%],y1[%],z1[%],v1[%],_borders |");
  gmic_help(___"(noargs)");
  gmic_help(_"Crop selected images from specified geometry.\n");
  gmic_help(_"(noargs) runs interactive mode (uses the instant window [0] if opened).\n");

  gmic_option("-autocrop","color1,color2,..","");
  gmic_help(_"Autocrop selected images using the specified background color.\n");

  gmic_option("-channels","{ [ind0] | v0[%] },_{ [ind1] | v1[%] }","");
  gmic_help(_"Select channels v0..v1 of selected images.\n");

  gmic_option("-slices","{ [ind0] | z0[%] },_{ [ind1] | z1[%] }","");
  gmic_help(_"Select slices z0..z1 of selected images.\n");

  gmic_option("-lines","{ [ind0] | y0[%] },_{ [ind1] | y1[%] }","");
  gmic_help(_"Select lines y0..y1 of selected images.\n");

  gmic_option("-columns","{ [ind0] | x0[%] },_{ [ind1] | x1[%] }","");
  gmic_help(_"Select columns x0..x1 of selected images.\n");

  gmic_option("-rotate","angle,_borders,_interpolation,_cx[%],_cy[%],_zoom","");
  gmic_help(_"Rotate selected images with a given angle.");
  gmic_help(_"'borders' can be { 0=zero | 1=nearest | 2=cyclic }.");
  gmic_help(_"'interpolation' can be { 0=none | 1=linear | 2=cubic }.");
  gmic_help(_"If center ('cx','cy') is specified, the rotation is done in-place.\n");

  gmic_option("-mirror","axis={x|y|z|c}","");
  gmic_help(_"Mirror selected images along specified axis.\n");

  gmic_option("-shift","sx[%],_sy[%],_sz[%],_sv[%],_borders","");
  gmic_help(_"Shift selected images by specified translation vector.\n");
  gmic_help(_"'borders' can be { 0=zero | 1=nearest | 2=cyclic }.\n");

  gmic_option("-transpose","","");
  gmic_help(_"Transpose selected images.\n");

  gmic_option("-invert","","");
  gmic_help(_"Compute inverse of the selected images, viewed as matrices.\n");

  gmic_option("-permute","permutation","");
  gmic_help(_"Permute selected image axes by specified permutation.\n");
  gmic_help(_"'permutation' is a combination of the character set {x|y|z|c},");
  gmic_help(_"for istance 'xycz', 'cxyz', ...\n");

  gmic_option("-unroll","axis={x|y|z|c}","");
  gmic_help(_"Unroll selected images along specified axis.\n");

  gmic_option("-split","axis={x|y|z|c},_parts>0 |","");
  gmic_help(___"patch_x>0,_patch_y>0,_patch_z>0,_patch_v>0,borders={0|1} |");
  gmic_help(___"value,_keep_splitting_values={0|1}");
  gmic_help(_"Split selected images along specified axis, patch or scalar value.");
  gmic_help(_"(eq. to '-s').\n");

  gmic_option("-append","axis={x|y|z|c},_alignement","");
  gmic_help(_"Append selected images along specified axis.");
  gmic_help(_"(eq. to '-a').\n");
  gmic_help(_"'alignement' can be { p=left | c=center | n=right }.\n");

  gmic_option("-warp","[indice],_is_relative={0|1},_interpolation={0|1},_borders,_nb_frames","");
  gmic_help(_"Warp selected image with specified displacement field.\n");
  gmic_help(_"'borders' can be { 0=zero | 1=nearest | 2=cyclic }.\n");

  gmic_section(" Image filtering\n"
               " ---------------\n");

  gmic_option("-blur","std>=0[%],_borders={0|1}","");
  gmic_help(_"Blur selected images by a quasi-gaussian recursive filter.\n");

  gmic_option("-bilateral","std_s>0[%],std_r>0","");
  gmic_help(_"Blur selected images by anisotropic bilateral filtering.\n");

  gmic_option("-denoise","std_s>=0,_std_p>=0,_patch_size>0,_lookup_size>0,_smoothness,_approx={0|1}","");
  gmic_help(_"Denoise selected images with a patch-averaging procedure.\n");

  gmic_option("-smooth","amplitude>=0,_sharpness>=0,_anisotropy,_alpha,_sigma,_dl>0,_da>0,_precision>0,","");
  gmic_help(___" interpolation,_fast_approx={0|1} |");
  gmic_help(___"nb_iters>=0,_sharpness>=0,_anisotropy,_alpha,_sigma,_dt>0,0 |");
  gmic_help(___"[indice],_amplitude>=0,_dl>0,_da>0,_precision>0,_interpolation,_fast_approx={0|1} |");
  gmic_help(___"[indice],_nb_iters>=0,_dt>0,0");
  gmic_help(_"Smooth selected images anisotropically using diffusion PDE's.\n");
  gmic_help(_"'_anisotropy' must be in [0,1].");
  gmic_help(_"'_interpolation' can be { 0=nearest, 1=linear, 2=runge-kutta }.\n");

  gmic_option("-edgetensors","sharpness>=0,_anisotropy,_alpha,_sigma,is_sqrt={0|1}","");
  gmic_help(_"Compute diffusion tensors for edge-preserving smoothing from selected images.\n");
  gmic_help(_"'_anisotropy' must be in [0,1].\n");

  gmic_option("-median","radius>=0","");
  gmic_help(_"Apply median filter of specified radius on selected images.\n");

  gmic_option("-sharpen","amplitude>=0 |","");
  gmic_help(___"amplitude>=0,1,_edge>=0,_alpha,_sigma");
  gmic_help(_"Sharpen selected images by inverse diffusion or shock filters methods.\n");

  gmic_option("-convolve","[indice],_borders={0|1}","");
  gmic_help(_"Convolve selected images by specified mask.\n");

  gmic_option("-correlate","[indice],_borders={0|1}","");
  gmic_help(_"Correlate selected images by specified mask.\n");

  gmic_option("-erode","size>=0' |","");
  gmic_help(___"size_x>=0,size_y>=0,_size_z>=0 |");
  gmic_help(___"[indice],_borders={0|1}");
  gmic_help(_"Erode selected images by specified mask.\n");

  gmic_option("-dilate","size>=0 |","");
  gmic_help(___"size_x>=0,size_y>=0,size_z>=0 |");
  gmic_help(___"[indice],_borders={0|1}");
  gmic_help(_"Dilate selected images by specified mask.\n");

  gmic_option("-inpaint","[indice]","");
  gmic_help(_"Inpaint selected images by specified mask.\n");

  gmic_option("-gradient","{x|y|z}..{x|y|z} |","");
  gmic_help(___"(no args)");
  gmic_help(_"Compute gradient components of selected images.\n");
  gmic_help(_"(no args) compute all significant components.\n");

  gmic_option("-hessian","{xx|xy|xz|yy|yz|zz}..{xx|xy|xz|yy|yz|zz} |","");
  gmic_help(___"(no args)");
  gmic_help(_"Compute hessian components of selected images.\n");
  gmic_help(_"(no args) compute all significant components.\n");

  gmic_option("-haar","scale>0","");
  gmic_help(_"Compute direct Haar multiscale wavelet transform of selected images.\n");

  gmic_option("-ihaar","scale>0","");
  gmic_help(_"Compute inverse Haar multiscale wavelet transform of selected images.\n");

  gmic_option("-fft","","");
  gmic_help(_"Compute direct Fourier transform of selected images.\n");

  gmic_option("-ifft","","");
  gmic_help(_"Compute inverse Fourier transform of selected images.\n");

  gmic_section(" Image creation and drawing\n"
               " --------------------------\n");

  gmic_option("-histogram","nb_levels>0[%],_valmin[%],_valmax[%]","");
  gmic_help(_"Compute histogram of selected images.\n");

  gmic_option("-distance","isovalue","");
  gmic_help(_"Compute unsigned distance functions to specified isovalue.\n");

  gmic_option("-eikonal","nb_iterations>=0,_band_size>=0","");
  gmic_help(_"Compute signed distance functions to 0, using Eikonal PDE.\n");

  gmic_option("-label","","");
  gmic_help(_"Label connected components of selected images.\n");

  gmic_option("-displacement","[indice],_smoothness>=0,_precision>0,_nbscales>=0,itermax>=0,","");
  gmic_help(___"is_backward={0|1}");
  gmic_help(_"Estimate displacement field between selected images and specified source.\n");
  gmic_help(_"If '_nbscales==0', the number of used scales is automatically estimated.\n");

  gmic_option("-sort","","");
  gmic_help(_"Sort values of selected images in increasing order.\n");

  gmic_option("-psnr","_maximum_value","");
  gmic_help(_"Compute PSNR values between selected images and store them in a matrix.\n");

  gmic_option("-point","x[%],y[%],_z[%],_opacity,_color1,..","");
  gmic_help(_"Set specified colored pixel on selected images.\n");

  gmic_option("-line","x0[%],y0[%],x1[%],y1[%],_opacity,_color1,..'","");
  gmic_help(_"Draw specified colored line on selected images.\n");

  gmic_option("-polygon","N,x1[%],y1[%],..,xN[%],yN[%],_opacity,_color1,..","");
  gmic_help(_"Draw specified colored N-polygon on selected images.\n");

  gmic_option("-spline","x0,y0,u0,v0,x1,y1,u1,v1,_opacity,_color1,..","");
  gmic_help(_"Draw specified colored spline curve on selected images.\n");

  gmic_option("-ellipse","x[%],y[%],r[%],R[%],_angle,_opacity,_color1,..","");
  gmic_help(_"Draw specified colored ellipse on selected images.\n");

  gmic_option("-text","text,_x[%],_y[%],_size>0,_opacity,_color1,..","");
  gmic_help(_"Draw specified colored text string on selected images.\n");

  gmic_option("-image","[indice],_x[%],_y[%],_z[%],_opacity,_[indice_mask]","");
  gmic_help(_"Draw specified sprite image on selected images.\n");

  gmic_option("-object3d","[indice],_x[%],_y[%],_z,_opacity,_is_zbuffer={0|1}","");
  gmic_help(_"Draw specified 3D object on selected images.\n");

  gmic_option("-plasma","alpha,_beta,_opacity","");
  gmic_help(_"Draw a random colored plasma on selected images.\n");

  gmic_option("-mandelbrot","z0r,z0i,z1r,z1i,_itermax>=0,_is_julia={0|1},_c0r,_c0i,_opacity","");
  gmic_help(_"Draw Mandelbrot/Julia fractals on selected images.\n");

  gmic_option("-quiver","[indice],_sampling>0,_factor,_quiver_type={0|1},_opacity,_color1,..","");
  gmic_help(_"Draw a 2D vector field on selected images.\n");

  gmic_option("-flood","x[%],_y[%],_z[%],_tolerance>=0,_opacity,_color1,..","");
  gmic_help(_"Flood-fill selected images using specified fill value and tolerance.\n");

  gmic_section(" List manipulation\n"
               " -----------------\n");

  gmic_option("-remove","","");
  gmic_help(_"Remove selected images from the list.");
  gmic_help(_"(eq. to '-rm').\n");

  gmic_option("-keep","","");
  gmic_help(_"Keep only selected images in the list.");
  gmic_help(_"(eq. to '-k').\n");

  gmic_option("-move","position","");
  gmic_help(_"Move selected images at specified position in the list.");
  gmic_help(_"(eq. to '-mv').\n");

  gmic_option("-reverse","","");
  gmic_help(_"Reverse position order of selected images.\n");

  gmic_option("-name","[indice] |","");
  gmic_help(___"name");
  gmic_help(_"Set name of selected images.\n");
  gmic_help(_"The specified name is either the one from another image or from the specified 'name'.\n");

  gmic_section(" 3D rendering\n"
               " ------------\n");

  gmic_option("-line3d","x0,y0,z0,x1,y1,z1","");
  gmic_help(_"Insert a 3D line object at the end of the list.\n");

  gmic_option("-triangle3d","x0,y0,z0,x1,y1,z1,x2,y2,z2","");
  gmic_help(_"Insert a 3D triangle object at the end of the list.\n");

  gmic_option("-quadrangle3d","x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3","");
  gmic_help(_"Insert a 3D quadrangle object at the end of the list.\n");

  gmic_option("-box3d","size |","");
  gmic_help(___"size_x,size_y,size_z");
  gmic_help(_"Insert a 3D box object at the end of the list.\n");

  gmic_option("-cone3d","radius,_size_z,_subdivisions>0","");
  gmic_help(_"Insert a 3D cone object at the end of the list.\n");

  gmic_option("-cylinder3d","radius,_height,_subdivisions>0","");
  gmic_help(_"Insert a 3D cylinder object at the end of the list.\n");

  gmic_option("-torus3d","radius1,_radius2,_subdivisions1>0,_subdivisions2>0","");
  gmic_help(_"Insert a 3D torus object at the end of the list.\n");

  gmic_option("-plane3d","size1_size2,_subdivisions1>0,_subdisivions2>0","");
  gmic_help(_"Insert a 3D plane object at the end of the list.\n");

  gmic_option("-sphere3d","radius,_recursions>=0","");
  gmic_help(_"Insert a 3D sphere object at the end of the list.\n");

  gmic_option("-elevation3d","z-factor |","");
  gmic_help(___"[indice] |");
  gmic_help(___"'formula',_x0,_y0,_x1,y1,_dx[%],_dy[%] |");
  gmic_help(___"(no args)");
  gmic_help(_"Compute 3D elevation of selected images with specified elevation map.\n");
  gmic_help(_"If a 'z-factor' is specified, the elevation maps are given by the norm of the");
  gmic_help(_"selected images themselves. Else, elevation values are taken from the image");
  gmic_help(_"[indice] or from the specified mathematical expression 'formula'.\n");

  gmic_option("-isoline3d","isovalue[%] |","");
  gmic_help(___"'formula',value,_x0,_y0,_x1,_y1,_dx>0[%],_dy>0[%]");
  gmic_help(_"Compute 3D isoline from selected images or from specified mathematical expression.\n");

  gmic_option("-isosurface3d","isovalue[%] |","");
  gmic_help(___"'formula',value,_x0,_y0,_z0,_x1,_y1,_z1,_dx>0[%],_dy>0[%],_dz>0[%]");
  gmic_help(_"Compute 3D isosurfaces from selected images or from specified mathematical expression.\n");

  gmic_option("-streamline3d","x,y,z,_L>=0,_dl>0,_interpolation,_is_backward={0|1},_is_oriented={0|1} |","");
  gmic_help(___"'formula',x,y,z,_L>=0,_dl>0,_interpolation,_is_backward={0|1},_is_oriented={0|1}");
  gmic_help(_"Compute 3D streamlines from selected vector fields or from specified mathematical expression.\n");
  gmic_help(_"'interpolation' can be { 0=nearest integer, 1=1st-order, 2=2nd-order, 3=4th-order }.\n");

  gmic_option("-add3d","tx,_ty,_tz |","");
  gmic_help(___"[indice] |");
  gmic_help(___"(noargs)");
  gmic_help(_"Shift selected 3D objects with specified translation vector, or merge them");
  gmic_help(_"with 3D object [indice], or merge selected 3D objects together.");
  gmic_help(_"(eq. to '-+3d').\n");

  gmic_option("-add3d","tx,_ty,_tz","");
  gmic_help(_"Shift selected 3D objects with the opposite of specified translation vector.");
  gmic_help(_"(eq. to '--3d').\n");

  gmic_option("-mul3d","factor |","");
  gmic_help(___"factor_x,factor_y,_factor_z");
  gmic_help(_"Scale selected 3D objects isotropically or anisotropically, with specified factors.");
  gmic_help(_"(eq. to '-*3d').\n");

  gmic_option("-div3d","factor |","");
  gmic_help(___"factor_x,factor_y,_factor_z");
  gmic_help(_"Scale selected 3D objects isotropically or anisotropically, with the inverse of");
  gmic_help(_"specified factors.");
  gmic_help(_"(eq. to '-/3d').\n");

  gmic_option("-center3d","","");
  gmic_help(_"Center positions of selected 3D objects.");
  gmic_help(_"(eq. to '-c3d').\n");

  gmic_option("-normalize3d","","");
  gmic_help(_"Normalize selected 3D objects to unit size.");
  gmic_help(_"(eq. to '-n3d').\n");

  gmic_option("-rotate3d","u,v,w,angle","");
  gmic_help(_"Rotate selected 3D objects around axis (u,v,w) with specified angle (in degree.).");
  gmic_help(_"(eq. to '-rot3d').\n");

  gmic_option("-color3d","R,G,B,_opacity","");
  gmic_help(_"Set color and opacity of selected 3D objects.");
  gmic_help(_"(eq. to '-col3d').\n");

  gmic_option("-opacity3d","opacity","");
  gmic_help(_"Set opacity of selected 3D objects.");
  gmic_help(_"(eq. to '-o3d').\n");

  gmic_option("-reverse3d","","");
  gmic_help(_"Invert primitive orientations of selected 3D objects.");
  gmic_help(_"(eq. to '-r3d').\n");

  gmic_option("-primitives3d","mode","");
  gmic_help(_"Set specified primitive mode for selected 3D objects.");
  gmic_help(_"(eq. to '-p3d').\n");
  gmic_help(_"'mode' can be { 0=points | 1=segments }.\n")

  gmic_option("-split3d","","");
  gmic_help(_"Split selected 3D objects data into 6 1-column data images :");
  gmic_help(_"'header', '(nb_vertices,nb_primitives)', 'vertices', 'primitives', 'colors', 'opacities' ].");
  gmic_help(_"(eq. to '-s3d').\n");

  gmic_option("-light3d","position_x,position_y,position_z","");
  gmic_help(_"Set position of the light for 3D rendering.");
  gmic_help(_"(eq. to '-l3d').\n");

  gmic_option("-focale3d","focale","");
  gmic_help(_"Set focale value for 3D rendering.");
  gmic_help(_"(eq. to '-f3d').\n");

  gmic_option("-specl3d","value","");
  gmic_help(_"Set amount of specular light for 3D rendering.");
  gmic_help(_"(eq. to '-sl3d').\n");

  gmic_option("-specs3d","value","");
  gmic_help(_"Set shininess of specular light for 3D rendering.");
  gmic_help(_"(eq. to '-ss3d').\n");

  gmic_option("-double3d","mode={0|1}","");
  gmic_help(_"Set/unset double-sided mode for 3D rendering.");
  gmic_help(_"(eq. to '-db3d').\n");

  gmic_option("-render3d","mode","");
  gmic_help(_"Set static 3D rendering mode.");
  gmic_help(_"(eq. to '-m3d').\n");
  gmic_help(_"'mode' can be { -1=bounding-box | 0=pointwise | 1=linear | 2=flat | 3=flat-shaded |");
  gmic_help(_"                 4=gouraud-shaded | 5=phong-shaded }.");


  gmic_option("-renderd3d","mode","");
  gmic_help(_"Set dynamic 3D rendering mode.");
  gmic_help(_"(eq. to '-md3d').\n");
  gmic_help(_"'mode' can be { -1=bounding-box | 0=pointwise | 1=linear | 2=flat | 3=flat-shaded |");
  gmic_help(_"                 4=gouraud-shaded | 5=phong-shaded }.");

  gmic_option("-background3d","R,_G,_B","");
  gmic_help(_"Define background color in 3D viewer.");
  gmic_help(_"(eq. to '-b3d').\n");

  gmic_section(" Program controls\n"
               " ----------------\n");

  gmic_option("-nop","","");
  gmic_help(_"Do nothing.\n");

  gmic_option("-skip","item","");
  gmic_help(_"Do nothing but skip specified item.\n");

  gmic_option("-return","","");
  gmic_help(_"Return from current command scope.\n");

  gmic_option("-exec","command","");
  gmic_help(_"Execute external command, using a system call.\n");

  gmic_option("-do","","");
  gmic_help(_"Start a 'do..while' code bloc.\n");

  gmic_option("-while","condition","");
  gmic_help(_"End a 'do..while' code bloc and go back to associated '-do'");
  gmic_help(_"if specified 'condition' is verified.\n");
  gmic_help(_"'condition' must be a number standing for { 0=false | other=true }.\n");

  gmic_option("-if","condition","");
  gmic_help(_"Start a 'if..elif..else..endif' code bloc and test");
  gmic_help(_"if specified 'condition' is verified.\n");
  gmic_help(_"'condition' must be a number standing for { 0=false | other=true }.\n");

  gmic_option("-elif","condition","");
  gmic_help(_"Start a 'elif..else..endif' code bloc if previous '-if' was not verified");
  gmic_help(_"and test if 'condition' is verified.\n");

  gmic_option("-else","","");
  gmic_help(_"Execute following commands if previous '-if' or '-elif' conditions failed.\n");

  gmic_option("-endif","","");
  gmic_help(_"End a 'if..elif..else..endif' code bloc.\n");

  gmic_option("-repeat","number","");
  gmic_help(_"Start 'number' iterations of a 'repeat..done' code loop.\n");

  gmic_option("-done","","");
  gmic_help(_"End a 'repeat..done' code loop, and go to associated '-repeat' if iterations remain.\n");

  gmic_option("-check","condition","");
  gmic_help(_"Check 'condition', and quit interpreter if not verified.\n");

  gmic_option("-quit","","");
  gmic_help(_"Quit interpreter.");
  gmic_help(_"(eq. to '-q').\n");

  gmic_option("-push","item","");
  gmic_help(_"Push 'item' on the global stack at selected positions.");
  gmic_help(_"(eq. to '-p').\n");
  gmic_help(_"Specified selection (if any) for this command stands for stack indices, not image indices.\n");

  gmic_option("-pop","","");
  gmic_help(_"Pop items from the global stack at selected positions.");
  gmic_help(_"(eq. to '-pp').\n");
  gmic_help(_"Specified selection (if any) for this command stands for stack indices, not image indices.\n");

  gmic_option("-local","","");
  gmic_help(_"Start a local environment with the selected images.");
  gmic_help(_"(eq. to '-l').\n");

  gmic_option("-endlocal","","");
  gmic_help(_"End the previous local environment.");
  gmic_help(_"(eq. to '-endl').\n");

  gmic_option("-patch","size_x,_size_y,_size_z,_size_v,borders","");
  gmic_help(_"Enable patch processing environment with the selected images.\n");

  gmic_option("-endpatch","","");
  gmic_help(_"End the previous patch processing environment.");
  gmic_help(_"(eq. to '-endp').\n");

  gmic_option("-print","","");
  gmic_help(_"Print informations on selected images.\n");

  gmic_option("-echo","message","");
  gmic_help(_"Output specified message on the standard output.");
  gmic_help(_"(eq. to '-e').\n");

  gmic_option("-error","message","");
  gmic_help(_"Print error message, and quit interpreter.\n");

  gmic_option("-warning","message","");
  gmic_help(_"Print warning message.\n");

  gmic_option("-progress","0<=value<=100 |","");
  gmic_help(___"-1");
  gmic_help(_"Set the progression state of the current process.\n");
  gmic_help(_"It is useful only when the G'MIC interpreter is used in an embedding application.\n");

  gmic_section(" Input/output\n"
               " ------------\n");

  gmic_option("-input","filename |","");
  gmic_help(___"{ width>0[%] | [indw] },{ _height>0[%] | [indh] },{ _depth>0[%] | [indd] },");
  gmic_help(___"{ _spectrum>0[%] | [inds] },_value1,.. |");
  gmic_help(___"[indice]x_nb_copies>0 |");
  gmic_help(___"(value1{,|;|/|^}value2{,|;|/|^}..)");
  gmic_help(_"Insert new image from a filename or from a copy of an existing image ['indice'],");
  gmic_help(_"or insert new image with specified values.");
  gmic_help(_"(eq. to '-i' | (no args)).\n");

  gmic_option("-output","filename,_format_specific_options","");
  gmic_help(_"Output selected images into one or several filenames.");
  gmic_help(_"(eq. to '-o').\n");

  gmic_option("-display","","");
  gmic_help(_"Display selected images with an interactive viewer.");
  gmic_help(_"(eq. to '-d').\n");

  gmic_option("-display3d","","");
  gmic_help(_"Display selected 3D objects with an interactive viewer.");
  gmic_help(_"(eq. to '-d3d').\n");

  gmic_option("-window","width>=-1,_height>=-1,_normalization,_fullscreen |","");
  gmic_help(___"(no args)");
  gmic_help(_"Display selected images into an instant window of specified size and");
  gmic_help(_"normalization type.");
  gmic_help(_"(eq. to '-w').\n");
  gmic_help(_"'normalization' can be { -1=keep same | 0=none | 1=always | 2=1st-time | 3=auto }.");
  gmic_help(_"'fullscreen' can be { 0=no, 1=yes }.");
  gmic_help(_"If 'width'=0, the instant window is closed. If width=-1, it is resized");
  gmic_help(_"to its current GUI-window size. You can also manage up to 10 different instant windows");
  gmic_help(_"by using the commands '-w0' (eq. to '-w'),'-w1',..,'-w9'.\n");

  gmic_option("-wait","delay","");
  gmic_help(___"(no args)");
  gmic_help(_"Wait a given delay or wait for any event from the instant window.\n");
  gmic_help(_"'delay' can be { <0=delay+flush |  0=event | >0=delay }.");
  gmic_help(_"Specified selection (if any) for this command stands for instant window indices, not image indices.\n");

  gmic_option("-plot","_plot_type,_vertex_type,_xmin,_xmax,_ymin,_ymax |","");
  gmic_help(___"'formula',_xmin,xmax,_ymin,_ymax,_resolution,_plot_type,_vertex_type");
  gmic_help(_"Display image or mathematical expression with an interactive viewer.");
  gmic_help(_"'plot_type' can be { 0=none | 1=lines | 2=splines | 3=bar }.");
  gmic_help(_"'vertex_type' can be { 0=none | 1=points | 2|3=crosses | 4|5=circles | 6|7=squares }.\n");

  gmic_option("-select","feature","");
  gmic_help(_"Interactively select a feature from selected images.\n");
  gmic_help(_"'feature' can be { 0=point | 1=segment | 2=rectangle | 3=ellipse }.");
  gmic_help(_"If the instant display window [0] is active, it is used for the selection.\n");

  // Print descriptions of user-defined custom commands.
  char line[256*1024] = { 0 }, name[4096] = { 0 }, args[4096] = { 0 }, desc[4096] = { 0 };
  bool first_description = true;
  for (int i = 1; i<argc-1; ++i) if (!std::strcmp("-m",argv[i]) || !std::strcmp("-command",argv[i])) {
    std::FILE *file = cimg::fopen(argv[i+1],"r");
    if (file) {
      int err = 0;
      while ((err=std::fscanf(file,"%262143[^\n] ",line)>=0)) {
        if (err) {
          name[0] = args[0] = desc[0] = 0;
          if (line[0]=='#' && std::sscanf(line,"#@gmic %4095[^:]:%4095[^:]:%4095[^:]",name,args,desc)>0) {
            if (first_description) {
              gmic_section(" Commands : User-defined custom commands\n"
                           " ---------------------------------------\n");
            }
            std::sprintf(line,"-%s",name); cimg::strpare(line); gmic_option(line,args,"");
            std::sprintf(line,_"%s\n",desc); cimg::strescape(line); gmic_help(line);
            first_description = false;
          }
        }
      }
    }
    cimg::fclose(file);
  }

  // Print descriptions of default commands.
  first_description = true;
  for (const char *data = data_gmic_def; *data; ) {
    if (*data=='\n') ++data;
    else {
      if (std::sscanf(data,"%262143[^\n]",line)>0) data+=std::strlen(line);
      name[0] = args[0] = desc[0] = 0;
      if (line[0]=='#' && std::sscanf(line,"#@gmic %4095[^:]:%4095[^:]:%4095[^:]",name,args,desc)>0) {
        if (first_description) {
          gmic_section(" Commands : Default custom commands\n"
                       " ----------------------------------\n");
        }
        std::sprintf(line,"-%s",name); cimg::strpare(line); gmic_option(line,args,"");
        std::sprintf(line,_"%s\n",desc); cimg::strescape(line); gmic_help(line);
        first_description = false;
      }
    }
  }

  gmic_section(" Viewers shortcuts\n"
               " -----------------\n");

  gmic_help(" The following shortcuts are available for the viewers of images or 3D objects :\n");
  gmic_help("   - CTRL+D : Increase window size.");
  gmic_help("   - CTRL+C : Decrease window size.");
  gmic_help("   - CTRL+R : Reset window to its initial size.");
  gmic_help("   - CTRL+F : Toggle fullscreen mode.");
  gmic_help("   - CTRL+S : Save current window snapshot.");
  gmic_help("   - CTRL+O : Save current instance of viewed image (or 3D object).\n");

  gmic_help(" Special shortcuts for the viewer of 2D/3D images are :\n");
  gmic_help("   - CTRL+P             : Play stack of 3D frames as a movie.");
  gmic_help("   - CTRL+(mousewheel)  : Zoom in/out.");
  gmic_help("   - SHIFT+(mousewheel) : Go left/right.");
  gmic_help("   - ALT+(mousewheel)   : Go up/down.");
  gmic_help("   - Numeric PAD        : Zoom in/out (+/-) and move zoomed region (digits).");
  gmic_help("   - BACKSPACE          : Reset zoom to its initial scale.\n");

  gmic_help(" Special shortcuts for the viewer of 3D objects are :\n");
  gmic_help("   - (mouse)+(left mouse button)   : Rotate object.");
  gmic_help("   - (mouse)+(right mouse button)  : Zoom object.");
  gmic_help("   - (mouse)+(middle mouse button) : Shift object.");
  gmic_help("   - (mousewheel)                  : Zoom in/out.");
  gmic_help("   - CTRL+F1 ... CTRL+F6           : Swap between rendering modes.");
  gmic_help("   - CTRL+Z                        : Enable/disable Z-buffer.");
  gmic_help("   - CTRL+A                        : Show/hide 3D axes.");
  gmic_help("   - CTRL+T                        : Enable/disable double-sided triangles.\n");

  gmic_section(" File options\n"
               " ------------\n");

  gmic_help(" 'gmic' is able to read/write most of the classical image file formats, including :\n");
  gmic_help("   - 2D grayscale/color images : PNG, JPEG, GIF, PNM, TIFF, BMP, ..");
  gmic_help("   - 3D volumetric images : DICOM, HDR, NII, PAN, CIMG, INR, ..");
  gmic_help("   - Video files : MPEG, AVI, MOV, OGG, FLV, ..");
  gmic_help("   - Generic data files : DLM, ASC, RAW, TXT.");
  gmic_help("   - 3D objects : OFF.\n");

  gmic_help(" Specific file format options :\n");
  gmic_help("   - For video files : you can read only sub-frames of the sequence with");
  gmic_help("     'video.ext,[first_frame[%][,last_frame[%][,step]]]'.");
  gmic_help("   - For RAW binary files : you must specify image dimensions with");
  gmic_help("     'file.raw,width[,height[,depth[,dim]]]]'.");
  gmic_help("   - For YUV files : you must specify the image dimensions and can read only sub-frames");
  gmic_help("     of the image sequence with");
  gmic_help("     'file.yuv,width,height[,first_frame[,last_frame[,step]]]'.");
  gmic_help("   - For JPEG files : you can specify the quality (in %) of an output jpeg file with");
  gmic_help("     'file.jpg,30%'.");
  gmic_help("   - If an input file has extension '.gmic', it is read as a G'MIC custom command file.\n");

  gmic_section(" Mathematical expressions\n"
               " ------------------------\n");

  gmic_help(" Some G'MIC commands may take a 'formula' as an argument. Actually, G'MIC has a simple");
  gmic_help(" but useful parser of mathematical expressions. This parser is able to understand the");
  gmic_help(" following functions, operators and variable names :\n");
  gmic_help("   - Variables 'x','y','z','c' returns the current pixel coordinates of an associated");
  gmic_help("     image, if any (else, their are equal to '0').");
  gmic_help("   - Variables 'w','h','d','s' return the associated image dimensions (if any).");
  gmic_help("   - Variable 'i' return the current pixel value of the associated image (if any),");
  gmic_help("     at coordinates (x,y,z,c).");
  gmic_help("   - Using 'i' as a function 'i(X,_Y,_Z,_C,_borders)' is possible. It returns");
  gmic_help("     the pixel value at another location (X,_Y,_Z,_C) of the associated image.");
  gmic_help("   - Variable '?' (or 'u') and 'g' returns random values following respectively");
  gmic_help("     uniform or gaussian distributions.");
  gmic_help("   - Using '?' (or 'u') as a function '?(min,max)' is possible. It returns");
  gmic_help("     a random number in specified value range [min,max].");
  gmic_help("   - Variables 'pi' and 'e' are replaced by their classical values (3.14.. and 2.71..).");
  gmic_help("   - Most of the classical arithmetic and logical operators can be used, as well as");
  gmic_help("     usual mathematical functions (same syntax as in C).");
  gmic_help("   - Function 'if(condition,action,action_else)' can be used for conditional tests.");
  gmic_help("   - User-defined variables can be assigned and re-used inside a mathematical");
  gmic_help("     expression, using the '=' operator. Variable names are case-sensitive.");
  gmic_help("   - Separators ';' can be used to separate multiple variable definitions.");
  gmic_help("   - When using separators ';', only the evaluation of the last expression is returned.\n");

  gmic_section(" Examples of use\n"
               " ---------------\n");
  gmic_help(" 'gmic' is a simple but quite complete interpreter of image processing commands");
  gmic_help(" and can be used in a wide variety of situations for various image processing tasks.");
  gmic_help(" Here are few examples of possible applications :\n");

  gmic_help("   - View images : ");
  gmic_help("       gmic file1.bmp file2.jpeg\n");

  gmic_help("   - Convert image file : ");
  gmic_help("       gmic input.bmp -o output.jpg\n");

  gmic_help("   - Create volumetric image from a movie sequence : ");
  gmic_help("       gmic input.mpg -a z -o output.hdr\n");

  gmic_help("   - Compute image gradient norm : ");
  gmic_help("       gmic input.bmp -gradient_norm\n");

  gmic_help("   - Denoise a color image :");
  gmic_help("       gmic image.jpg -denoise 30,10 -o denoised.jpg\n");

  gmic_help("   - Compose two images using overlay :");
  gmic_help("       gmic image1.jpg image2.jpg -compose_overlay -o composed.jpg\n");

  gmic_help("   - Evaluate a mathematical expression :");
  gmic_help("       gmic -e \"cos(pi/4)^2+sin(pi/4)^2={cos(pi/4)^2+sin(pi/4)^2}\"\n");

  gmic_help("   - Plot a 2D function :");
  gmic_help("       gmic 1000,1,1,2 -f \"X=3*(x-500)/500;X^2*sin(3*X^2)+if(c==0,u(0,-1),cos(X*10))\" -plot\n");

  gmic_help("   - Plot a 3D elevated function in random colors:");
  gmic_help("       gmic 128,128,1,3,\"?(0,255)\" -plasma 10,3 -blur 4 -sharpen 10000 \\");
  gmic_help("       128,128,1,1,\"X=(x-64)/6;Y=(y-64)/6;100*exp(-(X^2+Y^2)/30)*abs(cos(X)*sin(Y))\"\\");
  gmic_help("       -elevation3d[-2] [-1]\n");

  gmic_help("   - Plot the isosurface of a 3D volume :");
  gmic_help("       gmic -m3d 5 -md3d 5 -db3d 0 -isosurface3d \"'x^2+y^2+abs(z)^abs(4*cos(x*y*z*3))'\",3\n");

  gmic_help("   - Create a G'MIC 3D logo : ");
  gmic_help("       gmic 180,70,1,3 -text G\\'MIC,30,5,50,1,1 -blur 2 -n 0,100 --plasma 0.4 -+ \\");
  gmic_help("       -blur 1 -elevation3d -0.1 -md3d 4\n");

  gmic_help("   - Create a 3D ring of torii :");
  gmic_help("       gmic -repeat 20 -torus3d 15,2 -col3d[-1] \"{?(60,255)},{?(60,255)},{?(60,255)}\" \\");
  gmic_help("       -*3d[-1] 0.5,1 -if \"{@{>,-1}%2}\" -rot3d[-1] 0,1,0,90 -endif -+3d[-1] 70 -+3d \\");
  gmic_help("       -rot3d 0,0,1,18 -done -md3d 3 -m3d 5 -db3d 0\n");

  gmic_help("   - Create a vase from a 3D isosurface :");
  gmic_help("       gmic -md3d 4 -isosurface3d \"'x^2+2*abs(y/2)*sin(2*y)^2+z^2-3',0\" -sphere3d 1.5 --3d[-1] 0,5 \\");
  gmic_help("       -plane3d 15,15 -rot3d[-1] 1,0,0,90 -c3d[-1] -+3d[-1] 0,3.2 -col3d[-1] 180,150,255 \\");
  gmic_help("       -col3d[-2] 128,255,0 -col3d[-3] 255,128,0 -+3d\n");
  gmic_help(" See also the command file 'gmic_def.raw' for more examples of G'MIC commands.\n");

  gmic_help(" ** G'MIC comes with ABSOLUTELY NO WARRANTY; for details visit http://gmic.sourceforge.net **");

  return is_help_displayed;
}

//-----------------------
// Start main procedure.
//-----------------------
int main(int argc, char **argv) {

  // Display help if necessary.
  //---------------------------
  cimg::output(stdout);
  if (argc==1) {
    std::fprintf(cimg::output(),"[gmic] No options or data provided. Try '%s -h' for help.\n",cimg::basename(argv[0]));
    std::fflush(cimg::output());
    std::exit(0);
  }

  char name[1024] = { 0 };
  const char
    *const is_help1 = cimg_option("-h",(char*)0,0),
    *const is_help2 = cimg_option("-help",(char*)0,0),
    *const is_help3 = cimg_option("--help",(char*)0,0);

  if (is_help1 || is_help2 || is_help3) {
    const char
      *const is_help = is_help1?"-h":is_help2?"-help":"--help",
      *command = is_help1?is_help1:is_help2?is_help2:is_help3;
    if (!std::strcmp(is_help,command)) help(argc,argv);  // Display general help.
    else { // Display help only for a specified command.
      if (command[0]!='-') std::sprintf(name,"-%s",command);
      else if (command[1]!='-') std::strcpy(name,command);
      else std::strcpy(name,command+1);
      if (!help(argc,argv,name,true))
        std::fprintf(cimg::output(),"[gmic] Command '%s' has no description. Try '%s -h' for global help.\n\n",
                     name+1,cimg::basename(argv[0]));
    }
    std::exit(0);
  }

  // Launch G'MIC instance.
  //-----------------------
  cimg::output(stderr);
  CImgList<float> images;
  try { gmic(argc,argv,images); }
  catch (gmic_exception &e) {
    std::fprintf(cimg::output(),"\n[gmic] %s\n\n",e.message());
    if (*e.command()) {
      const char *_argv[2] = { "gmic", "-h" };
      if (*e.command()!='-') std::sprintf(name,"-%s",e.command()); else std::strcpy(name,e.command());
      help(2,_argv,name,false);
    }
    std::fflush(cimg::output());
    return -1;
  }
  return 0;
}
#endif

#endif // #ifdef cimg_plugin .. #else ..
