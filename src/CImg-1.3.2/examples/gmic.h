/*
  #
  #  File        : gmic.h
  #                ( C++ header file )
  #
  #  Description : GREYC's Magic Image Converter
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

#ifndef gmic_version
#define gmic_version 1328

// Define environement variables.
#ifndef cimg_debug
#define cimg_debug 1
#endif

#if defined(cimg_build)
#define cimg_plugin "examples/gmic.cpp"
#define cimg_location "../CImg.h"
#elif defined(gmic_build)
#define cimg_plugin "gmic.cpp"
#define cimg_location "./CImg.h"
#endif

// Define the structures to store images and image lists.
#if defined(cimg_build) || defined(gmic_build)
#include cimg_location
#else
#include <cstdio>
#include <cstring>

namespace cimg_library {
  template<typename T> struct CImg {
    unsigned int width;       // Number of image columns.
    unsigned int height;      // Number of image lines.
    unsigned int depth;       // Number of image slices.
    unsigned int spectrum;    // Number of image channels.
    bool is_shared;           // Tells if data buffer is shared by another structure (usually, false).
    T *data;                  // Pointer to first image pixel.
    ~CImg();
    CImg():width(0),height(0),depth(0),spectrum(0),is_shared(false),data(0) {}
    CImg<T>& assign(const unsigned int w, const unsigned int h=1, const unsigned int d=1, const unsigned int s=1);
  };

  template<typename T> struct CImgList {
    unsigned int size;           // Number of images in the list.
    unsigned int allocated_size; // Allocated items in the list (must be >size and a power of 2).
    CImg<T> *data;               // Pointer to first image of the list.
    ~CImgList();
    CImgList():size(0),allocated_size(0),data(0) {}
    CImgList<T>& assign(const unsigned int n);
  };
}
#endif
#define gmic_image cimg_library::CImg
#define gmic_list cimg_library::CImgList

// The lines below are necessary when using a non-standard compiler such as visualcpp6.
#ifdef cimg_use_visualcpp6
#define std
#endif
#ifdef min
#undef min
#undef max
#endif

// Define G'MIC exception class.
//------------------------------
struct gmic_exception {
  char _message[16384];
  char _command[1024];
  gmic_exception() { *_message = 0; *_command = 0; }
  gmic_exception(const char *const command, const char *const message) {
    if (command) std::strcpy(_command,command); else *_command = 0;
    if (message) std::strcpy(_message,message); else *_message = 0;
  }
  const char *message() const { return _message; }
  const char *command() const { return _command; }
};

// Define G'MIC interpreter class.
//--------------------------------
struct gmic {

  // Internal variables.
#if cimg_display!=0
  cimg_library::CImgDisplay instant_window[10];
#endif
  gmic_list<char> command_names, command_definitions, scope, stack;
  float focale3d, light3d_x, light3d_y, light3d_z, specular_light3d, specular_shine3d, _progress, *progress;
  bool is_released, is_debug, is_fullpath, is_begin, is_end, is_quit, is_double3d, patch_borders, check_elif;
  int verbosity_level, render3d, renderd3d;
  volatile int _cancel, *cancel;
  unsigned char background3d[3];
  unsigned int position, patch_w, patch_h, patch_d, patch_c;
  char *tmpstr;

  // Constructors - Destructors.
  gmic(const char *const command_line, const char *const custom_commands=0, const bool default_commands=true,
       float *const p_progress=0, int *const p_cancel=0);
  template<typename T> gmic(const int argc, const char *const *const argv, gmic_list<T>& images,
                            const char *const custom_commands=0, const bool default_commands=true,
                            float *const p_progress=0, int *const p_cancel=0);
  template<typename T> gmic(const char *const command_line, gmic_list<T>& images,
                            const char *const custom_commands=0, const bool default_commands=true,
                            float *const p_progress=0, int *const p_cancel=0);
  ~gmic();
  gmic& assign(const char *const custom_commands=0, const bool default_commands=true,
               float *const p_progress=0, int *const p_cancel=0);

  // Messages procedures.
  gmic_image<char> scope2string() const;
  gmic_list<char> command_line_to_CImgList(const char *const command_line) const;
  template<typename T>
  const gmic& error(const gmic_list<T>& list, const char *format, ...) const;
  const gmic& error(const char *format, ...) const;
  template<typename T>
  const gmic& _arg_error(const gmic_list<T>& list, const char *const command, const char *const argument) const;
  template<typename T>
  const gmic& warning(const gmic_list<T>& list, const char *format, ...) const;
  const gmic& warning(const char *format, ...) const;
  template<typename T>
  const gmic& debug(const gmic_list<T>& list, const char *format, ...) const;
  const gmic& debug(const char *format, ...) const;
  template<typename T>
  const gmic& print(const gmic_list<T>& list, const char *format, ...) const;
  const gmic& print(const char *format, ...) const;

  // Add custom G'MIC commands.
  gmic& add_commands(const char *const data_commands);
  gmic& add_commands(std::FILE *const file);

  // Return image selection from a string.
  gmic_image<unsigned int> selection2cimg(const char *const string, const unsigned int indice_max,
                                          const char *const command, const bool is_selection) const;

  // Return stringified version of image selection.
  char *selection2string(const gmic_image<unsigned int>& selection,
                         const gmic_list<char>& filenames,
                         const bool display_selection) const;

  // Display image data.
  template<typename T>
  bool display_images(const gmic_list<T>& images,
                      const gmic_list<char>& filenames,
                      const gmic_image<unsigned int>& selection,
                      const bool verbose) const;
  template<typename T>
  bool display_objects3d(const gmic_list<T>& images,
                         const gmic_list<char>& filenames,
                         const gmic_image<unsigned int>& selection,
                         const bool verbose) const;
  template<typename T>
  bool display_plots(const gmic_list<T>& images,
                     const gmic_list<char>& filenames,
                     const gmic_image<unsigned int>& selection,
                     const unsigned int plot_type, const unsigned int vertex_type,
                     const double xmin, const double xmax,
                     const double ymin, const double ymax,
                     const bool verbose) const;

  // Substitute '@' and '{}' expressions.
  template<typename T>
  bool substitute_item(const char *const source, char *const destination, const gmic_list<T>& images,
                       const gmic_list<unsigned int>& repeatdones) const;

  // Main parsing procedure.
  template<typename T>
  gmic& parse(const gmic_list<char>& command_line, unsigned int& position,
              gmic_list<T> &images, gmic_list<char> &filenames,
              gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
              gmic_list<unsigned int>& locals, const bool initial_call);
  gmic& parse_bool(const gmic_list<char>& command_line, unsigned int& position,
                   gmic_list<bool>& images, gmic_list<char> &filenames,
                   gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
                   gmic_list<unsigned int>& locals, const bool initial_call);
  gmic& parse_uchar(const gmic_list<char>& command_line, unsigned int& position,
                    gmic_list<unsigned char>& images, gmic_list<char> &filenames,
                    gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
                    gmic_list<unsigned int>& locals, const bool initial_call);
  gmic& parse_char(const gmic_list<char>& command_line, unsigned int& position,
                   gmic_list<char>& images, gmic_list<char> &filenames,
                   gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
                   gmic_list<unsigned int>& locals, const bool initial_call);
  gmic& parse_ushort(const gmic_list<char>& command_line, unsigned int& position,
                     gmic_list<unsigned short>& images, gmic_list<char> &filenames,
                     gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
                     gmic_list<unsigned int>& locals, const bool initial_call);
  gmic& parse_short(const gmic_list<char>& command_line, unsigned int& position,
                    gmic_list<short>& images, gmic_list<char> &filenames,
                    gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
                    gmic_list<unsigned int>& locals, const bool initial_call);
  gmic& parse_uint(const gmic_list<char>& command_line, unsigned int& position,
                   gmic_list<unsigned int>& images, gmic_list<char> &filenames,
                   gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
                   gmic_list<unsigned int>& locals, const bool initial_call);
  gmic& parse_int(const gmic_list<char>& command_line, unsigned int& position,
                  gmic_list<int>& images, gmic_list<char> &filenames,
                  gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
                  gmic_list<unsigned int>& locals, const bool initial_call);
  gmic& parse_float(const gmic_list<char>& command_line, unsigned int& position,
                    gmic_list<float>& images, gmic_list<char> &filenames,
                    gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
                    gmic_list<unsigned int>& locals, const bool initial_call);
  gmic& parse_double(const gmic_list<char>& command_line, unsigned int& position,
                     gmic_list<double>& images, gmic_list<char> &filenames,
                     gmic_list<unsigned int>& dowhiles, gmic_list<unsigned int>& repeatdones,
                     gmic_list<unsigned int>& locals, const bool initial_call);

}; // End of the 'gmic' class.

#endif

// Local Variables:
// mode: c++
// End:
