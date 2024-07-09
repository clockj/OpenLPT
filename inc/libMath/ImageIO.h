//
//	ImageIO.h 
//
//	Base class for various image formats
//	All image format types should inherit from Image 
//
//	Created by Shijie Zhong 07/06/2022 
//

#ifndef IMAGEIO_H
#define IMAGEIO_H


#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring> 

#include <tiff.h>
#include <tiffio.h>

#include "Matrix.h"
#include "STBCommons.h"


#define IMAGEIO_CHECK_CALL(call) \
    {if (call==false) {std::cerr << "ImageIO Error (line: " << __LINE__  << ")" << std::endl; throw error_io;}} 
#define IMAGEIO_CHECK_CALL_DEBUG(call) \
    {if (call==false) {std::cerr << "ImageIO Optional Call Fail (line: " << __LINE__ << ")" << std::endl;}}

#define TILE_MAX_WIDTH  (1 << 24)
#define TILE_MAX_HEIGHT (1 << 24)
#define MAX_TILE_SIZE   (1 << 30)
#define BITS_PER_BYTE    8

#ifndef __IPL_H__
   typedef unsigned char uchar;
   typedef unsigned short ushort;
#endif

struct ImageParam
{
    int n_row;
    int n_col;
    int bits_per_sample;
    int n_channel;
};

class ImageIO
{
private:
    // tiff image properties: used when saved image 
    int _n_row = 0;            // image height       
    int _n_col = 0;            // image width       
    int _bits_per_sample = 0;  // maximum intensity of the image [bit]
    int _n_channel = 0;        // number of channels

    int _is_tiled = 0;         // 1 for tiled; 0 for stripped
    int _tile_height0 = 0;     // tile height (length)
    int _tile_width0 = 0;      // tile width
    int _img_orientation = ORIENTATION_TOPLEFT; // image orientation

    int _img_id = -1;          // current image index
    std::vector<std::string> _img_path; // path to all images of one camera

    void init ();

public:
    ImageIO () {};
    ImageIO (std::string folder_path, std::string file_img_path);
    ImageIO (const ImageIO& img_io); // deep copy
    ~ImageIO () {};

    // Load path
    //  every filename should end with ';'
    // input: a main folder path, a file storing all the img paths 
    void loadImgPath (std::string folder_path, std::string file_img_path);

    // Load Image
    // input: id_img (image index)
    // output: intensity matrix
    Image loadImg (int img_id);

    // Save Image
    void saveImg (std::string save_path, Image const& image);

    // Set image info 
    void setImgParam (ImageParam const& img_param);

    // Get image info
    ImageParam getImgParam () const;
      
    // Get current image index
    int getCurrImgID () const {return _img_id;};

    // Get image path
    std::vector<std::string> getImgPath () const {return _img_path;};
};

#endif
