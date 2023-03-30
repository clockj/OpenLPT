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
#include <deque> 

#include <tiffio.h>

#include "Matrix.h"
#include "STBCommons.h"


class ImageIO
{
protected:
    // tiff image properties: used when saved image 
    int _n_row = 0;            
    int _n_col = 0;            
    int _bit_per_sample = 0;   // maximum intensity of the image [bit]
    int _sample_per_pixel = 0; // number of component per pixel
    int _strip_size = 0;       // one strip size, used in compressed tiff image  
    int _strip_max = 0;        // number of strips, used in compressed tiff image 

    int _img_id = -1;          // current image index
    std::deque<std::string> _img_path; // path to all images of one camera

public:
    ImageIO () {};
    ImageIO (const ImageIO& img_io); // deep copy
    ImageIO (int n_row, int n_col);
    ~ImageIO () {};

    // Load path
    //  every filename should end with ';'
    // input: a main folder path, a file storing all the img paths 
    void LoadImgPath (std::string folder_path, std::string file_img_path);

    // Load Image
    // input: id_img (image index)
    // output: intensity matrix
    Matrix<int> LoadImg (int img_id);

    // Save Image
    void SaveImage (std::string save_path, Matrix<int>& intensity_mtx);

    // Set image info 
    void SetImgNumOfRow (int n_row) {_n_row = n_row;};
    void SetImgNumOfCol (int n_col) {_n_col = n_col;};
    void SetBitPerSample (int bit_per_sample) {_bit_per_sample = bit_per_sample;};
    void SetSamplePerPixel (int sample_per_pixel) {_sample_per_pixel = sample_per_pixel;};
    void SetStripSize (int strip_size) {_strip_size = strip_size;};
    void SetStripMax  (int strip_max)  {_strip_max = strip_max;};   

    // Get image info
    int GetNumOfImg ()     {return _img_path.size();};
    int GetImgNumOfRow ()  {return _n_row;};
    int GetImgNumOfCol ()  {return _n_col;};
    int GetBitPerSample () {return _bit_per_sample;};
    int GetSamplePerPixel () {return _sample_per_pixel;};
    int GetStripSize () {return _strip_size;};
    int GetStripMax  () {return _strip_max;};   

    int GetCurrentImgID () {return _img_id;};
};

#endif
