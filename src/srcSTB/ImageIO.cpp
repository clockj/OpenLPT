#include "ImageIO.h"


ImageIO::ImageIO (const ImageIO& img)
    : _n_row(img._n_row), _n_col(img._n_col), _bit_per_sample(img._bit_per_sample), _sample_per_pixel(img._sample_per_pixel), _strip_size(img._strip_size), _strip_max(img._strip_max), _img_id(img._img_id), _img_path(img._img_path)
{}

ImageIO::ImageIO (int n_row, int n_col)
    : _n_row(n_row), _n_col(n_col)
{}


void ImageIO::LoadImgPath (std::string folder_path, std::string file_img_path)
{
    std::ifstream infile(folder_path + file_img_path, std::ios::in);
    std::string line;

    while (std::getline(infile, line)) {
        // std::cout << line << std::endl;
        // std::cout << '"' << line[line.size()-1] << '"' << std::endl;
        size_t comment_pos = line.find(';');
        if (comment_pos > 0) 
        {
            if (comment_pos < std::string::npos) 
            {
                line.erase(comment_pos);
            }
            _img_path.push_back(folder_path + line);
            // std::cout << _img_path[_img_path.size()-1] << std::endl;
        }
    }
    infile.close();

    if (_img_path.size() == 0)
    {
        std::cout << "There is no image path loaded." << std::endl;
    }
}


Matrix<int> ImageIO::LoadImg (int img_id)
{
    if (img_id >= _img_path.size())
    {
        std::cerr << "Image id: " << img_id 
                  << " is larger than total number of image: " 
                  << _img_path.size()
                  << std::endl;
        throw;
    }
    _img_id = img_id;
    std::string file = _img_path[img_id];
    
    TIFF* image;
    if ((image = TIFFOpen(file.c_str(), "r")) == NULL) 
    {
        throw std::invalid_argument("Could not open image!");
    }


    TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &_bit_per_sample);
    TIFFGetField(image, TIFFTAG_SAMPLESPERPIXEL, &_sample_per_pixel);
    // std::cout << "Bit_per_sample: "   << _bit_per_sample   << std::endl;
    // std::cout << "Sample_per_pixel: " << _sample_per_pixel << std::endl;
    
    // Find image size:
    //  _strip_size, _strip_max 
    tsize_t strip_size = TIFFStripSize(image);
    int strip_max = TIFFNumberOfStrips(image); 
    _strip_size = strip_size;
    _strip_max  = strip_max;

    // Find the number of rows and number of columns
    //  _n_col & _n_row are different from 
    //  stripMax & stripSize!
    if (TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &_n_col) == 0) 
    {
        throw std::invalid_argument("Tiff image does not define its width");
    }
    if (TIFFGetField(image, TIFFTAG_IMAGELENGTH, &_n_row) == 0) 
    {
        throw std::invalid_argument("Tiff image does not define its length");
    }
    // std::cout << "StripMax: "  << strip_max  << ", n_row: " << _n_row << std::endl; 
    // std::cout << "StripSize: " << strip_size << ", n_col: " << _n_col << std::endl;

    unsigned long buffer_size = strip_max * strip_size;
    unsigned long image_offset = 0;
    unsigned char* buffer;
    buffer = new unsigned char [buffer_size];
    
    for (int strip_count = 0; strip_count < strip_max; strip_count ++) 
    {
        long result = TIFFReadEncodedStrip(image, strip_count, buffer+image_offset, strip_size);
        if (result == -1) 
        {
            std::cerr << "Read error for tiff image" << std::endl;
            throw;
        }
        image_offset += result;
    }

    Matrix<int> intensity_mtx(_n_row, _n_col, -1);
    for (int iter_x = 0; iter_x < _n_row; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _n_col; iter_y ++)
        {
            intensity_mtx(iter_x, iter_y) 
            = static_cast<int>(buffer[_n_col * iter_x + iter_y]);
        }
    }

    TIFFClose(image);
    delete [] buffer;

    return intensity_mtx;
}


void ImageIO::SaveImage (std::string save_path, Matrix<int>& intensity_mtx)
{
    // if (_img_id < 0)
    // {
    //     std::cerr << "No image loaded yet!" 
    //               << "The image size is unknow."
    //               << std::endl;
    //     throw;
    // }

    int buffer_size = _strip_max * _strip_size;
    unsigned char* buffer;
    buffer = new unsigned char [buffer_size];

    for (int iter_x = 0; iter_x < _n_row; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _n_col; iter_y ++)
        {
            buffer[_n_col * iter_x + iter_y]
            = static_cast<unsigned char> (intensity_mtx(iter_x, iter_y));
        }
    }

    TIFF* image = TIFFOpen(save_path.c_str(), "w");

    TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, _bit_per_sample);
    TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, _sample_per_pixel);
    TIFFSetField(image, TIFFTAG_IMAGELENGTH, _n_row);
    TIFFSetField(image, TIFFTAG_IMAGEWIDTH, _n_col);
    TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, _n_row/_strip_max);
    TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK); // or PHOTOMETRIC_MINISWHITE

    // tsize_t strip_size = TIFFStripSize(image);
    // int strip_max = TIFFNumberOfStrips(image);
    // std::cout << "(strip_size, strip_max) = " 
    //           << "(" << strip_size << ", " 
    //           << strip_max << ")" 
    //           << std::endl;

    unsigned long image_offset = 0;
    for (int strip_count = 0; strip_count < _strip_max; strip_count ++) 
    {
        long result = TIFFWriteEncodedStrip(image, strip_count, buffer+image_offset, _strip_size);
        if (result == -1) 
        {
            std::cerr << "Write error for tiff image" << std::endl;
            std::cerr << "current strip_count: " << strip_count << std::endl;
            throw;
        }
        image_offset += result;
    }

    TIFFClose(image);

    delete[] buffer;
}
