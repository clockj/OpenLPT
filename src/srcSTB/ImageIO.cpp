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
    if (img_id >= int(_img_path.size()))
    {
        std::cerr << "Image id: " << img_id 
                  << " is larger than total number of image: " 
                  << _img_path.size()
                  << std::endl;
        throw error_size;
    }
    _img_id = img_id;
    std::string file = _img_path[img_id];
    
    TIFF* image;
    if ((image = TIFFOpen(file.c_str(), "r")) == NULL) 
    {
        std::cerr << "ImageIO::LoadImg: Could not open image!" << std::endl;
        throw error_io;
    }

    TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &_bit_per_sample);
    if (_bit_per_sample != 8)
    {
        std::cerr << "ImageIO::LoadImg: Could not open image with pixel bit = " 
                  << _bit_per_sample << std::endl;
        throw error_range;
    }

    // Find the number of rows and number of columns
    //  _n_col & _n_row are different from 
    //  stripMax & stripSize!
    if (TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &_n_col) == 0) 
    {
        std::cerr << "ImageIO::LoadImg: Tiff image does not define its width" 
                  << std::endl;
        throw error_io;
    }
    if (TIFFGetField(image, TIFFTAG_IMAGELENGTH, &_n_row) == 0) 
    {
        std::cerr << "ImageIO::LoadImg: Tiff image does not define its length"
                  << std::endl;
        throw error_io;
    }

    // int n_frame = TIFFNumberOfDirectories(image);
    // TIFFSetDirectory(image, 0);

    uint32* buffer = (uint32*) _TIFFmalloc(_n_row*_n_col * sizeof (uint32));
    TIFFReadRGBAImage(image, _n_col, _n_row, buffer, 0);

    Matrix<int> intensity_mtx(_n_row, _n_col, -1);
    for (int iter_x = 0; iter_x < _n_row; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _n_col; iter_y ++)
        {
            intensity_mtx(_n_row-1-iter_x, iter_y) = (int) TIFFGetG(buffer[iter_x*_n_col + iter_y]);
        }
    }

    _TIFFfree(buffer);
    TIFFClose(image);

    return intensity_mtx;
}


void ImageIO::SaveImage (std::string save_path, Matrix<int>& intensity_mtx)
{
    if (_img_id < 0)
    {
        std::cerr << "No image loaded yet!" 
                  << "The image size is unknow."
                  << std::endl;
        throw error_io;
    }


    uint8* buffer = (uint8*) _TIFFmalloc(_n_row*_n_col * sizeof (uint8));
    for (int iter_x = 0; iter_x < _n_row; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _n_col; iter_y ++)
        {
            buffer[_n_col*iter_x + iter_y] = static_cast<uint8> (intensity_mtx(iter_x, iter_y));
        }
    }

    TIFF* image = TIFFOpen(save_path.c_str(), "w");

    TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, _bit_per_sample);
    TIFFSetField(image, TIFFTAG_IMAGELENGTH, _n_row);
    TIFFSetField(image, TIFFTAG_IMAGEWIDTH, _n_col);
    TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

    for (int iter_x = 0; iter_x < _n_row; iter_x ++)
    {
        TIFFWriteScanline(image, &buffer[_n_row*iter_x], iter_x, 0);
    }

    _TIFFfree(buffer);
    TIFFClose(image);

}
