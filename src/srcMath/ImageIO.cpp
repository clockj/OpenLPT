#include "ImageIO.h"


ImageIO::ImageIO (const ImageIO& img)
    : _n_row(img._n_row), _n_col(img._n_col), _bits_per_sample(img._bits_per_sample), _n_channel(img._n_channel), _is_tiled(img._is_tiled), _tile_height0(img._tile_height0), _tile_width0(img._tile_width0), _img_orientation(img._img_orientation), _img_id(img._img_id), _img_path(img._img_path)
{}


void ImageIO::init ()
{
    _n_row = 0;
    _n_col = 0;
    _bits_per_sample = 0;
    _n_channel = 0;
    _is_tiled = 0;
    _tile_height0 = 0;
    _tile_width0 = 0;
    _img_orientation = ORIENTATION_TOPLEFT;
    _img_id = -1;
    _img_path.clear();
}


void ImageIO::loadImgPath (std::string folder_path, std::string file_img_path)
{
    std::ifstream infile(folder_path + file_img_path, std::ios::in);

    if (!infile.is_open())
    {
        std::cerr << "ImageIO::LoadImgPath: Could not open image path file!" << std::endl;
        throw error_io;
    }

    // Initialize image path
    init();

    std::string line;
    while (std::getline(infile, line)) 
    {
        _img_path.push_back(folder_path + line);
    }
    infile.close();

    if (_img_path.size() == 0)
    {
        std::cout << "There is no image path loaded." << std::endl;
    }
}


Image ImageIO::loadImg (int img_id)
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
    
    TIFF* tif;
    if ((tif = TIFFOpen(file.c_str(), "r")) == NULL) 
    {
        std::cerr << "ImageIO::LoadImg: Could not open image!" << std::endl;
        throw error_io;
    }

    // check is the image is colorful
    _n_channel = 1;

    #ifdef DEBUG
    IMAGEIO_CHECK_CALL_DEBUG(TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &_n_channel));
    #endif
    if (_n_channel != 1)
    {
        std::cerr << "ImageIO::LoadImg: current version not supported for colorful image! " 
                  << "Line: " << __LINE__
                  << std::endl;
        throw error_io;
    }

    // check image size
    IMAGEIO_CHECK_CALL(TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &_n_row));
    IMAGEIO_CHECK_CALL(TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &_n_col));
    IMAGEIO_CHECK_CALL((_n_row > 0 && _n_col > 0));
    
    // load image bits
    IMAGEIO_CHECK_CALL(TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &_bits_per_sample));
    IMAGEIO_CHECK_CALL((_bits_per_sample==8 || _bits_per_sample==16 || _bits_per_sample==32 || _bits_per_sample==64));

    // check is the image is tiled or stripped
    _is_tiled = TIFFIsTiled(tif) != 0;
    bool is_read_scanline = false;
    if (_is_tiled)
    {
        IMAGEIO_CHECK_CALL(TIFFGetField(tif, TIFFTAG_TILELENGTH, &_tile_height0));
        IMAGEIO_CHECK_CALL(TIFFGetField(tif, TIFFTAG_TILEWIDTH, &_tile_width0));
        IMAGEIO_CHECK_CALL((_tile_height0>0 && _tile_height0<=TILE_MAX_HEIGHT && _tile_width0>0 && _tile_width0<=TILE_MAX_WIDTH));

        // For debug
        // if (_tile_height0 > TILE_MAX_HEIGHT || _tile_width0 > TILE_MAX_WIDTH || _tile_height0 <= 0 || _tile_width0 <= 0)
        // {
        //     std::cerr << "ImageIO::LoadImg: Image tile size is out of range! " 
        //               << "_tile_height0: " << _tile_height0 << " "
        //               << "_tile_width0: " << _tile_width0 << " "
        //               << "Line: " << __LINE__
        //               << std::endl;
        //     throw error_io;
        // }
    }
    else 
    {
        _tile_width0 = _n_col;

        if (_bits_per_sample <= 32)
        {
            is_read_scanline = true;
            _tile_height0 = 1; // read each scanline
        }
        else 
        {
            _tile_height0 = _n_row;
        }
    }
    
    // create buffer to store image data
    const size_t buffer_bytes_per_row = (_n_channel*_tile_width0*_bits_per_sample + BITS_PER_BYTE - 1) / BITS_PER_BYTE;
    const size_t buffer_size = _tile_height0 * buffer_bytes_per_row;
    if (buffer_size > MAX_TILE_SIZE)
    {
        std::cerr << "ImageIO::LoadImg: Image buffer size is out of range! " 
                  << "buffer_size: " << buffer_size << " "
                  << "Line: " << __LINE__
                  << std::endl;
        throw error_io;
    }
    uchar* buffer = new uchar[buffer_size]; 
    
    // check image orientation
    _img_orientation = ORIENTATION_TOPLEFT;
    #ifdef DEBUG
    IMAGEIO_CHECK_CALL_DEBUG(TIFFGetField(tif, TIFFTAG_ORIENTATION, &_img_orientation));
    #endif
    bool vert_flip = _img_orientation == ORIENTATION_BOTLEFT || _img_orientation == ORIENTATION_BOTRIGHT || _img_orientation == ORIENTATION_LEFTBOT || _img_orientation == ORIENTATION_RIGHTBOT;

    // initialize data array to convert buffer
    const size_t data_num_per_row = _n_channel * _tile_width0 * _tile_height0;
    uint8*  buffer_8;
    uint16* buffer_16;
    uint32* buffer_32;
    uint64* buffer_64;
    switch (_bits_per_sample)
    {
    case 8:
        buffer_8 = new uint8[data_num_per_row];
        buffer_16 = nullptr;
        buffer_32 = nullptr;
        buffer_64 = nullptr;
        break;
    case 16:
        buffer_8 = nullptr;
        buffer_16 = new uint16[data_num_per_row];
        buffer_32 = nullptr;
        buffer_64 = nullptr;
        break;
    case 32:
        buffer_8 = nullptr;
        buffer_16 = nullptr;
        buffer_32 = new uint32[data_num_per_row];
        buffer_64 = nullptr;
        break;
    case 64:
        buffer_8 = nullptr;
        buffer_16 = nullptr;
        buffer_32 = nullptr;
        buffer_64 = new uint64[data_num_per_row];
        break;
    default:
        std::cerr << "ImageIO::LoadImg: Image bits per sample is out of range! " 
                  << "_bits_per_sample: " << _bits_per_sample << " "
                  << "Line: " << __LINE__
                  << std::endl;
        break;
    }

    // read image data
    Image image(_n_row, _n_col, -1);
    int tile_id = 0;
    for (int row = 0; row < _n_row; row += _tile_height0)
    {
        int tile_height = std::min(_tile_height0, _n_row - row);
        const int img_row = vert_flip ? _n_row-1-row : row; 

        for (int col = 0; col < _n_col; col += _tile_width0, tile_id ++)
        {
            int tile_width = std::min(_tile_width0, _n_col - col);

            if (_is_tiled)
            {
                IMAGEIO_CHECK_CALL((TIFFReadEncodedTile(tif, tile_id, buffer, buffer_size) >= 0));
            }
            else if (is_read_scanline)
            {
                IMAGEIO_CHECK_CALL((TIFFReadScanline(tif, (uint32*)buffer, row) >= 0));
            }
            else
            {
                IMAGEIO_CHECK_CALL((TIFFReadEncodedStrip(tif, tile_id, buffer, buffer_size) >= 0));
            }

            switch (_bits_per_sample)
            {
            case 8:
                std::memcpy(buffer_8, buffer, buffer_size);
                for (int i = 0; i < tile_height; i ++)
                {
                    for (int j = 0; j < tile_width; j ++)
                    {
                        image(img_row+i, col+j) = (double) buffer_8[i*tile_width+j];
                    }
                }
                break;
            case 16:
                std::memcpy(buffer_16, buffer, buffer_size);
                for (int i = 0; i < tile_height; i ++)
                {
                    for (int j = 0; j < tile_width; j ++)
                    {
                        image(img_row+i, col+j) = (double) buffer_16[i*tile_width+j];
                    }
                }
                break;
            case 32:
                std::memcpy(buffer_32, buffer, buffer_size);
                for (int i = 0; i < tile_height; i ++)
                {
                    for (int j = 0; j < tile_width; j ++)
                    {
                        image(img_row+i, col+j) = (double) buffer_32[i*tile_width+j];
                    }
                }
                break;
            case 64:
                std::memcpy(buffer_64, buffer, buffer_size);
                for (int i = 0; i < tile_height; i ++)
                {
                    for (int j = 0; j < tile_width; j ++)
                    {
                        image(img_row+i, col+j) = (double) buffer_64[i*tile_width+j];
                    }
                }
                break;
            default:
                break;
            } 
        }
    }
    
    switch (_bits_per_sample)
    {
    case 8:
        delete[] buffer_8;
        break;
    case 16:
        delete[] buffer_16;
        break;
    case 32:
        delete[] buffer_32;
        break;
    case 64:
        delete[] buffer_64;
        break;
    default:
        break;
    }

    _TIFFfree(buffer);
    TIFFClose(tif);

    return image;
}


void ImageIO::saveImg (std::string save_path, Image const& image)
{
    // check image size
    IMAGEIO_CHECK_CALL((_n_row>0 && _n_col>0 && _n_row==image.getDimRow() && _n_col==image.getDimCol()));

    TIFF* tif = TIFFOpen(save_path.c_str(), "w");

    switch (_bits_per_sample)
    {
    case 8:
        {
            uint8* buffer_8 = (uint8*) _TIFFmalloc(_n_row*_n_col * sizeof (uint8));
            for (int iter_x = 0; iter_x < _n_row; iter_x ++)
            {
                for (int iter_y = 0; iter_y < _n_col; iter_y ++)
                {
                    buffer_8[_n_col*iter_x + iter_y] = static_cast<uint8> (image(iter_x, iter_y));
                }
            }

            TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, _bits_per_sample);
            TIFFSetField(tif, TIFFTAG_IMAGELENGTH, _n_row);
            TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, _n_col);
            TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, _n_channel);
            TIFFSetField(tif, TIFFTAG_ORIENTATION, _img_orientation);
            TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

            for (int iter_x = 0; iter_x < _n_row; iter_x ++)
            {
                TIFFWriteScanline(tif, &buffer_8[_n_row*iter_x], iter_x, 0);
            }

            _TIFFfree(buffer_8);
            break;
        }
    
    case 16:
        {
            uint16* buffer_16 = (uint16*) _TIFFmalloc(_n_row*_n_col * sizeof (uint16));
            for (int iter_x = 0; iter_x < _n_row; iter_x ++)
            {
                for (int iter_y = 0; iter_y < _n_col; iter_y ++)
                {
                    buffer_16[_n_col*iter_x + iter_y] = static_cast<uint16> (image(iter_x, iter_y));
                }
            }

            TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, _bits_per_sample);
            TIFFSetField(tif, TIFFTAG_IMAGELENGTH, _n_row);
            TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, _n_col);
            TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, _n_channel);
            TIFFSetField(tif, TIFFTAG_ORIENTATION, _img_orientation);
            TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

            for (int iter_x = 0; iter_x < _n_row; iter_x ++)
            {
                TIFFWriteScanline(tif, &buffer_16[_n_row*iter_x], iter_x, 0);
            }

            _TIFFfree(buffer_16);
            break;
        }

    case 32:
        {
            uint32* buffer_32 = (uint32*) _TIFFmalloc(_n_row*_n_col * sizeof (uint32));
            for (int iter_x = 0; iter_x < _n_row; iter_x ++)
            {
                for (int iter_y = 0; iter_y < _n_col; iter_y ++)
                {
                    buffer_32[_n_col*iter_x + iter_y] = static_cast<uint32> (image(iter_x, iter_y));
                }
            }

            TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, _bits_per_sample);
            TIFFSetField(tif, TIFFTAG_IMAGELENGTH, _n_row);
            TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, _n_col);
            TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, _n_channel);
            TIFFSetField(tif, TIFFTAG_ORIENTATION, _img_orientation);
            TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

            for (int iter_x = 0; iter_x < _n_row; iter_x ++)
            {
                TIFFWriteScanline(tif, &buffer_32[_n_row*iter_x], iter_x, 0);
            }

            _TIFFfree(buffer_32);
            break;
        }

    case 64:
        {
            uint64* buffer_64 = (uint64*) _TIFFmalloc(_n_row*_n_col * sizeof (uint64));
            for (int iter_x = 0; iter_x < _n_row; iter_x ++)
            {
                for (int iter_y = 0; iter_y < _n_col; iter_y ++)
                {
                    buffer_64[_n_col*iter_x + iter_y] = static_cast<uint64> (image(iter_x, iter_y));
                }
            }

            TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, _bits_per_sample);
            TIFFSetField(tif, TIFFTAG_IMAGELENGTH, _n_row);
            TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, _n_col);
            TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, _n_channel);
            TIFFSetField(tif, TIFFTAG_ORIENTATION, _img_orientation);
            TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

            for (int iter_x = 0; iter_x < _n_row; iter_x ++)
            {
                TIFFWriteScanline(tif, &buffer_64[_n_row*iter_x], iter_x, 0);
            }

            _TIFFfree(buffer_64);
            break;
        }

    default:
        std::cerr << "ImageIO::saveImg: Image bits per sample is out of range! " 
                  << "_bits_per_sample: " << _bits_per_sample << " "
                  << "Line: " << __LINE__
                  << std::endl;
        break;
    }

    TIFFClose(tif);
}


void ImageIO::setImgParam (ImageParam const& img_param)
{
    _n_row = img_param.n_row;
    _n_col = img_param.n_col;
    _bits_per_sample = img_param.bits_per_sample;
    _n_channel = img_param.n_channel;
}


ImageParam ImageIO::getImgParam () const
{
    ImageParam img_param;
    img_param.n_row = _n_row;
    img_param.n_col = _n_col;
    img_param.bits_per_sample = _bits_per_sample;
    img_param.n_channel = _n_channel;

    return img_param;
}