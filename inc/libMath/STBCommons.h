#ifndef STBCOMMONS_H
#define STBCOMMONS_H

#include <iostream>
#include <string>

// calibration
#define SMALLNUMBER 1e-8
// #define SMALLNUMBER 2e-6
#define MAGSMALLNUMBER 1e-8

// intensity
#define INTSMALLNUMBER 1e-6

struct PixelRange 
{
    // left is closed, right is open 
    // [min, max)
    int _row_min = 0;
    int _row_max = 0;
    int _col_min = 0;
    int _col_max = 0;

    // Note: before using SetRowRange or SetColRange
    //       make sure it has been initialized!!!
    void SetRowRange (int row)
    {
        if (row > _row_max)
        {
            _row_max = row;
        }
        else if (row < _row_min)
        {
            _row_min = row;
        }
    }; // TODO: we can add a choice to increase search area

    // Note: before using SetRowRange or SetColRange
    //       make sure it has been initialized!!!
    void SetColRange (int col)
    {
        if (col > _col_max)
        {
            _col_max = col;
        }
        else if (col < _col_min)
        {
            _col_min = col;
        }
    }; // TODO: we can add a choice to increase search area

    // Note: before using SetRowRange or SetColRange
    //       make sure it has been initialized!!!
    void SetRange (int row, int col)
    {
        SetRowRange(row);
        SetColRange(col);
    };
    int GetNumOfRow ()
    {
        return _row_max - _row_min;
    }
    int GetNumOfCol()
    {
        return _col_max - _col_min;
    }
};

struct AxisLimit 
{
    double _x_min = 0;
    double _x_max = 0;
    double _y_min = 0;
    double _y_max = 0;
    double _z_min = 0;
    double _z_max = 0;

    void operator= (AxisLimit const& limit)
    {
        _x_min = limit._x_min;
        _x_max = limit._x_max;
        _y_min = limit._y_min;
        _y_max = limit._y_max;
        _z_min = limit._z_min;
        _z_max = limit._z_max;
    };

    bool check (double x, double y, double z)
    {
        if (x > _x_max || x < _x_min || 
            y > _y_max || y < _y_min || 
            z > _z_max || z < _z_min)
        {
            return false;
        }

        return true;
    }
};

enum ErrorTypeID
{
    error_size = 1,
    error_type,
    error_range,
    error_space,
    error_io,
    error_div0,
    error_parallel
};

enum ObjectTypeID
{
    type_tracer,
    type_bubble,
    type_filament
};

enum FrameTypeID
{
    PREV_FRAME,
    CURR_FRAME
};

#endif // !STBCOMMONS
