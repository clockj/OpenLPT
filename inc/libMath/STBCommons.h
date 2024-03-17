#ifndef STBCOMMONS_H
#define STBCOMMONS_H

#include <iostream>
#include <string>

#define SAVEPRECISION 8

#define SMALLNUMBER 1e-8
#define SQRTSMALLNUMBER 1e-6
#define MAGSMALLNUMBER 1e-8

// log small
#define LOGSMALLNUMBER 1e-4

// undistort
#define UNDISTORT_MAX_ITER 50
#define UNDISTORT_EPS 1e-5

struct PixelRange 
{
    // left is closed, right is open 
    // [min, max)
    int row_min = 0;
    int row_max = 0;
    int col_min = 0;
    int col_max = 0;

    // Note: before using SetRowRange or SetColRange
    //       make sure it has been initialized!!!
    void setRowRange (int row)
    {
        if (row > row_max)
        {
            row_max = row;
        }
        else if (row < row_min)
        {
            row_min = row;
        }
    }; // TODO: we can add a choice to increase search area

    // Note: before using SetRowRange or SetColRange
    //       make sure it has been initialized!!!
    void setColRange (int col)
    {
        if (col > col_max)
        {
            col_max = col;
        }
        else if (col < col_min)
        {
            col_min = col;
        }
    }; // TODO: we can add a choice to increase search area

    // Note: before using SetRowRange or SetColRange
    //       make sure it has been initialized!!!
    void setRange (int row, int col)
    {
        setRowRange(row);
        setColRange(col);
    };
    int getNumOfRow () const
    {
        return row_max - row_min;
    };
    int getNumOfCol() const
    {
        return col_max - col_min;
    };
};

struct AxisLimit 
{
    double x_min = 0;
    double x_max = 0;
    double y_min = 0;
    double y_max = 0;
    double z_min = 0;
    double z_max = 0;

    AxisLimit () {};

    AxisLimit (double x1, double x2, double y1, double y2, double z1, double z2) 
        : x_min(x1), x_max(x2), y_min(y1), y_max(y2), z_min(z1), z_max(z2) {};

    void operator= (AxisLimit const& limit)
    {
        x_min = limit.x_min;
        x_max = limit.x_max;
        y_min = limit.y_min;
        y_max = limit.y_max;
        z_min = limit.z_min;
        z_max = limit.z_max;
    };

    bool check (double x, double y, double z)
    {
        if (x > x_max || x < x_min || 
            y > y_max || y < y_min || 
            z > z_max || z < z_min)
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
