#ifndef OBJECTFINDER_HPP
#define OBJECTFINDER_HPP

#include "ObjectFinder.h"

void ObjectFinder2D::findTracer2D
(std::vector<Tracer2D>& tr2d_list, Image const& img, int max_intensity, int min_intensity)
{
    // Judge whether a given point is local maximum intensity or not
    //  skip first and last rows and columns   
    for (int row = 1; row < img.getDimRow()-1; row ++)
    {
        for (int col = 1; col < img.getDimCol()-1; col ++)
        {
            if (img(row, col) > max_intensity)
            {
                // the intensity is out of range 
                std::cerr << "ObjectFinder2D::findTracer2D error: "
                          << "The intensity at (row,col)=(" << row << "," << col << ") is "
                          << img(row, col) << ". " 
                          << "It is larger than max intensity " 
                          << max_intensity << std::endl;
                throw error_range;
            }

            if (img(row, col) >= min_intensity && myMATH::isLocalMax(img, row, col))
            {
                // use Gaussian distribution to find the maximum intensity
                //  => particle center position

                // Note: row => y 
                //       col => x 
                int x1 = (col - 1);
                int x2 =  col;     
                int x3 = (col + 1);
                int y1 = (row - 1);
                int y2 =  row;     
                int y3 = (row + 1);

                double ln_z1 = 0.0;
                double ln_z2 = 0.0;
                double ln_z3 = 0.0;

                // find the col value (coordinate: x)
                ln_z1 = img(y2, x1) < INTSMALLNUMBER ? std::log(INTSMALLNUMBER) : std::log(img(y2, x1));
                ln_z2 = img(y2, x2) < INTSMALLNUMBER ? std::log(INTSMALLNUMBER) : std::log(img(y2, x2));
                ln_z3 = img(y2, x3) < INTSMALLNUMBER ? std::log(INTSMALLNUMBER) : std::log(img(y2, x3));
    
                double xc = -0.5 * (  (ln_z1 * double((x2 * x2) - (x3 * x3))) 
                                    - (ln_z2 * double((x1 * x1) - (x3 * x3))) 
                                    + (ln_z3 * double((x1 * x1) - (x2 * x2))) ) 
                                 / (  (ln_z1 * double(x3 - x2)) 
                                    - (ln_z3 * double(x1 - x2)) 
                                    + (ln_z2 * double(x1 - x3)) );
                if (!std::isfinite(xc)) 
                {
                    continue;
                }

                // find the row value (coordinate: y)
                ln_z1 = img(y1, x2) < INTSMALLNUMBER ? std::log(INTSMALLNUMBER) : std::log(img(y1, x2));
                ln_z3 = img(y3, x2) < INTSMALLNUMBER ? std::log(INTSMALLNUMBER) : std::log(img(y3, x2));

                double yc = -0.5 * (  (ln_z1 * double((y2 * y2) - (y3 * y3))) 
                                    - (ln_z2 * double((y1 * y1) - (y3 * y3))) 
                                    + (ln_z3 * double((y1 * y1) - (y2 * y2))) ) 
                                 / (  (ln_z1 * double(y3 - y2)) 
                                    - (ln_z3 * double(y1 - y2)) 
                                    + (ln_z2 * double(y1 - y3)) );
                if (!std::isfinite(yc)) 
                {
                    continue;
                }

                tr2d_list.push_back(Tracer2D(Pt2D(xc, yc)));
            }
        }
    }
}

template<class T>
void ObjectFinder2D::findObject2D
(std::vector<T>& obj2d_list, Image const& img, std::vector<double> const& properties)
{
    if (typeid(T) == typeid(Tracer2D))
    {
        if (obj2d_list.size() > 0)
        {
            obj2d_list.clear();
        }

        findTracer2D(obj2d_list, img, properties[0], properties[1]);
    }
    else
    {
        std::cerr << "ObjectFinder: " 
                  << "class " << typeid(T).name()
                  << "is not included in ObjectFinder!" << std::endl;
        throw error_type;
    }
}

#endif 