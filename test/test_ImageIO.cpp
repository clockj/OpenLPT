#include<iostream>

#include "test.h"
#include "Matrix.h"
#include "ImageIO.h"

// test load function
bool test_function_1 ()
{
    ImageIO img_io;

    img_io.loadImgPath("../test/inputs/test_ImageIO/", "test_function_1.txt");
    Image img_8bits = img_io.loadImg(0);
    img_8bits.write("../test/results/test_ImageIO/test_8bits.csv");
    Image img_ans_8bits("../test/solutions/test_ImageIO/test_8bits.csv");
    if (img_8bits != img_ans_8bits)
    {
        std::cout << "img_8bits != img_ans_8bits" << std::endl;
        return false;
    }

    Image img_16bits = img_io.loadImg(1);
    img_16bits.write("../test/results/test_ImageIO/test_16bits.csv");
    Image img_ans_16bits("../test/solutions/test_ImageIO/test_16bits.csv");
    if (img_16bits != img_ans_16bits)
    {
        std::cout << "img_16bits != img_ans_16bits" << std::endl;
        return false;
    }

    return true;
}

// test save function
bool test_function_2 ()
{
    ImageIO img_io;
    ImageIO img_io_res;
    

    img_io.loadImgPath("../test/inputs/test_ImageIO/", "test_function_1.txt");
    img_io_res.loadImgPath("../test/", "inputs/test_ImageIO/test_function_2.txt");

    Image img_8bits = img_io.loadImg(0);
    img_io.saveImg("../test/results/test_ImageIO/test_8bits_save.tif", img_8bits);
    Image img_res_8bits = img_io_res.loadImg(0);
    if (img_8bits != img_res_8bits)
    {
        std::cout << "img_8bits != img_res_8bits" << std::endl;
        return false;
    }
    
    Image img_16bits = img_io.loadImg(1);
    img_io.saveImg("../test/results/test_ImageIO/test_16bits_save.tif", img_16bits);
    Image img_res_16bits = img_io_res.loadImg(1);
    if (img_16bits != img_res_16bits)
    {
        std::cout << "img_16bits != img_res_16bits" << std::endl;
        return false;
    }

    return true;
}

// test set/get parameter function
bool test_function_3 ()
{
    ImageIO img_io;
    ImageParam img_param;
    img_param.n_row = 100;
    img_param.n_col = 100;
    img_param.bits_per_sample = 8;
    img_param.n_channel = 3;
    img_io.setImgParam(img_param);
    ImageParam img_param_res = img_io.getImgParam();
    if (img_param.n_row != img_param_res.n_row || img_param.n_col != img_param_res.n_col || img_param.bits_per_sample != img_param_res.bits_per_sample || img_param.n_channel != img_param_res.n_channel)
    {
        std::cout << "img_param != img_param_res" << std::endl;
        return false;
    }
    return true;
}

int main ()
{
    fs::create_directory("../test/results/test_ImageIO/");

    IS_TRUE(test_function_1());
    IS_TRUE(test_function_2());
    IS_TRUE(test_function_3());

    return 0;
}