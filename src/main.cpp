#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <fenv.h>

#include <eigen3/Eigen/Dense>

#include "Fields/Field.hpp"
//#include "Mesh/Mesh.hpp"
//#include "Mesh/CellField.hpp"


int main(int argc, char** argv)
{
    //feenableexcept(FE_INVALID | FE_OVERFLOW);

    /*Eigen::ArrayXXf m(4,4);
    m <<  1, 2, 3, 4,
          5, 6, 7, 8,
          9,10,11,12,
          13,14,15,16;
    std::cout << "Block in the middle" << std::endl;
    std::cout << m.block<2,2>(1,0) << std::endl << std::endl;

    Eigen::ArrayXXf n(4,4);
    n <<  5, 6, 7, 8,
          5, 6, 7, 8,
          9,10,11,12,
          13,14,15,16;
    auto res = m+n;*/

    Eigen::Vector2d vec;

    Field<3,3> pole = Field<3,3>(12);
    Field<3,3> pole2 = Field<3,3>(12);

    pole[2] += 2;
    pole2[2] += 3;

    auto res = pole + pole2;

    Eigen::Vector2d u(-1,1), v(2,0);


    std::cout << u.array()*v.array() << std::endl << std::endl;


    return 0;
}