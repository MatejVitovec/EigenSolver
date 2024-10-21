#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <fenv.h>

#include <eigen3/Eigen/Dense>

#include "Fields/Field.hpp"
#include "Mesh/Mesh.hpp"
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

      /*Field<3,1> pole = Field<3,1>(12);
      Field<3,1> pole2 = Field<3,1>(12);

      Eigen::Array3d u(-1, 1, 4);

      pole[2] += u;
      auto res = pole + pole2;



      std::cout << res << std::endl;*/

      Mesh<3> mesh = Mesh<3>();

      mesh.loadGmsh2("../Meshes/GAMM.msh");

      //std::cout << mesh.nodes.getData()(0,0) << std::endl;

      return 0;
}