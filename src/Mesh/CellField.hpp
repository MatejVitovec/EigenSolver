#ifndef CELLFIELD_HPP
#define CELLFIELD_HPP

#include <vector>
#include <memory>
#include <string>
#include <eigen3/Eigen/Dense>

#include "Mesh.hpp"

template <int N>
class CellField
{
    typedef Eigen::Array<double, Eigen::Dynamic, N> DataArray;

    public:
        CellField() {}

        CellField(const Mesh& mesh_) : mesh(mesh_)
        {
            data = DataArray(mesh.cellSize);
        }

        /*operator DataArray&()
        {
            return data;
        }*/

        DataArray& getData()
        {
            return data;
        }



    protected:
        const Mesh& mesh;        
        DataArray data;
        
};

#endif // CELLFIELD_HPP