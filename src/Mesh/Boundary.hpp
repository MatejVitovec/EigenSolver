#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <vector>
#include <string>

class Boundary
{
    public:
        Boundary() : boundaryConditionName("empty") {};
        Boundary(std::string name) : boundaryConditionName(name) {};

        std::vector<int> indexToFaces;
        std::string boundaryConditionName;
};

#endif // BOUNDARY_HPP