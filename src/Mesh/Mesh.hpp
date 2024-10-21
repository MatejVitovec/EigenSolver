#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <memory>
#include <string>
#include <list>

#include "Boundary.hpp"
#include "../Fields/Field.hpp"

template <int DIM>
class Mesh
{
    public:
        enum CellType{GENERALCELL, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID};
        enum FaceType{GENERALFACE, TRIANGULAR, QUADRILATERAL};

        Mesh() {};

        //const std::vector<Vars<3>>& getNodeList() const;
        //const std::vector<Cell>& getCellList() const;
        //const std::vector<Face>& getFaceList() const;
        //const std::vector<Boundary>& getBoundaryList() const;
        //const std::vector<int>& getOwnerIndexList() const;
        //const std::vector<int>& getNeighborIndexList() const;

        //int getNodesSize() const;
        //int getFacesSize() const;
        //int getCellsSize() const;

        //void update();
        void loadGmsh2(std::string fileName);

        //void deleteBoundary(std::string boundaryConditionName);

        //size_t cellSize;

    private:
        void sortCells();
        void createFaces();

        //void updateCells();
        //void updateFaces();

        //bool checkFaces() const;

        Field<DIM,1> nodes;
        std::vector<std::vector<int>> cells;
        std::vector<std::vector<int>> faces;
        std::vector<Boundary> boundaries;
        std::vector<int> owners;
        std::vector<int> neighbors;

        std::vector<CellType> cellTypes;
        std::vector<FaceType> faceTypes;
        Field<1,1> CellVolumes;
        Field<DIM,1> CellCenters;
        Field<DIM,1> CellProjectedAreas;

        Field<DIM,1> FaceNormals;
        Field<DIM,1> FaceMidpoints;
        


        //load GMSH
        std::vector<std::string> readFile(std::string fileName);
        std::vector<std::vector<std::string>> parseBlockDataGmsh(const std::vector<std::string>& dataIn, std::string blockName);
        void createNodesGmsh(const std::vector<std::vector<std::string>>& nodesGmsh);
        void createCellsGmsh(const std::vector<std::vector<std::string>>& elementsGmsh);
        void createBoundariesGmsh(const std::vector<std::vector<std::string>>& physicalNamesGmsh, const std::vector<std::vector<std::string>>& elementsGmsh);

        //void stepDownEveryOverIndex(std::vector<int>& vec, int index) const;

        bool isFaceEqual(const std::vector<int>& face, const std::vector<int>& compFace) const;

};



template <int DIM>
std::vector<std::string> Mesh<DIM>::readFile(std::string fileName)
{
    std::ifstream stream;
    std::string line;
    std::vector<std::string> data;

    stream.open(fileName, std::ios_base::in);

    if (stream.is_open())
    {
        while (std::getline(stream, line))
        {
            line.pop_back();
            data.push_back(line);
        }
    }
    stream.close();

    return data;
}


template <int DIM>
std::vector<std::vector<std::string>> Mesh<DIM>::parseBlockDataGmsh(const std::vector<std::string>& dataIn, std::string blockName)
{
    std::vector<std::vector<std::string>> out;

    int lineIndex = 0;

    while (dataIn[lineIndex] != ("$" + blockName))
    {
        ++lineIndex;
    }
    ++lineIndex;

    while(dataIn[lineIndex] != ("$End" + blockName))
    {
        std::vector<std::string> rowData;

        std::string w = "";
        for (auto x : dataIn[lineIndex]) 
        {
            if (x == ' ')
            {
                rowData.push_back(w);
                w = "";
            }
            else
            {
                w = w + x;
            }
        }
        rowData.push_back(w);
        out.push_back(rowData);
        ++lineIndex;
    }

    return out;
}


template <int DIM>
void Mesh<DIM>::createNodesGmsh(const std::vector<std::vector<std::string>>& nodesGmsh)
{
    //int nodesSum = stoi(nodesGmsh[0][0]);

    nodes = Field<DIM,1>(nodesGmsh.size());

    for (int i = 1; i < nodesGmsh.size(); i++)
    {
        if(stoi(nodesGmsh[i][0]) == i)
        {
            for (size_t k = 0; k < DIM; k++)
            {
                nodes[i-1][k] = stof(nodesGmsh[i][k+1]);
            }
        }
        else
        {
            std::cout << "Chybejici Node, index:" << i << std::endl;
        }
    }
}


template <int DIM>
void Mesh<DIM>::createCellsGmsh(const std::vector<std::vector<std::string>>& elementsGmsh)
{
    int cellsSum = stoi(elementsGmsh[0][0]);

    cells.clear();

    for (int i = 1; i < elementsGmsh.size(); i++)
    {
        if(stoi(elementsGmsh[i][0]) == i)
        {
            int numOfTags = stoi(elementsGmsh[i][2]);

            if(DIM == 3)
            {
                switch (stoi(elementsGmsh[i][1]))
                {
                    case 4:
                        cells.push_back(std::vector<int>{stoi(elementsGmsh[i][3+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][4+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][5+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][6+numOfTags]) - 1});
                        cellTypes.push_back(TETRAHEDRON);
                        break;
                    case 5:
                        cells.push_back(std::vector<int>{stoi(elementsGmsh[i][3+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][4+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][5+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][6+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][7+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][8+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][9+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][10+numOfTags]) - 1});
                        cellTypes.push_back(HEXAHEDRON);
                        break;
                    case 6:
                        cells.push_back(std::vector<int>{stoi(elementsGmsh[i][3+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][4+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][5+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][6+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][7+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][8+numOfTags]) - 1});
                        cellTypes.push_back(PRISM);   
                        break;
                    case 7:
                        cells.push_back(std::vector<int>{stoi(elementsGmsh[i][3+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][4+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][5+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][6+numOfTags]) - 1,
                                                         stoi(elementsGmsh[i][7+numOfTags]) - 1});
                        cellTypes.push_back(PYRAMID);
                        break;
                    
                    default:
                        //std::cout << "Neplatný typ 3D bunky na pozici" << i << std::endl;
                        break;
                }
            }
            else if(DIM == 2)
            {
                
            }
            
        }
        else
        {
            std::cout << "Chybejici Cell, index:" << i << std::endl;
            cells.push_back(std::vector<int>({-1}));
        }
    }
}


template<int DIM>
void Mesh<DIM>::sortCells()
{


    Eigen::Array<double,DIM,1> zeroAxisForCellSort;
    //zeroAxisForCellSort(1) = -1.0;

    std::sort(cells.begin(), cells.end(), [zeroAxisForCellSort, this](std::vector<int> cellA, std::vector<int> cellB)
    {
        auto calcApproxCenter = [](std::vector<int> cell, const Field<DIM,1>& nodes)
        {
            Eigen::Array<double,DIM,1> cen;

            for (size_t i = 0; i < cell.size(); i++)
            {
                cen += nodes[cell[i]];
            }
            
            return cen/cell.size();
        };

        return (calcApproxCenter(cellA, nodes) - zeroAxisForCellSort).matrix().norm() < (calcApproxCenter(cellB, nodes) - zeroAxisForCellSort).matrix().norm();

        //return norm2(a.center - zeroAxisForCellSort) < norm2(b.center - zeroAxisForCellSort); 
    });
}


template <int DIM>
void Mesh<DIM>::createFaces()
{
    auto createFaces = [](std::vector<int> cell, CellType cellType, const Field<DIM,1>& nodes)
    {
        switch(cellType)
        {
            case TETRAHEDRON:
                return std::vector<std::vector<int>>{std::vector<int>{cell[0], cell[2], cell[1]},
                                                     std::vector<int>{cell[0], cell[1], cell[3]},
                                                     std::vector<int>{cell[0], cell[3], cell[2]},
                                                     std::vector<int>{cell[1], cell[2], cell[3]}};

            case HEXAHEDRON:
                return std::vector<std::vector<int>>{std::vector<int>{cell[0], cell[1], cell[5], cell[4]},
                                                     std::vector<int>{cell[2], cell[3], cell[7], cell[6]},
                                                     std::vector<int>{cell[0], cell[4], cell[7], cell[3]},
                                                     std::vector<int>{cell[4], cell[5], cell[6], cell[7]},
                                                     std::vector<int>{cell[1], cell[2], cell[6], cell[5]},
                                                     std::vector<int>{cell[0], cell[3], cell[2], cell[1]}};

            case PRISM:
                return std::vector<std::vector<int>>{std::vector<int>{cell[0], cell[2], cell[1]},
                                                     std::vector<int>{cell[3], cell[4], cell[5]},
                                                     std::vector<int>{cell[0], cell[1], cell[4], cell[3]},
                                                     std::vector<int>{cell[0], cell[3], cell[5], cell[2]},
                                                     std::vector<int>{cell[1], cell[2], cell[5], cell[4]}};


            case PYRAMID:
                return std::vector<std::vector<int>>{std::vector<int>{cell[0], cell[3], cell[2], cell[1]},
                                                     std::vector<int>{cell[0], cell[1], cell[4]},
                                                     std::vector<int>{cell[1], cell[2], cell[4]},
                                                     std::vector<int>{cell[2], cell[3], cell[4]},
                                                     std::vector<int>{cell[0], cell[4], cell[3]}};
            
            default:
                std::cout << "ERROR - Nedefinovaný typ buňky" << std::endl;
                return std::vector<std::vector<int>>{std::vector<int>{-1}};
        }
    };

    faces.clear();
    owners.clear();
    neighbors.clear();

    //vytvoreni sten + ownerList, neighborList
    std::list<int> nonNeighborList;
    for (int j = 0; j < cells.size(); j++)
    {
        std::vector<std::vector<int>> ownFaces = createFaces(cells[j], cellTypes[j], nodes);

        for (auto & ownFace : ownFaces)
        {
            bool existInList = false;

            std::list<int>::iterator it;
            for (it = nonNeighborList.begin(); it != nonNeighborList.end(); it++)
            {
                if (isFaceEqual(faces[*it], ownFace))
                {
                    existInList = true;
                    break;
                }                
            }

            if (existInList)
            {
                neighbors[*it] = j;
                nonNeighborList.erase(it);
            }
            else
            {
                nonNeighborList.push_front(faces.size());
                faces.push_back(ownFace);                
                owners.push_back(j);
                neighbors.push_back(-1);
            }
        }
    }
}


template <int DIM>
void Mesh<DIM>::createBoundariesGmsh(const std::vector<std::vector<std::string>>& physicalNamesGmsh, const std::vector<std::vector<std::string>>& elementsGmsh)
{
    boundaries.clear();

    std::vector<std::vector<int>> auxFaceList;
    auxFaceList.clear();
    std::vector<int> auxFacePhysicalGroupList;
    auxFacePhysicalGroupList.clear();

    for (int i = 1; i < elementsGmsh.size(); i++)
    {
        if(stoi(elementsGmsh[i][0]) == i)
        {
            int numOfTags = stoi(elementsGmsh[i][2]);

            if (DIM == 3)
            {
                switch (stoi(elementsGmsh[i][1]))
                {
                    case 2:
                        //TRIANGULAR
                        auxFaceList.push_back(std::vector<int>{stoi(elementsGmsh[i][3+numOfTags]) - 1,
                                                               stoi(elementsGmsh[i][4+numOfTags]) - 1,
                                                               stoi(elementsGmsh[i][5+numOfTags]) - 1});

                        auxFacePhysicalGroupList.push_back(stoi(elementsGmsh[i][3]));
                        break;
                        
                    case 3:
                        //QUADRILATERAL
                        auxFaceList.push_back(std::vector<int>{stoi(elementsGmsh[i][3+numOfTags]) - 1,
                                                               stoi(elementsGmsh[i][4+numOfTags]) - 1,
                                                               stoi(elementsGmsh[i][5+numOfTags]) - 1,
                                                               stoi(elementsGmsh[i][6+numOfTags]) - 1});

                        auxFacePhysicalGroupList.push_back(stoi(elementsGmsh[i][3]));
                        break;
                    
                    default:
                        //std::cout << "Neplatný typ 3D bunky na pozici" << i << std::endl;
                        break;
                }
            }
            else if(DIM == 2)
            {

            }
        }
        else
        {
            std::cout << "Chybejici Cell, index:" << i << std::endl;
        }
    }

    for (int j = 1; j < physicalNamesGmsh.size(); j++)
    {
        if(DIM == 3)
        {
            if(stoi(physicalNamesGmsh[j][0]) != 2) //only 2d phasicalNames
            {
                continue;
            }

            int physicalId = stoi(physicalNamesGmsh[j][1]);
            std::string physicalName = physicalNamesGmsh[j][2];
            physicalName.erase(0,1);
            physicalName.pop_back();

            Boundary auxBoundary = Boundary(physicalName);

            for (int j = 0; j < auxFacePhysicalGroupList.size(); j++)
            {
                if(auxFacePhysicalGroupList[j] == physicalId)
                {
                    for (int i = 0; i < faces.size(); i++)
                    {
                        if(isFaceEqual(faces[i], auxFaceList[j]))
                        {
                            auxBoundary.indexToFaces.push_back(i);
                        }
                    }                
                }
            }
            
            boundaries.push_back(auxBoundary);
        }
        else if(DIM == 2)
        {

        }
    }

    for(auto& boundary : boundaries)
    {
        std::sort(boundary.indexToFaces.begin(), boundary.indexToFaces.end());
    }
}


template <int DIM>
void Mesh<DIM>::loadGmsh2(std::string fileName)
{
    std::vector<std::string> stringData = readFile(fileName);

    createNodesGmsh(parseBlockDataGmsh(stringData, "Nodes"));
    createCellsGmsh(parseBlockDataGmsh(stringData, "Elements"));

    //sortCells();

    createFaces();
    createBoundariesGmsh(parseBlockDataGmsh(stringData, "PhysicalNames"), parseBlockDataGmsh(stringData, "Elements"));

    //update();

    //std::cout << "Mesh initialization complete" <<std::endl;
}



template <int DIM>
bool Mesh<DIM>::isFaceEqual(const std::vector<int>& face, const std::vector<int>& compFace) const
{
    bool isEqual = true;

    for (auto & compNodeIndex : compFace)
    {
        if(std::find(face.begin(), face.end(), compNodeIndex) == face.end())
        {
            isEqual = false;
        }
    }
        
    return isEqual;
};


#endif // MESH_HPP