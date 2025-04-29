//PIERANDREI MATTIAS	

#pragma once

#include <iostream>
#include "Eigen/Eigen"
#include "PolygonalMesh.hpp"

using namespace std;

namespace PolygonalLibrary
{
    /// import mesh
    /// return: true if the import was successful
    bool ImportMesh(PolygonalMesh& mesh);

    /// import Cell0Ds
    /// return: true if the import was successful
    bool ImportCell0Ds(PolygonalMesh& mesh);

    /// import Cell1Ds 
    /// return: true if the import was successful
    bool ImportCell1Ds(PolygonalMesh& mesh);

    /// import Cell2Ds 
    /// return: true if the import was successful
    bool ImportCell2Ds(PolygonalMesh& mesh);
	
	/// checking edges length
    /// return: true if are all not null
    bool ValidateEdges(const PolygonalMesh& mesh);
	
	
	/// checking areas
    /// return: true if are all non null
	bool CalculatePolygonArea(const PolygonalMesh& mesh) ;


    /// checking if for each marker there is an ID for each cell
    /// return: true if all marker are correctly stored
    bool CheckMarkers(const PolygonalMesh& mesh) ;


}
