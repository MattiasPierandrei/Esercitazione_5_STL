//PIERANDREI MATTIAS

#pragma once

#include <vector>
#include <map>
#include <array>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolygonalLibrary
{
struct PolygonalMesh
{
    unsigned int NumCell0Ds;
    unsigned int NumCell1Ds;
    unsigned int NumCell2Ds;

    // Cell0Ds (points)
    vector<unsigned int> Cell0DsId;     // ID points
    Eigen::MatrixXd Cell0DsCoordinates;      // Coordinates (X, Y) points
    map<unsigned int, list<unsigned int> > MarkerCell0Ds; // Marker points

    // Cell1Ds (edges)
    vector<unsigned int> Cell1DsId;      // ID edges
    Eigen::MatrixXi Cell1DsExtrema;           // Extremes (origin, end) edges
    map<unsigned int, list<unsigned int> > MarkerCell1Ds;     // Marker edges

    // Cell2Ds (polygons)
    vector<unsigned int> Cell2DsId;       // ID polygons
    vector<vector<unsigned int>> Cell2DsVertices;   // Vertices 
    vector<vector<unsigned int>> Cell2DsEdges;      // Edges
    map<unsigned int, list<unsigned int> > MarkerCell2Ds;   // Marker polygons
};
}