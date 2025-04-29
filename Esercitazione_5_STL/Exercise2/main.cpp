//PIERANDREI MATTIAS  Esercitazione_5_STL

#include <iostream>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

int main()
{
    // Initialize the mesh
    PolygonalMesh mesh;

    // Try to import the mesh
    if (!ImportMesh(mesh))
    {
        cerr << "Mesh import failed. Make sure the CSV files are present and formatted correctly." << endl;
        return 1;
    }


    // Validate that all edges have non-zero length
    if (!ValidateEdges(mesh))
    {
        cerr << "Edge validation failed!" << endl;
        return 1;
    }

    // Validate that all polygons have non-zero area
	if (!CalculatePolygonArea(mesh))
    {
        cerr << "Area validation failed!" << endl;
        return 1;
    }
		
	// Validate markers
    if (!CheckMarkers(mesh))
    {
        cerr << "Marker validation failed!" << endl;
        return 1;
    }
	
    cout << "Mesh successfully imported and validated!" << endl;
	cout << "- all markers are correctly stored" << endl;
	cout << "- each edge has non-zero length" << endl;
	cout << "- each polygon has a non-zero area" << endl;
	
	/*cout << "Number of points: " << mesh.Cell0DsCoordinates.cols() << endl;
	cout << "Number of segments: " << mesh.Cell1DsExtrema.cols() << endl;*/

    Gedim::UCDUtilities utilities;
    utilities.ExportPoints("./Output_Cell0Ds.inp", mesh.Cell0DsCoordinates);
    utilities.ExportSegments("./Output_Cell1Ds.inp", mesh.Cell0DsCoordinates, mesh.Cell1DsExtrema);
	
	
	//I have checked the imagine of the generated mesh and it is identical to the one shown in the provided images.
	
    return 0;
}
