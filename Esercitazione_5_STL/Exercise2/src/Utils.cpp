//PIERANDREI MATTIAS

#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

namespace PolygonalLibrary
{
	
	
/////////////////////////////////////////////////////////////////////////////////////////////
	
	
//import and check mesh
bool ImportMesh(PolygonalMesh& mesh)
{

    if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;

    return true;

}
	
	
/////////////////////////////////////////////////////////////////////////////////////////////
	
	
bool ImportCell0Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell0Ds.csv");
	
    if(file.fail()) {
        return false;
    }

    list<string> listLines;
	
	//reading file
    string line;
    while (getline(file, line))
        listLines.push_back(line);
    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell0Ds = listLines.size();
    if (mesh.NumCell0Ds == 0) {
        cerr << "No cell 0D" << endl;
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (const string& line : listLines)
    {
        string cleanedLine = line;
        
        // changing "," with "."
        std::replace(cleanedLine.begin(), cleanedLine.end(), ',', '.');

        // separating values by ";"
        stringstream ss(cleanedLine);
        string token;
        
        unsigned int id, marker;
        double x, y;

        // ID
        getline(ss, token, ';');
        id = std::stoi(token);  // in int
        // Marker
        getline(ss, token, ';');
        marker = std::stoi(token);  // in int
        // X
        getline(ss, token, ';');
        x = std::stod(token);  // in double
        // Y
        getline(ss, token, ';');
        y = std::stod(token);  // in double

        // check parsing error
        if (ss.fail()) {
            cerr << "Parsing error: " << line << endl;
            return false;
        }

        mesh.Cell0DsId.push_back(id);
        mesh.Cell0DsCoordinates(0, id) = x;
        mesh.Cell0DsCoordinates(1, id) = y;

        // Marker
        if(marker != 0)
		{
			auto it = mesh.MarkerCell0Ds.find(marker);
			if(it != mesh.MarkerCell0Ds.end() )
			{
				mesh.MarkerCell0Ds[marker].push_back(id);
			}
			else{
				mesh.MarkerCell0Ds.insert({marker, {id}});
			}
		}

    }

    return true;
}
	
	
/////////////////////////////////////////////////////////////////////////////////////////////
	
	
	
bool ImportCell1Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell1Ds.csv");

    if(file.fail()) {
        return false;
    }

    list<string> listLines;
    string line;
    
    // reading file
    while (getline(file, line)) {
        listLines.push_back(line);
    }
    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0) {
        cerr << "No cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

    for (const string& line : listLines) {
        string cleanedLine = line;
        
        // changing "," with "."
        std::replace(cleanedLine.begin(), cleanedLine.end(), ',', '.');
        
        // separating values by ";"
        stringstream ss(cleanedLine);
        string token;
        
        unsigned int id, marker;
        int vertex1, vertex2;

        // ID
        getline(ss, token, ';');
        id = std::stoi(token);  //in int

        // Marker
        getline(ss, token, ';');
        marker = std::stoi(token);  //in int

        // Vertices
        getline(ss, token, ';');
        vertex1 = std::stoi(token);  //in int
        getline(ss, token, ';');
        vertex2 = std::stoi(token);  //in int

        // check parsing error
        if (ss.fail()) {
            cerr << "Parsing error: " << line << endl;
            return false;
        }

        mesh.Cell1DsId.push_back(id);
        mesh.Cell1DsExtrema(0, id) = vertex1;
        mesh.Cell1DsExtrema(1, id) = vertex2;

        // Marker
        if (marker != 0) {
            auto it = mesh.MarkerCell1Ds.find(marker);
            if (it != mesh.MarkerCell1Ds.end()) {
                it->second.push_back(id);
            } else {
                mesh.MarkerCell1Ds.insert({marker, {id}});
            }
        }
    }

    return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////



bool ImportCell2Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell2Ds.csv");
    if (file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);
    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell2Ds = listLines.size();
    if (mesh.NumCell2Ds == 0)
        return false;

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (const string& line : listLines)
    {
        stringstream ss(line);
        string token;

        // ID
        getline(ss, token, ';');
        unsigned int id = stoi(token);
        // Marker
        getline(ss, token, ';');
        unsigned int marker = stoi(token);
		
        // vertices
        getline(ss, token, ';');
        unsigned int numVertices = stoi(token);
        vector<unsigned int> vertices(numVertices);
        for (unsigned int i = 0; i < numVertices; ++i) {
            getline(ss, token, ';');
            vertices[i] = stoi(token);
        }

        // edges
        getline(ss, token, ';');
        unsigned int numEdges = stoi(token);
        vector<unsigned int> edges(numEdges);
        for (unsigned int i = 0; i < numEdges; ++i) {
            getline(ss, token, ';');
            edges[i] = stoi(token);
        }


        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsEdges.push_back(edges);

        // marker
        if (marker != 0)
        {
            auto it = mesh.MarkerCell2Ds.find(marker);
            if (it != mesh.MarkerCell2Ds.end()) {
                it->second.push_back(id);
            } else {
                mesh.MarkerCell2Ds.insert({marker, {id}});
            }
        }
    }

    return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////



bool ValidateEdges(const PolygonalMesh& mesh) {
	
    for (unsigned int i = 0; i < mesh.NumCell1Ds; ++i) {
        // getting vertices
        int vertex1 = mesh.Cell1DsExtrema(0, i);
        int vertex2 = mesh.Cell1DsExtrema(1, i);

        // getting coordinates (from Cell0DsCoordinates)
        double x1 = mesh.Cell0DsCoordinates(0, vertex1);
        double y1 = mesh.Cell0DsCoordinates(1, vertex1);
        double x2 = mesh.Cell0DsCoordinates(0, vertex2);
        double y2 = mesh.Cell0DsCoordinates(1, vertex2);

        // length edge
        double length = std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));

        // check the length
        if (length == 0.0) {  
            std::cerr << "Error: edge ID " << mesh.Cell1DsId[i] << " has zero length." << endl;       
            return false;  
        }
    }

    return true;  
}




/////////////////////////////////////////////////////////////////////////////////////////////



bool CalculatePolygonArea(const PolygonalMesh& mesh) 
{
    // for every polygon
    for (unsigned int i = 0; i < mesh.NumCell2Ds; ++i) 
    {
        // getting vertices of currently polygon
        const vector<unsigned int>& vertices = mesh.Cell2DsVertices[i];

        double area = 0.0;
        size_t numVertices = vertices.size();

        // Gauss
        for (size_t j = 0; j < numVertices; ++j) 
        {
            unsigned int a = vertices[j];
            unsigned int b = vertices[(j + 1) % numVertices]; 

            // getting coordinates of a and b
            double x1 = mesh.Cell0DsCoordinates(0, a);
            double y1 = mesh.Cell0DsCoordinates(1, a);
            double x2 = mesh.Cell0DsCoordinates(0, b);
            double y2 = mesh.Cell0DsCoordinates(1, b);

            area += (x1 * y2) - (x2 * y1);
        }

        area = 0.5 * fabs(area);

        if (area == 0.0) 
        {
            cerr << "Polygon with zero area found at ID: " << mesh.Cell2DsId[i] << endl;
            return false;  
        }
    }

    return true;  //no zero area polygon
}




/////////////////////////////////////////////////////////////////////////////////////////////


bool CheckMarkers(const PolygonalMesh& mesh) {
    // Cells 0D
    for (const auto& [marker, ids] : mesh.MarkerCell0Ds) {
        if (marker == 0) continue;  

        for (unsigned int id : ids) {
            if (std::find(mesh.Cell0DsId.begin(), mesh.Cell0DsId.end(), id) == mesh.Cell0DsId.end()) {
                std::cerr << "Error: point with ID " << id << " not found" << std::endl;
                return false;
            }
        }
    }

    // Cells 1D
    for (const auto& [marker, ids] : mesh.MarkerCell1Ds) {
        if (marker == 0) continue;

        for (unsigned int id : ids) {
            if (std::find(mesh.Cell1DsId.begin(), mesh.Cell1DsId.end(), id) == mesh.Cell1DsId.end()) {
                std::cerr << "Error: edge with ID " << id << " not found" << std::endl;
                return false;
            }
        }
    }

    // Celle 2D
    for (const auto& [marker, ids] : mesh.MarkerCell2Ds) {
        if (marker == 0) continue;

        for (unsigned int id : ids) {
            if (std::find(mesh.Cell2DsId.begin(), mesh.Cell2DsId.end(), id) == mesh.Cell2DsId.end()) {
                std::cerr << "Error: polygon with ID " << id << " not found" << std::endl;
                return false;
            }
        }
    }

    return true;  // All non zero ID are okay
}



/////////////////////////////////////////////////////////////////////////////////////////////


}