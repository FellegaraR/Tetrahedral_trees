/*
    This file is part of the Tetrahedral Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Tetrahedral Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Tetrahedral Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Tetrahedral Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "reader.h"
#include <sstream>
#include <algorithm>
#include "geometry/geometry_wrapper.h"

bool Reader::read_mesh(Mesh& mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string line;
    getline(input, line);
    int tpos = line.find_first_of(' ');

    int num_vertices = atoi(line.substr(0, tpos).c_str());
    int num_tetrahedra = atoi(line.substr(tpos).c_str());

    if (num_vertices == 0 || num_tetrahedra == 0)
    {
        cerr << "This is not a valid .ts file: " << path << endl;
        return false;
    }

    mesh.reserve(num_vertices, num_tetrahedra);

    double x, y, z, field;
    //legge i vertici aggiustando il dominio..
    for (int i = 0; i < num_vertices; i++)
    {
        input >> x;
        input >> y;
        input >> z;
        input >> field;
        if (input.eof())
            break;

        Vertex v = Vertex(x, y, z, field);
        mesh.add_vertex(v);
        if (i == 0) {
            Box b = Box(v, v);
            mesh.set_domain(b);
        } else {
            mesh.get_domain().resize(v);
        }
    }

    int v[4];
    for (int i = 0; i < num_tetrahedra; i++)
    {
        for (int j = 0; j < 4; j++) {
            int index;
            input >> index;
            v[j] = index+1;
        }
        Tetrahedron t = Tetrahedron(v[0], v[1], v[2], v[3]);
        mesh.add_tetrahedron(t);
    }

    return true;
}

void Reader::read_queries(vector<Point>& points, string fileName)
{
    ifstream input(fileName.c_str());
    int size = 0;
    input >> size;
    points.reserve(size);
    double x = 0, y = 0, z = 0;
    while (input)
    {
        input >> x;
        input >> y;
        input >> z;
        if (input.eof())
            break;
        points.push_back(Point(x, y, z));
    }
}

void Reader::read_queries(vector<Box> &boxes, string fileName)
{
    ifstream input(fileName.c_str());
    int size = 0;
    input >> size;
    boxes.reserve(size);
    double x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0;
    while (input) {
        input >> x1;
        input >> y1;
        input >> z1;
        input >> x2;
        input >> y2;
        input >> z2;
        if (input.eof())
            break;
        Point min = Point(x1, y1, z1);
        Point max = Point(x2, y2, z2);
        Box b = Box(min, max);
        boxes.push_back(b);
    }
}

void Reader::read_leaf(Node_T* n, ifstream& input, vector<string>& tokens)
{
    string line;
    int numTetra = atoi(tokens.at(1).c_str());
    if(numTetra > 0)
    {
        getline(input, line);
        vector<string> tokens2;
        istringstream iss2(line);
        copy(istream_iterator<string > (iss2),
             istream_iterator<string > (),
             back_inserter<vector<string> >(tokens2));
        for (unsigned int i = 1; i < tokens2.size(); i++)
            n->add_tetrahedron(atoi(tokens2.at(i).c_str()));
    }
}

void Reader::read_leaf(Node_V *n, ifstream &input, vector<string>& tokens)
{
    string line;
    int numVertex = atoi(tokens.at(1).c_str());
    int numTetra = atoi(tokens.at(2).c_str());

    if (numVertex > 0)
    {
        getline(input, line);
        vector<string> tokens2;
        istringstream iss2(line);
        copy(istream_iterator<string > (iss2),
             istream_iterator<string > (),
             back_inserter<vector<string> >(tokens2));
        for (unsigned int i = 1; i < tokens2.size(); i++)
            n->add_vertex(atoi(tokens2.at(i).c_str()));
    }

    if(numTetra > 0)
    {
        getline(input, line);
        vector<string> tokens3;
        istringstream iss3(line);
        copy(istream_iterator<string > (iss3),
             istream_iterator<string > (),
             back_inserter<vector<string> >(tokens3));
        for (unsigned int i = 1; i < tokens3.size(); i++)
            n->add_tetrahedron(atoi(tokens3.at(i).c_str()));
    }
}
