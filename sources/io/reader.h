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

#ifndef _READER_H
#define	_READER_H

#include <vector>
#include <string>
#include <queue>
#include <cstdlib>
#include <iterator>
#include <fstream>
#include <iostream>
#include <sstream>

#include "basic_types/point.h"
#include "basic_types/box.h"
#include "basic_types/mesh.h"
#include "tetrahedral_trees/node_v.h"
#include "tetrahedral_trees/node_t.h"

using namespace std;
///A class that provides an interface for reading input-file and initializite the library structures
class Reader {
public:
    ///A public method that reads a file containing a tetrahedral mesh
    /*!
     * \param mesh a Mesh& argument, representing the mesh to initialize
     * \param path a string argument, representing the path to the mesh file
     * \return a boolean value, true if the file is correctly readed, false otherwise
     */
    static bool read_mesh(Mesh& mesh, string path);
    ///A public method that reads a file containing a list of points coordinate used into a point location
    /*!
     * \param points a vector<Point>& argument, representing the point list to initialize
     * \param fileName a string argument, representing the path to the points file
     */
    static void read_queries(vector<Point>& points, string fileName);
    ///A public method that reads a file containing a list of boxes used into a box query
    /*!
     * \param boxes a vector<Box>& argument, representing the box list to initialize
     * \param fileName a string argument, representing the path to the boxes file
     */
    static void read_queries(vector<Box>& boxes, string fileName);
    ///A public method that reads a file containing a tree
    /*!
     * \param tree a T& argument, representing the tree to initialize
     * \param n a N& argument, representing the root of the tree
     * \param fileName a string argument, representing the path to the mesh file
     * \return a boolean value, true if the file is correctly readed, false otherwise
     */
    template<class T, class N> static bool read_tree(T& tree, N& n, string fileName);
private:
    ///A constructor method
    Reader() {}
    ///A constructor method
    Reader(const Reader&) {}
    ///A destructor method
    virtual ~Reader() {}
    ///A private method that reads a node into the file and saves the information readed (Generic Version)
    /*!
     * \param n a N* argument, representing the empty node to be set
     * \param input an ifstream& argument, representing the stream to read
     * \param son_number an integer argument, representing the number of son every internal node has
     * \param mesh a Mesh& argument, representing the tetrahedral mesh
     */
    template<class N> static void read_node(N* n, ifstream& input, int son_number, Mesh& mesh);
    ///A protected method that reads a leaf into the file and saves the information readed
    /*!
     * \param n a Node_T* argument, representing the empty node to be set
     * \param input an ifstream& argument, representing the stream to read
     * \param tokens a vector<string>& argument, representing the a row of the file
     */
    static void read_leaf(Node_T* n, ifstream& input, vector<string>& tokens);
    ///A protected method that reads a leaf into the file and saves the information readed
    /*!
     * \param n a Node_V* argument, representing the empty node to be set
     * \param input an ifstream& argument, representing the stream to read
     * \param tokens a vector<string>& argument, representing the a row of the file
     */
    static void read_leaf(Node_V* n, ifstream& input, vector<string>& tokens);
};

template<class T, class N> bool Reader::read_tree(T& tree, N& root, string fileName)
{
    ifstream input(fileName.c_str());

    if (input.is_open() == false)
    {
        cerr << "Error in file " << fileName << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    queue<N*> coda;
    N* current;

    coda.push(&root);

    while (input.good())
    {
        current = coda.front();
        coda.pop();

        Reader::read_node(current,input,tree.get_decomposition().son_number(), tree.get_mesh());

        if (!current->is_leaf())
        {
            for (int i = 0; i < tree.get_decomposition().son_number(); i++)
            {
                N* s = new N();
                current->set_son(s, i);
                coda.push(current->get_son(i));
            }
        }

        if (input.eof() || coda.empty())
            break;
    }
    return true;
}

template<class N> void Reader::read_node(N* n, ifstream &input, int son_number, Mesh&)
{
    string line;
    vector<string> tokens;
    getline(input, line);
    istringstream iss(line);
    copy(istream_iterator<string > (iss), istream_iterator<string > (), back_inserter<vector<string> >(tokens));

    if (tokens.at(0) == "N")
    {
        n->init_sons(son_number);
    }
    else if (tokens.at(0) == "L")
    {
        Reader::read_leaf(n,input,tokens);
    }
}

#endif	/* _READER_H */

