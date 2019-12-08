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

#ifndef _WRITER_H
#define	_WRITER_H

#include <string>
#include <set>

#include <fstream>
#include <queue>
#include <iostream>

#include "tetrahedral_trees/tree.h"
#include "statistics/index_statistics.h"
#include "statistics/full_query_statistics.h"
#include "basic_types/box.h"

#include "tetrahedral_trees/node_v.h"
#include "tetrahedral_trees/node_t.h"

using namespace std;
///A class that provides an interface for writing to file or standard output some data structures or statistics
class Writer {
public:
    ///A public method that writes to file the tree structure
    /*!
     * \param fileName a string argument, representing the output path where the tree structure will be written
     * \param root a N& argument, representing the root of the tree to visit
     * \param division a D& argument, representing the space subdivision type used for the spatial index
     */
    template<class N, class D> static void write_tree(string fileName, N& root, D& division);
    ///A public method that writes to standard output the spatial index statistics
    /*!
     * \param indexStats an IndexStatistics& argument, representing the statistics to save
     */
    static void write_tree_stats(IndexStatistics& indexStats);
    ///A public method that writes to standard output the query statistics
    /*!
     * \param size an integer argument, representing the number of query executed
     * \param fullQueryStats an FullQueryStatistics& argument, representing the statistics to save
     * \param hit_ratio an integer representing the number of successfully answered queries
     */
    static void write_queries_stats(int size, FullQueryStatistics& fullQueryStats, int hit_ratio);
    ///A public method that writes to file a series of points that will be used as point location input
    /*!
     * \param points a set<Point>& argument, representing the points list to save
     * \param fileName a string argument, representing the file name
     */
    static void write_point_queries(set<Point> &points, string fileName);
    ///A public method that writes to file a series of boxes that will be used as box query input
    /*!
     * \param boxes a set<Box>& argument, representing the boxes list to save
     * \param fileName a string argument, representing the file name
     */
    static void write_box_queries(set<Box>& boxes, string fileName);

private:
    ///A constructor method
    Writer() {}
    ///A constructor method
    Writer(const Writer&) {}
    ///A destructor method
    virtual ~Writer() {}

    ///A private method that writes to an output stream a tree node information (Generic Version)
    /*!
     * \param output an ofstream& argument, representing the stream
     * \param n a N* argument, representing the node to save
     */
    static void write_node(ofstream& output, Node_T* n);
    ///A private method that writes to an output stream a tree node information (P_Node Version)
    /*!
     * \param output an ofstream& argument, represents the stream
     * \param n a P_Node* argument, represents the node to save
     */
    static void write_node(ofstream& output, Node_V *n);
};

template<class N, class D> void Writer::write_tree(string fileName, N& root, D& division)
{
    ofstream output(fileName.c_str());
    queue<N*> coda;
    N* visited;

    coda.push(&root);
    bool is_root = true;

    while (!coda.empty()) {
        visited = coda.front();
        coda.pop();

        string begin;

        if (visited->is_leaf())
            begin = "L";
        else
            begin = "N";

        if(!is_root)
            output << endl;
        output << begin << " ";

        is_root = false;

        Writer::write_node(output, visited);

        if (visited->is_leaf() == false)
            for (int i = 0; i < division.son_number(); i++)
                coda.push(visited->get_son(i));
    }
    output.close();
}

#endif	/* _WRITER_H */

