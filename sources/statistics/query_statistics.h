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

#ifndef _QUERYSTATISTICS_H
#define	_QUERYSTATISTICS_H

#include <list>
#include <vector>
#include <set>
#include <bm/bm.h>
#include "utilities/sorting.h"

using namespace std;
///A class representing a container used to store the statistics obtained from a query over a spatial index
class QueryStatistics{
public:
    ///A public variable representing the number of nodes visited during query
    int numNode;
    ///A public variable representing the number of leaves visited during query
    int numLeaf;
    ///A public variable representing the number of geometric tests executed during a query
    int numGeometricTest;
    ///A public variable representing the list of tetrahedra that has been checked during a box query    
    bm::bvector<> checkTetra;
    ///A public array, with an entry for each tetrahedron, containing the number of accesses a tetrahedron had during a query
    vector<int> access_per_tetra;
    ///A public variable representing the tetrahedra (without duplicates) satisfying a query
    vector<int> tetrahedra;

    ///A public variable encoding the number of successfull completely_contains in a box query
    /// i.e., the number of times a box query completely contains a leaf node
    int box_completely_contains_leaf_num;
    int box_completely_contains_bbox_num;
    int box_intersect_bbox_num;
    int box_no_intersect_bbox_num;
    int box_intersect_bbox_geom_tests_num;

    int avoided_tetra_geom_tests_num;

    int tetra_compl_cont_leaf_num;
    int tetra_compl_cont_bbox_num;

    bm::bvector<> avoid_to_check_tetra;

    ///A constructor method
    QueryStatistics()  { numNode=numLeaf=numGeometricTest=0; } // used for point locations
    ///A constructor method
    QueryStatistics(int num_t, int perc_res) // used for box queries
    {
        numNode=numLeaf=numGeometricTest=0;
        access_per_tetra = vector<int>(num_t,0);

        int reserving = num_t / perc_res;
        tetrahedra.reserve(reserving);

        box_completely_contains_leaf_num = box_completely_contains_bbox_num = 0;
        box_intersect_bbox_num = box_no_intersect_bbox_num = box_intersect_bbox_geom_tests_num = 0;

        tetra_compl_cont_leaf_num = tetra_compl_cont_bbox_num = 0;
        avoided_tetra_geom_tests_num = 0;
    }
    ///A destructor method
    virtual ~QueryStatistics()
    {
        numNode=numLeaf=numGeometricTest=0;
        checkTetra.reset();
        tetrahedra.clear();
        access_per_tetra.clear();

        box_completely_contains_leaf_num = box_completely_contains_bbox_num = 0;
        box_intersect_bbox_num = box_no_intersect_bbox_num = box_intersect_bbox_geom_tests_num = 0;

        tetra_compl_cont_leaf_num = tetra_compl_cont_bbox_num = 0;

        avoided_tetra_geom_tests_num = 0;

        avoid_to_check_tetra.reset();
    }

    /**
     * @brief A public procedure that resets the variables of this class (box and line queries wrapper)
     */
    inline void reset(bool )
    {
        numNode=numLeaf=numGeometricTest=0;
        checkTetra.reset();
        tetrahedra.clear();
        fill(access_per_tetra.begin(),access_per_tetra.end(),0);

        box_completely_contains_leaf_num = box_completely_contains_bbox_num = 0;
        box_intersect_bbox_num = box_no_intersect_bbox_num = box_intersect_bbox_geom_tests_num = 0;

        tetra_compl_cont_leaf_num = tetra_compl_cont_bbox_num = 0;

        avoided_tetra_geom_tests_num = 0;

        avoid_to_check_tetra.reset();
    }
    /**
     * @brief A public procedure that resets the variables of this class (point queries wrapper)
     */
    inline void reset()
    {
        numNode=numLeaf=numGeometricTest=0;
        tetrahedra.clear();
    }
};

#endif	/* _QUERYSTATISTICS_H */

