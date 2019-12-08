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

#ifndef _FULLQUERYSTATISTICS_H
#define	_FULLQUERYSTATISTICS_H

#include <limits>
using namespace std;

///A class storing the statistics obtained from a series of query
class FullQueryStatistics{
public:
    ///A constructor method
    FullQueryStatistics()
    {
        minTetra = std::numeric_limits<int>::max(); maxTetra=std::numeric_limits<int>::min();
        minNode = std::numeric_limits<int>::max(); maxNode=std::numeric_limits<int>::min();
        minLeaf = std::numeric_limits<int>::max(); maxLeaf=std::numeric_limits<int>::min();

        minGeometricTest = std::numeric_limits<int>::max(); maxGeometricTest=std::numeric_limits<int>::min();
        minUniqueTetraAccess = std::numeric_limits<int>::max(); maxUniqueTetraAccess=std::numeric_limits<int>::min();
        minMultipleTetraAccess = std::numeric_limits<int>::max(); maxMultipleTetraAccess=std::numeric_limits<int>::min();
        avgTetra = avgNode = avgLeaf = avgGeometricTest = avgMultipleTetraAccess = avgUniqueTetraAccess = 0.0;

        min_box_completely_contains_leaf_num = std::numeric_limits<int>::max();
        max_box_completely_contains_leaf_num = std::numeric_limits<int>::min();
        min_box_completely_contains_bbox_num = std::numeric_limits<int>::max();
        max_box_completely_contains_bbox_num = std::numeric_limits<int>::min();
        min_box_intersect_bbox_num = std::numeric_limits<int>::max();
        max_box_intersect_bbox_num = std::numeric_limits<int>::min();
        min_box_no_intersect_bbox_num = std::numeric_limits<int>::max();
        max_box_no_intersect_bbox_num = std::numeric_limits<int>::min();
        min_box_intersect_bbox_geom_tests_num = std::numeric_limits<int>::max();
        max_box_intersect_bbox_geom_tests_num = std::numeric_limits<int>::min();
        avg_box_completely_contains_leaf_num = avg_box_completely_contains_bbox_num = 0.0;
        avg_box_intersect_bbox_num = avg_box_no_intersect_bbox_num = avg_box_intersect_bbox_geom_tests_num = 0.0;

        min_tetra_compl_cont_leaf_num = std::numeric_limits<int>::max();
        max_tetra_compl_cont_leaf_num = std::numeric_limits<int>::min();
        min_tetra_compl_cont_bbox_num = std::numeric_limits<int>::max();
        max_tetra_compl_cont_bbox_num = std::numeric_limits<int>::min();
        avg_tetra_compl_cont_leaf_num = avg_tetra_compl_cont_bbox_num = 0.0;

        min_avoided_tetra_geom_tests_num = std::numeric_limits<int>::max();
        max_avoided_tetra_geom_tests_num = std::numeric_limits<int>::min();
        avg_avoided_tetra_geom_tests_num = 0.0;
    }

    ///A public variable representing the minimum number of tetrahedra found during query
    int minTetra;
    ///A public variable representing the maximum number of tetrahedra found during query
    int maxTetra;
    ///A public variable representing the minimum number of node visited during query
    int minNode;
    ///A public variable representing the maximum number of node visited during query
    int maxNode;
    ///A public variable representing the minimum number of leaf node visited during query
    int minLeaf;
    ///A public variable representing the minimum number of leaf node visited during query
    int maxLeaf;
    ///A public variable representing the minimum number of geometric test executed during query
    int minGeometricTest;
    ///A public variable representing the maximum number of geometric test executed during query
    int maxGeometricTest;
    ///A public variable representing the average number of tetrahedra found during query
    double avgTetra;
    ///A public variable representing the average number of node visited during query
    double avgNode;
    ///A public variable representing the average number of leaf node visited during query
    double avgLeaf;
    ///A public variable representing the average number of geometric test executed during query
    double avgGeometricTest;
    ///A public variable representing the minimum number of tetrahedra that have been uniquely accessed during a query
    int minUniqueTetraAccess;
    ///A public variable representing the maximum number of tetrahedra that have been uniquely accessed during a query
    int maxUniqueTetraAccess;
    ///A public variable representing the average number of tetrahedra that have been uniquely accessed during a query
    double avgUniqueTetraAccess;
    ///A public variable representing the minimum number of tetrahedra that have been accessed several times during a query
    int minMultipleTetraAccess;
    ///A public variable representing the maximum number of tetrahedra that have been accessed several times during a query
    int maxMultipleTetraAccess;
    ///A public variable representing the average number of tetrahedra that have been accessed several times during a query
    double avgMultipleTetraAccess;

    int min_box_completely_contains_leaf_num;
    double avg_box_completely_contains_leaf_num;
    int max_box_completely_contains_leaf_num;

    int min_box_completely_contains_bbox_num;
    double avg_box_completely_contains_bbox_num;
    int max_box_completely_contains_bbox_num;

    int min_box_intersect_bbox_num;
    double avg_box_intersect_bbox_num;
    int max_box_intersect_bbox_num;

    int min_box_no_intersect_bbox_num;
    double avg_box_no_intersect_bbox_num;
    int max_box_no_intersect_bbox_num;

    int min_box_intersect_bbox_geom_tests_num;
    double avg_box_intersect_bbox_geom_tests_num;
    int max_box_intersect_bbox_geom_tests_num;

    int min_tetra_compl_cont_leaf_num;
    double avg_tetra_compl_cont_leaf_num;
    int max_tetra_compl_cont_leaf_num;

    int min_tetra_compl_cont_bbox_num;
    double avg_tetra_compl_cont_bbox_num;
    int max_tetra_compl_cont_bbox_num;

    int min_avoided_tetra_geom_tests_num;
    double avg_avoided_tetra_geom_tests_num;
    int max_avoided_tetra_geom_tests_num;

};

#endif	/* _FULLQUERYSTATISTICS_H */

