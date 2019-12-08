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

#ifndef _INDEXSTATISTICS_H
#define	_INDEXSTATISTICS_H

#include <vector>

using namespace std;

///A class representing a container used to store the statistics obtained from a spatial index
class IndexStatistics{
public:
    ///A constructor method
    IndexStatistics()
    {
        numNode=0;
        numFullLeaf=numEmptyLeaf=0;

        minTreeDepth=-1;
        avgTreeDepth=maxTreeDepth=0;
        minVertexInFullLeaf=-1;
        avgVertexInFullLeaf=maxVertexInFullLeaf=0;

        min_partially_indexed_tetra=-1;
        avg_partially_indexed_tetra=max_partially_indexed_tetra=0;

        min_overlapping_tetra=-1;
        avg_overlapping_tetra=max_overlapping_tetra=0;

        min_completely_indexed_tetra=-1;
        avg_completely_indexed_tetra=max_completely_indexed_tetra=0;

        numTin1Leaf = numTin2Leaf = numTin3Leaf = numTin4Leaf = numTinMoreLeaf = 0;
        min_leaves_for_tetra=-1;
        avg_leaves_for_tetra=max_leaves_for_tetra=0;
        avg_weighted_leaves_for_tetra = 0;

        t_list_length = 0;
        real_t_list_length = 0;
    }

    ///A public variable representing the number of tree nodes
    int numNode;
    ///A public variable representing the number of leaf containing significant information
    int numFullLeaf;
    ///A public variable representing the number of empty leaf
    int numEmptyLeaf;
    ///A public variable representing the minimum tree depth
    int minTreeDepth;
    ///A public variable representing the maximum tree depth
    int maxTreeDepth;
    ///A public variable representing the average tree depth
    double avgTreeDepth;
    ///A public variable representing the minimum number of vertex per full leaf node
    int minVertexInFullLeaf;
    ///A public variable representing the m aximum number of vertex per full leaf node
    int maxVertexInFullLeaf;
    ///A public variable representing the average number of vertex per full leaf node
    double avgVertexInFullLeaf;
    ///A public variable representing the minimum number of tetrahedra partially indexed by a leaf node
    int min_partially_indexed_tetra;
    ///A public variable representing the maximum number of tetrahedra partially indexed by a leaf node
    int max_partially_indexed_tetra;
    ///A public variable representing the average number of tetrahedra partially indexed by a leaf node
    double avg_partially_indexed_tetra;
    ///A public variable representing the minimum number of tetrahedra overlapping-only a leaf node
    int min_overlapping_tetra;
    ///A public variable representing the maximum number of tetrahedra overlapping-only a leaf node
    int max_overlapping_tetra;
    ///A public variable representing the average number of tetrahedra overlapping-only a leaf node
    double avg_overlapping_tetra;
    ///A public variable representing the minimum number of tetrahedra completely indexed by a leaf node
    int min_completely_indexed_tetra;
    ///A public variable representing the maximum number of tetrahedra completely indexed by a leaf node
    int max_completely_indexed_tetra;
    ///A public variable representing the average number of tetrahedra completely indexed by a leaf node
    double avg_completely_indexed_tetra;
    ///A public array, with an entry for each tetrahedron, containing the number of leaf indexing each tetrahedron
    vector<int> num_leaves_for_tetra;
    ///A public variable representing the number of tetrahedra indexed in exactly one leaf
    int numTin1Leaf;
    ///A public variable representing the number of tetrahedra indexed in exactly two leaves
    int numTin2Leaf;
    ///A public variable representing the number of tetrahedra indexed in exactly three leaves
    int numTin3Leaf;
    ///A public variable representing the number of tetrahedra indexed in exactly four leaves
    int numTin4Leaf;
    ///A public variable representing the number of tetrahedra indexed in more than four leaves
    int numTinMoreLeaf;
    ///A public variable representing the minimum number of leaves indexing a single tetrahedron
    int min_leaves_for_tetra;
    ///A public variable representing the maximum number of leaves indexing a single tetrahedron
    int max_leaves_for_tetra;
    ///A public variable representing the average number of leaves indexing a single tetrahedron (chi)
    double avg_leaves_for_tetra;
    ///A public variable representing the average weighted number of leaves indexing a single tetrahedron (chi_w)
    double avg_weighted_leaves_for_tetra;
    ///A public variable representing the summation of the compressed tetrahedra arrays
    int t_list_length;
    ///A public variable representing the summation of the un-compressed tetrahedra arrays
    int real_t_list_length;
};

#endif	/* _INDEXSTATISTICS_H */
