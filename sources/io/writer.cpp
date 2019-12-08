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

#include "tetrahedral_trees/node.h"
#include "basic_types/mesh.h"
#include "writer.h"

void Writer::write_node(ofstream& output, Node_T *n)
{
    if (n->is_leaf())
    {
        output << n->get_real_t_array_size();
        if (n->get_real_t_array_size() > 0)
        {
            output << endl << "  T ";
            for(RunIterator runIt = n->t_array_begin_iterator(), runEnd = n->t_array_end_iterator(); runIt != runEnd; ++runIt)
                output << *runIt << " ";
        }
    }
}


void Writer::write_node(ofstream &output, Node_V *n)
{
    if (n->is_leaf())
    {
        output << n->get_real_v_array_size() << " " << n->get_real_t_array_size();

        if (n->get_real_v_array_size() > 0)
        {
            output << endl << "  V ";
            for(RunIterator runIt = n->v_array_begin_iterator(), runEnd = n->v_array_end_iterator(); runIt != runEnd; ++runIt)
                output << *runIt << " ";
        }

        if (n->get_real_t_array_size() > 0)
        {
            output << endl << "  T ";
            for(RunIterator runIt = n->t_array_begin_iterator(), runEnd = n->t_array_end_iterator(); runIt != runEnd; ++runIt)
                output << *runIt << " ";
        }
    }
}

void Writer::write_tree_stats(IndexStatistics& indexStats)
{
    cout << indexStats.numNode << " ";
    cout << indexStats.numFullLeaf << " ";    
    cout << indexStats.numEmptyLeaf << " ";
    cout << indexStats.minTreeDepth << " ";
    cout << indexStats.avgTreeDepth << " ";
    cout << indexStats.maxTreeDepth << " ";
    cout << indexStats.avgVertexInFullLeaf << " ";
    cout << indexStats.avg_completely_indexed_tetra << " ";
    cout << indexStats.avg_partially_indexed_tetra << " ";
    cout << indexStats.avg_overlapping_tetra << " "; // different from zero only for PM criterion
    cout << indexStats.avg_leaves_for_tetra << " ";
    cout << indexStats.avg_weighted_leaves_for_tetra << " ";
    cout << indexStats.max_leaves_for_tetra << " ";

    cout << (indexStats.numTin1Leaf*100)/static_cast<double>(indexStats.num_leaves_for_tetra.size()) << " ";
    cout << (indexStats.numTin2Leaf*100)/static_cast<double>(indexStats.num_leaves_for_tetra.size()) << " ";
    cout << (indexStats.numTin3Leaf*100)/static_cast<double>(indexStats.num_leaves_for_tetra.size()) << " ";
    cout << (indexStats.numTin4Leaf*100)/static_cast<double>(indexStats.num_leaves_for_tetra.size()) << " ";
    cout << (indexStats.numTinMoreLeaf*100)/static_cast<double>(indexStats.num_leaves_for_tetra.size()) << " ";
    cout << endl;
    // for debug only (not for tables)
    cout << "internal_tetra_per_leaf " << indexStats.min_completely_indexed_tetra << " " << indexStats.avg_completely_indexed_tetra << " " << indexStats.max_completely_indexed_tetra << endl;
    cout << "partial_tetra_per_leaf " << indexStats.min_partially_indexed_tetra << " " << indexStats.avg_partially_indexed_tetra << " " << indexStats.max_partially_indexed_tetra << endl;
    cout << "overlapping_tetra_per_leaf " << indexStats.min_overlapping_tetra << " " << indexStats.avg_overlapping_tetra << " " << indexStats.max_overlapping_tetra << endl;
    cout << "leaf_per_tetra " << indexStats.min_leaves_for_tetra << " " << indexStats.avg_leaves_for_tetra << " " << indexStats.max_leaves_for_tetra << " " << endl;
    cout << "chi_star " << indexStats.avg_weighted_leaves_for_tetra << endl;
    cout << "t_list_length " << indexStats.t_list_length << endl;
    cout << "real_t_list_length " << indexStats.real_t_list_length << endl;
    return;
}

void Writer::write_queries_stats(int size, FullQueryStatistics& fullQueryStats, int hit_ratio)
{    
    cerr << "==query_stats=="<<endl;

    cerr << "nodes_visited: ";
    cerr << fullQueryStats.minNode << " ";
    cerr << fullQueryStats.avgNode / static_cast<double>(size) << " ";
    cerr << fullQueryStats.maxNode << endl;
    cerr << "leaves_visited: ";
    cerr << fullQueryStats.minLeaf << " ";
    cerr << fullQueryStats.avgLeaf / static_cast<double>(size) << " ";
    cerr << fullQueryStats.maxLeaf << endl;

    cerr << "tetra_num: ";
    cerr << fullQueryStats.minTetra << " ";
    cerr << fullQueryStats.avgTetra / static_cast<double>(size) << " ";
    cerr << fullQueryStats.maxTetra << endl;
    cerr << "tetra_compl_cont_leaf_num: ";
    cerr << fullQueryStats.min_tetra_compl_cont_leaf_num << " ";
    cerr << fullQueryStats.avg_tetra_compl_cont_leaf_num / static_cast<double>(size) << " ";
    cerr << fullQueryStats.max_tetra_compl_cont_leaf_num << endl;
    cerr << "tetra_compl_cont_bbox_num: ";
    cerr << fullQueryStats.min_tetra_compl_cont_bbox_num << " ";
    cerr << fullQueryStats.avg_tetra_compl_cont_bbox_num / static_cast<double>(size) << " ";
    cerr << fullQueryStats.max_tetra_compl_cont_bbox_num << endl;

    //should be printed only for box queries
    if(fullQueryStats.maxMultipleTetraAccess > fullQueryStats.minMultipleTetraAccess)
    {
        cerr << "unique_tetra_access: ";
        cerr << fullQueryStats.minUniqueTetraAccess << " ";
        cerr << fullQueryStats.avgUniqueTetraAccess / static_cast<double>(size) << " ";
        cerr << fullQueryStats.maxUniqueTetraAccess << endl;
        cerr << "multiple_tetra_access: ";
        cerr << fullQueryStats.minMultipleTetraAccess << " ";
        cerr << fullQueryStats.avgMultipleTetraAccess / static_cast<double>(size) << " ";
        cerr << fullQueryStats.maxMultipleTetraAccess << endl;
    }        

    cerr << "box_completely_contains_leaf_num: ";
    cerr << fullQueryStats.min_box_completely_contains_leaf_num << " ";
    cerr << fullQueryStats.avg_box_completely_contains_leaf_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.max_box_completely_contains_leaf_num << endl;
    cerr << "box_completely_contains_bbox_num: ";
    cerr << fullQueryStats.min_box_completely_contains_bbox_num << " ";
    cerr << fullQueryStats.avg_box_completely_contains_bbox_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.max_box_completely_contains_bbox_num << endl;
    cerr << "box_intersect_bbox_num: ";
    cerr << fullQueryStats.min_box_intersect_bbox_num << " ";
    cerr << fullQueryStats.avg_box_intersect_bbox_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.max_box_intersect_bbox_num << endl;
    cerr << "box_no_intersect_bbox_num: ";
    cerr << fullQueryStats.min_box_no_intersect_bbox_num << " ";
    cerr << fullQueryStats.avg_box_no_intersect_bbox_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.max_box_no_intersect_bbox_num << endl;

    cerr << "geometric_tests_executed: ";
    cerr << fullQueryStats.minGeometricTest << " ";
    cerr << fullQueryStats.avgGeometricTest / static_cast<double>(size) << " ";
    cerr << fullQueryStats.maxGeometricTest << endl;

    cerr << "box_intersect_bbox_geom_tests_num: ";
    cerr << fullQueryStats.min_box_intersect_bbox_geom_tests_num << " ";
    cerr << fullQueryStats.avg_box_intersect_bbox_geom_tests_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.max_box_intersect_bbox_geom_tests_num << endl;

    cerr << "avoided_geometric_tests: ";
    cerr << fullQueryStats.min_avoided_tetra_geom_tests_num << " ";
    cerr << fullQueryStats.avg_avoided_tetra_geom_tests_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.max_avoided_tetra_geom_tests_num << endl;

    cerr << "hit_ratio: " << hit_ratio << endl;

    cerr << "compact_stats: ";
    cerr << fullQueryStats.avgNode / static_cast<double>(size) << " ";
    cerr << fullQueryStats.avgLeaf / static_cast<double>(size) << " ";
    cerr << fullQueryStats.avgTetra / static_cast<double>(size) << " ";
    cerr << fullQueryStats.avg_tetra_compl_cont_leaf_num / static_cast<double>(size) << " ";
    cerr << fullQueryStats.avg_tetra_compl_cont_bbox_num / static_cast<double>(size) << " ";
//    if(fullQueryStats.maxMultipleTetraAccess > fullQueryStats.minMultipleTetraAccess)
//    {
//        cerr << fullQueryStats.avgUniqueTetraAccess / static_cast<double>(size) << " ";
//        cerr << fullQueryStats.avgMultipleTetraAccess / static_cast<double>(size) << " ";
//    }
    cerr << fullQueryStats.avg_box_completely_contains_leaf_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.avg_box_completely_contains_bbox_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.avg_box_intersect_bbox_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.avg_box_no_intersect_bbox_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.avgGeometricTest / static_cast<double>(size) << " ";
    cerr << fullQueryStats.avg_box_intersect_bbox_geom_tests_num / static_cast<double>(size)<< " ";
    cerr << fullQueryStats.avg_avoided_tetra_geom_tests_num / static_cast<double>(size)<< " ";
    cerr << endl;

}

void Writer::write_point_queries(set<Point>& points, string fileName)
{
    ofstream output(fileName.c_str());
    output << points.size() << endl;
    for(set<Point>::iterator it=points.begin(); it!=points.end(); ++it)
    {
        const Point &p = *it;
        output << p.get_x() << " " << p.get_y() << " " << p.get_z() << endl;
    }
    output.close();
}

void Writer::write_box_queries(set<Box> &boxes, string fileName)
{
    ofstream output(fileName.c_str());
    output << boxes.size() << endl;
    for(set<Box>::iterator it=boxes.begin(); it!=boxes.end(); ++it)
    {
        output << *it << endl;
    }
    output.close();
}

