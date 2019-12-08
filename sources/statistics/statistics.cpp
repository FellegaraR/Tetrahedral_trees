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

#include "statistics.h"

void Statistics::compute_leaf_statistics(Node_T &n, Box& dom, Mesh& mesh, bool)
{
    int num_t_completely = 0;
    int num_t_partially = 0;
    int num_t_overlapping = 0;

    this->indexStats.t_list_length += n.get_t_array_size();
    this->indexStats.real_t_list_length += n.get_real_t_array_size();

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& runIt = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*runIt);
        if(n.completely_indexes_tetrahedron_vertices_dom(tet,dom,mesh))
            num_t_completely++;
        else if(n.indexes_tetrahedron_vertices_dom(tet,dom,mesh))
            num_t_partially++;
        else
            num_t_overlapping++;
        this->indexStats.num_leaves_for_tetra[*runIt-1]++;
    }

    if((num_t_completely + num_t_overlapping + num_t_partially) > 0)
    {
        this->indexStats.numFullLeaf++;
        set_leaf_tetrahedra_stats(num_t_completely,num_t_partially,num_t_overlapping);
    }
    else
    {
        this->indexStats.numEmptyLeaf++;
    }
}

void Statistics::compute_leaf_statistics(Node_V &n, Box& dom, Mesh& mesh, bool reindex)
{
    int num_t_completely = 0;
    int num_t_partially = 0;
    int num_t_overlapping = 0;

    this->indexStats.t_list_length += n.get_t_array_size();
    this->indexStats.real_t_list_length += n.get_real_t_array_size();

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& runIt = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*runIt);
        if((reindex && n.completely_indexes_tetrahedron_vertices(tet)) || n.completely_indexes_tetrahedron_vertices_dom(tet,dom,mesh))
            num_t_completely++;
        else if((reindex && n.indexes_tetrahedron_vertices(tet)) || n.indexes_tetrahedron_vertices_dom(tet,dom,mesh))
            num_t_partially++;
        else
            num_t_overlapping++;
        this->indexStats.num_leaves_for_tetra[*runIt-1]++;
    }


    if((num_t_completely + num_t_overlapping + num_t_partially) > 0)
    {
        this->indexStats.numFullLeaf++;
        set_leaf_vertices_stats(n.get_real_v_array_size());
        set_leaf_tetrahedra_stats(num_t_completely,num_t_partially,num_t_overlapping);
    }
    else
    {        
        this->indexStats.numEmptyLeaf++;
    }
}

void Statistics::set_leaf_vertices_stats(int num_vertex)
{
    if(this->indexStats.minVertexInFullLeaf==-1 || this->indexStats.minVertexInFullLeaf > num_vertex)
        this->indexStats.minVertexInFullLeaf = num_vertex;
    if(this->indexStats.maxVertexInFullLeaf < num_vertex)
        this->indexStats.maxVertexInFullLeaf = num_vertex;
    this->indexStats.avgVertexInFullLeaf += num_vertex;
}

void Statistics::set_leaf_tetrahedra_stats(int num_t_completely, int num_t_partially, int num_t_overlapping)
{
    if(this->indexStats.min_completely_indexed_tetra==-1 || this->indexStats.min_completely_indexed_tetra > num_t_completely)
        this->indexStats.min_completely_indexed_tetra =  num_t_completely;
    if(this->indexStats.max_completely_indexed_tetra < num_t_completely)
        this->indexStats.max_completely_indexed_tetra = num_t_completely;
    this->indexStats.avg_completely_indexed_tetra += num_t_completely;

    if(this->indexStats.min_partially_indexed_tetra==-1 || this->indexStats.min_partially_indexed_tetra > num_t_partially)
        this->indexStats.min_partially_indexed_tetra = num_t_partially;
    if(this->indexStats.max_partially_indexed_tetra < num_t_partially)
        this->indexStats.max_partially_indexed_tetra = num_t_partially;
    this->indexStats.avg_partially_indexed_tetra += num_t_partially;

    if(this->indexStats.min_overlapping_tetra==-1 || this->indexStats.min_overlapping_tetra > num_t_overlapping)
        this->indexStats.min_overlapping_tetra = num_t_overlapping;
    if(this->indexStats.max_overlapping_tetra < num_t_overlapping)
        this->indexStats.max_overlapping_tetra = num_t_overlapping;
    this->indexStats.avg_overlapping_tetra += num_t_overlapping;
}

void Statistics::calc_remaining_index_statistics()
{
    for(int_vect_iter iter=this->indexStats.num_leaves_for_tetra.begin(); iter!=this->indexStats.num_leaves_for_tetra.end(); ++iter)
    {
        if(*iter == 1)
            this->indexStats.numTin1Leaf++;
        else if(*iter == 2)
            this->indexStats.numTin2Leaf++;
        else if(*iter == 3)
            this->indexStats.numTin3Leaf++;
        else if(*iter == 4)
            this->indexStats.numTin4Leaf++;
        else
            this->indexStats.numTinMoreLeaf++;

        if(this->indexStats.min_leaves_for_tetra == -1 || this->indexStats.min_leaves_for_tetra > *iter)
            this->indexStats.min_leaves_for_tetra = *iter;
        if(this->indexStats.max_leaves_for_tetra < *iter)
            this->indexStats.max_leaves_for_tetra = *iter;

        this->indexStats.avg_leaves_for_tetra += *iter;

        if(*iter!=1)
            this->indexStats.avg_weighted_leaves_for_tetra += *iter;
    }

    if(this->indexStats.numNode > 0)
        this->indexStats.avgTreeDepth /= (this->indexStats.numEmptyLeaf + this->indexStats.numFullLeaf);
    if(this->indexStats.numFullLeaf > 0)
    {
        this->indexStats.avgVertexInFullLeaf /= this->indexStats.numFullLeaf;
        this->indexStats.avg_partially_indexed_tetra /= this->indexStats.numFullLeaf;

        if(this->indexStats.avg_overlapping_tetra != 0)
            this->indexStats.avg_overlapping_tetra /= this->indexStats.numFullLeaf;
        if(this->indexStats.avg_completely_indexed_tetra != 0)
            this->indexStats.avg_completely_indexed_tetra /= this->indexStats.numFullLeaf;
    }
    if(this->indexStats.num_leaves_for_tetra.size() > 0)
    {
        this->indexStats.avg_leaves_for_tetra /= this->indexStats.num_leaves_for_tetra.size();
        this->indexStats.avg_weighted_leaves_for_tetra /= (this->indexStats.numTin2Leaf+this->indexStats.numTin3Leaf+this->indexStats.numTin4Leaf+this->indexStats.numTinMoreLeaf);
    }
}

void Statistics::check_inconsistencies()
{
    if(this->indexStats.numEmptyLeaf == 0)
        this->indexStats.min_overlapping_tetra = 0;
    if(this->indexStats.numFullLeaf == 0)
    {
        this->indexStats.minVertexInFullLeaf = 0;
        this->indexStats.min_partially_indexed_tetra = 0;
    }
    if(this->indexStats.numNode == 0)
        this->indexStats.minTreeDepth = 0;
    if(this->indexStats.num_leaves_for_tetra.size() == 0)
        this->indexStats.avg_leaves_for_tetra = 0;

    if(this->indexStats.minVertexInFullLeaf==-1)
        this->indexStats.minVertexInFullLeaf=0;
    if(this->indexStats.min_overlapping_tetra==-1)
        this->indexStats.min_overlapping_tetra=0;

    return;
}

int Statistics::compute_queries_statistics(QueryStatistics &qS)
{
    int hit_ratio=0;

    if(qS.tetrahedra.size() > 0)
        hit_ratio++;

    //tetrahedra statistics
    if(this->fullQueryStats.minTetra > static_cast<int>(qS.tetrahedra.size()))
        this->fullQueryStats.minTetra = static_cast<int>(qS.tetrahedra.size());
    if(this->fullQueryStats.maxTetra < static_cast<int>(qS.tetrahedra.size()))
        this->fullQueryStats.maxTetra = static_cast<int>(qS.tetrahedra.size());
    this->fullQueryStats.avgTetra += static_cast<int>(qS.tetrahedra.size());

    if(this->fullQueryStats.min_tetra_compl_cont_leaf_num > qS.tetra_compl_cont_leaf_num)
        this->fullQueryStats.min_tetra_compl_cont_leaf_num = qS.tetra_compl_cont_leaf_num;
    if(this->fullQueryStats.max_tetra_compl_cont_leaf_num < qS.tetra_compl_cont_leaf_num)
        this->fullQueryStats.max_tetra_compl_cont_leaf_num = qS.tetra_compl_cont_leaf_num;
    this->fullQueryStats.avg_tetra_compl_cont_leaf_num += qS.tetra_compl_cont_leaf_num;

    if(this->fullQueryStats.min_tetra_compl_cont_bbox_num > qS.tetra_compl_cont_bbox_num)
        this->fullQueryStats.min_tetra_compl_cont_bbox_num = qS.tetra_compl_cont_bbox_num;
    if(this->fullQueryStats.max_tetra_compl_cont_bbox_num < qS.tetra_compl_cont_bbox_num)
        this->fullQueryStats.max_tetra_compl_cont_bbox_num = qS.tetra_compl_cont_bbox_num;
    this->fullQueryStats.avg_tetra_compl_cont_bbox_num += qS.tetra_compl_cont_bbox_num;

    //nodes statistics
    if(this->fullQueryStats.minNode > qS.numNode)
        this->fullQueryStats.minNode = qS.numNode;
    if(this->fullQueryStats.maxNode < qS.numNode)
        this->fullQueryStats.maxNode = qS.numNode;
    this->fullQueryStats.avgNode += qS.numNode;

    //leaves statistics
    if(this->fullQueryStats.minLeaf > qS.numLeaf)
        this->fullQueryStats.minLeaf = qS.numLeaf;
    if(this->fullQueryStats.maxLeaf < qS.numLeaf)
        this->fullQueryStats.maxLeaf = qS.numLeaf;
    this->fullQueryStats.avgLeaf += qS.numLeaf;

    //geometric tests statistics
    if(this->fullQueryStats.minGeometricTest > qS.numGeometricTest)
        this->fullQueryStats.minGeometricTest = qS.numGeometricTest;
    if(this->fullQueryStats.maxGeometricTest < qS.numGeometricTest)
        this->fullQueryStats.maxGeometricTest = qS.numGeometricTest;
    this->fullQueryStats.avgGeometricTest += qS.numGeometricTest;

    //local analysis of access-per-tetra
    int multiple = 0;
    int mult_counter = 0;
    int unique = 0;
    for(int_vect_iter iter=qS.access_per_tetra.begin(); iter!=qS.access_per_tetra.end(); ++iter)
    {
        if(*iter==1)
            unique++;
        else
        {
            multiple += *iter;
            mult_counter++;
        }
    }

    //unique tetrahedra access statistics
    if(this->fullQueryStats.minUniqueTetraAccess > unique)
        this->fullQueryStats.minUniqueTetraAccess = unique;
    if(this->fullQueryStats.maxUniqueTetraAccess < unique)
        this->fullQueryStats.maxUniqueTetraAccess = unique;
    this->fullQueryStats.avgUniqueTetraAccess += unique;

    //multiple tetrahedra access statistics
    if(multiple != 0)
    {
        if(this->fullQueryStats.minMultipleTetraAccess > multiple)
            this->fullQueryStats.minMultipleTetraAccess = multiple;
        if(this->fullQueryStats.maxMultipleTetraAccess < multiple)
            this->fullQueryStats.maxMultipleTetraAccess = multiple;
        this->fullQueryStats.avgMultipleTetraAccess += multiple;
    }

    // in the following of this function are computed statistics about the effectiveness of bounding box strategy
    //(1) box completely contains leaf node
    if(this->fullQueryStats.min_box_completely_contains_leaf_num > qS.box_completely_contains_leaf_num)
        this->fullQueryStats.min_box_completely_contains_leaf_num = qS.box_completely_contains_leaf_num;
    if(this->fullQueryStats.max_box_completely_contains_leaf_num < qS.box_completely_contains_leaf_num)
        this->fullQueryStats.max_box_completely_contains_leaf_num = qS.box_completely_contains_leaf_num;
    this->fullQueryStats.avg_box_completely_contains_leaf_num += qS.box_completely_contains_leaf_num;
    //(2) box completely contains a run bounding box
    if(this->fullQueryStats.min_box_completely_contains_bbox_num > qS.box_completely_contains_bbox_num)
        this->fullQueryStats.min_box_completely_contains_bbox_num = qS.box_completely_contains_bbox_num;
    if(this->fullQueryStats.max_box_completely_contains_bbox_num < qS.box_completely_contains_bbox_num)
        this->fullQueryStats.max_box_completely_contains_bbox_num = qS.box_completely_contains_bbox_num;
    this->fullQueryStats.avg_box_completely_contains_bbox_num += qS.box_completely_contains_bbox_num;
    //(3) box intersects a run bounding box
    if(this->fullQueryStats.min_box_intersect_bbox_num > qS.box_intersect_bbox_num)
        this->fullQueryStats.min_box_intersect_bbox_num = qS.box_intersect_bbox_num;
    if(this->fullQueryStats.max_box_intersect_bbox_num < qS.box_intersect_bbox_num)
        this->fullQueryStats.max_box_intersect_bbox_num = qS.box_intersect_bbox_num;
    this->fullQueryStats.avg_box_intersect_bbox_num += qS.box_intersect_bbox_num;
    //(4) box does not intersect a run bounding box
    if(this->fullQueryStats.min_box_no_intersect_bbox_num > qS.box_no_intersect_bbox_num)
        this->fullQueryStats.min_box_no_intersect_bbox_num = qS.box_no_intersect_bbox_num;
    if(this->fullQueryStats.max_box_no_intersect_bbox_num < qS.box_no_intersect_bbox_num)
        this->fullQueryStats.max_box_no_intersect_bbox_num = qS.box_no_intersect_bbox_num;
    this->fullQueryStats.avg_box_no_intersect_bbox_num += qS.box_no_intersect_bbox_num;
    //(4) number of tetra-in-box tests executed for the intersect case
    if(this->fullQueryStats.min_box_intersect_bbox_geom_tests_num > qS.box_intersect_bbox_geom_tests_num)
        this->fullQueryStats.min_box_intersect_bbox_geom_tests_num = qS.box_intersect_bbox_geom_tests_num;
    if(this->fullQueryStats.max_box_intersect_bbox_geom_tests_num < qS.box_intersect_bbox_geom_tests_num)
        this->fullQueryStats.max_box_intersect_bbox_geom_tests_num = qS.box_intersect_bbox_geom_tests_num;
    this->fullQueryStats.avg_box_intersect_bbox_geom_tests_num += qS.box_intersect_bbox_geom_tests_num;
    //(5) number of avoided tetra-in-xxx tests thanks to bounding boxes strategy
    if(this->fullQueryStats.min_avoided_tetra_geom_tests_num > qS.avoided_tetra_geom_tests_num)
        this->fullQueryStats.min_avoided_tetra_geom_tests_num = qS.avoided_tetra_geom_tests_num;
    if(this->fullQueryStats.max_avoided_tetra_geom_tests_num < qS.avoided_tetra_geom_tests_num)
        this->fullQueryStats.max_avoided_tetra_geom_tests_num = qS.avoided_tetra_geom_tests_num;
    this->fullQueryStats.avg_avoided_tetra_geom_tests_num += qS.avoided_tetra_geom_tests_num;
    //

    return hit_ratio;
}
