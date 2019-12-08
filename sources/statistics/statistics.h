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

#ifndef _STATISTICS_H
#define	_STATISTICS_H

#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <sstream>

#include "basic_types/mesh.h"
#include "statistics/query_statistics.h"
#include "statistics/index_statistics.h"
#include "statistics/full_query_statistics.h"
#include "tetrahedral_trees/node_v.h"
#include "tetrahedral_trees/node_t.h"
#include "geometry/geometry_wrapper.h"
#include "io/writer.h"
#include "io/reader.h"

using namespace std;
///A class providing an interface for computing statistics over a spatial index and a several queries executed on it
class Statistics {
public:
    ///A constructor method
    Statistics()
    {
        this->indexStats = IndexStatistics();
        this->fullQueryStats = FullQueryStatistics();
    }
    ///A copy-constructor method
    Statistics(const Statistics& orig)
    {
        this->indexStats = orig.indexStats;
        this->fullQueryStats = orig.fullQueryStats;
    }
    ///A destructor method
    virtual ~Statistics() {}
    ///A public method that computes index statistics over a given tree
    /*!
     * This method prints the results on standard output and it doesn't write the statistics into a file
     *
     * \param tree a T& argument, representing the tree where the statistics are executed
     * \param reindex a boolean saying if the spatial coherence has been exploited or not on the index
     */
    template<class T> void get_index_statistics(T& tree, bool reindex);
    ///A public method that updates the local query statistics container with the statistics gathered after a query
    /*!
     * \param qS a QueryStatistics& argument containing the results of a query
     */
    int compute_queries_statistics(QueryStatistics &qS);
    /**
     * @brief A public procedure that returns the local variable containing the queries statistics
     * @return the statical summary of queries
     */
    inline FullQueryStatistics& get_query_statistics() { return this->fullQueryStats; }

private:
    ///A private variable in which the index statistics are saved
    IndexStatistics indexStats;
    ///A private variable in which the queries statistics are saved
    FullQueryStatistics fullQueryStats;
    ///A private method that simulates a tree visit and computes the index statistics
    /*!
     * \param n a N& argument, representing the current node to visit
     * \param dom a Box& argument, representing the domain of the current node
     * \param level an integer representing the node level in the hierarchy
     * \param mesh a Mesh& argument representing the tetrahedral mesh
     * \param division a D& argument representing the adopted spatial subdivision
     * \param reindex a boolean saying if the spatial coherence has been exploited or not on the index
     */
    template<class N,class D> void visit_tree(N& n, Box& dom, int level, Mesh& mesh, D& division, bool reindex);
    ///A private method that initializes the vector used for the index statistics operation
    /*!
     * \param mesh a Mesh& argument, represents the tetrahedral mesh
     */
    inline void init_vector(Mesh& mesh) { this->indexStats.num_leaves_for_tetra.assign(mesh.get_num_tetrahedra(),0); }
    ///A private method that computes some index statistics that it is not possible to compute during the visit
    void calc_remaining_index_statistics();
    ///A private method that checks common inconsistencies of index statistics operation
    void check_inconsistencies();

    /**
     * \brief A private procedure that computes the index statistics in a leaf (wrapper for T-Ttrees and RT-Ttrees)
     *
     * \param n a N& argument, representing the current node to visit
     * \param dom a Box& argument, representing the domain of the current node
     * \param mesh a Mesh& argument representing the tetrahedral mesh
     * \param reindex a boolean saying if the spatial coherence has been exploited or not on the index
     */
    void compute_leaf_statistics(Node_T& n,Box& dom,Mesh& mesh,bool reindex);
    /**
     * \brief A private procedure that computes the index statistics in a leaf (wrapper for PT-Ttrees)
     *
     * \param n a N& argument, representing the current node to visit
     * \param dom a Box& argument, representing the domain of the current node
     * \param mesh a Mesh& argument representing the tetrahedral mesh
     * \param reindex a boolean saying if the spatial coherence has been exploited or not on the index
     */
    void compute_leaf_statistics(Node_V& n,Box& dom,Mesh& mesh,bool reindex);
    /**
     * @brief A private procedure that the statistics connected to the vertices indexed by a leaf
     * @param num_vertex an integer referring to the number of vertices indexed by a leaf
     */
    void set_leaf_vertices_stats(int num_vertex);
    /**
     * @brief A private procedure that the statistics connected to the vertices indexed by a leaf
     * @param num_t_completely an integer containing the number of tetrahedra completely indexed by the leaf
     * @param num_t_partially an integer containing the number of tetrahedra partially indexed by the leaf
     * @param num_t_overlapping an integer containing the number of tetrahedra simply crossing (i.e., no indexed vertices) by the leaf
     */
    void set_leaf_tetrahedra_stats(int num_t_completely, int num_t_partially, int num_t_overlapping);
};

template<class T> void Statistics::get_index_statistics(T& tree, bool reindex)
{
    init_vector(tree.get_mesh());
    visit_tree(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_decomposition(),reindex);
    calc_remaining_index_statistics();
    check_inconsistencies();
    Writer::write_tree_stats(this->indexStats);
    return;
}

template<class N,class D> void Statistics::visit_tree(N &n, Box& dom, int level, Mesh& mesh, D& division, bool reindex)
{
    this->indexStats.numNode++;

    if(n.is_leaf())
    {
        if(this->indexStats.minTreeDepth==-1 || this->indexStats.minTreeDepth > level)
            this->indexStats.minTreeDepth = level;
        if(this->indexStats.maxTreeDepth < level)
            this->indexStats.maxTreeDepth = level;
        this->indexStats.avgTreeDepth += level;

        this->compute_leaf_statistics(n,dom,mesh,reindex);
    }
    else
    {
        int son_level = level + 1;

        for(int i=0;i<division.son_number();i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            this->visit_tree(*n.get_son(i),son_dom,son_level,mesh,division,reindex);
        }
    }
    return;
}

#endif	/* _STATISTICS_H */

