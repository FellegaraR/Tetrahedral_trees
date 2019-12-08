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

#ifndef TOPOLOGICALQUERIES_H
#define TOPOLOGICALQUERIES_H

#include <set>
#include <map>
#include <boost/dynamic_bitset.hpp>

#include "basic_types/vertex.h"
#include "basic_types/tetrahedron.h"
#include "basic_types/mesh.h"
#include "basic_types/box.h"
#include "tetrahedral_trees/node_v.h"
#include "tetrahedral_trees/node_t.h"
#include "geometry/geometry_wrapper.h"
#include "geometry/geometry_distortion.h"

using namespace std;

/**
 * @brief The Topological_Queries class provides an interface for executing topological queries on the Tetrahedral trees
 * NOTA: of this class are documented only the public procedures.
 */
class Topological_Queries
{
public:
    /**
     * @brief A constructor method
     */
    Topological_Queries() { }
    /**
     * @brief A copy-constructor method
     */
    Topological_Queries(const Topological_Queries&) {}
    /**
     * @brief A destructor method
     */
    virtual ~Topological_Queries() {}
    ///A public method that excutes windowed VTop queries, reading the boxes from file
    /*!
     * This method prints the results on standard output
     *
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     * \param query_path a string argument, representing the file path of the query input
     * \param reindexed a boolean, true if the index and the mesh are spatially reordered
     */
    template<class N, class D> void windowed_VT(N &n, Box &dom, Mesh &mesh, D &division, string query_path, bool reindexed);
    ///A public method that excutes windowed curvature queries, reading the boxes from file
    /*!
     * This method prints the results on standard output
     *
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     * \param query_path a string argument, representing the file path of the query input
     * \param reindexed a boolean, true if the index and the mesh are spatially reordered
     */
    template<class N, class D> void windowed_Distortion(N &n, Box &dom, Mesh &mesh, D &division, string query_path, bool reindexed);
    ///A public method that excutes windowed TT queries, reading the boxes from file
    /*!
     * This method prints the results on standard output
     *
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     * \param query_path a string argument, representing the file path of the query input
     */
    template<class N, class D> void windowed_TT(N &n, Box &dom, Mesh &mesh, D &division, string query_path);
    ///A public method that excutes linearized TT queries, reading the boxes from file
    /*!
     * This method prints the results on standard output
     *
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     * \param query_path a string argument, representing the file path of the query input
     */
    template<class N, class D> void linearized_TT(N &n, Box &dom, Mesh &mesh, D &division, string query_path);

    template<class N, class D> void batched_VT(N &n, Box &dom, Mesh &mesh, D &division, bool reindex);
    template<class N, class D> void batched_TT(N &n, Mesh &mesh, D &division);

private:
    // windowed VT - auxiliary functions
    template<class D> void windowed_VT(Node_T &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &vt);
    template<class N, class D> void windowed_VT_no_reindex(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &vt);
    template<class D> void windowed_VT(Node_V &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &vt);
    void windowed_VT_Leaf(Node_T& n, Box &dom, Box &b, Mesh& mesh, map<int,vector<int> > &vt);
    void windowed_VT_Leaf(Node_V& n, Box &b, Mesh& mesh, map<int,vector<int> > &vt);
    template<class N> void windowed_VT_Leaf_no_reindex(N& n, Box &dom, Box &b, Mesh& mesh, map<int,vector<int> > &vt);
    void update_resulting_VT(int v, int t, map<int,vector<int> > &vt);
    // windowed distortion - auxiliary functions
    template<class D> void windowed_Distortion(Node_T &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,double> &dist);
    template<class N, class D> void windowed_Distortion_no_reindex(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,double> &dist);
    template<class D> void windowed_Distortion(Node_V &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,double> &dist);
    void windowed_Distortion_Leaf(Node_T& n, Box &dom, Box &b, Mesh& mesh, map<int,double> &dist);
    template<class N> void windowed_Distortion_Leaf_no_reindex(N& n, Box &dom, Box &b, Mesh& mesh, map<int,double> &dist);
    void windowed_Distortion_Leaf(Node_V &n, Box &b, Mesh& mesh, map<int,double> &dist);
    void update_resulting_distortion(int v, double d, map<int,double> &dist);
    void finalize_Distortion_Leaf(int v_start, vector<vector<int> > &all_vt, vector<double> &local_distortion, boost::dynamic_bitset<> &isVBorder, Mesh& mesh, map<int,double> &dist);
    // windowed TT - auxiliary functions
    template<class N, class D> void windowed_TT(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &tt, bm::bvector<> &checkTetra);
    template<class N> void windowed_TT_Leaf_test(N& n, Box &b, Mesh& mesh, map<int,vector<int> > &tt, bm::bvector<> &checkTetra);
    template<class N> void windowed_TT_Leaf_add(N& n, Mesh& mesh, map<int,vector<int> > &tt, bm::bvector<> &checkTetra);
    void add_faces(int t_id, vector<triangle_tetrahedron_tuple> &faces, Mesh &mesh, map<int,vector<int> >::const_iterator &iter, map<int,vector<int> > &tt);
    void pair_adjacent_tetrahedra(vector<triangle_tetrahedron_tuple> &faces, Mesh &mesh, map<int, vector<int> > &tt);
    void update_resulting_TT(int pos, int t1, int t2, map<int,vector<int> > &tt);
    void init_TT_entry(int t1, map<int,vector<int> > &tt);
    // linearized TT - auxiliary functions
    template<class N, class D> void linearized_TT(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &tt, bm::bvector<> &checkTetra);
    template<class N> void linearized_TT_Leaf(N& n, Box &b, Mesh& mesh, map<int,vector<int> > &tt, bm::bvector<> &checkTetra);
    // windowed and linearized TT auxiliary function
    void finalize_TT_Leaf(vector<triangle_tetrahedron_tuple> &faces, map<int,vector<int> > &tt, Mesh &mesh);

    template<class N, class D> void batched_VT_visit(N &n, Box &dom, int level, Mesh &mesh, D &division, bool stats, int &max_entries);
    template<class N, class D> void batched_VT_no_reindex(N &n, Box &dom, int level, Mesh &mesh, D &division, bool stats, int &max_entries);
    void batched_VT_leaf(Node_T &n, Box &dom, Mesh &mesh, bool stats, int &max_entries);
    void batched_VT_leaf(Node_V &n, Box &, Mesh &mesh, bool stats, int &max_entries);
    void batched_VT_no_reindex_leaf(Node_T &n, Box &dom, Mesh &mesh, bool stats, int &max_entries);
    void batched_VT_no_reindex_leaf(Node_V &n, Box &dom, Mesh &mesh, bool stats, int &max_entries);

    template<class N, class D> void batched_TT_visit(N &n, Mesh &mesh, D &division, vector<vector<int> > &tt, bool stats, int &max_entries);
    template<class N> void batched_TT_leaf(N &n, Mesh &mesh, vector<vector<int> > &tt, bool stats, int &max_entries);
};

#include "topological_queries_windowed.h"
#include "topological_queries_batched.h"

#endif // TOPOLOGICALQUERIES_H
