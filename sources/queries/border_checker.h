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

#ifndef BORDERCHECKER_H
#define BORDERCHECKER_H

#include <iostream>
#include <algorithm>
#include <map>

#include "utilities/sorting.h"
#include "tetrahedral_trees/node_v.h"
#include "tetrahedral_trees/node_t.h"

/**
 * @brief The Border_Checker class exploits the borders of a tetrahedral mesh indexed by a Tetrahedral tree
 * The mesh borders are needed when computing the windowed 3D Curvature.
 */
class Border_Checker
{
public:
    Border_Checker() {}
    /**
     * @brief A public procedure that visits recursively a tree and exploit the mesh borders.
     * This version is compatible with T-Ttrees and RT-Ttrees on which the spatial coherence has been exploited
     * as well as on Tetrahedral trees on which the spatial coherence is not exploited.
     *
     * @param n a N& parameter representing the current node
     * @param dom a Box& parameter representing the current domain of node n
     * @param level an integer argument representing the level of n in the hierarchy
     * @param mesh a Mesh& parameter representing the indexed tetrahedral mesh
     * @param division a D& parameter representing the spatial subdivision of the tree
     */
    template<class D> void calc_mesh_borders(Node_T &n, Box &dom, int level, Mesh &mesh, D &division);
    /**
     * @brief A public procedure that visits recursively a tree and exploit the mesh borders.
     * This version requires a tree on which the spatial coherence has been exploited, and is compatible with P-Ttrees.
     *
     * @param n a P_Node& parameter representing the current node
     * @param dom a Box& parameter representing the current domain of node n
     * @param level an integer argument representing the level of n in the hierarchy
     * @param mesh a Mesh& parameter representing the indexed tetrahedral mesh
     * @param division a D& parameter representing the spatial subdivision of the tree
     */
    template<class D> void calc_mesh_borders(Node_V &n, Box &dom, int level, Mesh &mesh, D &division);

private:
    /**
     * @brief A private procedure that exploits the mesh border in a leaf block.
     * This version requires a tree on which the spatial coherence has been exploited, and is compatible with P-Ttrees and PT-Ttrees.
     *
     * @param n a N& parameter representing the current node
     * @param mesh a Mesh& parameter representing the indexed tetrahedral mesh
     */
    void calc_mesh_borders(Node_V& n, Mesh& mesh);
    /**
     * @brief A private procedure that exploits the mesh border in a leaf block.
     * This version is compatible with T-Ttrees and RT-Ttrees on which the spatial coherence has been exploited
     * as well as on Tetrahedral trees on which the spatial coherence is not exploited.
     *
     * @param n a N& parameter representing the current node
     * @param dom a Box& parameter representing the current domain of node n
     * @param mesh a Mesh& parameter representing the indexed tetrahedral mesh
     */
    void calc_mesh_borders(Node_T& n, Box &dom, Mesh& mesh);
    /**
     * @brief A private procedure that sets the border of the sub-mesh indexed by the leaf block
     *
     * @param faces an array containing the tuples composed by the three indices forming a face and the tetrahedron in its co-boundary
     * @param mesh a Mesh& parameter representing the indexed tetrahedral mesh
     * @return true if we change the borders, false otherwise
     */
    bool set_mesh_borders(vector<triangle_tetrahedron_tuple> faces, Mesh& mesh);
    /**
     * @brief A private procedure that extracts the three triangular faces incident in a vertex in the boundary of a tetrahedron
     *
     * @param t a Tetrahedron& parameter, representing the current tetrahedron
     * @param t_id an integer parameter, representing the position index of t
     * @param v_pos an integer parameter representing the vertex position index in the boundary of t
     * @param faces an array containing the faces incident in the vertex
     */
    void get_incident_triangles(Tetrahedron &t, int t_id, int v_pos, vector<triangle_tetrahedron_tuple> &faces);
};

template<class D> void Border_Checker::calc_mesh_borders(Node_T &n, Box &dom, int level, Mesh &mesh, D &division)
{
    if (n.is_leaf())
    {
        this->calc_mesh_borders(n,dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->calc_mesh_borders(*n.get_son(i),son_dom,son_level,mesh,division);
        }
    }
}

template<class D> void Border_Checker::calc_mesh_borders(Node_V &n, Box &dom, int level, Mesh &mesh, D &division)
{
    if (n.is_leaf())
    {
        this->calc_mesh_borders(n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            this->calc_mesh_borders(*n.get_son(i),dom,level,mesh,division);
        }
    }
}

#endif // BORDERCHECKER_H
