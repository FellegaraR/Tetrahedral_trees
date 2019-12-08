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

#ifndef GEOMETRY_WRAPPER_H
#define GEOMETRY_WRAPPER_H

#include "basic_types/mesh.h"
#include "geometry.h"
/**
 * @brief The Geometry_Wrapper class provides an interface for executing geometric tests for generating trees and answering queries
 */
class Geometry_Wrapper : public Geometry
{
public:
    /**
     * @brief A public static method that computes the centroid of a tetrahedron
     *
     * @param t_id an integer representing the position index of the tetrahedron
     * @param p a Point& that will contains the centroid coordinates
     * @param mesh a Mesh&, the tetrahedral mesh
     */
    static void get_tetrahedron_centroid(int t_id, Point& p, Mesh &mesh);
    /**
     * @brief A public static method that computes the point-in-tetrahedron geometric test
     *
     * @param t_id an integer representing the position index of the tetrahedron
     * @param p a Point& representing the point to test
     * @param mesh a Mesh&, the tetrahedral mesh
     * @return true if the point is contained in the tetrahedron, false otherwise
     */
    static bool point_in_tetra(int t_id, Point& point, Mesh &mesh);
    /**
     * @brief A public static method that computes the tetrahedron-in-box geometric test
     * NOTA: this procedure is used during the generation process of a tree.
     * It considers a tetrahedron internal even if it has only a vertex inside the box
     *
     * @param t_id an integer representing the position index of the tetrahedron
     * @param box a Box& representing the domain of a block
     * @param mesh a Mesh&, the tetrahedral mesh
     * @return true if exists an intersection between t_id and box, false otherwise
     */
    static bool tetra_in_box_build(int t_id, Box& box, Mesh& mesh);
    /**
     * @brief A public static method that computes the tetrahedron-in-box geometric test
     * NOTA: this procedure is used in box queries.
     * It considers all the faces of the box open as we want to output the tetrahedra
     * with a real intersection with the box and not only a vertex adjacent to one of box faces.
     *
     * @param t_id an integer representing the position index of the tetrahedron
     * @param box a Box& representing the domain of a block
     * @param mesh a Mesh&, the tetrahedral mesh
     * @return true if exists a real intersection between t_id and box, false otherwise
     */
    static bool tetra_in_box(int t_id, Box& box, Mesh& mesh);
    /**
     * @brief A public static method that computes the line-in-box geometric tests
     * NOTA: the procedure is used to check if a line intersects the domain of a box node in the hierarchy
     *
     * @param v1 a Point& argument, representing the first extreme of the line
     * @param v2 a Point& argument, representing the second extreme of the line
     * @param box a Box& representing the domain of a block
     * @return true if the line intersects the block, false otherwise
     */
    static bool line_in_box(const Point &v1, const Point &v2, Box& box); //only for line in leaf_domain test
    /**
     * @brief A public static method that computes the line-in-box geometric tests
     * NOTA: the procedure is used to check if a line intersects the bounding box of a run of tetrahedra
     *
     * @param v1 a Point& argument, representing the first extreme of the line
     * @param v2 a Point& argument, representing the second extreme of the line
     * @param bb a Box& representing the domain of a block
     * @return true if the line intersects the bounding box, false otherwise
     */
    static bool line_in_bounding_box(const Point &v1, const Point &v2, Box& bb); //only for line in run bounding box test
    /**
     * @brief A public static method that computes the line-in-tetrahedron geometric tests
     *
     * @param v1 a Point& argument, representing the first extreme of the line
     * @param v2 a Point& argument, representing the second extreme of the line
     * @param t_id an integer representing the position index of the tetrahedron
     * @param mesh a Mesh&, the tetrahedral mesh
     * @return true if the line intersects the tetrahedron, false otherwise
     */
    static bool line_in_tetra(const Point& v1, const Point& v2, int t_id, Mesh &mesh); // same algorithm without distance computation
    /**
     * @brief A public static method that reorder the triangular faces of the mesh tetrahedra
     *
     * @param mesh a Mesh&, the tetrahedral mesh
     */
    static void set_faces_ordering(Mesh &mesh);

private:
    //used in line_in_tetra
    static void ordered_TF(Tetrahedron &t, int pos, vector<int> &f);
    //used in set_faces_ordering
    static void set_face_orientation(Tetrahedron &tet, Mesh &mesh);
    static int four_point_turn_wrapper(const Point &v0, const Point &v1, const Point &v2, const Point &op);
};

#endif // GEOMETRY_WRAPPER_H
