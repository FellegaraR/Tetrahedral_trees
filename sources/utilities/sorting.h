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

#ifndef _SORTING_H
#define	_SORTING_H

#include "basic_types/mesh.h"
#include "sorting_structure.h"
#include <vector>
#include <set>
#include <list>

///A protected method that creates a list of ordered pairs of vertices-incident tetrahedra.
/*!
 * The ordered list has common vertex between tetrahedra in adiacent position of the list
 *
 * \param vert_vec an vector<vertex_tetrahedron_pair>& argument, an array that will contain the sorted array of tuples
 * \param tetrahedra a const vector<int>& argument, the list of tetrahedra from which extract the pairs
 * \param mesh a Mesh& argument representing the tetrahedral mesh
 */
void sorting_vertices(vector<vertex_tetrahedron_pair> &vert_vec, const vector<int> &tetrahedra, Mesh &mesh);

/**
 * @brief A procedure that sort an array of triangle_tetrahedron tuples
 * @param faces
 */
void sorting_faces(vector<triangle_tetrahedron_tuple>& faces);


#endif	/* _SORTING_H */
