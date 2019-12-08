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

#include "sorting.h"
#include <algorithm>

void sorting_vertices(vector<vertex_tetrahedron_pair> &vert_vec, const vector<int> &tetrahedra, Mesh &mesh)
{
    int k=0;
    //create the array of pairs vertex-tetrahedron
    for(unsigned int t=0;t<tetrahedra.size();t++)
    {
        Tetrahedron &tet = mesh.get_tetrahedron(tetrahedra[t]);
        for(int i=0;i<tet.vertices_num();i++)
        {
            vert_vec[k].v = tet.TV(i);
            vert_vec[k].t = tetrahedra[t];
            k++;
        }
    }
    //order the array
    std::sort(vert_vec.begin(),vert_vec.end());
}

void sorting_faces(vector< triangle_tetrahedron_tuple >& faces) { std::sort(faces.begin(),faces.end()); }
