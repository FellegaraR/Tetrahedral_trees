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

#include "node_t.h"

void Node_T::get_v_range(int &v_start, int &v_end, Box& dom, Mesh& mesh)
{
    v_start = v_end = -1;

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& runIt = itPair.first;
        Tetrahedron& t = mesh.get_tetrahedron(*runIt);
        for(int v=0; v<t.vertices_num(); v++)
        {
            int v_id = t.TV(v);
            if((v_start == -1 && v_end == -1) || (v_id < v_start || v_id >= v_end))
                //I check only the vertices outside the current range
            {
                if(dom.contains(mesh.get_vertex(v_id),mesh.get_domain().get_max()))
                {
                    if(v_start == -1 || v_start > v_id)
                        v_start = v_id;
                    if(v_end == -1 || v_end <= v_id)
                        v_end = v_id+1;// v_end is +1 because it represents the first index outside the leaf
                }
            }
        }
    }
}
