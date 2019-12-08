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

#include "border_checker.h"

void Border_Checker::calc_mesh_borders(Node_V &n, Mesh &mesh)
{
    if(n.get_v_array_size() == 0)
        return;

    vector< vector<triangle_tetrahedron_tuple> > all_faces;
    vector<triangle_tetrahedron_tuple> faces;
    all_faces.assign(n.get_v_end()-n.get_v_start(),faces);

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& t = mesh.get_tetrahedron(*tet_id);
        for(int j=0; j<t.vertices_num(); j++)
        {
            int real_index = t.TV(j);
            if(n.indexes_vertex(real_index))
            {
                //we insert the three triangular faces incident in vertex v_pos
                get_incident_triangles(t,*tet_id,j,all_faces[real_index-n.get_v_start()]);
            }
        }
    }

    for(vector< vector<triangle_tetrahedron_tuple> >::iterator iter=all_faces.begin(); iter!=all_faces.end(); ++iter)
    {
        if(iter->size() > 0)
        {
            this->set_mesh_borders(*iter,mesh);
        }
    }
}

void Border_Checker::calc_mesh_borders(Node_T &n, Box &dom, Mesh& mesh)
{
    map< int, vector<triangle_tetrahedron_tuple> > all_faces;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& t = mesh.get_tetrahedron(*tet_id);
        for(int j=0; j<t.vertices_num(); j++)
        {
            int real_index = t.TV(j);
            if(dom.contains(mesh.get_vertex(real_index),mesh.get_domain().get_max()))
            {
                //we insert the three triangular faces incident in vertex v_pos
                vector<triangle_tetrahedron_tuple> faces;
                get_incident_triangles(t,*tet_id,j,faces);

                map< int, vector<triangle_tetrahedron_tuple> >::iterator iter = all_faces.find(real_index);
                if(iter == all_faces.end())
                {
                    all_faces.insert(make_pair(real_index,faces));
                }
                else
                {
                    for(vector<triangle_tetrahedron_tuple>::iterator f=faces.begin(); f!=faces.end(); ++f)
                        iter->second.push_back(*f);
                }
            }
        }
    }

    for(map< int, vector<triangle_tetrahedron_tuple> >::iterator iter=all_faces.begin(); iter!=all_faces.end(); ++iter)
    {
        if(iter->second.size() > 0)
        {
            this->set_mesh_borders(iter->second,mesh);
        }
    }
}

bool Border_Checker::set_mesh_borders(vector<triangle_tetrahedron_tuple> faces, Mesh &mesh)
{
    bool borderChange = false;

    sorting_faces(faces);

    unsigned int j=0;
    while(j<faces.size())
    {
        if(j+1<faces.size())
        {
            if(faces[j] != faces[j+1])
            {
                Tetrahedron& tet = mesh.get_tetrahedron(faces[j].t);

                for(int v=0;v<tet.vertices_num();v++)
                {
                    int v_ind = tet.TV(v);
                    if(faces[j].has_not(v_ind))
                    {
                        if(tet.is_border_face(v)) //identified a triangular face on the mesh border
                            borderChange = true;
                        tet.setTV(v,-v_ind); //we flag with a negative index the vertex opposite to the face on the border
                        break;
                    }
                }

                j++;

            }
            else if(faces[j] == faces[j+1])
            {
                j+=2;
            }
        }
        else
        {
            Tetrahedron& tet = mesh.get_tetrahedron(faces[j].t);

            for(int v=0;v<tet.vertices_num();v++)
            {
                int v_ind = tet.TV(v);
                if(faces[j].has_not(v_ind))
                {
                    if(tet.is_border_face(v))
                        borderChange = true;
                    tet.setTV(v,-v_ind);
                    break;
                }
            }

            j++;

        }

        if(j>=faces.size())
            break;
    }

    return borderChange;
}

void Border_Checker::get_incident_triangles(Tetrahedron &t, int t_id, int v_pos, vector<triangle_tetrahedron_tuple> &faces)
{
    //we insert the three triangular faces incident in vertex v_pos
    triangle_tetrahedron_tuple new_item;

    for(int i=1; i<t.vertices_num(); i++)
    {
        for(int j=i+1; j<t.vertices_num(); j++)
        {
            new_item.sort_and_set(abs(t.TV(v_pos)),abs(t.TV((v_pos+i)%t.vertices_num())),abs(t.TV((v_pos+j)%t.vertices_num())),t_id);
            faces.push_back(new_item);
        }
    }
}
