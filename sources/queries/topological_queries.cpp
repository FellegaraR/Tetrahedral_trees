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

#include "topological_queries.h"

void Topological_Queries::windowed_VT_Leaf(Node_T &n, Box &dom, Box &b, Mesh& mesh, map<int,vector<int> > &vt)
{
    int v_start;
    int v_end;

    n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    vector<vector<int> > local_vt;  // local smaller structure... in the end inserted into the global map..
    local_vt.assign(v_end-v_start,vector<int>());

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<tet.vertices_num(); v++)
        {
            int real_v_index = tet.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (n.indexes_vertex(v_start,v_end,real_v_index) && b.contains_with_all_closed_faces(mesh.get_vertex(real_v_index)))
                local_vt[real_v_index-v_start].push_back(*tet_id);
        }
    }

    // finally put the local VT into the global associative array
    for(unsigned i=0; i<local_vt.size(); i++)
    {
        if(local_vt[i].size() > 0)
        {
            int real_v_index = i+v_start;
            vt.insert(make_pair(real_v_index,local_vt[i]));
        }
    }
}

void Topological_Queries::windowed_VT_Leaf(Node_V& n, Box &b, Mesh& mesh, map<int,vector<int> > &vt)
{
    if(n.get_v_array_size() == 0)
        return;

    int v_start = n.get_v_start();
    int v_end = n.get_v_end();

//    if(v_start == v_end) //no internal vertices..
//        return;

    vector<vector<int> > local_vt;  // local smaller structure... in the end inserted into the global map..
    local_vt.assign(v_end-v_start,vector<int>());

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<tet.vertices_num(); v++)
        {
            int real_v_index = tet.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (n.indexes_vertex(real_v_index) && b.contains_with_all_closed_faces(mesh.get_vertex(real_v_index)))
                local_vt[real_v_index-v_start].push_back(*tet_id);
        }
    }

    // finally put the local VT into the global associative array
    for(unsigned i=0; i<local_vt.size(); i++)
    {
        if(local_vt[i].size() > 0)
        {
            int real_v_index = i+v_start;
            vt.insert(make_pair(real_v_index,local_vt[i]));
        }
    }
}

void Topological_Queries::update_resulting_VT(int v, int t, map<int,vector<int> > &vt)
{
    map<int,vector<int> >::iterator iter = vt.find(v);
    if(iter == vt.end())
    {
        // to insert into the map
        vector<int> local;
        local.push_back(t);
        vt.insert(make_pair(v,local));
    }
    else
    {
        // add the tetra to the current
        iter->second.push_back(t);
    }
}

void Topological_Queries::windowed_Distortion_Leaf(Node_T& n, Box &dom, Box &b, Mesh& mesh, map<int,double> &dist)
{
    int v_start;
    int v_end;

    n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range

    if(v_start == v_end) //no internal vertices..
        return;

    vector<vector<int> > all_vt;
    all_vt.assign(v_end-v_start,vector<int>());

    vector<double> local_distortion;  // local smaller structure... in the end inserted into the global map..
    local_distortion.assign(v_end-v_start,0);

    boost::dynamic_bitset<> isVBorder(v_end - v_start);

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<tet.vertices_num(); v++)
        {
            int real_v_index = tet.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (n.indexes_vertex(v_start,v_end,real_v_index) && b.contains_with_all_closed_faces(mesh.get_vertex(real_v_index)))
            {
                //add the current tetrahedron to the vt of the internal vertex
                all_vt[real_v_index-v_start].push_back(*tet_id);
                //for now the distortions array is global but we can transform it locally (this is a debug choice)
                double partial_distortion = Geometry_Distortion::get_trihedral_angle(tet,real_v_index,mesh);
                local_distortion[real_v_index-v_start] += partial_distortion;
                //se non so se il vertice è di bordo.. allora provo a capirlo
                if(!isVBorder[real_v_index - v_start])
                {
                    //controllo se almeno uno degli altri tre è negativo e nel caso metto il vertice corrente come di bordo
                    for(int j=1; j<tet.vertices_num(); j++)
                    {
                        if(tet.is_border_face((j+v)%tet.vertices_num()))
                        {
                            isVBorder[real_v_index - v_start] = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    finalize_Distortion_Leaf(v_start,all_vt,local_distortion,isVBorder,mesh,dist);
}

void Topological_Queries::windowed_Distortion_Leaf(Node_V &n, Box &b, Mesh& mesh, map<int,double> &dist)
{
    if(n.get_v_array_size() == 0)
        return;

    int v_start = n.get_v_start();
    int v_end = n.get_v_end();

//    if(v_start == v_end) //no internal vertices..
//        return;

    vector<vector<int> > all_vt;
    all_vt.assign(v_end-v_start,vector<int>());

    vector<double> local_distortion;  // local smaller structure... in the end inserted into the global map..
    local_distortion.assign(v_end-v_start,0);

    boost::dynamic_bitset<> isVBorder(v_end - v_start);

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<tet.vertices_num(); v++)
        {
            int real_v_index = tet.TV(v);
            //if a vertex has the partial vt != from 0 then must be into the search box...
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (n.indexes_vertex(real_v_index) && b.contains_with_all_closed_faces(mesh.get_vertex(real_v_index)))
            {
                //add the current tetrahedron to the vt of the internal vertex
                all_vt[real_v_index-v_start].push_back(*tet_id);
                //for now the distortions array is global but we can transform it locally (this is a debug choice)
                double partial_distortion = Geometry_Distortion::get_trihedral_angle(tet,real_v_index,mesh);
                local_distortion[real_v_index-v_start] += partial_distortion;
                //se non so se il vertice è di bordo.. allora provo a capirlo
                if(!isVBorder[real_v_index - v_start])
                {
                    //controllo se almeno uno degli altri tre è negativo e nel caso metto il vertice corrente come di bordo
                    for(int j=1; j<4; j++)
                    {
                        if(tet.is_border_face((j+v)%tet.vertices_num()) < 0)
                        {
                            isVBorder[real_v_index - v_start] = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    finalize_Distortion_Leaf(v_start,all_vt,local_distortion,isVBorder,mesh,dist);
}

void Topological_Queries::update_resulting_distortion(int v, double d, map<int,double> &dist)
{
    map<int,double>::iterator iter = dist.find(v);
    if(iter == dist.end())
    {
        // to insert into the map
        dist.insert(make_pair(v,d));
    }
    else
    {
        // add the tetra to the current
        iter->second = iter->second + d;
    }
}

void Topological_Queries::finalize_Distortion_Leaf(int v_start, vector<vector<int> > &all_vt, vector<double> &local_distortion, boost::dynamic_bitset<> &isVBorder, Mesh& mesh, map<int,double> &dist)
{
    for(unsigned i=0; i<all_vt.size(); i++)
    {
        if(all_vt[i].size() > 0) //I have an internal vertex that is inside the search box
        {
            int real_v_index = i+v_start;

            if(isVBorder[i])
            {
                local_distortion[i] = - local_distortion[i];
                for(int_vect_iter t = all_vt[i].begin(); t != all_vt[i].end(); ++t)
                {
                    Tetrahedron& tet = mesh.get_tetrahedron(*t);
                    double partial_distortion = Geometry_Distortion::get_trihedral_angle_3D(tet,real_v_index,mesh);
                    local_distortion[i] += partial_distortion;
                }
            }
            else
            {
                local_distortion[i] = (4*PI - local_distortion[i]);
            }

            dist.insert(make_pair(real_v_index,local_distortion[i]));
        }
    }
}

void Topological_Queries::add_faces(int t_id, vector<triangle_tetrahedron_tuple> &faces, Mesh &mesh, map<int,vector<int> >::const_iterator &iter, map<int, vector<int> > &tt)
{
    Tetrahedron& tet = mesh.get_tetrahedron(t_id);

    triangle_tetrahedron_tuple face;
    if(iter==tt.end()) // no entries in tt map
    {
        for(int v=0; v<4; v++)
        {
            tet.face_tuple(v,face,t_id);
            faces.push_back(face);
        }
    }
    else
    {
        const vector<int> &partial_tt = iter->second;
        for(unsigned i=0; i<partial_tt.size(); i++)
        {
            if(partial_tt[i]==-1) //the adj is unset
            {
                tet.face_tuple(i,face,t_id);
                faces.push_back(face);
            }
        }
    }
}

void Topological_Queries::pair_adjacent_tetrahedra(vector<triangle_tetrahedron_tuple> &faces, Mesh &, map<int, vector<int> > &tt)
{
    unsigned j=0;
    while(j<faces.size())
    {
        if(j+1<faces.size())
        {
            if(faces[j] == faces[j+1])
            {
                update_resulting_TT(faces[j].f_pos,faces[j].t,faces[j+1].t,tt);
                update_resulting_TT(faces[j+1].f_pos,faces[j+1].t,faces[j].t,tt);
                j+=2;
            }
            else
            {
                j++;
            }
        }
        else
        {
            j++;
        }

        if(j>=faces.size())
            break;
    }
}

void Topological_Queries::update_resulting_TT(int pos, int t1, int t2, map<int,vector<int> > &tt)
{
    map<int,vector<int> >::iterator iter = tt.find(t1);
    if(iter == tt.end())
    {
        cout<<"[update_resulting_TT] something wrong goes here..."<<endl;
        int a; cin>>a;
    }
    else
    {
        // add the tetra to the current
        iter->second[pos] = t2;
    }
}

void Topological_Queries::init_TT_entry(int t1, map<int,vector<int> > &tt)
{
    vector<int> local;
    local.assign(4,-1);
    tt.insert(make_pair(t1,local));
}

void Topological_Queries::finalize_TT_Leaf(vector<triangle_tetrahedron_tuple> &faces, map<int,vector<int> > &tt, Mesh &mesh)
{
    // (4) order the faces array
    sorting_faces(faces);
    // (5) set the adjacencies on these faces
    pair_adjacent_tetrahedra(faces,mesh,tt);
}

void Topological_Queries::batched_VT_leaf(Node_V &n, Box &, Mesh &mesh, bool stats, int &max_entries)
{
    // here we have a reindexed index thus, if there are no vertices indexed the array size is zero
    if(n.get_v_array_size() == 0)
        return;

    int v_start = n.get_v_start();
    int v_end = n.get_v_end();
//    if(v_start == v_end) //no internal vertices..
//        return;

    vector<vector<int> > local_vt;  // local smaller structure... in the end inserted into the global map..
    local_vt.assign(v_end-v_start,vector<int>());

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<4; v++)
        {
            int real_v_index = tet.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (n.indexes_vertex(real_v_index))
                local_vt[real_v_index-v_start].push_back(*tet_id);
        }
    }

    if(stats)
    {
        int entries = 0;

        for(vector<vector<int> >::iterator it = local_vt.begin(); it != local_vt.end(); ++it)
        {
            entries += it->size();
        }

        if(max_entries < entries)
            max_entries = entries;
    }
}

void Topological_Queries::batched_VT_leaf(Node_T &n, Box &dom, Mesh &mesh, bool stats, int &max_entries)
{
    int v_start;
    int v_end;

    n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    vector<vector<int> > local_vt;  // local smaller structure... in the end inserted into the global map..
    local_vt.assign(v_end-v_start,vector<int>());

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<tet.vertices_num(); v++)
        {
            int real_v_index = tet.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (n.indexes_vertex(v_start,v_end,real_v_index))
                local_vt[real_v_index-v_start].push_back(*tet_id);
        }
    }


    if(stats)
    {
        int entries = 0;

        for(vector<vector<int> >::iterator it = local_vt.begin(); it != local_vt.end(); ++it)
        {
            entries += it->size();
        }

        if(max_entries < entries)
            max_entries = entries;
    }
}

void Topological_Queries::batched_VT_no_reindex_leaf(Node_V &n, Box &dom, Mesh &mesh, bool stats, int &max_entries)
{
    if(n.get_v_array_size() == 0)
        return; // no vertices.. skip the current leaf block

    map<int,vector<int> > local_vt;  // local smaller structure... in the end inserted into the global map..

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<tet.vertices_num(); v++)
        {
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (dom.contains(mesh.get_vertex(tet.TV(v)),mesh.get_domain().get_max()))
                update_resulting_VT(tet.TV(v),*tet_id,local_vt);
        }
    }

    if(stats)
    {
        int entries = 0;

        for(map<int,vector<int> >::iterator it = local_vt.begin(); it != local_vt.end(); ++it)
        {
            entries += it->second.size();
        }

        if(max_entries < entries)
            max_entries = entries;
    }
}

void Topological_Queries::batched_VT_no_reindex_leaf(Node_T &n, Box &dom, Mesh &mesh, bool stats, int &max_entries)
{
    map<int,vector<int> > local_vt;  // local smaller structure... in the end inserted into the global map..

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<tet.vertices_num(); v++)
        {
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (dom.contains(mesh.get_vertex(tet.TV(v)),mesh.get_domain().get_max()))
                update_resulting_VT(tet.TV(v),*tet_id,local_vt);
        }
    }

    if(stats)
    {
        int entries = 0;

        for(map<int,vector<int> >::iterator it = local_vt.begin(); it != local_vt.end(); ++it)
        {
            entries += it->second.size();
        }

        if(max_entries < entries)
            max_entries = entries;
    }
}
