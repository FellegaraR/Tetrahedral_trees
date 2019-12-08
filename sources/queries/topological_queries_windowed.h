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
#include "border_checker.h"
#include "utilities/timer.h"
#include "io/reader.h"
#include "utilities/sorting.h"
#include "io/writer.h"

template<class N, class D> void Topological_Queries::windowed_VT(N &n, Box &dom, Mesh &mesh, D &division, string query_path, bool reindexed)
{
    vector<Box> boxes;
    Reader::read_queries(boxes,query_path);

    map<int,vector<int> > results;

    Timer time;
    double tot_time = 0;


    for(unsigned j=0;j<boxes.size();j++)
    {
        time.start();
        if(reindexed)
            windowed_VT(n,dom,0,boxes[j],mesh,division,results);
        else
            windowed_VT_no_reindex(n,dom,0,boxes[j],mesh,division,results);
        time.stop();
        tot_time += time.get_elapsed_time();

        //debug print
        cout<<"for box "<<j<<" vertices found: "<<results.size()<<endl;
        results.clear();
    }    
    cerr<<"extracting windowed VT "<<tot_time<<endl;
}

template<class D> void Topological_Queries::windowed_VT(Node_T &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &vt)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        this->windowed_VT_Leaf(n,dom,b,mesh,vt);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_VT(*n.get_son(i), son_dom, son_level, b, mesh, division, vt);
        }
    }
}

template<class N, class D> void Topological_Queries::windowed_VT_no_reindex(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &vt)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        this->windowed_VT_Leaf_no_reindex(n,dom,b,mesh,vt);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_VT_no_reindex(*n.get_son(i), son_dom, son_level, b, mesh, division, vt);
        }
    }
}

template<class D> void Topological_Queries::windowed_VT(Node_V &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &vt)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        // here we do not need the son dom to understand if a vertex is inside the leaf ==>> we use the vertices range
        this->windowed_VT_Leaf(n,b,mesh,vt);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_VT(*n.get_son(i), son_dom, son_level, b, mesh, division, vt);
        }
    }
}

template<class N> void Topological_Queries::windowed_VT_Leaf_no_reindex(N& n, Box &dom, Box &b, Mesh& mesh, map<int,vector<int> > &vt)
{
    map<int,vector<int> > local_vt;  // local smaller structure... in the end inserted into the global map..

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<tet.vertices_num(); v++)
        {
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (dom.contains(mesh.get_vertex(tet.TV(v)),mesh.get_domain().get_max()) &&
                    b.contains_with_all_closed_faces(mesh.get_vertex(tet.TV(v))))
                update_resulting_VT(tet.TV(v),*tet_id,local_vt);
        }
    }

    vt.insert(local_vt.begin(),local_vt.end());
}

template<class N, class D> void Topological_Queries::windowed_Distortion(N &n, Box &dom, Mesh &mesh, D &division, string query_path, bool reindexed)
{
    vector<Box> boxes;
    Reader::read_queries(boxes,query_path);

    Border_Checker checker = Border_Checker();
    Timer time;
    time.start();
    checker.calc_mesh_borders(n,dom,0,mesh,division);
    time.stop();
    time.print_elapsed_time("updating borders ");

    map<int,double> results;
    double tot_time = 0;

    for(unsigned j=0;j<boxes.size();j++)
    {
        time.start();
        if(reindexed)
            windowed_Distortion(n,dom,0,boxes[j],mesh,division,results);
        else
            windowed_Distortion_no_reindex(n,dom,0,boxes[j],mesh,division,results);
        time.stop();
        tot_time += time.get_elapsed_time();

        //debug print
        cout<<"for box "<<j<<" vertices found: "<<results.size()<<endl;
        results.clear();
    }
    cerr<<"extracting windowed distortion "<<tot_time<<endl;
}

template<class D> void Topological_Queries::windowed_Distortion(Node_T &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,double> &dist)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        this->windowed_Distortion_Leaf(n,dom,b,mesh,dist);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_Distortion(*n.get_son(i), son_dom, son_level, b, mesh, division, dist);
        }
    }
}

template<class N, class D> void Topological_Queries::windowed_Distortion_no_reindex(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,double> &dist)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        this->windowed_Distortion_Leaf_no_reindex(n,dom,b,mesh,dist);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_Distortion_no_reindex(*n.get_son(i), son_dom, son_level, b, mesh, division, dist);
        }
    }
}

template<class D> void Topological_Queries::windowed_Distortion(Node_V &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,double> &dist)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        this->windowed_Distortion_Leaf(n,b,mesh,dist);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_Distortion(*n.get_son(i), son_dom, son_level, b, mesh, division, dist);
        }
    }
}

template<class N> void Topological_Queries::windowed_Distortion_Leaf_no_reindex(N& n, Box &dom, Box &b, Mesh& mesh, map<int,double> &dist)
{
    map<int,vector<int> > vt;
    set<int> border_v;
    map<int,double> local_dist; // local smaller structure... in the end inserted into the global map..

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
        for(int v=0; v<tet.vertices_num(); v++)
        {
            int real_v_index = tet.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (dom.contains(mesh.get_vertex(real_v_index),mesh.get_domain().get_max()) &&
                    b.contains_with_all_closed_faces(mesh.get_vertex(real_v_index)))
            {
                update_resulting_VT(real_v_index,*tet_id,vt);

                double partial_distortion = Geometry_Distortion::get_trihedral_angle(tet,real_v_index,mesh);
                update_resulting_distortion(real_v_index,partial_distortion,local_dist);

                set<int>::iterator iter = border_v.find(real_v_index);
                if(iter == border_v.end())
                {
                    //controllo se almeno uno degli altri tre è negativo e nel caso metto il vertice corrente come di bordo
                    for(int j=1; j<tet.vertices_num(); j++)
                    {
                        if(tet.is_border_face((j+v)%tet.vertices_num()))
                        {
                            border_v.insert(real_v_index);
                            break;
                        }
                    }
                }
            }
        }
    }

    for(map<int,vector<int> >::iterator iter=vt.begin(); iter!=vt.end(); ++iter)
    {
        int real_v_index = iter->first;

        map<int,double>::iterator it_v = local_dist.find(real_v_index);
        set<int>::iterator it = border_v.find(real_v_index);
        if(it == border_v.end())
        {
            it_v->second = - it_v->second;
            for(int_vect_iter t = iter->second.begin(); t != iter->second.end(); ++t)
            {
                Tetrahedron& tet = mesh.get_tetrahedron(*t);
                double partial_distortion = Geometry_Distortion::get_trihedral_angle_3D(tet,real_v_index,mesh);
                it_v->second += partial_distortion;
            }
        }
        else
        {
            it_v->second = ( 4*PI - it_v->second );
        }
    }

    dist.insert(local_dist.begin(),local_dist.end());
}

template<class N, class D> void Topological_Queries::windowed_TT(N &n, Box &dom, Mesh &mesh, D &division, string query_path)
{
    vector<Box> boxes;
    Reader::read_queries(boxes,query_path);

    map<int,vector<int> > results;
    Timer time;
    double tot_time = 0.0;

    bm::bvector<> checkTetra;

    for(unsigned j=0;j<boxes.size();j++)
    {
        time.start();
        windowed_TT(n,dom,0,boxes[j],mesh,division,results,checkTetra);
        time.stop();
        tot_time += time.get_elapsed_time();

        //debug print
        cout<<"for box "<<j<<" tetrahedra found: "<<results.size()<<endl;
        results.clear();

        checkTetra.reset();
    }
    cerr<<"extracting windowed TT "<<tot_time<<endl;
}

template<class N, class D> void Topological_Queries::windowed_TT(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &tt, bm::bvector<> &checkTetra)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        if(b.completely_contains(dom))
            this->windowed_TT_Leaf_add(n,mesh,tt,checkTetra);
        else
            this->windowed_TT_Leaf_test(n,b,mesh,tt,checkTetra);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_TT(*n.get_son(i), son_dom, son_level, b, mesh, division, tt, checkTetra);
        }
    }
}

template<class N> void Topological_Queries::windowed_TT_Leaf_test(N& n, Box &b, Mesh& mesh, map<int,vector<int> > &tt, bm::bvector<> &checkTetra)
{
    vector<triangle_tetrahedron_tuple> faces;
    Box bb;
    pair<int,int> run;

    for(vector<int>::iterator it=n.get_t_array_begin(); it!=n.get_t_array_end(); ++it)
    {
        if(n.get_run_bounding_box(it,bb,mesh,run))
        {
            if(b.completely_contains(bb))
            {
                for(int t_id=run.first; t_id<=run.second; t_id++)
                {
                    map<int,vector<int> >::const_iterator entry = tt.find(t_id);

                    checkTetra[t_id]=true;

                    //if the run is completely contained.. simply add..
                    if(entry == tt.end()) // first time for the current tetrahedron
                        init_TT_entry(t_id,tt);
                    add_faces(t_id,faces,mesh,entry,tt);
                }
            }
            else if(b.intersects(bb))
            {
                for(int t_id=run.first; t_id<=run.second; t_id++)
                {
                    map<int,vector<int> >::const_iterator entry = tt.find(t_id);

                    //if I have an entry into the result or I have an intersection with the box
                    if(entry != tt.end() || (!checkTetra[t_id] && Geometry_Wrapper::tetra_in_box(t_id,b,mesh)))
                    {

                        if(entry == tt.end()) // first time for the current tetrahedron
                            init_TT_entry(t_id,tt);
                        add_faces(t_id,faces,mesh,entry,tt);
                    }

                    checkTetra[t_id]=true;
                }
            }
        }
        else
        {
            map<int,vector<int> >::const_iterator entry = tt.find(*it);

            //if I have an entry into the result or I have an intersection with the box
            if(entry != tt.end() || (!checkTetra[*it] && Geometry_Wrapper::tetra_in_box(*it,b,mesh)))
            {

                if(entry == tt.end()) // first time for the current tetrahedron
                    init_TT_entry(*it,tt);
                add_faces(*it,faces,mesh,entry,tt);
            }

            checkTetra[*it]=true;
        }
    }
    finalize_TT_Leaf(faces,tt,mesh);
}

template<class N> void Topological_Queries::windowed_TT_Leaf_add(N& n, Mesh& mesh, map<int,vector<int> > &tt, bm::bvector<> &checkTetra)
{
    vector<triangle_tetrahedron_tuple> faces;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;

        map<int,vector<int> >::const_iterator entry = tt.find(*tet_id);

        checkTetra[*tet_id]=true;

        //if I have an entry into the result or I have an intersection with the box
        if(entry == tt.end()) // first time for the current tetrahedron
            init_TT_entry(*tet_id,tt);
        add_faces(*tet_id,faces,mesh,entry,tt);
    }
    finalize_TT_Leaf(faces,tt,mesh);
}

template<class N, class D> void Topological_Queries::linearized_TT(N &n, Box &dom, Mesh &mesh, D &division, string query_path)
{
    vector<Box> boxes;
    Reader::read_queries(boxes,query_path);

    map<int,vector<int> > results;
    Timer time;
    double tot_time = 0.0;

    bm::bvector<> checkTetra;

    for(unsigned j=0;j<boxes.size();j++)
    {
        time.start();
        linearized_TT(n,dom,0,boxes[j],mesh,division,results,checkTetra);
        tot_time += time.get_elapsed_time();

        //debug print
        cout<<"for box "<<j<<" tetrahedra found: "<<results.size()<<endl;
        results.clear();
        checkTetra.reset();
    }
    cerr<<"extracting linearized TT "<<tot_time<<endl;
}

template<class N, class D> void Topological_Queries::linearized_TT(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, map<int,vector<int> > &tt, bm::bvector<> &checkTetra)
{
    if(!Geometry_Wrapper::line_in_box(b.get_min(),b.get_max(),dom))
        return;

    if (n.is_leaf())
    {
        this->linearized_TT_Leaf(n,b,mesh,tt,checkTetra);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->linearized_TT(*n.get_son(i), son_dom, son_level, b, mesh, division, tt,checkTetra);
        }
    }
}

template<class N> void Topological_Queries::linearized_TT_Leaf(N& n, Box &b, Mesh& mesh, map<int,vector<int> > &tt, bm::bvector<> &checkTetra)
{
    vector<triangle_tetrahedron_tuple> faces;

    Box bb;
    pair<int,int> run;

    for(vector<int>::iterator it=n.get_t_array_begin(); it!=n.get_t_array_end(); ++it)
    {
        if(n.get_run_bounding_box(it,bb,mesh,run))
        {
            if(Geometry_Wrapper::line_in_bounding_box(b.get_min(),b.get_max(),bb))
            {
                for(int t_id=run.first; t_id<=run.second; t_id++)
                {
                    map<int,vector<int> >::const_iterator entry = tt.find(t_id);

                    //if I have an entry into the result or I have an intersection with the box
                    if(entry != tt.end() || (!checkTetra[t_id] && Geometry_Wrapper::line_in_tetra(b.get_min(),b.get_max(),t_id,mesh)))
                    {
                        if(entry == tt.end()) // first time for the current tetrahedron
                            init_TT_entry(t_id,tt);
                        add_faces(t_id,faces,mesh,entry,tt);
                    }

                    checkTetra[t_id] = true;
                }
            }
        }
        else
        {
            map<int,vector<int> >::const_iterator entry = tt.find(*it);

            //if I have an entry into the result or I have an intersection with the box
            if(entry != tt.end() || (!checkTetra[*it] && Geometry_Wrapper::line_in_tetra(b.get_min(),b.get_max(),*it,mesh)))
            {
                if(entry == tt.end()) // first time for the current tetrahedron
                    init_TT_entry(*it,tt);
                add_faces(*it,faces,mesh,entry,tt);
            }

            checkTetra[*it] = true;
        }
    }
    finalize_TT_Leaf(faces,tt,mesh);
}
