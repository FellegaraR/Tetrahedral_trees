#ifndef TOPOLOGICAL_QUERIES_BATCHED
#define TOPOLOGICAL_QUERIES_BATCHED

#include "topological_queries.h"

template<class N, class D> void Topological_Queries::batched_VT(N &n, Box &dom, Mesh &mesh, D &division, bool reindex)
{
//    cout<<"batched_VT"<<endl;

    Timer time;
    int max_entities = 0;

    time.start();
    if(reindex)
        this->batched_VT_visit(n,dom,0,mesh,division,false,max_entities);
    else
        this->batched_VT_no_reindex(n,dom,0,mesh,division,false,max_entities);
    time.stop();
    time.print_elapsed_time("[TIME] extracting bactched VT: ");

    if(reindex)
        this->batched_VT_visit(n,dom,0,mesh,division,true,max_entities);
    else
        this->batched_VT_no_reindex(n,dom,0,mesh,division,true,max_entities);
    cerr<<"[STATS] maximum number of entities: "<<max_entities<<endl;
}

template<class N, class D> void Topological_Queries::batched_VT_visit(N &n, Box &dom, int level, Mesh &mesh, D &division, bool stats, int &max_entries)
{
//    cout<<"ci arrivo"<<endl;
//    cout<<n<<" "<<dom<<endl;
    if (n.is_leaf())
    {        
//        cout<<"batched_VT_leaf"<<endl;
        this->batched_VT_leaf(n,dom,mesh,stats,max_entries);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            N& son = *n.get_son(i);
            this->batched_VT_visit(son, son_dom, son_level, mesh, division, stats, max_entries);
        }
    }
}

template<class N, class D> void Topological_Queries::batched_VT_no_reindex(N &n, Box &dom, int level, Mesh &mesh, D &division, bool stats, int &max_entries)
{
//    cout<<"batched_VT_no_reindex"<<endl;
    if (n.is_leaf())
    {
        this->batched_VT_no_reindex_leaf(n,dom,mesh,stats,max_entries);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->batched_VT_no_reindex(*n.get_son(i), son_dom, son_level, mesh, division, stats, max_entries);
        }
    }
}

template<class N, class D> void Topological_Queries::batched_TT(N &n, Mesh &mesh, D &division)
{
    int max_entities = 0;

    vector<vector<int> > tt;
    vector<int> tmp;
    tmp.assign(4,-1);
    tt.assign(mesh.get_num_tetrahedra(),tmp);

    Timer time;
    time.start();
    this->batched_TT_visit(n,mesh,division,tt,false,max_entities);
    time.stop();
    time.print_elapsed_time("[TIME] extracting bactched TT: ");

    tt.clear();
    tmp.assign(4,-1);
    tt.assign(mesh.get_num_tetrahedra(),tmp);

    this->batched_TT_visit(n,mesh,division,tt,true,max_entities);
    cerr<<"[STATS] maximum number of faces: "<<max_entities<<endl;
}

template<class N, class D> void Topological_Queries::batched_TT_visit(N &n, Mesh &mesh, D &division, vector<vector<int> > &tt, bool stats, int &max_entries)
{
    if (n.is_leaf())
    {
        this->batched_TT_leaf(n,mesh,tt,stats,max_entries);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            this->batched_TT_visit(*n.get_son(i), mesh, division, tt, stats, max_entries);
        }
    }
}

template<class N> void Topological_Queries::batched_TT_leaf(N &n, Mesh &mesh, vector<vector<int> > &tt, bool stats, int &max_entries)
{
    vector<triangle_tetrahedron_tuple> faces;
    triangle_tetrahedron_tuple face;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);

        for(int v=0; v<tet.vertices_num(); v++)
        {
            if(tt[*tet_id-1][v]==-1) // the entry is not initialized
            {
                tet.face_tuple(v,face,*tet_id);
                faces.push_back(face);
            }
        }
    }

    // (*) order the faces array
    sorting_faces(faces);
    // (*) set the adjacencies on these faces
    unsigned j=0;
    while(j<faces.size())
    {
        if(j+1<faces.size())
        {
            if(faces[j] == faces[j+1])
            {
                tt[faces[j].t -1][faces[j].f_pos] = faces[j+1].t;
                tt[faces[j+1].t -1][faces[j+1].f_pos] = faces[j].t;
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

    if(stats)
    {
        if(max_entries < (int)faces.size())
            max_entries = (int)faces.size();
    }
}

#endif // TOPOLOGICAL_QUERIES_BATCHED

