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

#ifndef REINDEXER_H
#define REINDEXER_H

#include <vector>
#include <set>
#include <iostream>
#include <map>

#include "basic_types/mesh.h"
#include "p_tree.h"
#include "pt_tree.h"
#include "node_t.h"

/**
 * @brief A class that exploits the spatial coherence of the indexed enties of a Tetrahedral tree and resort accordingly these in the indexed mesh and Tetrahedral tree representation
 * This class also compresses the tree representation using the Compressed encoding.
 */
class Reindexer
{
public:
    Reindexer() { this->indices_counter=1; }
    /**
     * @brief A public method that exploit the spatial coherence of vertices and tetrahedra and compresses their representation within the tree
     *
     * @param tree a T& argument representing the tree to reindex and compress
     */
    template<class T> void reindex_tree_and_mesh(T& tree); //for RT-T and T-T trees
    /**
     * @brief A public method that exploit the spatial coherence of vertices and tetrahedra and compresses their representation within the tree
     *
     * @param tree a P_Tree& argument representing the tree to reindex and compress
     */
    template<class D> void reindex_tree_and_mesh(P_Tree<D>& tree);
    /**
     * @brief A public method that exploit the spatial coherence of vertices and tetrahedra and compresses their representation within the tree
     *
     * @param tree a PT_Tree& argument representing the tree to reindex and compress
     */
    template<class D> void reindex_tree_and_mesh(PT_Tree<D>& tree);

private:
    // FOR VERTICES
    /**
     * @brief A private method that exploits the spatial coherence of the vertices, compresses their representation in the tree and resort the corresponding array of the mesh
     * The procedure is compatible for T-Ttrees and RT-Ttrees
     *
     * @param n a Node_V& argument, represents the node
     * @param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * @param division a D& argument, representing the tree subdivision
     * @param mesh a Mesh& argument, the tetrahedral mesh
     */
    template<class D> void reindex_vertices(Node_T& n, Box& domain, int level, D& division, Mesh& mesh);
    /**
     * @brief A private method that exploits the spatial coherence of the vertices, compresses their representation in the tree and resort the corresponding array of the mesh
     * The procedure is compatible for P-Ttrees and PT-Ttrees
     *
     * @param n a Node_V& argument, represents the node
     * @param division a D& argument, representing the tree subdivision
     */
    template<class D> void reindex_vertices(Node_V& n, D& division);
    /**
     * @brief A private method that resort the vertices array of the mesh
     *
     * @param mesh a Mesh& argument, the tetrahedral mesh
     */
    void update_mesh_vertices(Mesh& mesh);

    // FOR TETRAHEDRA
    /**
     * @brief A private method that exploits the spatial coherence of the tetrahedra by extracting the tetrahedra-leaves association
     * The procedure is compatible for T-Ttrees and RT-Ttrees
     *
     * @param n a Node_V& argument, represents the node
     * @param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * @param division a D& argument, representing the tree subdivision
     * @param mesh a Mesh& argument, the tetrahedral mesh
     */
    template<class D> void extract_tetra_leaves_association(Node_T& n, Box& dom, int level, D& division, Mesh& mesh);
    /**
     * @brief A private method that exploits the spatial coherence of the tetrahedra by extracting the tetrahedra-leaves association
     * The procedure is compatible for P-Ttrees and PT-Ttrees
     *
     * @param n a Node_V& argument, represents the node
     * @param division a D& argument, representing the tree subdivision
     * @param mesh a Mesh& argument, the tetrahedral mesh
     */
    template<class D> void extract_tetra_leaves_association(Node_V& n, D& division, Mesh& mesh);
    /**
     * @brief A privater method that extracts the dual association leaves-tetrahedra
     */
    void extract_leaves_tetra_association(Mesh &mesh);

    /**
     * @brief A private method that exploits the spatial coherence of the tetrahedra and compresses their representation in the tree
     *
     * @param n a N& argument, represents the node
     * @param division a D& argument, representing the tree subdivision
     * @param mesh a Mesh& argument, the tetrahedral mesh
     */
    template<class N,class D> void reindex_tetrahedra(N& n, D& division, Mesh& mesh);
    /**
     * @brief A private method that compresses the tetrahedra array of a leaf block
     *
     * @param n a N& argument, represents the node
     * @param division a D& argument, representing the tree subdivision
     * @param mesh a Mesh& argument, the tetrahedral mesh
     */
    template<class N> void compress_t_array(N& n,int_vect &new_t_list);
    /**
     * @brief A private method that resort the tetrahedra array of the mesh
     *
     * @param mesh a Mesh& argument, the tetrahedral mesh
     */
    void update_mesh_tetrahedra(Mesh& mesh);
    /**
     * @brief A private method that resets the auxiliary variables of the class
     *
     */
    inline void reset()
    {
        this->coherent_indices.clear();
        this->indices_counter=1;
    }

    ///A private vector containing the coherent position indices of vertices/tetrahedra
    int_vect coherent_indices;
    ///A private variable representing the counter of position indices
    int indices_counter;
    ///A private nested vector of pairs representing the tetrahedra-leaves association
    vector<vector<pair<int,int> > > tetra_leaves_association;
};

template<class T> void Reindexer::reindex_tree_and_mesh(T& tree)
{ 
    coherent_indices.assign(tree.get_mesh().get_num_vertices(),-1);
    reindex_vertices(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_decomposition(),tree.get_mesh());
    update_mesh_vertices(tree.get_mesh());
    reset();

    coherent_indices.assign(tree.get_mesh().get_num_tetrahedra(),-1);
    tetra_leaves_association.assign(tree.get_mesh().get_num_tetrahedra(),vector<pair<int,int> >());
    extract_tetra_leaves_association(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_decomposition(),tree.get_mesh());
    extract_leaves_tetra_association(tree.get_mesh());

    reindex_tetrahedra(tree.get_root(),tree.get_decomposition(),tree.get_mesh());
    update_mesh_tetrahedra(tree.get_mesh());

    reset();
    return;
}

template<class D> void Reindexer::reindex_tree_and_mesh(P_Tree<D> &tree)
{
    coherent_indices.assign(tree.get_mesh().get_num_vertices(),-1);
    reindex_vertices(tree.get_root(),tree.get_decomposition());
    update_mesh_vertices(tree.get_mesh());
    reset();

    coherent_indices.assign(tree.get_mesh().get_num_tetrahedra(),-1);
    tetra_leaves_association.assign(tree.get_mesh().get_num_tetrahedra(),vector<pair<int,int> >());
    extract_tetra_leaves_association(tree.get_root(),tree.get_decomposition(),tree.get_mesh());
    extract_leaves_tetra_association(tree.get_mesh());

    reindex_tetrahedra(tree.get_root(),tree.get_decomposition(),tree.get_mesh());
    update_mesh_tetrahedra(tree.get_mesh());

    reset();
    return;
}

template<class D> void Reindexer::reindex_tree_and_mesh(PT_Tree<D> &tree)
{
    coherent_indices.assign(tree.get_mesh().get_num_vertices(),-1);
    reindex_vertices(tree.get_root(),tree.get_decomposition());
    update_mesh_vertices(tree.get_mesh());
    reset();

    coherent_indices.assign(tree.get_mesh().get_num_tetrahedra(),-1);
    tetra_leaves_association.assign(tree.get_mesh().get_num_tetrahedra(),vector<pair<int,int> >());
    extract_tetra_leaves_association(tree.get_root(),tree.get_decomposition(),tree.get_mesh());
    extract_leaves_tetra_association(tree.get_mesh());

    reindex_tetrahedra(tree.get_root(),tree.get_decomposition(),tree.get_mesh());
    update_mesh_tetrahedra(tree.get_mesh());

    reset();
    return;
}

template<class D> void Reindexer::reindex_vertices(Node_T& n, Box &domain, int level, D& division, Mesh &mesh)
{
    if (n.is_leaf())
    {
        set<int> contained_vertices;
        //we recollect the vertices geometrically contained by the leaf block
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& tet_id = itPair.first;
            Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
            for(int j=0;j<tet.vertices_num();j++)
            {
                if(domain.contains(mesh.get_vertex(tet.TV(j)),mesh.get_domain().get_max()))
                {
                    int real_v_index = tet.TV(j);
                    contained_vertices.insert(real_v_index);
                }
            }
        }
        //set the coherent position indices for those vertices
        for(set<int>::iterator it=contained_vertices.begin(); it!=contained_vertices.end();it++)
        {
            coherent_indices[*it-1] = indices_counter;
            indices_counter++;
        }
    }
    else
    {        
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                Box son_dom = division.compute_domain(domain,level,i);
                int son_level = level +1;
                reindex_vertices(*n.get_son(i),son_dom,son_level,division,mesh);
            }
        }
    }
}

template<class D> void Reindexer::reindex_vertices(Node_V &n, D& division)
{
    if (n.is_leaf())
    {
        //set the coherent position indices for the indexed vertices
        if(n.get_real_v_array_size()>0)
        {
            int start = indices_counter;
            for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
            {
                RunIterator const& v_id = itPair.first;
                coherent_indices[*v_id-1] = indices_counter;
                indices_counter++;
            }
            int end = indices_counter;
            n.clear_v_array();
            n.set_v_range(start,end);
//            cout<<n.get_v_start()<<" "<<n.get_v_end()<<endl;
//            int a; cin>>a;
        }        
    }
    else
    {        
        int start = indices_counter;
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                reindex_vertices(*n.get_son(i),division);
            }
        }
        int end = indices_counter;
        n.set_v_range(start,end);
    }
}

template<class N,class D> void Reindexer::reindex_tetrahedra(N& n, D& division, Mesh &mesh)
{
    if (n.is_leaf())
    {
//        cout<<n<<endl;
//        int_vect tmp = n.get_t_array();
//        cout<<"before: ";
//        cout<<"size: "<<tmp.size()<<endl;
//        for (unsigned i=0; i < tmp.size(); i++)
//            cout << tmp[i] << " ";
//        cout << endl;

        int_vect new_t_list;
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& tet_id = itPair.first;
            new_t_list.push_back(coherent_indices[*tet_id-1]);
        }
        n.clear_t_array();

        if(new_t_list.size()>0)
        {
            compress_t_array(n,new_t_list);
        }

//        tmp = n.get_t_array();
//        cout<<"after: ";
//        cout<<"size: "<<tmp.size()<<endl;
//        for (unsigned i=0; i < tmp.size(); i++)
//            cout << tmp[i] << " ";
//        cout << endl;
//        int a; cin>>a;
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                reindex_tetrahedra(*n.get_son(i),division,mesh);
            }
        }
    }
}

template<class N> void Reindexer::compress_t_array(N& n, int_vect &new_t_list)
{
    sort(new_t_list.begin(),new_t_list.end());

    int count=0;
    int start_t_id = new_t_list[0];

    if(new_t_list.size()==1)
    {
       n.add_tetrahedron(start_t_id);
       return;
    }
    //otherwise visit the t_list

    //now obtain the new encoding
    for(unsigned i=0; i<new_t_list.size(); i++)
    {
        if((i+1<new_t_list.size()) && (new_t_list[i]+1) == new_t_list[i+1]) //I have a consecutive range in the t_list
        {
            count++;
        }
        else //found a possible new range
        {
            if(count > 1) //more than two consecutive tetrahedra
            {
                //the range should be from start_t_id to start_t_id + count
                //example: 12 - 13 - 14 - 15
                //is encoded: -12 - 3
                //the first index is implicit
                n.add_tetrahedron(-start_t_id);
                n.add_tetrahedron(count);
            }
            else //less or equal to two
            {
                n.add_tetrahedron(start_t_id);
                if(count==1)
                    n.add_tetrahedron(start_t_id+count);
            }
            //re-init
            count = 0;
            start_t_id = new_t_list[i+1];
        }
    }
}

template<class D> void Reindexer::extract_tetra_leaves_association(Node_T &n, Box& dom, int level, D& division, Mesh &mesh)
{
    if (n.is_leaf())
    {
        pair<int,int> leaf;
        //get the range of internal vertices
        n.get_v_range(leaf.first, leaf.second,dom,mesh);

        //I only check a leaf if it contains at least a vertex
        if(leaf.first == -1 && leaf.second == -1)
            return;

        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& tet_id = itPair.first;
            Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);
            if(n.indexes_tetrahedron_vertices_dom(tet,dom,mesh))
            {
                tetra_leaves_association[*tet_id-1].push_back(leaf);
            }
        }
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                Box son_dom = division.compute_domain(dom,level,i);
                int son_level = level +1;
                extract_tetra_leaves_association(*n.get_son(i),son_dom,son_level,division,mesh);
            }
        }
    }
}

template<class D> void Reindexer::extract_tetra_leaves_association(Node_V &n, D& division, Mesh &mesh)
{
    if (n.is_leaf())
    {
        //I only check a leaf if it contains at least a vertex
        if(n.get_v_array_size() == 0)
            return;

        pair<int,int> leaf = make_pair(n.get_v_start(),n.get_v_end());
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& tet_id = itPair.first;
            Tetrahedron& tet = mesh.get_tetrahedron(*tet_id);

            if(n.indexes_tetrahedron_vertices(tet))
            {
                tetra_leaves_association[*tet_id-1].push_back(leaf);
            }
        }
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                extract_tetra_leaves_association(*n.get_son(i),division,mesh);
            }
        }
    }
}

#endif // REINDEXER_H
