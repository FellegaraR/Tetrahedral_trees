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

#ifndef T_TREE_H
#define	T_TREE_H

#include "tree.h"
#include "node_t.h"

#include <iostream>
#include <algorithm>

///An inner-class, implementing Tree, that represents a tree which uses the pm criterion, without vertices
template<class D>
class T_Tree : public Tree<Node_T,D>{
public:
    ///A constructor method
    T_Tree(const T_Tree& orig) : Tree<Node_T,D>(orig) { this->tetrahedra_threshold = orig.tetrahedra_threshold; }
     ///A constructor method
    /*!
     * \param maxT an integer argument, represents the maximum number of tetrahedra admitted for a tree node
     */
    T_Tree(int maxT);
    ///A destructor method
    virtual ~T_Tree() {}
    ///A public method that returns the maximum number of tetrahedra admitted for a node
    /*!
     * \return an integer value representing the maximum number of tetrahedra admitted
     */
    inline int get_tetrahedra_threshold() { return this->tetrahedra_threshold; }
    ///A public method that sets the maximum number of tetrahedra admitted for a node
    /*!
     * \param maxT an integer argument, represents the maximum number of tetrahedra
     */
    inline void set_tetrahedra_threshold(int maxT) { this->tetrahedra_threshold = maxT; }
    ///A public method that builds the tree
    void build_tree();
private:
    ///A private variable representing the maximum number of tetrahedra admitted for a node
    int tetrahedra_threshold;
    ///A private method that adds a tetrahedron to the tree structure
    /*!
     * This method checks if a tetrahedron has a proper intersection with the node, and then insert the tetrahedron into
     * the corresponding array if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     * If the maximum number of tetrahedra for a node is exceeded, a splitting operation is activated only if
     * the number of tetrahedra exceeedes the threshold and happens not all tetrahedra are incident in a common vertex
     *
     * \param n a Node_T& argument, represents the node in which we try to insert the tetrahedron
     * \param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param t an integer argument, represents the tetrahedron to insert
     */
    void add_tetrahedron(Node_T& n, Box& domain, int level, int t);
    ///A protected method that split a node, creating the sons node, following the current division type
    /*!
     * This method also reinsert all the tetrahedron, which are in the splitted node, into the sons node
     * \param n a Node_T& argument, represents the node to split
     * \param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     */
    void split(Node_T& n, Box& domain, int level);
    ///A public method that checks if a node contains the maximum number of tetrahedra admitted
    /*!
     * \param n a Node_T& argument, represents the node to check
     * \param mesh a Mesh& argument, representing the tetrahedral mesh
     * \return a boolean value, true if the limit is exceeded, false otherwise
     */
    bool is_full(Node_T &n, Mesh &mesh);
};

template<class D> T_Tree<D>::T_Tree(int maxT)
{
    this->tetrahedra_threshold = maxT;
    this->mesh = Mesh();
    this->root = Node_T();
}

template<class D> void T_Tree<D>::build_tree()
{
    for(int i=1;i<=this->mesh.get_num_tetrahedra();i++)
    {
        this->add_tetrahedron(this->root,this->mesh.get_domain(),0,i);
    }
}

template<class D> void T_Tree<D>::add_tetrahedron(Node_T& n, Box& domain, int level, int t)
{
    if (!Geometry_Wrapper::tetra_in_box_build(t,domain,this->mesh)) return;

    if(n.is_leaf())
    {
        n.add_tetrahedron(t);
        if(is_full(n,this->get_mesh()))
            this->split(n,domain,level);
    }
    else
    {
        for(int i=0;i<this->decomposition.son_number();i++)
        {
            Box son_dom = this->decomposition.compute_domain(domain,level,i);
            int son_level = level+1;
            this->add_tetrahedron(*n.get_son(i),son_dom,son_level,t);
        }
    }
}

template<class D> void T_Tree<D>::split(Node_T& n, Box& domain, int level)
{
    n.init_sons(this->decomposition.son_number());
    for(int i=0;i<this->decomposition.son_number();i++)
    {
        Node_T* s = new Node_T();
        n.set_son(s,i);
    }

    for(RunIterator runIt = n.t_array_begin_iterator(), runEnd = n.t_array_end_iterator(); runIt != runEnd; ++runIt)
        this->add_tetrahedron(n,domain,level,*runIt);

    n.clear_t_array();
}

template<class D> bool T_Tree<D>::is_full(Node_T &n, Mesh &mesh)
{
    int t_size = n.get_t_array_size();
    if(t_size > this->tetrahedra_threshold)
    {
        //check if the tetrahedra are all incident in a common vertex
        vector<vertex_tetrahedron_pair> vert_vec;
        vert_vec.assign(t_size*4,vertex_tetrahedron_pair());
        sorting_vertices(vert_vec,n.get_t_array(),mesh);
        int count = 1;
        for(int i=0;i<(t_size*4)-1;i++)
        {
            if(vert_vec[i].v == vert_vec[i+1].v){
                count++;
                if(count > this->tetrahedra_threshold && count==t_size)
                {
                    //if we find this vertex, then we have no convinience at splitting..
                    return false;
                }
            }
            else{
                count = 1;
            }
        }
        return true;
    }
    return false;
}


#endif	/* T_TREE_H */

