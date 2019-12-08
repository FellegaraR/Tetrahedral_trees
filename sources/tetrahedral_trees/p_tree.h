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

#ifndef P_TREE_H
#define	P_TREE_H

#include "tree.h"
#include "node_v.h"
#include <geometry/geometry_wrapper.h>

///An inner-class, implementing Tree, that represents a tree which uses the Point-based criterion
template<class D> class P_Tree : public Tree<Node_V,D>
{
public:
    ///A constructor method
    /*!
     * \param vertices_per_leaf an integer argument, represents the maximum number of vertices admitted for a tree node
     */
    P_Tree(int vertices_per_leaf);
    ///A constructor method
    P_Tree(const P_Tree& orig) : Tree<Node_V,D>(orig) { this->vertices_threshold = orig.vertices_threshold; }
    ///A destructor method
    virtual ~P_Tree() {}
    ///A public method that sets the maximum number of vertices admitted for a node
    /*!
     * \param maxV an integer argument, represents the maximum number of vertices
     */
    inline void set_vertices_threshold(int maxV) { this->vertices_threshold = maxV; }
    ///A public method that returns the maximum number of vertices admitted for a node
    /*!
     * \return an integer value representing the maximum number of vertices admitted
     */
    inline int get_vertices_threshold() { return this->vertices_threshold; }
    ///A public method that builds the tree
    /*!
     * This method before inserts all the mesh vertices and then all the tetrahedra.
     * The insertion of tetrahedra does not change the hierarchy.
     */
    void build_tree();

private:
    ///A private variable representing the maximum number of vertices admitted for a node
    int vertices_threshold;
    ///A private method that adds a tetrahedron to the tree structure
    /*!
     * This method checks if a vertex is contained by a node, and then insert the vertex into
     * the node list if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     * If the maximum number of vertices for a node is exceeded, a splitting operation is activated.
     *
     * \param n a Node_V& argument, represents the node in which we try to insert the vertex
     * \param domain a Box& argument, represents the node domain
     * \param level an integer argument representing the level of n in the hierarchy
     * \param v an integer argument, represents the vertex to insert
     */
    void add_vertex(Node_V& n, Box& domain, int level, int v);
    ///A private method that adds a tetrahedron to the tree structure
    /*!
     * This method checks if a tetrahedron has a proper intersection with the node, and then insert the tetrahedron into
     * the node list if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     *
     * \param n a Node_V& argument, represents the node in which we try to insert the tetrahedron
     * \param domain a Box& argument, represents the node domain
     * \param level an integer argument representing the level of n in the hierarchy
     * \param t an integer argument, represents the tetrahedron to insert
     */
    void add_tetrahedron(Node_V& n, Box& domain, int level, int t);
    ///A protected method that split a node, creating the sons node, following the current division type
    /*!
     * This method also reinsert all the vertices and the tetrahedron, which are in the splitted node, into the sons node
     * \param n a Node_V& argument, represents the node to split
     * \param domain a Box& argument, represents the node domain
     * \param level an integer argument representing the level of n in the hierarchy
     */
    void split(Node_V& n, Box& domain, int level);
    ///A public method that checks if a node has raised the maximum block capacity
    /*!
     * \param n a Node_V& argument, represents the node to check
     * \return a boolean value, true if the limit is exceeded, false otherwise
     */
    inline bool is_full(Node_V &n) { return (n.get_v_array_size() > this->vertices_threshold); }
};

template<class D> P_Tree<D>::P_Tree(int vertices_per_leaf)
{
    this->vertices_threshold = vertices_per_leaf;
    this->mesh = Mesh();
    this->root = Node_V();
}

template<class D> void P_Tree<D>::build_tree()
{
//    cout<<"Tree_Domain: "<<this->mesh.get_domain()<<endl;
    for(int i=1;i<=this->mesh.get_num_vertices(); i++)
    {
//        cout<<"P: "<<this->mesh.get_vertex(i)<<endl;
        this->add_vertex(this->root,this->mesh.get_domain(),0,i);
//        int a; cin>>a;
    }    
    for(int i=1;i<=this->mesh.get_num_tetrahedra();i++)
    {
//        cout<<"T: "<<this->mesh.get_tetrahedron(i)<<endl;
        this->add_tetrahedron(this->root,this->mesh.get_domain(),0,i);
//        int a; cin>>a;
    }
}

template<class D> void P_Tree<D>::add_vertex(Node_V& n, Box& domain, int level, int v)
{
    if (n.is_leaf())
    {
        n.add_vertex(v);
        if(is_full(n))
            this->split(n,domain,level);
    }
    else
    {        
        for (int i = 0; i < this->decomposition.son_number(); i++)
        {
            Box son_dom = this->decomposition.compute_domain(domain,level,i);
            int son_level = level +1;            
            if (son_dom.contains(this->mesh.get_vertex(v),this->mesh.get_domain().get_max()))
            {
                this->add_vertex(*n.get_son(i),son_dom,son_level,v);
                break;
            }
        }
    }
}

template<class D> void P_Tree<D>::add_tetrahedron(Node_V& n, Box& domain, int level, int t)
{
    if (!Geometry_Wrapper::tetra_in_box_build(t,domain,this->mesh)) return;

    if(n.is_leaf())
    {
        n.add_tetrahedron(t);
    }
    else
    {
        for(int i=0;i<this->decomposition.son_number();i++)
        {
            Box son_dom = this->decomposition.compute_domain(domain,level,i);
            int son_level = level +1;
            this->add_tetrahedron(*n.get_son(i),son_dom,son_level,t);
        }
    }
}

template<class D> void P_Tree<D>::split(Node_V& n, Box& domain, int level)
{
    n.init_sons(this->decomposition.son_number());
    //initilize the son nodes
    for(int i=0;i<this->decomposition.son_number();i++)
    {
        Node_V* s = new Node_V();
        n.set_son(s,i);
    }

    //re-insert the vertices
    for(RunIterator runIt = n.v_array_begin_iterator(), runEnd = n.v_array_end_iterator(); runIt != runEnd; ++runIt)
        this->add_vertex(n,domain,level,*runIt);

    //re-insert the tetrahedra
    for(RunIterator runIt = n.t_array_begin_iterator(), runEnd = n.t_array_end_iterator(); runIt != runEnd; ++runIt)
        this->add_tetrahedron(n,domain,level,*runIt);

    //delete the arrays of the node
    n.clear_v_array();
    n.clear_t_array();
}

#endif	/* P_TREE_H */

