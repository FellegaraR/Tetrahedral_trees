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

#ifndef RT_TREE_H
#define	RT_TREE_H

#include "tree.h"
#include "node_t.h"
#include <geometry/geometry_wrapper.h>

///An inner-class, implementing Tree, that represents a tree which uses the RT-T criterion
template<class D>
class RT_Tree : public Tree<Node_T,D>{
public:
    ///A constructor method
    /*!
     * \param tetrahedra_per_leaf an integer argument, represents the maximum number of tetrahedra admitted for a tree node
     */
    RT_Tree(int tetrahedra_per_leaf);
    ///A constructor method
    RT_Tree(const RT_Tree& orig) : Tree<Node_T,D>(orig) { this->tetrahedra_threshold = orig.tetrahedra_threshold; }
    ///A destructor method
    virtual ~RT_Tree() {}
    ///A public method that sets the maximum number of tetrahedra admitted for a node
    /*!
     * \param maxT an integer argument, represents the maximum number of tetrahedra
     */
    inline void set_tetrahedra_threshold(int maxT) { this->tetrahedra_threshold = maxT; }
    ///A public method that returns the maximum number of tetrahedra admitted for a node
    /*!
     * \return an integer value representing the maximum number of tetrahedra admitted
     */
    inline int get_tetrahedra_threshold() { return this->tetrahedra_threshold; }
    ///A public method that builds the tree
    /**
     * Only the tetrahedra are inserted into the hierarchy.
     */
    void build_tree();

private:
    ///A private variable representing the maximum number of tetrahedra admitted for a node
    int tetrahedra_threshold;
    ///A private method that adds a tetrahedron to the tree structure
    /*!
     * This method checks if a tetrahedron has a proper intersection with the node, and then insert the tetrahedron into
     * the node list if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     * If the maximum number of tetrahedra for a node is exceeded, a splitting operation is activated only once.
     *
     * \param n a Node_T& argument, represents the node in which we try to insert the tetrahedron
     * \param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param t an interger argument, represents the tetrahedron to insert
     */
    void add_tetrahedron(Node_T& n, Box& domain, int level, int t);
    ///A private method that reinsert a tetrahedron to the tree structure only once
    /*!
     * This method checks if a tetrahedron has a proper intersection with the node, and then insert the tetrahedron into
     * the node list if the node is a leaf. The block threshold is not checked in this method, thus, the every time
     * a split operation is called only on the splitted parent node.
     *
     * \param n a Node_T& argument, represents the node in which we try to insert the tetrahedron
     * \param domain a Box& argument, represents the node domain
     * \param t an integer argument, represents the tetrahedron to insert
     */
    void reinsert_tetrahedron_once(Node_T& n, Box& domain, int t);
    ///A protected method that split a node, creating the sons node, following the current division type
    /*!
     * This method also reinsert all the tetrahedron, which are in the splitted node, into the sons node.
     * But unlike the other criteria it calls the addTetrahedron_once method, which doesn't cause a recursive insertion into the sons node.
     *
     * \param n a Node_T& argument, represents the node to split
     * \param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     */
    void split(Node_T& n, Box& domain, int level);

    ///A public method that checks if a node contains the maximum number of tetrahedra admitted
    /*!
     * \param n a Node_T& argument, represents the node to check
     * \return a boolean value, true if the limit is exceeded, false otherwise
     */
    inline bool is_full(Node_T &n) { return (n.get_t_array_size() > this->tetrahedra_threshold); }
};

template<class D> RT_Tree<D>::RT_Tree(int tetrahedra_per_leaf)
{
    this->tetrahedra_threshold = tetrahedra_per_leaf;
    this->mesh = Mesh();
    this->root = Node_T();
}

template<class D> void RT_Tree<D>::build_tree()
{
    for (int i = 1; i <= this->mesh.get_num_tetrahedra(); i++)
    {
        this->add_tetrahedron(this->root, this->mesh.get_domain(), 0, i);
    }
}

template<class D> void RT_Tree<D>::add_tetrahedron(Node_T& n, Box& domain, int level, int t)
{
    if (!Geometry_Wrapper::tetra_in_box_build(t,domain,this->mesh)) return;

    if (n.is_leaf())
    {
        n.add_tetrahedron(t);
        if (is_full(n))
        {
            this->split(n, domain, level);
        }
    }
    else
    {
        for (int i = 0; i<this->decomposition.son_number(); i++)
        {
            Box son_dom = this->decomposition.compute_domain(domain, level, i);
            int son_level = level +1;
            this->add_tetrahedron(*n.get_son(i), son_dom, son_level, t);
        }
    }
}

template<class D> void RT_Tree<D>::reinsert_tetrahedron_once(Node_T& n, Box& domain, int t)
{
    if (!Geometry_Wrapper::tetra_in_box_build(t,domain,this->mesh)) return;

    if (n.is_leaf())
    {
        n.add_tetrahedron(t);
    }
}

template<class D> void RT_Tree<D>::split(Node_T& n, Box& domain, int level)
{
    n.init_sons(this->decomposition.son_number());

    for (int i = 0; i<this->decomposition.son_number(); i++)
    {
        Node_T* s = new Node_T();
        n.set_son(s,i);
    }

    // we reinsert the tetrahedra in the son nodes only once
    // thus, without splitting any further the space
    for (int j = 0; j<this->decomposition.son_number(); j++)
    {
        Box son_dom = this->decomposition.compute_domain(domain,level, j);
        for(RunIterator runIt = n.t_array_begin_iterator(), runEnd = n.t_array_end_iterator(); runIt != runEnd; ++runIt)
        {
            this->reinsert_tetrahedron_once(*n.get_son(j) , son_dom, *runIt);
        }
    }

    n.clear_t_array();
}



#endif	/* RT_TREE_H */

