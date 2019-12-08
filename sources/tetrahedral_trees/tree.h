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

#ifndef TREE_H
#define	TREE_H

#include "basic_types/mesh.h"

///A super-class not instantiable representing a generic tree
template<class N, class D> class Tree
{
public:
    ///A public method that returns the mesh associated to the tree
    /*!
     * \return a Mesh& variable, representing the mesh
     */
    inline Mesh& get_mesh() { return this->mesh; }
    ///A public method that returns the root node of the tree
    /*!
     * \return a N& variable, representing the root
     */
    inline N& get_root() { return this->root; }
    ///A public method that returns the division type associated to the tree
    /*!
     * \return a D& variable, representing the division type of the tree
     */
    inline D& get_decomposition() { return this->decomposition; }
    ///A public pure virtual method, implemented by the heirs class, that builds the tree
    virtual void build_tree()=0;
    
protected:
    ///A constructor method
    Tree() {}
    ///A constructor method
    Tree(const Tree& orig)
    {
        this->decomposition = orig.decomposition;
        this->mesh = orig.mesh;
        this->root = orig.root;
    }
    ///A destructor method
    virtual ~Tree() {}

    ///A protected variable representing the mesh associated to the tree
    Mesh mesh;
    ///A protected variable representing the root node of the tree
    N root;
    ///A protected variable representing the division type of the tree
    D decomposition;

    ///A protected pure virtual method, implemented by the heirs class, that split a node, creating the sons node, following the current division type
    /*!
     * \param n a N& argument, represents the node to split
     * \param domain a Box& argument, represents the node domain
     * \param level an integer argument, representing the node level in the hierarchy
     */
    virtual void split(N& n, Box& domain, int level)=0;

private:

};

#endif	/* TREE_H */

