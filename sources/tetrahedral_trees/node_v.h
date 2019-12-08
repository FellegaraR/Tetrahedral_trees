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

#ifndef NODE_V_H
#define NODE_V_H

#include "tetrahedral_trees/node.h"

/**
 * @brief The Node_V class extends the Node class, encoding also a list of vertices.
 * This class is used by the P-Ttrees and PT-Ttrees
 */
class Node_V : public Node<Node_V>
{
public:
    ///A constructor
    Node_V() : Node<Node_V>() { }
    ///A copy-constructor
    Node_V(const Node_V& orig) : Node<Node_V>(orig) { this->vertices = orig.vertices; }
    ///A destructor
    virtual ~Node_V() {}
    ///A public method that add a vertex index to the corresponding node list
    /*!
     * \param ind an integer argument, representing the vertex index
     */
    inline void add_vertex(int ind) { this->vertices.push_back(ind); }
    ///A public method that free space occupied by the two lists
    inline void clear_v_array() { this->vertices.clear(); }
    /**
     * @brief A public method that set the vertices range after exploiting the vertices spatial coherence
     *
     * @param start an integer with the first vertex position index of the node
     * @param end an integer with the first vertex position index outside the node
     */
    inline void set_v_range(int start, int end) { vertices.push_back(-start); vertices.push_back(end-start-1); }
    ///
    inline int get_v_start() const { return abs(vertices[0]); }
    ///
    inline int get_v_end() const { return (abs(vertices[0])+vertices[1])+1; }
    //NOTA: to use only after the index is built/loaded from file
    /**
     * @brief A public method that checks if a vertex is indexed by the node
     *
     * @param v_id an integer representing the position index of the vertex
     * @return true if v_id is indexed, false otherwise
     */
    inline bool indexes_vertex(int v_id) { return (v_id >= get_v_start() && v_id < get_v_end()); }
    /**
     * @brief A public method that checks if a tetrahedron is indexed by the node
     * The method checks if at least one of the vertices of the tetrahedron is indexed by the node
     * NOTA: only range comparisons are executed
     *
     * @param t a Tetrahedron& representing the tetrahedron
     * @return true if at least one vertex is indexed, false otherwise
     */
    inline bool indexes_tetrahedron_vertices(Tetrahedron &t)
    {
        if(this->get_v_array_size()==0) // no vertices
            return false;
        for(int i=0; i<t.vertices_num(); i++)
            if(this->indexes_vertex(t.TV(i)))
                return true;
        return false;
    }
    /**
     * @brief A public method that checks if a tetrahedron is completely indexed by a node
     * The method checks if all the vertices of the tetrahedron are indexed by the node
     * NOTA: only range comparisons are executed
     *
     * @param t a Tetrahedron& representing the tetrahedron
     * @return true if all the vertices are indexed, false otherwise
     */
    inline bool completely_indexes_tetrahedron_vertices(Tetrahedron &t)
    {
        if(this->get_v_array_size()==0) // no vertices
            return false;
        for(int i=0; i<t.vertices_num(); i++)
            if(!this->indexes_vertex(t.TV(i)))
                return false;
        return true;
    }
    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */
    friend std::ostream& operator<<(std::ostream& out, const Node_V& p)
    {
        if(p.is_leaf())
        {
            if(p.get_v_array_size() == 2) // it should be reindexed
                out <<"Leaf["<< p.get_v_start() << " " << p.get_v_end()<<"]";
            else
            {
                out << "Leaf[s->"<< p.get_real_v_array_size() << "]";
            }
        }
        else
            out <<"Node["<< p.get_v_start() << " " << p.get_v_end()<<"]";
        return out;
    }
    /**
     * @brief A public method that returns a pair of iterators to the vertices array
     *
     * @return RunIteratorPair
     */
    inline RunIteratorPair make_v_array_iterator_pair() { return run_iterator<int>::make_run_iterator_pair(vertices); }
    /**
     * @brief A public method that returns the iterator to the begin of the vertices array
     *
     * @return RunIterator
     */
    inline RunIterator v_array_begin_iterator() { return run_iterator<int>(vertices.begin(),vertices.end()); }
    /**
     * @brief A public method that returns the iterator to the end of the vertices array
     *
     * @return RunIterator
     */
    inline RunIterator v_array_end_iterator() { return run_iterator<int>(vertices.end()); }
    /**
     * @brief A public method that returns the number of indexed vertices
     *
     * @return int
     */
    inline int get_real_v_array_size() const { return run_iterator<int>(vertices.begin(),vertices.end()).elementCountFast(vertices); }
    /**
     * @brief A public method that returns the size of the vertices array
     *
     * @return int
     */
    inline int get_v_array_size() const { return this->vertices.size(); }

protected:    
   int_vect vertices;
};

#endif // NODE_V_H
