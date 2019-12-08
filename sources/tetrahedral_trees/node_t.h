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

#ifndef NODE_T_H
#define NODE_T_H

#include "tetrahedral_trees/node.h"

/**
 * @brief The Node_T class extends the Node class, providing an instantiable definition of it.
 * This class is used by the T-Ttrees and RT-Ttrees
 */
class Node_T : public Node<Node_T>
{
public:
    ///A constructor
    Node_T() : Node<Node_T>() { }
    ///A copy-constructor
    Node_T(const Node_T& orig) : Node<Node_T>(orig) { }
    ///A destructor
    virtual ~Node_T() {}
    friend std::ostream& operator<<(std::ostream& out, const Node_T& p)
    {
        if(p.is_leaf())
            out <<"Leaf["<< p.get_real_t_array_size() <<"]";
        else
            out <<"Node";
        return out;
    }

    /**
     * @brief A public method that returns the range of vertices completely contained into the leaf
     * NOTA: this method only works after the reordering of vertices array
     *
     * @param v_start an integer that it will contains the first vertex indexed by the leaf
     * @param v_end an integer that it will contains the first vertex outside the leaf
     * @param dom a Box& representing the node domain
     * @param mesh a Mesh& representing the tetrahedral mesh
     */
    void get_v_range(int &v_start, int &v_end, Box& dom, Mesh& mesh);
    /**
     * @brief A public method that checks if a vertex is indexed by the current node
     * NOTA: this method is a wrapper thought for T-Ttrees and RT-Ttrees that do not explicitly encode the vertices
     * within each node. For P-Ttrees and PT-Ttrees a wrapper function has been defined
     *
     * @param v_start an integer that it will contains the first vertex indexed by the leaf
     * @param v_end an integer that it will contains the first vertex outside the leaf
     * @param v_id an integer representing the position index of a vertex
     * @return true if v_id is into the leaf, false otherwise
     */
    inline bool indexes_vertex(int v_start, int v_end, int v_id) { return (v_id >= v_start && v_id < v_end); }

protected:
};

#endif // NODE_T_H
