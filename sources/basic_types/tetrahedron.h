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

#ifndef _TETRAHEDRON_H
#define	_TETRAHEDRON_H

using namespace std;

#include "utilities/sorting_structure.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

///A class representing a tetrahedron of the mesh
class Tetrahedron {
public:
    ///A constructor method
    Tetrahedron() {}
    /**
     * @brief A copy-constructor method
     * @param orig
     */
    Tetrahedron(const Tetrahedron& orig)
    {
        for (int i = 0; i < orig.vertices_num(); i++)
        {
            this->vertices[i] = orig.vertices[i];
        }
    }
    ///A constructor method
    /*!
     * \param v1 a int argument, represents the first tetrahedron vertex
     * \param v2 a int argument, represents the second tetrahedron vertex
     * \param v3 a int argument, represents the third tetrahedron vertex
     * \param v4 a int argument, represents the fourth tetrahedron vertex
     */
    Tetrahedron(int v1, int v2, int v3, int v4) { this->set(v1,v2,v3,v4); }
    ///A destructor method
    virtual ~Tetrahedron() {}
    /**
     * @brief A public method that sets the current tetrahedron
     *
     * \param v1 a int argument, represents the first tetrahedron vertex
     * \param v2 a int argument, represents the second tetrahedron vertex
     * \param v3 a int argument, represents the third tetrahedron vertex
     * \param v4 a int argument, represents the fourth tetrahedron vertex
     */
    inline void set(int v1, int v2, int v3, int v4)
    {
        this->vertices[0] = v1;
        this->vertices[1] = v2;
        this->vertices[2] = v3;
        this->vertices[3] = v4;
    }
    ///A public method that returns the vertex at the pos-th position in the boundary array
    /*!
     * \param pos an integer argument, represents the vertex position into the list
     * \return an integer value, representing the position index of the vertex
     */
    inline int TV(int pos) const {  return abs(this->vertices[pos]); }
    /**
     * @brief A public procedure that updates the index of a vertex in the boundary array
     * @param pos an integer argument, represents the vertex position into the list
     * @param newId an integer representing the new vertex in the boundary
     */
    inline void setTV(int pos, int newId) { this->vertices[pos] = newId; }
    /**
     * @brief A public procedure that returns an edge in the boundary of the tetrahedron
     *
     * @param pos an integer representing the edge position in the boundary
     * @param e an integer vector that it is set with the sorted edge extrema
     */
    void TE(int pos, vector<int>& e);
    /**
     * @brief A public procedure that returns a triangular face in the boundary of the tetrahedron
     *
     * @param pos an integer representing the face position in the boundary
     * @param f an integer vector that it is set with the sorted face vertices
     */
    inline void TF(int pos, vector<int> &f)
    {
        f.assign(3,0);
        for(int i=0; i<3; i++)
            f[i] = this->TV((pos+i+1)%4);
        sort(f.begin(),f.end());
    }
    /**
     * @brief A public procedure that returns the tuple composed by the three indices forming the face and the tetrahedron id
     *
     * @param pos an integer representing the face position in the boundary
     * @param f a triangle_tetrahedron_tuple& that represents the tuple
     * @param t_id an integer representing the tetrahedron index
     */
    inline void face_tuple(int pos, triangle_tetrahedron_tuple &f, int t_id) const { f.sort_and_set(this->TV((pos+1)%4),this->TV((pos+2)%4),this->TV((pos+3)%4),t_id,pos); }
    /**
     * @brief A public procedure that checks if the tetrahedron has an input vertex_tetrahedron_struct
     *
     * @param v an integer representing the vertex index to search
     * @return true if the tetrahedron has the vertex in its boundary, false otherwise
     */
    inline bool has_vertex(int v) const
    {
        for(int i=0;i<vertices_num();i++)
        {
            if(v == TV(i))
            {
                return true;
            }
        }
        return false;
    }
    /**
     * @brief A public method that checks if a triangular face is on the mesh borders
     *
     * NOTA: a vertex index is make negative if the opposite triangular fase is on the mesh borders
     *       (see Border_Checker class for additional informations)
     *
     * @param pos an integer representing the edge position in the boundary
     * @return true if the edge is on the mesh borders, false otherwise
     */
    inline bool is_border_face(int pos) { return (this->vertices[pos] < 0); }

    /**
     * @brief operator ==
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator== (const Tetrahedron &p, const Tetrahedron &q)
    {
        bool b[p.vertices_num()];
        b[0] = false; b[1] = false; b[2] = false; b[3] = false;
        for(int i=0;i<p.vertices_num();i++)
        {
            for(int j=0;j<q.vertices_num();j++)
            {
                if(!b[j] && p.TV(i)==q.TV(j))
                {
                    b[j] = true;
                    break;
                }
            }
        }
        return b[0] && b[1] && b[2] && b[3];
    }
    /**
     * @brief operator !=
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator!= (const Tetrahedron &p, const Tetrahedron &q) { return !(p==q); }
    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Tetrahedron& p)
    {
        out <<"T[" << p.vertices[0] << " " << p.vertices[1] << " " << p.vertices[2] << " " << p.vertices[3] << "]";
        return out;
    }
    /**
     * @brief vertices_num
     * @return
     */
    inline int vertices_num() const { return 4; }

private:
    ///A private variable representing the array of vertices in the boundary
    int vertices[4];
};

#endif	/* _TETRAHEDRON_H */

