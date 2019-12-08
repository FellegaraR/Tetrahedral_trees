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

#ifndef _MESH_H
#define	_MESH_H

#include <vector>

#include "vertex.h"
#include "tetrahedron.h"
#include "box.h"

using namespace std;
///A class representing a tetrahedral mesh
class Mesh
{
public:
    ///A constructor method
    Mesh()
    {
        domain = Box();
        vertices = vector<Vertex>();
        tetrahedra = vector<Tetrahedron>();
    }
    ///A copy-constructor method
    Mesh(const Mesh& orig)
    {
        this->tetrahedra = orig.tetrahedra;
        this->vertices = orig.vertices;
        this->domain = orig.domain;
    }
    ///A destructor method
    virtual ~Mesh()
    {
        vertices.clear();
        tetrahedra.clear();
    }
    ///A public method that returns the vertex at the i-th position in the mesh list
    /*!
     * \param id an integer argument, representing the position in the list
     * \return a Vertex&, the vertex at the id-th position in the list
     */
    inline Vertex& get_vertex(int id) { return this->vertices[id-1]; }
    ///A public method that returns the tetrahedron at the i-th position in the mesh list
    /*!
     * \param id an integer argument, representing the position in the list
     * \return a Tetrahedron&, the tetrahedron at the id-th position in the list
     */
    inline Tetrahedron& get_tetrahedron(int id) { return this->tetrahedra[id-1]; }
    ///A public method that returns the mesh domain
    /*!
     * \return a Box&, the mesh domain
     */
    inline Box& get_domain() { return this->domain; }
    ///A public method that returns the number of mesh vertices
    /*!
     * \return an integer, representing the number of vertices
     */
    inline int get_num_vertices() { return this->vertices.size(); }
    ///A public method that returns the number of mesh tetrahedra
    /*!
     * \return an integer, representing the number of tetrahedra
     */
    inline int get_num_tetrahedra() { return this->tetrahedra.size(); }
    ///A public method that sets the mesh domain
    /*!
     * \param d a Box& argument, representing the domain to set
     */
    inline void set_domain(Box& d) { this->domain = d; }
    ///A public method that adds a vertex to the vertices list
    /*!
     * \param v a Vertex& argument, representing the vertex to add
     */
    inline void add_vertex(Vertex& v) { this->vertices.push_back(v); }
    ///A public method that adds a tetrahedron to the tetrahedra list
    /*!
     * \param t a Tetrahedron& argument, representing the tetrahedron to add
     */
    inline void add_tetrahedron(Tetrahedron& t) { this->tetrahedra.push_back(t); }
    ///A public method that initializes the space needed by the vertices and tetrahedra arrays
    /*!
     * \param numV an integer, represents the number of mesh vertices
     * \param numT an integer, represents the number of mesh tetrahedra
     */
    inline void reserve(int numV, int numT)
    {
        this->vertices.reserve(numV);
        this->tetrahedra.reserve(numT);
    }
    ///A public method that initializes the space needed by the vertices array
    /*!
     * \param numV an integer, represents the number of mesh vertices
     */
    inline void reserve_vertices_space(int numV) { this->vertices.reserve(numV); }
    ///A public method that resets the vertices array
    inline void reset_vertices() { this->vertices.clear(); }
    ///A public method that initializes the space needed by the tetrahedra array
    /*!
     * \param numT an integer, represents the number of mesh tetrahedra
     */
    inline void reserve_tetrahedra_space(int numT) { this->tetrahedra.reserve(numT); }
    ///A public method that resets the tetrahedra array
    inline void reset_tetrahedra() { this->tetrahedra.clear(); }

private:
    ///A private varible representing the mesh domain
    Box domain;
    ///A private varible representing the vertices array of the mesh
    vector<Vertex> vertices;
    ///A private varible representing the tetrahedra array of the mesh
    vector<Tetrahedron> tetrahedra;
};

#endif	/* _MESH_H */

