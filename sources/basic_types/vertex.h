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

#ifndef _VERTEX_H
#define	_VERTEX_H

#include "point.h"
using namespace std;

///An inner-class, extending Point, representing a vertex in a tetrahedral mesh
class Vertex : public Point
{
public:
    ///A constructor method
    Vertex() : Point() { this->field_value = 0; }
    ///A copy-constructor method
    Vertex(const Vertex& orig) : Point(orig) { this->field_value=orig.field_value; }
    ///A constructor method
    /*!
     * \param x a double argument, representing the x coordinate
     * \param y a double argument, representing the y coordinate
     * \param z a double argument, representing the z coordinate
     * \param field a double argument, representing the vertex field
     */
    Vertex(double x, double y, double z, double field) : Point(x,y,z)
    {
       this->field_value = field;
    }
    ///A destructor method
    virtual ~Vertex() { }
    ///A public method that returns the vertex field value
    /*!
     * \return an integer, representing the field value
     */
    inline double get_field() { return field_value; }
    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Vertex& p)
    {
        out <<"[" << p.coords[0] << " " << p.coords[1] << " " << p.coords[2] << " " << p.field_value << "]";
        return out;
    }
    //this function returns the norm of vector vec-v
    /**
     * @brief A public procedure that computes the norm value of the vector this-v
     *
     * @param v a Vertex& representing the other side of the vector
     * @return the norm value
     */
    inline double norm(Vertex& v)
    {
        double xdist = v.get_x()-this->coords[0];
        double ydist = v.get_y()-this->coords[1];
        double zdist = v.get_z()-this->coords[2];
        double fdist = v.get_field()-this->field_value;

        return sqrt(xdist*xdist + ydist*ydist + zdist*zdist + fdist*fdist);
    }
    //this function returns the scalar products between vectors v1-vec and v2-vec
    /**
     * @brief A public procedure returning the scalar product of the vectors v1-this and v2-this
     *
     * @param v1 a Vertex& representing a vertex
     * @param v2 a Vertex& representing a vertex
     * @return the scalar product value
     */
    inline double scalar_product(Vertex& v1,Vertex& v2)
    {
        return(((v1.get_x()-this->coords[0])*(v2.get_x()-this->coords[0]))
               +((v1.get_y()-this->coords[1])*(v2.get_y()-this->coords[1]))
               +((v1.get_z()-this->coords[2])*(v2.get_z()-this->coords[2]))
               +((v1.get_field()-this->field_value)*(v2.get_field()-this->field_value)));
    }
private:
    ///A private variable that represents the vertex field value
    double field_value;
};

#endif	/* _VERTEX_H */

