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

#ifndef _POINT_H
#define	_POINT_H

#include <iostream>
#include <cmath>
using namespace std;

///A class representing a tridimensional point in the space
class Point {
public:
    ///A constructor method
    Point()
    {
        this->coords[0]=0;
        this->coords[1]=0;
        this->coords[2]=0;
    }
    ///A copy-constructor method
    Point(const Point& orig)
    {
        this->coords[0] = orig.coords[0];
        this->coords[1] = orig.coords[1];
        this->coords[2] = orig.coords[2];
    }
    ///A constructor method
    /*!
     * \param x a double argument, representing the x coordinate
     * \param y a double argument, representing the y coordinate
     * \param z a double argument, representing the z coordinate
     */
    Point(double x, double y, double z)
    {
        this->coords[0]=x;
        this->coords[1]=y;
        this->coords[2]=z;
    }
    ///A destructor method
    virtual ~Point() {}
    /**
     * @brief operator ==
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator== (const Point& p, const Point &q) { return ((p.coords[0] == q.coords[0]) && (p.coords[1] == q.coords[1]) && (p.coords[2] == q.coords[2])); }
    /**
     * @brief operator !=
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator!= (const Point& p, const Point &q) { return !(p == q); }
    /**
     * @brief operator <
     * @param s
     * @return
     */
    inline bool operator<(const Point& s) const
    {
        return ((this->coords[0] < s.coords[0]) ||
                (this->coords[0] == s.coords[0] && this->coords[1] < s.coords[1]) ||
                (this->coords[0] == s.coords[0] && this->coords[1] == s.coords[1] && this->coords[2] < s.coords[2]));
    }
    /**
     * @brief operator >
     * @param s
     * @return
     */
    inline bool operator>(const Point& s) const
    {
        return ((this->coords[0] > s.coords[0]) ||
                (this->coords[0] == s.coords[0] && this->coords[1] > s.coords[1]) ||
                (this->coords[0] == s.coords[0] && this->coords[1] == s.coords[1] && this->coords[2] > s.coords[2]));
    }
    /**
     * @brief operator +
     * @param s
     * @return
     */
    inline Point operator+(const Point &s) const { return Point(this->coords[0]+s.get_x(),this->coords[1]+s.get_y(),this->coords[2]+s.get_z()); }
    /**
     * @brief operator -
     * @param s
     * @return
     */
    inline Point operator-(const Point &s) const { return Point(this->coords[0]-s.get_x(),this->coords[1]-s.get_y(),this->coords[2]-s.get_z()); }
    /**
     * @brief operator *
     * @param f
     * @return
     */
    inline Point operator*(const double &f) const { return Point(this->coords[0]*f,this->coords[1]*f,this->coords[2]*f); }
    ///
    inline friend std::ostream& operator<<(std::ostream& out, const Point& p)
    {
        out << p.coords[0] << " " << p.coords[1] << " " << p.coords[2] ;
        return out;
    }
    ///A public method that returns the x coordinate
    /*!
     * \return a double representing the x coordinate
     */
    inline double get_x() const { return this->coords[0]; }
    ///A public method that returns the y coordinate
    /*!
     * \return a double representing the y coordinate
     */
    inline double get_y() const { return this->coords[1]; }
    ///A public method that returns the z coordinate
    /*!
     * \return a double representing the z coordinate
     */
    inline double get_z() const { return this->coords[2]; }
    /**
     * @brief A public procedure returning the coordinate at a given position
     *
     * @param pos an integer representing the coordinate position in the point array
     * @return a double representing the coordinate
     */
    inline double get_c(int pos) const { return this->coords[pos]; }
    /**
     * @brief A public procedure that initializes a coordinate of the point
     *
     * @param pos an integer representing the coordinate position in the point array
     * @param c a double representing the coordinate value
     */
    inline void set_c(int pos, double c) { this->coords[pos] = c; }
    /**
     * @brief A public method that initializes all the coordinates of a point
     * @param x a double representing the value on the x-axis
     * @param y a double representing the value on the y-axis
     * @param z a double representing the value on the z-axis
     */
    inline void set(double x, double y, double z)
    {
        this->coords[0]=x;
        this->coords[1]=y;
        this->coords[2]=z;
    }

    /**
     * @brief A public method that initializes all the coordinates of a point
     * @param p a Point& argument containing the the coordinates value to set
     */
    inline void set(Point &p)
    {
        set(p.get_x(),p.get_y(),p.get_z());
    }
    /**
     * @brief A public procedure that computes the norm value of the vector this-v
     *
     * @param v a Point& representing the other side of the vector
     * @return the norm value
     */
    inline double norm_3D(Point& v) { return(sqrt(((v.get_x()-coords[0])*(v.get_x()-coords[0]))
                                                  +((v.get_y()-coords[1])*(v.get_y()-coords[1]))
                                                  +((v.get_z()-coords[2])*(v.get_z()-coords[2])))); }
    /**
     * @brief A public procedure that computes the norm value of th vector composed by the point
     *
     * @return the norm value
     */
    inline double norm_3D() { return sqrt(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]); }

    //this function returns the scalar products between vectors v1-vec and v2-vec
    /**
     * @brief A public procedure returning the cross product of the vectors v1-this and v2-this
     *
     * @param v1 a Point&, a point
     * @param v2 a Point&, a point
     * @return the cross product value
     */
    inline double cross_3D(Point& v1,Point& v2) { return(((v1.get_x()-coords[0])*(v2.get_x()-coords[0]))
                                                           +((v1.get_y()-coords[1])*(v2.get_y()-coords[1]))
                                                           +((v1.get_z()-coords[2])*(v2.get_z()-coords[2]))); }
    /**
     * @brief A public procedure returning the cross product of the vectors v1 and this
     *
     * @param v1 a Vertex& representing a vertex
     * @return the cross product value
     */
    inline Point cross_3D(Point &v1) {
        return Point(coords[1]*v1.get_z() - v1.get_y()*coords[2], coords[2]*v1.get_x() - v1.get_z()*coords[0], coords[0]*v1.get_y() - v1.get_x()*coords[1]);
    }
    /**
     * @brief A public procedure returning the dot product between two vertices
     *
     * @param v1 a Point&, a point
     * @return the dot product value
     */
    inline double dot_3D(const Point &v1) const { return (coords[0]*v1.get_x() + coords[1]*v1.get_y() + coords[2]*v1.get_z()); }
    /**
     * @brief A public procedure that computes the distance between two points
     * @param v a Point&, a point
     * @return the distance value
     */
    inline double distance_3D(const Point& v) const
    {
        double xdist = this->coords[0]-v.get_x();
        double ydist = this->coords[1]-v.get_y();
        double zdist = this->coords[2]-v.get_z();
        return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
    }
    /**
     * @brief get_dimension
     * @return
     */
    inline int get_dimension() const { return 3; }


protected:
    /// A protected array representing the point coordinates
    double coords[3];
};

#endif	/* _POINT_H */

