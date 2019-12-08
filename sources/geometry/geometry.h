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

#ifndef _GEOMETRY_H
#define	_GEOMETRY_H

#include <cmath>
#include <iostream>
#include <stdio.h>
using namespace std;

// ------------ defines and macros for tetra_in_box --------------------------
#ifndef PI
#define PI 3.14159265358979323846
#endif
/* ------------------------------------------------------------------------ */
/*                    Tolerance used in computations                        */
/* ------------------------------------------------------------------------ */
#define ZERO (10E-14)
#define Coincide(a,b) (fabs((a)-(b))<=ZERO)
/* ------------------------------------------------------------------------ */
/*                     Basic arithmetic functions                           */
/* ------------------------------------------------------------------------ */
#define Square(x) ((x)*(x))
/*
Calculate determinant  |a b|
                       |c d|
*/
#define Det2D(a,b,c,d)  ( (a*d)-(b*c) )
/*
Calculate determinant  |a1 a2 a3|
                       |b1 b2 b3|
                       |c1 c2 c3|
*/
#define Det3D(a1,a2,a3,b1,b2,b3,c1,c2,c3)	\
    ( a1*Det2D(b2,b3,c2,c3) - a2*Det2D(b1,b3,c1,c3) + a3*Det2D(b1,b2,c1,c2) )

/*
Calculate determinant |a1 a2 a3 a4|
                      |b1 b2 b3 b4|
                      |c1 c2 c3 c4|
                      |d1 d2 d3 d4|
*/
#define Det4D(a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,d1,d2,d3,d4)	\
    ( a1*Det3D(b2,b3,b4,c2,c3,c4,d2,d3,d4) - \
    a2*Det3D(b1,b3,b4,c1,c3,c4,d1,d3,d4) + \
    a3*Det3D(b1,b2,b4,c1,c2,c4,d1,d2,d4) - \
    a4*Det3D(b1,b2,b3,c1,c2,c3,d1,d2,d3))
/*
 *  Turns
 */
#define LEFT_TURN -1
#define UP_TURN -1
#define NO_TURN 0
#define RIGHT_TURN 1
#define DOWN_TURN 1

class Geometry
{
protected:
    // ------------ Main functions -------------
    // this version of tetra_in_box execute the following tests:
    // - if all vertices are on the same side of the box, no intersection
    // - consider all the faces of the box open (if a vertex is adjacent to one of these face is considered external)
    // - check if one of the vertex is inside the tetrahedron (PointInTetra_strict)
    // - check if the box center is inside the tetrahedron (PointInTetra_strict)
    // - check if one triangular facet intersects box (ClipTriangle3D_strict)
    // - check if some triangular face of the tetrahedron is coplanar to a squared face of the box and the rest of
    //   tetrahedron lies on the interior side of the box defined by such squared face (ClipTriangle2D_strict)
    //This version is suited for box query, while for building we have to do external tests if we use different assumption on the closeness of box faces
    static int tetra_in_box_strict(double * minF, double * maxF, double **c);

    // this version of tetra_in_box execute the following tests:
    // - consider all the faces of the box closed (if a vertex is adjacent to one of these face is considered internal)
    // - check if one of the vertex is inside the tetrahedron (PointInTetra)
    // - check if one triangular facet intersects box (ClipTriangle3D)
    static int tetra_in_box(double * minF, double * maxF, double **c);

    // --------------- Auxiliaries for tetra_in_box -----------------

    static int PointInTriangle2D (double x, double y,
                           double x1, double y1,
                           double x2, double y2,
                           double x3, double y3);

    /*
    Return 1 if the triangle is at least partially inside the box.
    */
    static int ClipTriangle3D(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                       double x[3], double y[3], double z[3]); /* triangle */

    /*
    Return 1 if the triangle is at least partially inside the box.
    Version returning false if the triangle is tangent to the box.
    */
    static int ClipTriangle3D_strict(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                              double x[3], double y[3], double z[3]); /* triangle */

    /*
    Return 1 if the triangle is at least partially inside the box.
    Version returning false if the triangle is tangent to the box.
    The flag means:
    0 - result is 1 iff there is proper intersection = the (2D) interior of the triangle
        intersects the (3D) interior of the box
    1 - result is 1 also if the (2D) interior of the triangle intersects the (2D) interior
        of one of the faces of the box, and this face corresponds to the plane of minX,
        minY, or minZ (i.e., one of the faces considered as closed in the spatial index)
    The further three flags mean
    1 - result is 1 also if the (2D) interior of the triangle intersects the (2D) interior
        of the face of the box corresponding to the specified max coordinate
        (this is to be used in the case of boxes adjacent to the max x,y,z of the entire domain
    */
    static int ClipTriangle3D_strict(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                              double x[3], double y[3], double z[3], /* triangle */
                              int flag, /* flag if minX, minY minZ faces are considered closed */
                              int flagMaxX, int flagMaxY, int flagMaxZ); /* flag if other faces are closed */

    /*
    Return 1 if the segment is at least partially inside the box.
    */
    static int ClipLine3D (double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                    double x1, double y1, double z1, double x2, double y2, double z2); /* line */
    //esclude incidenza in un vertice o lato della box
    static int ClipLine3D_strict (double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                           double x1, double y1, double z1, double x2, double y2, double z2); /* line */
    //considera chiusi solo i tre edge incidenti nel punto minimo della box
    static int ClipLine3D_middle (double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                           double x1, double y1, double z1, double x2, double y2, double z2); /* line */
    static int ClipLine3D_middle (double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                           double x1, double y1, double z1, double x2, double y2, double z2, /* line */
                           bool flag_minX, bool flag_maxX, bool flag_minY, bool flag_maxY, bool flag_minZ, bool flag_maxZ);

    /*
    Check if edge (x1,y1,z1) - (x2,y2,z2) intersects triangle
    whose x,y,z vertex coordinates are in the arrays
    */
    static int EdgeIntersectTriangle_strict(double x1, double y1, double z1, double x2, double y2, double z2, /* edge */
                                     double x[3], double y[3], double z[3]); /* triangle */

    /*
    Restrict the admissible interval [u1,u2] by intersecting it with the
    half-line, solution of inequality u*p <= q.
    Return 1 if the resulting interval is not empty, 0 otherwise.
    */
    static int ClipTest3D (double p, double q, double * u1, double * u2);
    static int ClipTest3D_strict (double p, double q, double * u1, double * u2);

    /*
    Here the tetrahedron is given as three arrays x[4], y[4], z[4] containing
    the coordinates of its four vertices
    */
    static int PointInTetra (double xp, double yp, double zp, double * v1, double * v2, double * v3, double * v4);
    static int PointInTetra_strict(double xp, double yp, double zp, double * v1, double * v2, double * v3, double * v4);


    /*
    Calculate sign of determinant  |a b|
                                   |c d|
    */
    static int DetSign2D (double a, double b, double c, double d);
    /*
    Calculate sign of determinant |a1 a2 a3|
                                  |b1 b2 b3|
                                  |c1 c2 c3|
    */
    static int DetSign3D(double a1, double a2, double a3,
                  double b1, double b2, double b3,
                  double c1, double c2, double c3);
    /*
    Calculate sign of determinant |a1 a2 a3 a4|
                                  |b1 b2 b3 b4|
                                  |c1 c2 c3 c4|
                                  |d1 d2 d3 d4|
    */
    static int DetSign4D(double a1, double a2, double a3, double a4,
                  double b1, double b2, double b3, double b4,
                  double c1, double c2, double c3, double c4,
                  double d1, double d2, double d3, double d4 );

    /* ------------------------------------------------------------------------ */
    /*                                 Turns                                    */
    /* ------------------------------------------------------------------------ */

    static int FourPointTurn(double x, double y, double z, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);
    static int PointTurn2D(double x,double y,double x1,double y1,double x2,double y2)	{ return Geometry::DetSign2D((x)-(x1), (y)-(y1), (x2)-(x1), (y2)-(y1)); }

    /* ------------------------------------------------------------------------ */
    /*                       Intersection test w.r.t. a box (2D)                */
    /* ------------------------------------------------------------------------ */
    /*
    Restrict the admissible interval [u1,u2] by intersecting it with the
    half-line, solution of inequality u*p <= q.
    Return 1 if the resulting interval is not empty, 0 otherwise.
     */
    static int ClipTest2D_strict(double p, double q, double * u1, double * u2);
    /*
    Return 1 if the segment is at least partially inside the box.
     */
    static int ClipLine2D_strict(double minX, double minY, double maxX, double maxY, /* box */
                          double x1, double y1, double x2, double y2); /* line */
    /*
    Return 1 iff edge (x1,y1) (x2,y2) overlaps edge having x=x0 and y01<=y<y02
     */
    static int OverlapXSegment(double x1, double y1, double x2, double y2,
                        double x0, double y01, double y02);
    /*
    Return 1 if the triangle is at least partially inside the box.
     */
    static int ClipTriangle2D_strict(double minX, double minY, double maxX, double maxY, /* box */
                              double x[3], double y[3]); /* triangle */
};



#endif	/* _GEOMETRY_H */

