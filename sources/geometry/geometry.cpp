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

#include "geometry.h"

// Begin definitions for tetra_in_box

/******************************************************************************/

int Geometry::DetSign2D(double a, double b, double c, double d)
{
    double t1, t2;
    t1 = (a * d);
    t2 = (b * c);
    if (t1 > (t2 + ZERO)) return 1;
    if (t2 > (t1 + ZERO)) return -1;
    return 0;
}

/*
Calculate sign of determinant |a1 a2 a3|
                              |b1 b2 b3|
                              |c1 c2 c3|
 */
int Geometry::DetSign3D(double a1, double a2, double a3,
              double b1, double b2, double b3,
              double c1, double c2, double c3)
{
    double d = Det3D(a1, a2, a3, b1, b2, b3, c1, c2, c3);
    if (fabs(d) <= ZERO) return 0;
    return ( (d > 0.0) ? 1 : -1);
}

/*
Calculate sign of determinant |a1 a2 a3 a4|
                              |b1 b2 b3 b4|
                              |c1 c2 c3 c4|
                                                          |d1 d2 d3 d4|
 */
int Geometry::DetSign4D(double a1, double a2, double a3, double a4,
              double b1, double b2, double b3, double b4,
              double c1, double c2, double c3, double c4,
              double d1, double d2, double d3, double d4)
{
    double d = Det4D(a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4);
    if (fabs(d) <= ZERO) return 0;
    return ( (d > 0.0) ? 1 : -1);
}

/******************************************/
int Geometry::FourPointTurn(double x, double y, double z, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
    int d = DetSign3D((x1)-(x), (y1)-(y), (z1)-(z),
                      (x2)-(x1), (y2)-(y1), (z2)-(z1),
                      (x3)-(x1), (y3)-(y1), (z3)-(z1));
    return d;
}
/******************************************/

/******************************************************************************/

int Geometry::PointInTriangle2D(double x, double y,
                      double x1, double y1,
                      double x2, double y2,
                      double x3, double y3)
{
    if ((Geometry::PointTurn2D(x, y, x1, y1, x2, y2) == LEFT_TURN) &&
            (Geometry::PointTurn2D(x, y, x2, y2, x3, y3) == LEFT_TURN) &&
            (Geometry::PointTurn2D(x, y, x3, y3, x1, y1) == LEFT_TURN))
        return 1;
    if ((PointTurn2D(x, y, x1, y1, x2, y2) == RIGHT_TURN) &&
            (PointTurn2D(x, y, x2, y2, x3, y3) == RIGHT_TURN) &&
            (PointTurn2D(x, y, x3, y3, x1, y1) == RIGHT_TURN))
        return 1;
    return 0;
}

/* ------------------------------------------------------------------------ */
/*                       Intersection test w.r.t. a box                     */
/* ------------------------------------------------------------------------ */

/*
WITH A TWO-DIMENSIONAL BOX
Implementation based on the line clipping algorithm by Liang-Barsky.
The basic idea is the following:
Consider the parametric equation of the segment to be clipped:

x(u) = x1 + u*(x2-x1)             with  0<=u<=1
y(u) = y1 + u*(y2-y1)
The segment intersects the box if there exists u, 0<=u<=1, such that
minX <= x(u) <= maxX   and   minY <= y(u) <= maxY
or, equivalently
minX-x1 <= u*(x2-x1) <= maxX-x1   and   minY-y1 <= u*(y2-y1) <= maxY-y1

We have four inequalities that give four conditions on u.
Check if the four conditions are mutually consistent and
consistent with condition 0<=u<=1.
 */

/*
Restrict the admissible interval [u1,u2] by intersecting it with the
half-line, solution of inequality u*p <= q.
Return 1 if the resulting interval is not empty, 0 otherwise.
 */
int Geometry::ClipTest2D_strict(double p, double q, double * u1, double * u2)
{
    double r;

    if (p < 0.0) {
        r = q / p;
        if (r >= (*u2)) return 0;
        else if (r > (*u1)) (*u1) = r;
    } else {
        if (p > 0.0) {
            r = q / p;
            if (r <= (*u1)) return 0;
            else if (r < (*u2)) (*u2) = r;
        } else {
            /* p==0.0 line parallel to clipping edge */
            if (q <= 0.0) return 0;
        }
    }
    return 1;
}

/*
Return 1 if the segment is at least partially inside the box.
 */
int Geometry::ClipLine2D_strict(double minX, double minY, double maxX, double maxY, /* box */
                      double x1, double y1, double x2, double y2) /* line */
{
    double u1 = 0.0, u2 = 1.0; /* admissible interval, initially all [0,1] */
    double dx = x2 - x1, dy = y2 - y1;
    if ( ClipTest2D_strict(-dx, x1 - minX, &u1, &u2) &&
         ClipTest2D_strict(dx, maxX - x1, &u1, &u2) &&
         ClipTest2D_strict(-dy, y1 - minY, &u1, &u2) &&
         ClipTest2D_strict(dy, maxY - y1, &u1, &u2) )
    {
        return 1;
    }
    return 0;
}

/*
Return 1 iff edge (x1,y1) (x2,y2) overlaps edge having x=x0 and y01<=y<y02
 */
int Geometry::OverlapXSegment(double x1, double y1, double x2, double y2,
                    double x0, double y01, double y02)
{
    if ((x1 != x0) || (x2 != x0)) return 0; /* edge does not lie on x=x0 */
    if ((y1 <= y01) && (y2 <= y01)) return 0; /* edge <= y01 */
    if ((y1 >= y02) && (y2 >= y02)) return 0; /* edge >= y02 */
    return 1;
}

/*
Return 1 if the triangle is at least partially inside the box.
 */
int Geometry::ClipTriangle2D_strict(double minX, double minY, double maxX, double maxY, /* box */
                          double x[3], double y[3]) /* triangle */
{
    int i;
    /* if all vertices are on the same side of the box, no intersection */
    //  cout<<"controllo che tutti i vertici siano da uno stesso lato del triangolo"<<endl;
    if (x[0] <= minX && x[1] <= minX && x[2] <= minX) return 0;
    if (x[0] >= maxX && x[1] >= maxX && x[2] >= maxX) return 0;
    if (y[0] <= minY && y[1] <= minY && y[2] <= minY) return 0;
    if (y[0] >= maxY && y[1] >= maxY && y[2] >= maxY) return 0;


    /* if a vertex is inside the box, then the triangle intersects the box */
    for (i = 0; i < 3; i++) {
        if ((x[i] < maxX) && (x[i] > minX) && (y[i] < maxY) && (y[i] > minY)) {
            return 1;
        }
    }
    /* if an edge is at least partially inside the box, then the triangle
     intersects the box */
    for (i = 0; i < 3; i++) {
        if (ClipLine2D_strict(minX, minY, maxX, maxY,
                              x[i], y[i], x[(i + 1) % 3], y[(i + 1) % 3])) {
            return 1;
        }
    }

    /* none of the three triangle edges intersects the box */
    /* check if the triangle completely contains the box
     by applying the point-in-triangle test on the center of the box */
    if (PointInTriangle2D(0.5 * (minX + maxX), 0.5 * (minY + maxY),
                          x[0], y[0], x[1], y[1], x[2], y[2])) {
        return 1;
    }

    /* if one triangle edge is alinged with, and contains one box edge,
     then if the triangle is extended on the interior side of the box
     w.r.t. to that box edge, then there is intersection */
    for (i = 0; i < 3; i++)
    {
        if (OverlapXSegment(x[i], y[i], x[(i + 1) % 3], y[(i + 1) % 3], minX, minY, maxY)
                && (x[(i + 2) % 3] > minX)) /* third triangle vertex on the right */
            return 1;
        if (OverlapXSegment(x[i], y[i], x[(i + 1) % 3], y[(i + 1) % 3], maxX, minY, maxY)
                && (x[(i + 2) % 3] < maxX)) /* third triangle vertex on the left */
            return 1;
        if (OverlapXSegment(y[i], x[i], y[(i + 1) % 3], x[(i + 1) % 3], minY, minX, maxX)
                && (y[(i + 2) % 3] > minY)) /* third triangle vertex below */
            return 1;
        if (OverlapXSegment(y[i], x[i], y[(i + 1) % 3], x[(i + 1) % 3], maxY, minX, maxX)
                && (y[(i + 2) % 3] < maxY)) /* third triangle vertex below */
            return 1;
    }
    return 0;
}

/******************************************************************************/

/*
Here the tetrahedron is given as three arrays x[4], y[4], z[4] containing
the coordinates of its four vertices
 */

int Geometry::PointInTetra(double xp, double yp, double zp,
                 double * v1, double * v2, double * v3, double * v4)
{
    if (xp == v1[0] && yp == v1[1] && zp == v1[2]) return 1;
    if (xp == v2[0] && yp == v2[1] && zp == v2[2]) return 1;
    if (xp == v3[0] && yp == v3[1] && zp == v3[2]) return 1;
    if (xp == v4[0] && yp == v4[1] && zp == v4[2]) return 1;

    //  static double * v[4] = {v1, v2, v3, v4};
    int d1, d2, d3, d4, orientation;
    orientation = DetSign4D(v1[0], v1[1], v1[2], 1,
                            v2[0], v2[1], v2[2], 1,
                            v3[0], v3[1], v3[2], 1,
                            v4[0], v4[1], v4[2], 1);

    d1 = DetSign4D(xp, yp, zp, 1,
                   v2[0], v2[1], v2[2], 1,
                   v3[0], v3[1], v3[2], 1,
                   v4[0], v4[1], v4[2], 1);

    if (d1 != orientation && d1 != 0) return 0;

    d2 = DetSign4D(v1[0], v1[1], v1[2], 1,
                   xp, yp, zp, 1,
                   v3[0], v3[1], v3[2], 1,
                   v4[0], v4[1], v4[2], 1);

    if (d2 != orientation && d2 != 0) return 0;

    d3 = DetSign4D(v1[0], v1[1], v1[2], 1,
                   v2[0], v2[1], v2[2], 1,
                   xp, yp, zp, 1,
                   v4[0], v4[1], v4[2], 1);

    if (d3 != orientation && d3 != 0) return 0;

    d4 = DetSign4D(v1[0], v1[1], v1[2], 1,
                   v2[0], v2[1], v2[2], 1,
                   v3[0], v3[1], v3[2], 1,
                   xp, yp, zp, 1);

    if (d4 != orientation && d4 != 0) return 0;

    return 1;
}

int Geometry::PointInTetra_strict(double xp, double yp, double zp, double * v1, double * v2, double * v3, double * v4)
{
    double d1, d2, d3, d4, d;
    d = DetSign4D(v1[0], v1[1], v1[2], 1,
                  v2[0], v2[1], v2[2], 1,
                  v3[0], v3[1], v3[2], 1,
                  v4[0], v4[1], v4[2], 1);

    d1 = DetSign4D(xp, yp, zp, 1,
                   v2[0], v2[1], v2[2], 1,
                   v3[0], v3[1], v3[2], 1,
                   v4[0], v4[1], v4[2], 1);

    if (d1 != d) return 0;

    d2 = DetSign4D(v1[0], v1[1], v1[2], 1,
                   xp, yp, zp, 1,
                   v3[0], v3[1], v3[2], 1,
                   v4[0], v4[1], v4[2], 1);

    if (d2 != d) return 0;

    d3 = DetSign4D(v1[0], v1[1], v1[2], 1,
                   v2[0], v2[1], v2[2], 1,
                   xp, yp, zp, 1,
                   v4[0], v4[1], v4[2], 1);

    if (d3 != d) return 0;

    d4 = DetSign4D(v1[0], v1[1], v1[2], 1,
                   v2[0], v2[1], v2[2], 1,
                   v3[0], v3[1], v3[2], 1,
                   xp, yp, zp, 1);

    if (d4 != d) return 0;
    return 1;
}

/******************************************************************************/

/*
Restrict the admissible interval [u1,u2] by intersecting it with the
half-line, solution of inequality u*p <= q.
Return 1 if the resulting interval is not empty, 0 otherwise.
 */
int Geometry::ClipTest3D(double p, double q, double * u1, double * u2)
{
    double r;
    if (p < 0.0)
    {
        r = q / p;
        if (r > (*u2)) return 0;
        else if (r > (*u1)) (*u1) = r;
    }
    else
    {
        if (p > 0.0)
        {
            r = q / p;
            if (r < (*u1)) return 0;
            else if (r < (*u2)) (*u2) = r;
        }
        else
        {
            /* p==0.0 line parallel to clipping edge */
            if (q < 0.0) return 0;
        }
    }
    return 1;
}

int Geometry::ClipTest3D_strict(double p, double q, double * u1, double * u2)
{
    double r;
    if (p < 0.0)
    {
        r = q / p;
        if (r >= (*u2)) return 0;
        else if (r > (*u1)) (*u1) = r;
    }
    else
    {
        if (p > 0.0) {
            r = q / p;
            if (r <= (*u1)) return 0;
            else if (r < (*u2)) (*u2) = r;
        }
        else
        {
            if (q <= 0.0) return 0;
        }
    }
    return 1;
}

/******************************************************************************/

/*
Return 1 if the segment is at least partially inside the box.
 */
int Geometry::ClipLine3D(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
               double x1, double y1, double z1, double x2, double y2, double z2) /* line */
{
    double u1 = 0.0, u2 = 1.0; /* admissible interval, initially all [0,1] */
    double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;

    if (ClipTest3D(-dx, x1 - minX, &u1, &u2) &&
            ClipTest3D(dx, maxX - x1, &u1, &u2) &&
            ClipTest3D(-dy, y1 - minY, &u1, &u2) &&
            ClipTest3D(dy, maxY - y1, &u1, &u2) &&
            ClipTest3D(-dz, z1 - minZ, &u1, &u2) &&
            ClipTest3D(dz, maxZ - z1, &u1, &u2))
    {
        return 1;
    }
    return 0;
}

int Geometry::ClipLine3D_strict(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                      double x1, double y1, double z1, double x2, double y2, double z2) /* line */
{
    double u1 = 0.0, u2 = 1.0; /* admissible interval, initially all [0,1] */
    double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;

    if (ClipTest3D_strict(-dx, x1 - minX, &u1, &u2) &&
            ClipTest3D_strict(dx, maxX - x1, &u1, &u2) &&
            ClipTest3D_strict(-dy, y1 - minY, &u1, &u2) &&
            ClipTest3D_strict(dy, maxY - y1, &u1, &u2) &&
            ClipTest3D_strict(-dz, z1 - minZ, &u1, &u2) &&
            ClipTest3D_strict(dz, maxZ - z1, &u1, &u2))
    {
        return 1;
    }
    return 0;
}

int Geometry::ClipLine3D_middle(double minX, double minY, double minZ,
                      double maxX, double maxY, double maxZ, /* box */
                      double x1, double y1, double z1,
                      double x2, double y2, double z2) /* line */
{
    double u1 = 0.0, u2 = 1.0; /* admissible interval, initially all [0,1] */
    double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;

    if (ClipTest3D(-dx, x1 - minX, &u1, &u2) &&
            ClipTest3D_strict(dx, maxX - x1, &u1, &u2) &&
            ClipTest3D(-dy, y1 - minY, &u1, &u2) &&
            ClipTest3D_strict(dy, maxY - y1, &u1, &u2) &&
            ClipTest3D(-dz, z1 - minZ, &u1, &u2) &&
            ClipTest3D_strict(dz, maxZ - z1, &u1, &u2))
    {
        return 1;
    }

    return 0;
}

int Geometry::ClipLine3D_middle (double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                       double x1, double y1, double z1, double x2, double y2, double z2, /* line */
                       bool flag_minX, bool flag_maxX, bool flag_minY, bool flag_maxY, bool flag_minZ, bool flag_maxZ)
{
    double u1 = 0.0, u2 = 1.0; /* admissible interval, initially all [0,1] */
    double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;

    if(flag_minX)
    {
         if(!ClipTest3D(-dx, x1 - minX, &u1, &u2))
             return 0;
    }
    if(flag_maxX)
    {
        if(!ClipTest3D_strict(dx, maxX - x1, &u1, &u2))
            return 0;
    }

    if(flag_minY)
    {
        if(!ClipTest3D(-dy, y1 - minY, &u1, &u2))
            return 0;
    }
    if(flag_maxY)
    {
        if(!ClipTest3D_strict(dy, maxY - y1, &u1, &u2))
            return 0;
    }

    if(flag_minZ)
    {
        if(!ClipTest3D(-dz, z1 - minZ, &u1, &u2))
            return 0;
    }
    if(flag_maxZ)
    {
        if(!ClipTest3D_strict(dz, maxZ - z1, &u1, &u2))
            return 0;
    }
    return 1;
}

/******************************************************************************/

/*
Return 1 if the triangle is at least partially inside the box.
 */
int Geometry::ClipTriangle3D(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                   double x[3], double y[3], double z[3]) /* triangle */
{
    int i;

    /* if a vertex is inside the box, then the triangle intersects the box */
    for (i = 0; i < 3; i++)
    {
        if ((x[i] < maxX) && (x[i] > minX) && (y[i] < maxY) && (y[i] > minY) &&
                (z[i] < maxZ) && (z[i] > minZ))
        {
            return 1;
        }
    }
    /* if an edge is at least partially inside the box, then the triangle
     intersects the box */
    for (i = 0; i < 3; i++)
    {
        if (ClipLine3D(minX, minY, minZ, maxX, maxY, maxZ,
                       x[i], y[i], z[i],
                       x[(i + 1) % 3], y[(i + 1) % 3], z[(i + 1) % 3]))
        {
            return 1;
        }
    }
    /* none of the three triangle edges intersects the box */
    /* check if the triangle completely contains the box
     by applying the point-in-triangle test on any vertex of the box */
    if ((x[0] < maxX) && (x[0] > minX) &&
            (y[0] < maxY) && (y[0] > minY) &&
            (z[0] < maxZ) && (z[0] > minZ))
    {
        return 1;
    }
    return 0;
}

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
int Geometry::ClipTriangle3D_strict(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                          double x[3], double y[3], double z[3], /* triangle */
                          int flag, /* flag if minX, minY minZ faces are considered closed */
                          int flagMaxX, int flagMaxY, int flagMaxZ) /* flag if other faces are closed */
{
    int i;
    int s1, s2;
    /* if a vertex is inside the box, then the triangle intersects the box */
    for (i = 0; i < 3; i++)
    {
        if ((x[i] < maxX) && (x[i] > minX) && (y[i] < maxY) && (y[i] > minY) &&
                (z[i] < maxZ) && (z[i] > minZ))
        {
            return 1;
        }
    }
    /* if an edge is at least partially inside the box, then the triangle
     intersects the box */
    for (i = 0; i < 3; i++) {
        if (ClipLine3D_strict(minX, minY, minZ, maxX, maxY, maxZ,
                              x[i], y[i], z[i],
                              x[(i + 1) % 3], y[(i + 1) % 3], z[(i + 1) % 3]))
        {
            return 1;
        }
    }
    /* none of the three triangle edges intersects the box */
    /* check if the triangle completely cuts the box
     by checking if one of the edges of the box has its endpoints on opposite
     sides with respect to the triangle */
    s1 = DetSign4D(minX, minY, minZ, 1,
                   x[0], y[0], z[0], 1,
                   x[1], y[1], z[1], 1,
                   x[2], y[2], z[2], 1);
    s2 = DetSign4D(minX, minY, maxZ, 1,
                   x[0], y[0], z[0], 1,
                   x[1], y[1], z[1], 1,
                   x[2], y[2], z[2], 1);
    if ( (s1!=0) && (s1==-s2) )
        if (PointInTriangle2D(minX,minY, x[0],y[0], x[1],y[1], x[2],y[2]))
            return 1;
    s1 = DetSign4D(minX, minY, minZ, 1,
                   x[0], y[0], z[0], 1,
                   x[1], y[1], z[1], 1,
                   x[2], y[2], z[2], 1);
    s2 = DetSign4D(maxX, minY, minZ, 1,
                   x[0], y[0], z[0], 1,
                   x[1], y[1], z[1], 1,
                   x[2], y[2], z[2], 1);
    if ( (s1!=0) && (s1==-s2) )
        if (PointInTriangle2D(minY,minZ, y[0],z[0], y[1],z[1], y[2],z[2]))
            return 1;
    s1 = DetSign4D(minX, minY, minZ, 1,
                   x[0], y[0], z[0], 1,
                   x[1], y[1], z[1], 1,
                   x[2], y[2], z[2], 1);
    s2 = DetSign4D(minX, maxY, minZ, 1,
                   x[0], y[0], z[0], 1,
                   x[1], y[1], z[1], 1,
                   x[2], y[2], z[2], 1);
    if ( (s1!=0) && (s1==-s2) )
        if (PointInTriangle2D(minX,minZ, x[0],z[0], x[1],z[1], x[2],z[2]))
            return 1;

    /* return 1 also if the triangle is coplanar to one of the closed faces of the box
       and they properly overlap */
    if (flag==1)
    {
        // faccia minX
        if ( (x[0]==minX) && (x[1]==minX) && (x[2]==minX) )
            if (ClipTriangle2D_strict(minY,minZ, maxY,maxZ, y, z)) return 1;
        // faccia minY
        if ( (y[0]==minY) && (y[1]==minY) && (y[2]==minY) )
            if (ClipTriangle2D_strict(minX,minZ, maxX,maxZ, x, z)) return 1;
        // faccia minZ
        if ( (z[0]==minZ) && (z[1]==minZ) && (z[2]==minZ) )
            if (ClipTriangle2D_strict(minX,minY, maxX,maxY, x, y)) return 1;
    }

    /* return 1 also if the triangle is coplanar the indicated open face of the box
       and they properly  */
    if (flagMaxX==1)
    {
        // faccia maxX
        if ( (x[0]==maxX) && (x[1]==maxX) && (x[2]==maxX) )
            if (ClipTriangle2D_strict(minY,minZ, maxY,maxZ, y, z)) return 1;
    }
    if (flagMaxY==1)
    {
        // faccia maxY
        if ( (y[0]==maxY) && (y[1]==maxY) && (y[2]==maxY) )
            if (ClipTriangle2D_strict(minX,minZ, maxX,maxZ, x, z)) return 1;
    }
    if (flagMaxZ==1)
    {
        // faccia maxZ
        if ( (z[0]==maxZ) && (z[1]==maxZ) && (z[2]==maxZ) )
            if (ClipTriangle2D_strict(minX,minY, maxX,maxY, x, y)) return 1;
    }

    return 0;
}

/*
Check if edge (x1,y1,z1) - (x2,y2,z2) intersects triangle
whose x,y,z vertex coordinates are in the arrays
 */
int Geometry::EdgeIntersectTriangle_strict(double x1, double y1, double z1, double x2, double y2, double z2, /* edge */
                                 double x[3], double y[3], double z[3]) /* triangle */
{
    int turn1, turn2;
    turn1 = FourPointTurn(x1, y1, z1,
                          x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2]);
    turn2 = FourPointTurn(x2, y2, z2,
                          x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2]);
    /* if edge vertices on opposite sides w.r.t. triangle then intersection */
    if ((turn1 == -turn2) && (turn1 != NO_TURN)) return 1;
    /* else no intersection */
    return 0;
}


int Geometry::ClipTriangle3D_strict(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, /* box */
                          double x[3], double y[3], double z[3]) /* triangle */
{
    int i;

    /* if all vertices are on the same side of the box , no intersection */
    if ((x[0] <= minX) && (x[1] <= minX) && (x[2] <= minX)) return 0;
    if ((y[0] <= minY) && (y[1] <= minY) && (y[2] <= minY)) return 0;
    if ((z[0] <= minZ) && (z[1] <= minZ) && (z[2] <= minZ)) return 0;
    if ((x[0] >= maxX) && (x[1] >= maxX) && (x[2] >= maxX)) return 0;
    if ((y[0] >= maxY) && (y[1] >= maxY) && (y[2] >= maxY)) return 0;
    if ((z[0] >= maxZ) && (z[1] >= maxZ) && (z[2] >= maxZ)) return 0;

    /* if a vertex is inside the box, then the triangle intersects the box */
    for (i = 0; i < 3; i++)
    {
        if ((x[i] < maxX) && (x[i] > minX) && (y[i] < maxY) && (y[i] > minY) &&
                (z[i] < maxZ) && (z[i] > minZ))
        {
            return 1;
        }
    }
    /* if an edge is at least partially inside the box, then the triangle
     intersects the box */

    for (i = 0; i < 3; i++)
    {
        if (ClipLine3D_strict(minX, minY, minZ, maxX, maxY, maxZ,
                              x[i], y[i], z[i],
                              x[(i + 1) % 3], y[(i + 1) % 3], z[(i + 1) % 3]))
        {
            return 1;
        }
    }

    /* none of the three triangle edges intersects the box */

    /* triangle may have its edges on the box faces and its interior
     totally inside the box: test middle point of triangle */    
    double centerX = (x[0] + x[1] + x[2]) / 3.0;
    double centerY = (y[0] + y[1] + y[2]) / 3.0;
    double centerZ = (z[0] + z[1] + z[2]) / 3.0;
    if ((centerX > minX) && (centerX < maxX) &&
            (centerY > minY) && (centerY < maxY) &&
            (centerZ > minZ) && (centerZ < maxZ))
        return 1;

    /* box may intersect properly the interior of the triangle */
    /* check if some box edge intersects the triangle */
    /* check edges parallel to x axis */
    if (EdgeIntersectTriangle_strict(minX, minY, minZ, maxX, minY, minZ, x, y, z)
            && PointInTriangle2D(minY, minZ, y[0], z[0], y[1], z[1], y[2], z[2]))
        return 1;
    if (EdgeIntersectTriangle_strict(minX, maxY, minZ, maxX, maxY, minZ, x, y, z)
            && PointInTriangle2D(maxY, minZ, y[0], z[0], y[1], z[1], y[2], z[2]))
        return 1;
    if (EdgeIntersectTriangle_strict(minX, maxY, maxZ, maxX, maxY, maxZ, x, y, z)
            && PointInTriangle2D(maxY, maxZ, y[0], z[0], y[1], z[1], y[2], z[2]))
        return 1;
    if (EdgeIntersectTriangle_strict(minX, minY, maxZ, maxX, minY, maxZ, x, y, z)
            && PointInTriangle2D(minY, maxZ, y[0], z[0], y[1], z[1], y[2], z[2]))
        return 1;

    /* check edges parallel to y axis */
    if (EdgeIntersectTriangle_strict(minX, minY, minZ, minX, maxY, minZ, x, y, z)
            && PointInTriangle2D(minX, minZ, x[0], z[0], x[1], z[1], x[2], z[2]))
        return 1;
    if (EdgeIntersectTriangle_strict(minX, minY, maxZ, minX, maxY, maxZ, x, y, z)
            && PointInTriangle2D(minX, maxZ, x[0], z[0], x[1], z[1], x[2], z[2]))
        return 1;
    if (EdgeIntersectTriangle_strict(maxX, minY, minZ, maxX, maxY, minZ, x, y, z)
            && PointInTriangle2D(maxX, minZ, x[0], z[0], x[1], z[1], x[2], z[2]))
        return 1;
    if (EdgeIntersectTriangle_strict(maxX, minY, maxZ, maxX, maxY, maxZ, x, y, z)
            && PointInTriangle2D(maxX, maxZ, x[0], z[0], x[1], z[1], x[2], z[2]))
        return 1;

    /* check edges parallel to z axis */
    if (EdgeIntersectTriangle_strict(minX, minY, minZ, minX, minY, maxZ, x, y, z)
            && PointInTriangle2D(minX, minY, x[0], y[0], x[1], y[1], x[2], y[2]))
        return 1;
    if (EdgeIntersectTriangle_strict(minX, maxY, minZ, minX, maxY, maxZ, x, y, z)
            && PointInTriangle2D(minX, maxY, x[0], y[0], x[1], y[1], x[2], y[2]))
        return 1;
    if (EdgeIntersectTriangle_strict(maxX, minY, minZ, maxX, minY, maxZ, x, y, z)
            && PointInTriangle2D(maxX, minY, x[0], y[0], x[1], y[1], x[2], y[2]))
        return 1;
    if (EdgeIntersectTriangle_strict(maxX, maxY, minZ, maxX, maxY, maxZ, x, y, z)
            && PointInTriangle2D(maxX, maxY, x[0], y[0], x[1], y[1], x[2], y[2]))
        return 1;

    /* otherwise no intersection */
    return 0;
}

// this version of tetra_in_box execute the following tests:
// - if all vertices are on the same side of the box, no intersection
// - consider all the faces of the box open (if a vertex is adjacent to one of these face is considered external)
// - check if one of the vertex is inside the tetrahedron (PointInTetra_strict)
// - check if the box center is inside the tetrahedron (PointInTetra_strict)
// - check if one triangular facet intersects box (ClipTriangle3D_strict)
// - check if some triangular face of the tetrahedron is coplanar to a squared face of the box and the rest of
//   tetrahedron lies on the interior side of the box defined by such squared face (ClipTriangle2D_strict)
//This version is suited for box query, while for building we have to do external tests if we use different assumption on the closeness of box faces
int Geometry::tetra_in_box_strict(double * minF, double * maxF, double **c)
{
    /* if all vertices are on the same side of the box, no intersection */
    for (int j = 0; j < 3; j++)
    {
        if ((c[0][j] <= minF[j]) && (c[1][j] <= minF[j]) &&
                (c[2][j] <= minF[j]) && (c[3][j] <= minF[j]))
        {
            return 0;
        }
        if ((c[0][j] >= maxF[j]) && (c[1][j] >= maxF[j]) &&
                (c[2][j] >= maxF[j]) && (c[3][j] >= maxF[j]))
        {
            return 0;
        }
    }

    int i, j;

    /* check if one tetrahedron vertex is inside box */
    for (i = 0; i < 4; i++)
    {
        if ((minF[0] < c[i][0]) && (c[i][0] < maxF[0]) &&
                (minF[1] < c[i][1]) && (c[i][1] < maxF[1]) &&
                (minF[2] < c[i][2]) && (c[i][2] < maxF[2]))
        {
            return 1;
        }
    }

    /* check if one box vertex is inside tetrahedron */
    if (PointInTetra_strict(minF[0], minF[1], minF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra_strict(minF[0], minF[1], maxF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra_strict(minF[0], maxF[1], minF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra_strict(maxF[0], minF[1], minF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra_strict(maxF[0], maxF[1], maxF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra_strict(maxF[0], maxF[1], minF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra_strict(maxF[0], minF[1], maxF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra_strict(minF[0], maxF[1], maxF[2], c[0], c[1], c[2], c[3]))
    {
        return 1;
    }

    /* check if box center is inside tetrahedron */
    if (PointInTetra_strict(0.5 * (minF[0] + maxF[0]),
                            0.5 * (minF[1] + maxF[1]),
                            0.5 * (minF[2] + maxF[2]), c[0], c[1], c[2], c[3]))
    {
        return 1;
    }

    /* check if one triangular facet intersects box */
    double* x = new double[3];
    double* y = new double[3];
    double* z = new double[3];

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 3; j++)
        {
            x[j] = c[(i + j) % 4][0];
            y[j] = c[(i + j) % 4][1];
            z[j] = c[(i + j) % 4][2];
        }
        if (ClipTriangle3D_strict(minF[0], minF[1], minF[2],
                                  maxF[0], maxF[1], maxF[2], x, y, z))
        {
            //DEALLOCAZIONE
            delete x;
            delete y;
            delete z;
            //DEALLOCAZIONE
            return 1;
        }
    }

    //DEALLOCAZIONE
    delete x;
    delete y;
    delete z;
    //DEALLOCAZIONE

    /* check if some triangular face of the tetrahedron is coplanar to a
     squared face of the box and the rest of tetrahedron lies on the interior
     side of the box defined by such squared face. */
    double* xt = new double[3];
    double* yt = new double[3];
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 3; j++)
        {
            /* if all three triangle vertices are on the plane with min j-coordinate */
            if ((c[i][j] == minF[j]) && (c[(i + 1) % 4][j] == minF[j]) && (c[(i + 2) % 4][j] == minF[j]))
            {
                /* and if triangular and squared face intersect on plane */
                xt[0] = c[i][(j + 1) % 3];
                xt[1] = c[(i + 1) % 4][(j + 1) % 3];
                xt[2] = c[(i + 2) % 4][(j + 1) % 3];

                yt[0] = c[i][(j + 2) % 3];
                yt[1] = c[(i + 1) % 4][(j + 2) % 3];
                yt[2] = c[(i + 2) % 4][(j + 2) % 3];

                if (ClipTriangle2D_strict(minF[(j + 1) % 3], minF[(j + 2) % 3],
                                          maxF[(j + 1) % 3], maxF[(j + 2) % 3], xt, yt))
                { /* if fourth vertex has larger j-coordinate, then they intersect */
                    if (c[(i + 3) % 4][j] > minF[j])
                    {
                        //DEALLOCAZIONE
                        delete xt;
                        delete yt;
                        //DEALLOCAZIONE
                        return 1;
                    }
                }

            }
            /* if all three triangle vertices are on the plane with maxn j-coordinate */
            if ((c[i][j] == maxF[j]) && (c[(i + 1) % 4][j] == maxF[j]) && (c[(i + 2) % 4][j] == maxF[j]))
            {
                /* and if triangular and squared face intersect on plane */
                xt[0] = c[i][(j + 1) % 3];
                xt[1] = c[(i + 1) % 4][(j + 1) % 3];
                xt[2] = c[(i + 2) % 4][(j + 1) % 3];

                yt[0] = c[i][(j + 2) % 3];
                yt[1] = c[(i + 1) % 4][(j + 2) % 3];
                yt[2] = c[(i + 2) % 4][(j + 2) % 3];

                if (ClipTriangle2D_strict(minF[(j + 1) % 3], minF[(j + 2) % 3],
                                          maxF[(j + 1) % 3], maxF[(j + 2) % 3], xt, yt))
                { /* if fourth vertex has smaller j-coordinate, then they intersect */
                    if (c[(i + 3) % 4][j] < maxF[j])
                    {
                        //DEALLOCAZIONE
                        delete xt;
                        delete yt;
                        //DEALLOCAZIONE
                        return 1;
                    }
                }

            }
        }
    }
    //DEALLOCAZIONE
    delete xt;
    delete yt;
    //DEALLOCAZIONE
    /* box and tetrahedron do not intersect each other */
    return 0;
}


/******************************************************************************/

/* intersezione box e tetraedro */
// this version of tetra_in_box execute the following tests:
// - consider all the faces of the box closed (if a vertex is adjacent to one of these face is considered internal)
// - check if one of the vertex is inside the tetrahedron (PointInTetra)
// - check if one triangular facet intersects box (ClipTriangle3D)
int Geometry::tetra_in_box(double * minF, double * maxF, double **c) {

    int i, j;


    /* check if one tetrahedron vertex is inside box */
    for (i = 0; i < 4; i++)
    {
        if ((minF[0] <= c[i][0]) && (c[i][0] <= maxF[0]) &&
                (minF[1] <= c[i][1]) && (c[i][1] <= maxF[1]) &&
                (minF[2] <= c[i][2]) && (c[i][2] <= maxF[2]))
        {
            return 1;
        }
    }

    /* check if one box vertex is inside tetrahedron */
    if (PointInTetra(minF[0], minF[1], minF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra(minF[0], minF[1], maxF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra(minF[0], maxF[1], minF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra(maxF[0], minF[1], minF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra(maxF[0], maxF[1], maxF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra(maxF[0], maxF[1], minF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra(maxF[0], minF[1], maxF[2], c[0], c[1], c[2], c[3]) ||
            PointInTetra(minF[0], maxF[1], maxF[2], c[0], c[1], c[2], c[3]))
    {
        return 1;
    }

    /* check if one triangular facet intersects box */
    double* x = new double[3]; /* vertex coordinates */
    double* y = new double[3];
    double* z = new double[3];

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 3; j++)
        {
            x[j] = c[(i + j) % 4][0];
            y[j] = c[(i + j) % 4][1];
            z[j] = c[(i + j) % 4][2];
        }
        if (ClipTriangle3D(minF[0], minF[1], minF[2],
                           maxF[0], maxF[1], maxF[2], x, y, z))
        {
            //DEALLOCAZIONE
            delete x;
            delete y;
            delete z;
            //DEALLOCAZIONE
            return 1;
        }
    }

    //DEALLOCAZIONE
    delete x;
    delete y;
    delete z;
    //DEALLOCAZIONE

    /* box and tetrahedron do not intersect each other */
    return 0;
}

// End definitions for tetra_in_box
