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

#include "geometry_distortion.h"

// ------------------------------------------------------ (1) --------------------------------------------------------- //
// ------------------------------------------------------ (1) --------------------------------------------------------- //

double Geometry_Distortion::get_trihedral_angle(Tetrahedron &t, int v, Mesh& mesh)
{
    int other_vert[3];
    double    prodscalv1vv2=0, prodscalv1vv3=0, prodscalv2vv3=0;
    double    normavv1=0,normavv2=0,normavv3=0;

    int j = 0;
    for(int pos = 0; pos < t.vertices_num(); pos++) /// before <= ???? (changed to <)
    {
        if( t.TV(pos) != v )
        {
            other_vert[j] = t.TV(pos);
            j++;
        }
    }

    //vertices coordinates
    Vertex& v0 = mesh.get_vertex(v);
    Vertex& v1 = mesh.get_vertex(other_vert[1]);
    Vertex& v2 = mesh.get_vertex(other_vert[2]);
    Vertex& v3 = mesh.get_vertex(other_vert[0]);

    //computing all the scalar product
    prodscalv1vv2 = v0.scalar_product(v1,v2);
    prodscalv1vv3 = v0.scalar_product(v1,v3);
    prodscalv2vv3 = v0.scalar_product(v2,v3);

    //computing all norm
    normavv1 = v0.norm(v1);
    normavv2 = v0.norm(v2);
    normavv3 = v0.norm(v3);

    return computeTrihedralAngle(prodscalv1vv2,prodscalv1vv3,prodscalv2vv3,normavv1,normavv2,normavv3);
}

double Geometry_Distortion::get_trihedral_angle_3D(Tetrahedron &t, int v, Mesh &mesh)
{
    int other_vert[3];
    double prodscalv1vv2=0, prodscalv1vv3=0, prodscalv2vv3=0;
    double normavv1=0,normavv2=0,normavv3=0;

    int j = 0;
    for(int pos = 0; pos < 4; pos++) /// before <= ???? (changed to <)
    {
        if( t.TV(pos) != v )
        {
            other_vert[j] = t.TV(pos);
            j++;
        }
    }

    //vertices coordinates
    Vertex& v0 = mesh.get_vertex(v);
    Vertex& v1 = mesh.get_vertex(other_vert[1]);
    Vertex& v2 = mesh.get_vertex(other_vert[2]);
    Vertex& v3 = mesh.get_vertex(other_vert[0]);

    //computing all the scalar product
    prodscalv1vv2 = v0.cross_3D(v1,v2);
    prodscalv1vv3 = v0.cross_3D(v1,v3);
    prodscalv2vv3 = v0.cross_3D(v2,v3);

    //computing all norm
    normavv1 = v0.norm_3D(v1);
    normavv2 = v0.norm_3D(v2);
    normavv3 = v0.norm_3D(v3);

    return computeTrihedralAngle(prodscalv1vv2,prodscalv1vv3,prodscalv2vv3,normavv1,normavv2,normavv3);
}

double Geometry_Distortion::computeTrihedralAngle(double prodscalv1vv2, double prodscalv1vv3, double prodscalv2vv3, double normavv1, double normavv2, double normavv3)
{
    double cosalpha=0, cosbeta=0, cosgamma=0;
    double senalpha=0, senbeta=0, sengamma=0;
    //dihedral angles
    double A=0,B=0,C=0;

    //computing angles at v
    cosalpha = getCos(prodscalv2vv3,normavv2,normavv3);
    senalpha = getSin(cosalpha);

    cosbeta  = getCos(prodscalv1vv3,normavv1,normavv3);
    senbeta = getSin(cosbeta);

    cosgamma = getCos(prodscalv1vv2,normavv1,normavv2);
    sengamma = getSin(cosgamma);

    A = getDihedralAngle(cosalpha,cosbeta,cosgamma,senbeta,sengamma);
    B = getDihedralAngle(cosbeta,cosalpha,cosgamma,senalpha,sengamma);
    C = getDihedralAngle(cosgamma,cosalpha,cosbeta,senalpha,senbeta);

    return( A + B + C - PI );
}

double Geometry_Distortion::getDihedralAngle(double cos1, double cos2, double cos3, double sen1, double sen2)
{
    double ang = ( cos1 - ( cos2  * cos3 ) ) / ( sen1  * sen2 );
    ang = acos(ang);
    return ang;
}

double Geometry_Distortion::getSin(double cos)
{
    double x = 1.0 - ( cos * cos );
    double sin = sqrt(x);
    return sin;
}

double Geometry_Distortion::getCos(double prod12, double norm1, double norm2)
{
    double cos = prod12 / ( norm1 * norm2 );
    return cos;
}
