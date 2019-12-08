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

#ifndef GEOMETRY_APPS_H
#define GEOMETRY_APPS_H

#include "basic_types/mesh.h"
#define PI 3.14159265358979323846

class Geometry_Distortion
{
public:
    /**
     * @brief A public method that computes the trihedral angle of tetrahedron t in vertex v (4D with field value)
     *
     * @param t a Tetrahedron& argument, representing the tetrahedron
     * @param v an integer representing the position index of the vertex in the boudary array of t
     * @param mesh a Mesh& argument, representing the tetrahedral mesh
     * @return the trihedral angle value
     */
    static double get_trihedral_angle(Tetrahedron& t, int v, Mesh& mesh);
    /**
     * @brief A public method that computes the trihedral angle of tetrahedron t in vertex v (3D)
     *
     * @param t a Tetrahedron& argument, representing the tetrahedron
     * @param v an integer representing the position index of the vertex in the boudary array of t
     * @param mesh a Mesh& argument, representing the tetrahedral mesh
     * @return the trihedral angle value
     */
    static double get_trihedral_angle_3D(Tetrahedron& t, int v,Mesh& mesh);

private:
    static double computeTrihedralAngle(double prodscalv1vv2, double prodscalv1vv3, double prodscalv2vv3, double normavv1, double normavv2, double normavv3);
    static double getDihedralAngle(double cos1, double cos2, double cos3, double sen1, double sen2);
    static double getSin(double cos);
    static double getCos(double prod12, double norm1, double norm2);
};

#endif // GEOMETRY_APPS_H
