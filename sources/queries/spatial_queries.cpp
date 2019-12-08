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

#include "spatial_queries.h"

bool Spatial_Queries::atomic_point_in_tetra_test(int tet_id, Point& p, QueryStatistics& qS, Mesh& mesh)
{
    qS.numGeometricTest++;
    if(Geometry_Wrapper::point_in_tetra(tet_id,p,mesh))
    {
        qS.tetrahedra.push_back(tet_id);
        return true;
    }
    return false;
}

void Spatial_Queries::atomic_tetra_in_box_test(int tet_id, Box &b, QueryStatistics& qS, Mesh& mesh, bool get_stats)
{
    if(get_stats)
        qS.access_per_tetra[tet_id]++;

    if(!qS.checkTetra[tet_id])
    {
        qS.checkTetra[tet_id]=true;

        if(get_stats)
            qS.numGeometricTest++;

        if (Geometry_Wrapper::tetra_in_box(tet_id,b,mesh))
        {
            qS.tetrahedra.push_back(tet_id);
        }
    }
}

void Spatial_Queries::atomic_line_in_tetra_test(int tet_id, Box &b, QueryStatistics& qS, Mesh& mesh, bool get_stats)
{
    if(get_stats)
        qS.access_per_tetra[tet_id]++;

    if(!qS.checkTetra[tet_id])
    {
        qS.checkTetra[tet_id]=true;

        if(get_stats)
            qS.numGeometricTest++;

        if (Geometry_Wrapper::line_in_tetra(b.get_min(),b.get_max(),tet_id,mesh))
        {
            qS.tetrahedra.push_back(tet_id);
        }
    }
}
