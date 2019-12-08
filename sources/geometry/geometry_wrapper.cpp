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

#include "geometry_wrapper.h"
#include <boost/dynamic_bitset.hpp>

void Geometry_Wrapper::get_tetrahedron_centroid(int t_id, Point& p, Mesh &mesh)
{
    Tetrahedron &tet = mesh.get_tetrahedron(t_id);

    Vertex& v0 = mesh.get_vertex(tet.TV(0));
    Vertex& v1 = mesh.get_vertex(tet.TV(1));
    Vertex& v2 = mesh.get_vertex(tet.TV(2));
    Vertex& v3 = mesh.get_vertex(tet.TV(3));

    for(int i=0; i<3; i++)
        p.set_c(i,(v0.get_c(i) + v1.get_c(i) + v2.get_c(i) + v3.get_c(i)) / 4.0);
}

bool Geometry_Wrapper::point_in_tetra(int t_id, Point& point, Mesh &mesh)
{
    Tetrahedron &tet = mesh.get_tetrahedron(t_id);
    double **c;
    c = new double* [4];
    for (int i = 0; i < tet.vertices_num(); i++)
    {
        c[i] = new double[3];
        Vertex &v = mesh.get_vertex(tet.TV(i));
        c[i][0] = v.get_x();
        c[i][1] = v.get_y();
        c[i][2] = v.get_z();
    }

    bool ret = false;

    if(PointInTetra(point.get_x(),point.get_y(),point.get_z(),c[0],c[1],c[2],c[3]))
        ret = true;

    for (int i = 0; i < 4; i++)
        delete c[i];
    delete c;

    return ret;
}

bool Geometry_Wrapper::tetra_in_box_build(int t_id, Box& box, Mesh& mesh)
{
    Tetrahedron &t = mesh.get_tetrahedron(t_id);

    //I consider internal a tetrahedron that at least has a vertex inside the box (also only a vertex)
    for(int v=0; v<t.vertices_num(); v++)
    {
        if(box.contains(mesh.get_vertex(t.TV(v)),mesh.get_domain().get_max()))
            return true;
    }

    return Geometry_Wrapper::tetra_in_box(t_id,box,mesh);
}

bool Geometry_Wrapper::tetra_in_box(int t_id, Box& box, Mesh& mesh)
{
    Tetrahedron &t = mesh.get_tetrahedron(t_id);

    // The library used only supports doubles, not doubles
    double* minf = new double[3];
    minf[0] = box.get_min().get_x(); minf[1] = box.get_min().get_y(); minf[2] = box.get_min().get_z();
    double* maxf = new double[3];
    maxf[0] = box.get_max().get_x(); maxf[1] = box.get_max().get_y(); maxf[2] = box.get_max().get_z();

    double **c;
    c = new double* [4];
    for (int i = 0; i < t.vertices_num(); i++)
    {
        c[i] = new double[3];
        Vertex& v = mesh.get_vertex(t.TV(i));

        c[i][0] = v.get_x();
        c[i][1] = v.get_y();
        c[i][2] = v.get_z();
    }

    bool ret = false;

    if(tetra_in_box_strict(minf, maxf, c))
        ret = true;

    //DEALLOCAZIONE
    delete minf;
    delete maxf;
    for (int i = 0; i < t.vertices_num(); i++) {
        delete c[i];
    }
    delete c;
    //DEALLOCAZIONE

    return ret;
}

bool Geometry_Wrapper::line_in_box(const Point& v1, const Point& v2, Box& box)
{
    return ClipLine3D_middle(box.get_min().get_x(),box.get_min().get_y(),box.get_min().get_z(),
                         box.get_max().get_x(),box.get_max().get_y(),box.get_max().get_z(),
                         v1.get_x(),v1.get_y(),v1.get_z(),
                         v2.get_x(),v2.get_y(),v2.get_z());
}

bool Geometry_Wrapper::line_in_bounding_box(const Point &v1, const Point &v2, Box& bb)
{
    return ClipLine3D(bb.get_min().get_x(),bb.get_min().get_y(),bb.get_min().get_z(),
                             bb.get_max().get_x(),bb.get_max().get_y(),bb.get_max().get_z(),
                             v1.get_x(),v1.get_y(),v1.get_z(),
                             v2.get_x(),v2.get_y(),v2.get_z());
}

bool Geometry_Wrapper::line_in_tetra(const Point& v1, const Point& v2, int t_id, Mesh &mesh)
{
    Tetrahedron &tet = mesh.get_tetrahedron(t_id);
    Point d = v2 - v1;
    double tfirst = 0.0;
    double tlast = 1.0;

    vector<int> f;
    f.assign(3,0);

    for(int i=0; i<tet.vertices_num(); i++)
    {
        Geometry_Wrapper::ordered_TF(tet,i,f);
        Point sub_ba = mesh.get_vertex(f[1]) - mesh.get_vertex(f[0]);
        Point sub_ca = mesh.get_vertex(f[2]) - mesh.get_vertex(f[0]);
        Point n = sub_ba.cross_3D(sub_ca);

        Point sub_v1a = v1 - mesh.get_vertex(f[0]);
        double N = - sub_v1a.dot_3D(n);
        double D = d.dot_3D(n);

        if(D==0)// then S is parallel to the current face
        {
            if(N<0) //then v1 is outside the current face
                return false;
        }
        else
        {
            double t = N / D;

            if (D < 0) //then segment S is entering OMEGA across face Fi
            {
                tfirst = max(tfirst,t);
                if (tfirst > tlast) //then segment S enters OMEGA after leaving
                        return false; //since S cannot intersect OMEGA
            }
            else if (D > 0) //then segment S is leaving OMEGA across face Fi
            {
                tlast = min(tlast,t);
                if (tlast < tfirst) //then segment S leaves OMEGA before entering
                    return false; //since S cannot intersect OMEGA
            }
        }
    }
    return true;
}

void Geometry_Wrapper::ordered_TF(Tetrahedron &t, int pos, vector<int> &f)
{
    switch(pos)
    {
    case(0):
        f[0] = t.TV(0);
        f[1] = t.TV(1);
        f[2] = t.TV(2);
        break;
    case(1):
        f[0] = t.TV(1);
        f[1] = t.TV(3);
        f[2] = t.TV(2);
        break;
    case(2):
        f[0] = t.TV(3);
        f[1] = t.TV(0);
        f[2] = t.TV(2);
        break;
    case(3):
        f[0] = t.TV(1);
        f[1] = t.TV(0);
        f[2] = t.TV(3);
        break;
    default:
        cerr<<"[ordered_TF] wrong face position. must be from 0 to 3."<<endl;
        int a; cin>>a;
        break;
    }
}

void Geometry_Wrapper::set_faces_ordering(Mesh &mesh)
{
    for(int i=1; i<=mesh.get_num_tetrahedra(); i++)
        Geometry_Wrapper::set_face_orientation(mesh.get_tetrahedron(i),mesh);
}
void Geometry_Wrapper::set_face_orientation(Tetrahedron &tet, Mesh &mesh)
{
    int turn;
    int new_0, new_1, new_2, new_3;

    turn = Geometry_Wrapper::four_point_turn_wrapper(mesh.get_vertex(tet.TV(0)), mesh.get_vertex(tet.TV(1)), mesh.get_vertex(tet.TV(2)), mesh.get_vertex(tet.TV(3)));
    if(turn == RIGHT_TURN)
    {
        //nothing to do, already coherent
        return;
    }

    turn = Geometry_Wrapper::four_point_turn_wrapper(mesh.get_vertex(tet.TV(1)), mesh.get_vertex(tet.TV(0)), mesh.get_vertex(tet.TV(2)), mesh.get_vertex(tet.TV(3)));
    if(turn == RIGHT_TURN)
    {
        new_0 = tet.TV(1);
        new_1 = tet.TV(0);
        new_2 = tet.TV(2);
        new_3 = tet.TV(3);
        tet.set(new_0,new_1,new_2,new_3);
        return;
    }

    turn = Geometry_Wrapper::four_point_turn_wrapper(mesh.get_vertex(tet.TV(2)), mesh.get_vertex(tet.TV(1)), mesh.get_vertex(tet.TV(0)), mesh.get_vertex(tet.TV(3)));
    if(turn == RIGHT_TURN)
    {
        new_0 = tet.TV(2);
        new_1 = tet.TV(1);
        new_2 = tet.TV(0);
        new_3 = tet.TV(3);
        tet.set(new_0,new_1,new_2,new_3);
        return;
    }

    cout<<"not found a face at the right of the first one"<<endl;
    cout<<tet<<endl;
    cout<<mesh.get_vertex(tet.TV(0))<<endl;
    cout<<mesh.get_vertex(tet.TV(1))<<endl;
    cout<<mesh.get_vertex(tet.TV(2))<<endl;
    cout<<mesh.get_vertex(tet.TV(3))<<endl;
    cout<<Geometry_Wrapper::four_point_turn_wrapper(mesh.get_vertex(tet.TV(0)), mesh.get_vertex(tet.TV(1)), mesh.get_vertex(tet.TV(2)), mesh.get_vertex(tet.TV(3)))<<" ";
    cout<<Geometry_Wrapper::four_point_turn_wrapper(mesh.get_vertex(tet.TV(1)), mesh.get_vertex(tet.TV(0)), mesh.get_vertex(tet.TV(2)), mesh.get_vertex(tet.TV(3)))<<" ";
    cout<<Geometry_Wrapper::four_point_turn_wrapper(mesh.get_vertex(tet.TV(2)), mesh.get_vertex(tet.TV(1)), mesh.get_vertex(tet.TV(0)), mesh.get_vertex(tet.TV(3)))<<endl;
}

int Geometry_Wrapper::four_point_turn_wrapper(const Point &v0, const Point &v1, const Point &v2, const Point &op)
{
    return FourPointTurn(op.get_x(), op.get_y(), op.get_z(), v0.get_x(), v0.get_y(), v0.get_z(), v1.get_x(), v1.get_y(), v1.get_z(), v2.get_x(), v2.get_y(), v2.get_z());
}
