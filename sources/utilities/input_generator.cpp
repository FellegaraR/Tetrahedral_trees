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

#include "io/writer.h"
#include <cstdlib>
#include <sstream>
#include "input_generator.h"

void Input_Generator::generate_random_point_inputs(Box &region, unsigned num_entries, string output)
{
    set<Point> points;
    Point p;
    while(points.size()<num_entries)
    {
        Input_Generator::generate_random_point(p,region);
        pair<set<Point>::iterator,bool> ret = points.insert(p);
        if(ret.second)
        {
            cout<<points.size()<<" "<<flush;
        }

    }
    cout<<endl;
    Writer::write_point_queries(points,output + "_point.pqin");
}

void Input_Generator::generate_near_point_inputs(Box& region, unsigned num_entries, Mesh &mesh, string output)
{
    srand((unsigned)time(NULL));
    set<Point> points;
    set<int> t_ids;
    Point centroid = Point();
    int rand_t_id;
    while(points.size()<num_entries)
    {
        rand_t_id = rand() % mesh.get_num_tetrahedra();
        pair<set<int>::iterator,bool> ret1 = t_ids.insert(rand_t_id);
        if(ret1.second)
        {
            Geometry_Wrapper::get_tetrahedron_centroid(rand_t_id,centroid,mesh);
            if(!region.contains_with_all_closed_faces(centroid))
                continue;
            pair<set<Point>::iterator,bool> ret = points.insert(centroid);
            if(ret.second)
            {
                cout<<points.size()<<" "<<flush;
            }
        }
    }
    cout<<endl;
    Writer::write_point_queries(points,output + "_point.pqin");
}

void Input_Generator::generate_random_point(Point &p, Box& region)
{
    srand(time(NULL));
    double c;
    for(int i=0; i<3; i++)
    {
        c = (double)(((double)rand())/((double)RAND_MAX)*(region.get_max().get_c(i)-region.get_min().get_c(i))+region.get_min().get_c(i));
        p.set_c(i,c);
    }
}

void Input_Generator::generate_random_versor(Point &p)
{
    srand(time(NULL));
    double c;
    for(int i=0; i<3; i++)
    {
        c = (double)(((double)rand())/((double)RAND_MAX));
        p.set_c(i,c);
    }
}

void Input_Generator::generate_random_box_inputs(Box& region, double ratio, unsigned num_entries, string output)
{
    double edge = region.get_diagonal()*ratio;

    set<Box> boxes;
    stringstream ost;
    ost << ratio;

    Input_Generator::generate_random_boxes(region,boxes,num_entries,edge);

    Writer::write_box_queries(boxes,output + "_box_"+ost.str()+".bqin");
}

void Input_Generator::generate_near_box_inputs(Box& region, double ratio, unsigned num_entries, Mesh& mesh, string output)
{
    double edge = region.get_diagonal()*ratio;

    set<Box> boxes;
    stringstream ost;
    ost << ratio;

    Input_Generator::generate_near_boxes(region,boxes,num_entries,edge,mesh);

    Writer::write_box_queries(boxes,output + "_box_"+ost.str()+".bqin");
}

void Input_Generator::generate_random_line_inputs(Box& region, double ratio, unsigned num_entries, string output)
{
    double edge = region.get_diagonal()*ratio;

    set<Box> boxes;
    stringstream ost;
    ost << ratio;

    Input_Generator::generate_random_boxes(region,boxes,num_entries,edge);

    Writer::write_box_queries(boxes,output + "_line_"+ost.str()+".lqin");
}

void Input_Generator::generate_near_line_inputs(Box& region, double ratio, unsigned num_entries, Mesh& mesh, string output)
{
    double edge = region.get_diagonal()*ratio;

    set<Box> boxes;
    stringstream ost;
    ost << ratio;

    Input_Generator::generate_near_boxes(region,boxes,num_entries,edge,mesh);

    Writer::write_box_queries(boxes,output + "_line_"+ost.str()+".lqin");
}

void Input_Generator::generate_random_boxes(Box &region, set<Box> &boxes, unsigned num_entries, double edge)
{
    Point min;
    Point max;

    while(boxes.size()<num_entries)
    {
        Input_Generator::generate_random_point(min,region);
        max.set((min.get_x() + edge),(min.get_y() + edge),(min.get_z() + edge));

        //check that max point is into the domain
        if(!region.contains_with_all_closed_faces(max))
            continue;

        Box b = Box(min,max);
        pair<set<Box>::iterator,bool> ret = boxes.insert(b);
        if(ret.second)
        {
            cout<<boxes.size()<<" "<<flush;
        }
    }
    cout<<endl;
}

void Input_Generator::generate_near_boxes(Box &region, set<Box> &boxes, unsigned num_entries, double edge, Mesh& mesh)
{
    srand((unsigned)time(NULL));
    set<int> t_ids;
    Point centroid = Point();
    Point max;
    int rand_t_id;

    while(boxes.size()<num_entries)
    {
        rand_t_id = rand() % mesh.get_num_tetrahedra();
        pair<set<int>::iterator,bool> ret1 = t_ids.insert(rand_t_id);
        if(ret1.second)
        {
            Geometry_Wrapper::get_tetrahedron_centroid(rand_t_id,centroid,mesh);
            max.set((centroid.get_x() + edge),(centroid.get_y() + edge),(centroid.get_z() + edge));
            //check that max point is into the domain
            if(!region.contains_with_all_closed_faces(max))
                continue;

            Box b = Box(centroid,max);
            pair<set<Box>::iterator,bool> ret = boxes.insert(b);
            if(ret.second)
            {
                cout<<boxes.size()<<" "<<flush;
            }
        }
    }
    cout<<endl;
}

void Input_Generator::generate_random_lines(Box &region, set<Box> &boxes, unsigned num_entries, double edge)
{
    Point min;
    Point versor;
    Point max;

    while(boxes.size()<num_entries)
    {
        Input_Generator::generate_random_point(min,region);
        Input_Generator::generate_random_versor(versor);

        max.set((min.get_x() + versor.get_x()*edge),(min.get_y() + versor.get_y()*edge),(min.get_z() + versor.get_z()*edge));

        //check that max point is into the domain
        if(!region.contains_with_all_closed_faces(max))
            continue;

        Box b = Box(min,max);
        pair<set<Box>::iterator,bool> ret = boxes.insert(b);
        if(ret.second)
        {
            cout<<boxes.size()<<" "<<flush;
        }
    }
    cout<<endl;
}

void Input_Generator::generate_near_lines(Box &region, set<Box> &boxes, unsigned num_entries, double edge, Mesh& mesh)
{
    srand((unsigned)time(NULL));
    set<int> t_ids;
    Point centroid = Point();
    Point versor;
    Point max;
    int rand_t_id;

    while(boxes.size()<num_entries)
    {
        rand_t_id = rand() % mesh.get_num_tetrahedra();
        pair<set<int>::iterator,bool> ret1 = t_ids.insert(rand_t_id);
        if(ret1.second)
        {
            Geometry_Wrapper::get_tetrahedron_centroid(rand_t_id,centroid,mesh);
            Input_Generator::generate_random_versor(versor);

            max.set((centroid.get_x() + versor.get_x()*edge),(centroid.get_y() + versor.get_y()*edge),(centroid.get_z() + versor.get_z()*edge));
            //check that max point is into the domain
            if(!region.contains_with_all_closed_faces(max))
                continue;

            Box b = Box(centroid,max);
            pair<set<Box>::iterator,bool> ret = boxes.insert(b);
            if(ret.second)
            {
                cout<<boxes.size()<<" "<<flush;
            }
        }
    }
    cout<<endl;
}
