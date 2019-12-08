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

#include "ok_subdivision.h"
#include <boost/dynamic_bitset.hpp>

Box OK_Subdivision::compute_domain(Box& parent_dom, int, int child_ind)
{
//    cout<<"p_dom: "<<parent_dom<<endl;
//    cout<<"child_ind: "<<child_ind<<endl;
    boost::dynamic_bitset<> son_id_bits(3,child_ind); // 3D
//    cout<<"son_id_bits: "<<son_id_bits<<endl;
    Box sonDom_tmp = Box();

//    Point &min = sonDom_tmp.get_min();
//    Point &max = sonDom_tmp.get_max();

    Point &p_min = parent_dom.get_min();
    Point &p_max = parent_dom.get_max();

//    for(int i=0; i<min.get_dimension(); i++) // 3D
//    {
//        if(son_id_bits[i])
//        {
//            min.set_c(i,p_min.get_c(i)+(p_max.get_c(i)-p_min.get_c(i))/2.0);
//            max.set_c(i,p_max.get_c(i));
//        }
//        else
//        {
//            min.set_c(i,p_min.get_c(i));
//            max.set_c(i,p_min.get_c(i)+(p_max.get_c(i)-p_min.get_c(i))/2.0);
//        }
//    }
//    cout<<sonDom_tmp<<endl;

//    return sonDom;
    Box sonDom = Box();
    double xmin, ymin, zmin, xmax, ymax, zmax;

    /// this method has been pushed-back from 2016 implementation
    /// as we were visiting differently some of the children (o and 3, and 4 and 7)
    /// (the commented implementation above is the one used in Stellar tree and Terrain trees..)
    /// we keep this to be coherent with the experiments we executed back then in 2016..
    /// as different visit ordering lead to different reindexing
    /// and, later, different geom-tests executed for a query
    //setto le x minime e massime della box
    if(child_ind == 0 || child_ind == 1 || child_ind == 4 || child_ind == 5)
    {
        xmin = p_min.get_x()+(p_max.get_x()-p_min.get_x())/2.0;
        xmax = p_max.get_x();
    }
    else // child_ind = 2 / 3 / 6 / 7
    {
        xmin = p_min.get_x();
        xmax = p_min.get_x()+(p_max.get_x()-p_min.get_x())/2.0;
    }

    //setto le y minime e massime della box
    if(child_ind == 0 || child_ind == 2 || child_ind == 4 || child_ind == 6)
    {
        ymin = p_min.get_y()+(p_max.get_y()-p_min.get_y())/2.0;
        ymax = p_max.get_y();
    }
    else // child_ind = 1 / 3 / 5 / 7
    {
        ymin = p_min.get_y();
        ymax = p_min.get_y()+(p_max.get_y()-p_min.get_y())/2.0;
    }

    //setto le z minime e massime della box
    if(child_ind == 0 || child_ind == 1 || child_ind == 2 || child_ind == 3)
    {
        zmin = p_min.get_z();
        zmax = p_min.get_z()+(p_max.get_z()-p_min.get_z())/2.0;
    }
    else // child_ind = 4 / 5 / 6 / 7
    {
        zmin = p_min.get_z()+(p_max.get_z()-p_min.get_z())/2.0;
        zmax = p_max.get_z();
    }

    sonDom.set_min(xmin,ymin,zmin);
    sonDom.set_max(xmax,ymax,zmax);

//    cout<<sonDom<<endl;

//    int a; cin>>a;
    return sonDom;
}
