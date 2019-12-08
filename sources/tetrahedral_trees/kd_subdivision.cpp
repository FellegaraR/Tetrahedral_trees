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

#include "kd_subdivision.h"

Box KD_Subdivision::compute_domain(Box& parent_dom, int level, int child_ind)
{
    int coord_to_change = level % 3; // 3D

    Box sonDom = parent_dom;

    if(child_ind == 1)
    {
        sonDom.get_min().set_c(coord_to_change,parent_dom.get_min().get_c(coord_to_change)+(parent_dom.get_max().get_c(coord_to_change)-parent_dom.get_min().get_c(coord_to_change))/2.0);
    }
    else if(child_ind == 0)
    {
        sonDom.get_max().set_c(coord_to_change,parent_dom.get_min().get_c(coord_to_change)+(parent_dom.get_max().get_c(coord_to_change)-parent_dom.get_min().get_c(coord_to_change))/2.0);
    }

    return sonDom;
}
