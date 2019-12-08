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


#ifndef BASIC_STRUCTURE
#define BASIC_STRUCTURE

#include <vector>
#include <set>
#include <queue>
#include <map>

typedef std::queue<int> int_queue;
typedef std::vector<int> int_vect;
typedef int_vect::iterator int_vect_iter;
typedef int_vect::const_iterator int_vect_const_iter;

typedef std::set<int> int_set;
typedef int_set::iterator int_set_iter;
typedef int_set::const_iterator int_set_const_iter;

typedef int_vect VT;
typedef int_set VV;
typedef int_vect ET;

typedef std::vector<VT> leaf_VT;
typedef std::vector<VV> leaf_VV;
typedef std::map<int_vect,ET> leaf_ET;

#endif // BASIC_STRUCTURE

