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

#ifndef SORTING_STRUCTURE_H
#define SORTING_STRUCTURE_H

#include <algorithm>

///A container used to store couple of vertex and tetrahedron indexes
struct vertex_tetrahedron_pair
{
    ///The vertex index
    int v;
    ///The tetrahedron index
    int t;

    inline vertex_tetrahedron_pair() { v = t = 0; }
    bool operator < (const vertex_tetrahedron_pair& p) const { return (v < p.v); }
};

/// for border checker
///A container used to store quadruplet of vertex, vertex, vertex and tetrahedron indexes
struct triangle_tetrahedron_tuple
{
    ///The first vertex index
    int v1;
    ///The second vertex index
    int v2;
    ///
    int v3;
    ///The tetrahedron index
    int t;
    /// face position in t
    short f_pos;

    triangle_tetrahedron_tuple() { v1 = v2 = v3 = t = f_pos = 0; }

    bool operator < (const triangle_tetrahedron_tuple& p) const { return (v1 < p.v1) || (v1 == p.v1 && v2 < p.v2) || (v1 == p.v1 && v2 == p.v2 && v3 < p.v3); }
    inline bool operator==(const triangle_tetrahedron_tuple& s) const { return (v1 == s.v1 && v2 == s.v2 && v3 == s.v3); }
    inline bool operator!=(const triangle_tetrahedron_tuple& s) const { return !(*this==s); }
    //
    inline void sort_and_set(int vid1, int vid2, int vid3, int tid)
    {
        int minV = min(vid1,min(vid2,vid3));
        int maxV = max(vid1,max(vid2,vid3));
        int medV=-1;

        if(minV < vid1 && maxV > vid1) medV = vid1;
        else if(minV < vid2 && maxV > vid2) medV = vid2;
        else
            medV = vid3;

        v1 = minV;
        v2 = medV;
        v3 = maxV;
        t = tid;
    }
    inline void sort_and_set(int vid1, int vid2, int vid3, int tid, short f_p)
    {
        sort_and_set(vid1, vid2, vid3, tid);
        f_pos = f_p;
    }
    inline bool has_not(int v_ind) { return (v_ind != v1 && v_ind != v2 && v_ind != v3); }
};

#endif // SORTING_STRUCTURE_H
