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

#ifndef _KDDIVISION_H
#define	_KDDIVISION_H

#include "subdivision.h"

///An class, that extends Subdivision, defining the KD subdivision
class KD_Subdivision : public Subdivision
{
public:
    ///A constructor method
    KD_Subdivision() {}
    ///A copy-constructor method
    KD_Subdivision(const KD_Subdivision& orig) : Subdivision(orig){}
    ///A destructor method
    virtual ~KD_Subdivision() {}

    ///Public method that returns the number of son nodes, in this case the costant value 2
    /*!
     * \return an integer value representing the number of sons
     */
    inline int son_number() { return 2; }
    ///Public method that computes the box domain of a son node.
    /*!
     * In this case the domain is computed considering the parent node level.
     * The subdivision splits alternately the three coordinate axes.
     *
     * \param parent_dom a Box& argument, representing the node domain
     * \param level an integer argument representing the level of n in the hierarchy
     * \param child_ind a integer argument, representing the son index position in the sub-tree
     * \return a Box value, representing the domain of the son node
     */
    Box compute_domain(Box& parent_dom, int level, int child_ind);
private:

};

#endif	/* _KDDIVISION_H */

