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

#ifndef _SUBDIVISION_H
#define	_SUBDIVISION_H

#include "basic_types/box.h"

///A super-class, not instantiable, representing a generic spatial subdivision of a tree
class Subdivision
{
public:
    ///A destructor method
    virtual ~Subdivision() {}
    ///Public pure virtual method that returns the number of son nodes
    /*!
     * \return an integer value representing the number of sons
     */
    virtual int son_number()=0;
    ///Public pure virtual method that computes the box domain of a node
    /*!
     * \param parent_dom a Box& argument, representing the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param child_ind a integer argument, representing the index position of the son in the sub-tree
     * \return a Box value, representing the domain of the son
     */
    virtual Box compute_domain(Box& parent_dom, int level, int child_ind)=0;

protected:
    ///A constructor method
    Subdivision() {}
    ///A copy-constructor method
    Subdivision(const Subdivision&) {}

private:

};

#endif	/* _SUBDIVISION_H */

