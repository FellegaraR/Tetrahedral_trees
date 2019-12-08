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

#ifndef _NODE_H
#define	_NODE_H

#include <vector>
#include <cstddef>
#include <set>
#include <bm/bm.h>
#include "basic_types/box.h"
#include "basic_types/mesh.h"
#include "run_iterator.h"

/**
 * @brief A super-class, not instantiable, that represents a generic node of the tree with associated an array of tetrahedra
 */
template<class N> class Node
{
public:
    ///A destructor method
    virtual ~Node() {}
    ///A public method that sets a node son
    /*!
     * \param n a N* argument, represents the son to set
     * \param pos an integer, represents the son position in the array
     */
    inline void set_son(N* n, int pos) { this->sons[pos] = n; }
    ///A public method that checks if the node is a leaf node
    /*!
     * \return a boolean, true if the node is a leaf, false otherwise
     */
    inline bool is_leaf() const { return !(this->sons!=NULL); }
    ///A public method that returns a node son
    /*!
     * \param i an integer argument, represents the son position into the list
     * \return a Node*, representing the son at the i-th position
     */
    inline N* get_son(int i) { return this->sons[i]; }
    /**
     * @brief A public method that initializes the sons array
     * @param son_number an integer containing the sons number
     */
    inline void init_sons(int son_number) { this->sons = new N*[son_number]; }
    ///A public method that adds a tetrahedron index to the array
    /*!
     * \param ind an integer argument, representing the tetrahedron index
     */
    inline void add_tetrahedron(int ind) { this->tetrahedra.push_back(ind); }    

    ///A public method that returns the run_iterator pair to navigate the tetrahedra array
    inline RunIteratorPair make_t_array_iterator_pair() { return run_iterator<int>::make_run_iterator_pair(tetrahedra); }
    ///A public method that returns the begin run_iterator to navigate the tetrahedra array
    inline RunIterator t_array_begin_iterator() { return run_iterator<int>(tetrahedra.begin(),tetrahedra.end()); }
    ///A public method that returns the end run_iterator to navigate the tetrahedra array
    inline RunIterator t_array_end_iterator() { return run_iterator<int>(tetrahedra.end()); }

    /**
     * @brief A public method that returns the number of indexed tetrahedra
     * NOTA: this method expand the runs and returns the real number of top d-cells indexed in the node
     *
     * @return int
     */
    inline int get_real_t_array_size() const { return run_iterator<int>(tetrahedra.begin(),tetrahedra.end()).elementCountFast(tetrahedra); }
    /**
     * @brief A public method that returns the size of the tetrahedral array
     *
     * @return int
     */
    inline int get_t_array_size() { return this->tetrahedra.size(); }
    /**
     * @brief A public method returning the tetrahedra array
     *
     * @return
     */
    inline int_vect get_t_array() const { return this->tetrahedra; }
    /**
     * @brief A public method that clears the space used by the tetrahedra array
     */
    inline void clear_t_array() { tetrahedra.clear(); }
    ///A public method that return the begin iterator of the tetrahedra array for explicitly unroll the runs of tetrahedra
    inline int_vect_iter get_t_array_begin() { return this->tetrahedra.begin(); }
    ///A public method that return the end iterator of the tetrahedra array for explicitly unroll the runs of tetrahedra
    inline int_vect_iter get_t_array_end() { return this->tetrahedra.end(); }

    // geometric procedures //
    /**
     * @brief A public method that checks if all the four vertices of a tetrahedron are indexed by the node
     *
     * @param t a Tetrahedron& argument
     * @param domain a Box& representing the node domain
     * @param mesh a Mesh& representing the tetrahedral mesh
     * @return true if all the four vertices are indexed, false otherwise
     */
    inline bool completely_indexes_tetrahedron_vertices_dom(Tetrahedron &t, Box& domain, Mesh& mesh)
    {
        for(int v=0; v<t.vertices_num(); v++)
            if(!domain.contains(mesh.get_vertex(t.TV(v)),mesh.get_domain().get_max()))
                return false;
        return true;
    }
    /**
     * @brief A public method that checks if at least one vertex of a tetrahedron is indexed by the node
     *
     * @param t a Tetrahedron& argument
     * @param domain a Box& representing the node domain
     * @param mesh a Mesh& representing the tetrahedral mesh
     * @return true if at least one vertex is indexed, false otherwise
     */
    inline bool indexes_tetrahedron_vertices_dom(Tetrahedron &t, Box& domain, Mesh& mesh)
    {
        for(int v=0; v<t.vertices_num(); v++)
            if(domain.contains(mesh.get_vertex(t.TV(v)),mesh.get_domain().get_max()))
                return true;
        return false;
    }
    /**
     * @brief A public method that computes the bounding box of a run
     * @param id an iterator to the current array entry
     * @param bb a Box& argument, that is set with the run bounding box (if a run is found)
     * @param mesh a Mesh& representing the tetrahedral mesh
     * @param run a pair that will contains the run, if any
     * @return true if a run has been encounter, false otherwise
     */
    bool get_run_bounding_box(int_vect_iter &id, Box& bb, Mesh &mesh, pair<int,int> &run);

protected:    
    ///A constructor method
    Node() { this->sons = NULL; }
    ///A copy-constructor method
    Node(const Node& orig)
    {
        this->sons = orig.sons;
        this->tetrahedra = orig.tetrahedra;
    }
    ///A protected variable representing the list of node sons
    N** sons;
    ///A private variable representing the list containing the tetrahedra indexed by the node
    int_vect tetrahedra;
};

template<class N> bool Node<N>::get_run_bounding_box(int_vect_iter &id, Box& bb, Mesh &mesh, pair<int,int> &run)
{
    if(*id<0) //I have a run
    {
        run.first = abs(*id);
        ++id;
        run.second = run.first + *id;

        double min_p[3]={std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max()};
        double max_p[3]={-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max()};

//        // init phase
//        int t_id=run.first;
//        Tetrahedron &t_first = mesh.get_tetrahedron(t_id);
//        min_p[0] = mesh.get_vertex(t_first.TV(0)).get_x();
//        min_p[1] = mesh.get_vertex(t_first.TV(0)).get_y();
//        min_p[2] = mesh.get_vertex(t_first.TV(0)).get_z();
//        max_p[0] = mesh.get_vertex(t_first.TV(0)).get_x();
//        max_p[1] = mesh.get_vertex(t_first.TV(0)).get_y();
//        max_p[2] = mesh.get_vertex(t_first.TV(0)).get_z();

//        for(int i=1; i<t_first.vertices_num(); i++)
//        {
//            Vertex &v = mesh.get_vertex(t_first.TV(i));
//            for(int j=0;j<v.get_dimension();j++)
//            {
//                if(v.get_c(j) < min_p[j])
//                    min_p[j] = v.get_c(j);
//                else if(v.get_c(j) > max_p[j])
//                    max_p[j] = v.get_c(j);
//            }
//        }
//        t_id++;
        for(int t_id=run.first; t_id<=run.second; t_id++)
        {
            Tetrahedron &tet = mesh.get_tetrahedron(t_id);
            for(int i=0; i<tet.vertices_num(); i++)
            {
                Vertex &v = mesh.get_vertex(tet.TV(i));
                for(int j=0;j<v.get_dimension();j++)
                {
                    if(v.get_c(j) < min_p[j])
                        min_p[j] = v.get_c(j);
                    if(v.get_c(j) > max_p[j])
                        max_p[j] = v.get_c(j);
                }
            }
        }
        //save the computed bounding box
        bb.set_min(min_p[0],min_p[1],min_p[2]);
        bb.set_max(max_p[0],max_p[1],max_p[2]);
        return true;
    }
    else
    {
        return false;
    }
}

#endif	/* _NODE_H */

