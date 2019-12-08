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

#ifndef SPATIAL_QUERIES_H
#define SPATIAL_QUERIES_H

#include "statistics/statistics.h"
#include "utilities/timer.h"

/**
 * @brief The Spatial_Queries class provides an interface for executing spatial queries on the Tetrahedral trees
 */
class Spatial_Queries
{
public:
    Spatial_Queries() {}

    ///A public method that excutes point locations, reading the points from file
    /*!
     * This method prints the results on standard output
     *
     * \param tree a T& argument, represents the tree where the statistics are executed
     * \param query_path a string argument, representing the file path of the query input
     * \param stats a Statistics& argument, representing the object for computing the associated statistics
     */
    template<class T> void exec_point_locations(T& tree, string query_path, Statistics &stats);
    ///A public method that excutes box queries, reading the boxes from file
    /*!
     * This method prints the results on standard output
     *
     * \param tree a T& argument, represents the tree where the statistics are executed
     * \param query_path a string argument, representing the file path of the query input
     * \param stats a Statistics& argument, representing the object for computing the associated statistics
     */
    template<class T> void exec_box_queries(T& tree, string query_path, Statistics &stats);
    ///A public method that excutes line queries, reading the lines from file
    /*!
     * This method prints the results on standard output
     *
     * \param tree a T& argument, represents the tree where the statistics are executed
     * \param query_path a string argument, representing the file path of the query input
     * \param stats a Statistics& argument, representing the object for computing the associated statistics
     */
    template<class T> void exec_line_queries(T& tree, string query_path, Statistics &stats);

private:
    ///A private method that executes a single point location on a Tetrahedral tree
    /*!
     * \param n a N& argument, representing the current node to visit
     * \param dom a Box& argument, representing the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param p a Point& argument, representing the point query
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision
     */
    template<class N, class D> void exec_point_query(N& n, Box& dom, int level, Point& p, QueryStatistics& qS, Mesh& mesh, D& division);
    ///A private method that executes a single box query on a Tetrahedral tree
    /*!
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param b a Box& argument, representing the box query
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     */
    template<class N, class D> void exec_box_query(N& n, Box& dom, int level, Box& b, QueryStatistics& qS, Mesh& mesh, D& division, bool get_stats);
    ///A private method that executes a single line query on a Tetrahedral tree
    /*!
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param b a Box& argument, representing the line query (the line is represented by the minimum and maximum points of the box)
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     */
    template<class N, class D> void exec_line_query(N& n, Box &dom, int level, Box& b, QueryStatistics& qS, Mesh& mesh, D& division, bool get_stats); // the line segment is represented by b

    ///A private method that executes a point location in a leaf block
    /*!
     * \param n a N& argument, representing the current leaf
     * \param p a Point& argument, representing the point query
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     */
    template<class N> void exec_point_query_leaf(N& n, Point& p, QueryStatistics& qS, Mesh& mesh);
    /**
     * @brief A private method executing a point-in-tetra test on a tetrahedron
     * @param tet_id an integer representing the current tetrahedron to test
     * @param p a Point& argument, representing the point query
     * @param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * @param mesh a Mesh& argument, representing the current mesh
     * @return true if an intersection exists, false otherwise
     */
    bool atomic_point_in_tetra_test(int tet_id, Point& p, QueryStatistics& qS, Mesh& mesh);
    ///A private method that executes a box query in a leaf
    /*!
     * \param n a N& argument, representing the actual leaf
     * \param b a Box& argument, representing the box query
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     * \param get_stats a boolean, true if statistics must be computed, false otherwise
     */
    template<class N> void exec_box_query_leaf_test(N& n, Box& b, QueryStatistics& qS, Mesh& mesh, bool get_stats);
    /**
     * @brief A private method executing a tetra-in-box test on a tetrahedron
     *
     * @param tet_id an integer representing the current tetrahedron to test
     * @param b a Box& argument, representing the box query
     * @param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * @param mesh a Mesh& argument, representing the current mesh
     * @param get_stats a boolean, true if statistics must be computed, false otherwise
     */
    void atomic_tetra_in_box_test(int tet_id, Box& b, QueryStatistics& qS, Mesh& mesh, bool get_stats);
    /**
     * @brief A private method that adds the tetrahedra in a leaf to the result set
     * NOTA: this procedures simply add all the tetrahedra as the domain of the leaf is completely contained by the box QueryStatistics
     *
     * @param n a N& argument, representing the actual leaf
     * @param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * @param get_stats a boolean, true if statistics must be computed, false otherwise
     */
    template<class N> void add_tetrahedra_to_box_query_result(N& n, QueryStatistics& qS, bool get_stats);

    ///A private method that executes a line query in a leaf
    /*!
     * \param n a N& argument, representing the actual leaf
     * \param b a Box& argument, representing the box query
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     * \param get_stats a boolean, true if statistics must be computed, false otherwise
     */
    template<class N> void exec_line_query_leaf(N& n, Box &b, QueryStatistics& qS, Mesh& mesh, bool get_stats);
    /**
     * @brief A private method executing a line-in-tetra test on a tetrahedron
     *
     * @param tet_id an integer representing the current tetrahedron to test
     * @param b a Box& argument, representing the line query
     * @param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * @param mesh a Mesh& argument, representing the current mesh
     * @param get_stats a boolean, true if statistics must be computed, false otherwise
     */
    void atomic_line_in_tetra_test(int tet_id, Box &b, QueryStatistics& qS, Mesh& mesh, bool get_stats);
};

template<class T> void Spatial_Queries::exec_point_locations(T& tree, string query_path, Statistics &stats)
{
    QueryStatistics qS = QueryStatistics();
    vector<Point> points;
    Reader::read_queries(points,query_path);

    Timer time;
    double tot_time = 0;
    int hit_ratio = 0;

    for(unsigned i=0;i<points.size();i++)
    {
        time.start();
        this->exec_point_query(tree.get_root(),tree.get_mesh().get_domain(),0,points[i],qS, tree.get_mesh(),tree.get_decomposition());
        time.stop();
        tot_time += time.get_elapsed_time();

        //debug print
        if(qS.tetrahedra.size()>0)
            cout<<"found tetra for point "<<i<<endl;
        else
            cout<<"nothing found for point "<<i<<endl;

        hit_ratio += stats.compute_queries_statistics(qS);
        qS.reset();
    }
    cerr<<"[TIME] exec point locations "<<tot_time<<endl;

    Writer::write_queries_stats(points.size(),stats.get_query_statistics(),hit_ratio);
    points.clear();
}

template<class T> void Spatial_Queries::exec_box_queries(T& tree, string query_path, Statistics &stats)
{
    QueryStatistics qS = QueryStatistics(tree.get_mesh().get_num_tetrahedra(),4);

    vector<Box> boxes;
    Reader::read_queries(boxes,query_path);

    Timer time;
    double tot_time = 0;
    int hit_ratio = 0;

    for(unsigned j=0;j<boxes.size();j++)
    {
//        cout<<"B: "<<boxes[j]<<endl;
        // exec for timings
        time.start();
        this->exec_box_query(tree.get_root(),tree.get_mesh().get_domain(),0,boxes[j],qS, tree.get_mesh(),tree.get_decomposition(),false);
        time.stop();
        tot_time += time.get_elapsed_time();

        // exec again for stats
        qS.reset(false);
        this->exec_box_query(tree.get_root(),tree.get_mesh().get_domain(),0,boxes[j],qS, tree.get_mesh(),tree.get_decomposition(),true);

        //debug print
        cout<<qS.tetrahedra.size()<<" intersect box "<<j<<endl;
//        int a; cin>>a;

        hit_ratio += stats.compute_queries_statistics(qS);
        qS.reset(true);
    }
    cerr<<"[TIME] exec box queries "<<tot_time<<endl;

    Writer::write_queries_stats(boxes.size(),stats.get_query_statistics(),hit_ratio);
    boxes.clear();
}

template<class T> void Spatial_Queries::exec_line_queries(T& tree, string query_path, Statistics &stats)
{
    QueryStatistics qS = QueryStatistics(tree.get_mesh().get_num_tetrahedra(),8);

    vector<Box> boxes;
    Reader::read_queries(boxes,query_path);

    Timer time;
    double tot_time = 0;
    int hit_ratio = 0;

    for(unsigned j=0;j<boxes.size();j++)
    {
        // exec for timings
        time.start();
        this->exec_line_query(tree.get_root(),tree.get_mesh().get_domain(),0,boxes[j],/*line_length,*/qS, tree.get_mesh(),tree.get_decomposition(),false);
        std::sort(qS.tetrahedra.begin(),qS.tetrahedra.end());
        vector<int>::iterator last_pos = std::unique(qS.tetrahedra.begin(),qS.tetrahedra.end());
        qS.tetrahedra.resize(std::distance(qS.tetrahedra.begin(),last_pos));
        time.stop();
        tot_time += time.get_elapsed_time();

        qS.reset(false);
        this->exec_line_query(tree.get_root(),tree.get_mesh().get_domain(),0,boxes[j],qS, tree.get_mesh(),tree.get_decomposition(),true);
        std::sort(qS.tetrahedra.begin(),qS.tetrahedra.end());
        std::unique(qS.tetrahedra.begin(),qS.tetrahedra.end());
        qS.tetrahedra.resize(std::distance(qS.tetrahedra.begin(),last_pos));
        hit_ratio += stats.compute_queries_statistics(qS);

        //debug print
        cout<<qS.tetrahedra.size()<<" intersect line "<<j<<" "<<boxes[j]<<endl;
        qS.reset(true);
    }
    cerr<<"[TIME] exec line queries "<<tot_time<<endl;
    cerr<<"avg geom test: "<<stats.get_query_statistics().avgGeometricTest/(double)hit_ratio<<endl;

    Writer::write_queries_stats(boxes.size(),stats.get_query_statistics(),hit_ratio);
    boxes.clear();
}

template<class N, class D> void Spatial_Queries::exec_point_query(N &n, Box &dom, int level, Point &p, QueryStatistics &qS, Mesh &mesh, D &division)
{
    qS.numNode++;

    if (n.is_leaf())
    {
        qS.numLeaf++;
        this->exec_point_query_leaf(n,p,qS,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            if(son_dom.contains(p,mesh.get_domain().get_max()))
            {
                this->exec_point_query(*n.get_son(i),son_dom,son_level,p,qS,mesh,division);
                break;
            }
        }
    }
}

template<class N> void Spatial_Queries::exec_point_query_leaf(N &n, Point &p, QueryStatistics &qS, Mesh &mesh)
{
    Box bb;
    pair<int,int> run;

    for(vector<int>::iterator it=n.get_t_array_begin(); it!=n.get_t_array_end(); ++it)
    {
        if(n.get_run_bounding_box(it,bb,mesh,run))
        {
            if(bb.contains(p,mesh.get_domain().get_max()))
            {
                for(int t_id=run.first; t_id<=run.second; t_id++)
                {
                    if(atomic_point_in_tetra_test(t_id,p,qS,mesh))
                        return;
                }
            }
        }
        else
        {
            if(atomic_point_in_tetra_test(*it,p,qS,mesh))
                return;
        }
    }
}

template<class N, class D> void Spatial_Queries::exec_box_query(N &n, Box &dom, int level, Box &b, QueryStatistics &qS, Mesh &mesh, D &division, bool get_stats)
{
    if(get_stats)
        qS.numNode++;

    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
//        qS.tetrahedra.clear(); /// debug
//        cerr<<"L: "<<dom<<endl;
//        cerr<<"L-Tnum: "<<n.get_real_t_array_size()<<endl;
//        cerr<<"search-box: "<<b<<endl;
        if(get_stats)
            qS.numLeaf++;
        if(b.completely_contains(dom))
        {
//            cerr<<"box completely_contains dom"<<endl;
            if(get_stats)
                qS.box_completely_contains_leaf_num++;
            this->add_tetrahedra_to_box_query_result(n,qS,get_stats);
        }
        else
            this->exec_box_query_leaf_test(n,b,qS,mesh,get_stats);

//        ///////////////////////////// debug
//        int_vect intersecting = qS.tetrahedra;
//        qS.tetrahedra.clear();

//        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
//        {
//            RunIterator const& tet_id = itPair.first;
//            atomic_tetra_in_box_test(*tet_id,b,qS,mesh,get_stats);
//        }

//        if(intersecting.size() != qS.tetrahedra.size())
//        {
//            cerr<<"L: "<<dom<<endl;
//            cerr<<"L-Tnum: "<<n.get_real_t_array_size()<<endl;
//            cerr<<"search-box: "<<b<<endl;

//            cerr<<"found tetra (bf): "<<qS.tetrahedra.size()<<endl;
//            cerr<<"found tetra (run): "<<intersecting.size()<<endl;

//            int_vect diff(qS.tetrahedra.size());
//            int_vect_iter it = set_difference(qS.tetrahedra.begin(),qS.tetrahedra.end(),intersecting.begin(),intersecting.end(),diff.begin());
//            diff.resize(it-diff.begin());

//            cout<<"missing tetrahedra: ";
//            for (size_t i=0; i<diff.size(); i++)
//            {
//                cout<<diff[i]<<" ";
//            }
//            cout<<endl;
//            qS.tetrahedra.clear();
//            this->exec_box_query_leaf_test(n,b,qS,mesh,true);
//            int a; cin>>a;
//        }
//        ///////////////////////////// debug
//        cerr<<qS.tetrahedra.size()<<endl;
//        int a; cin>>a;
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->exec_box_query(*n.get_son(i), son_dom, son_level, b, qS, mesh,division, get_stats);
        }
    }
}

template<class N> void Spatial_Queries::exec_box_query_leaf_test(N &n, Box &b, QueryStatistics &qS, Mesh &mesh, bool get_stats)
{
    Box bb;
    pair<int,int> run;

    for(vector<int>::iterator it=n.get_t_array_begin(); it!=n.get_t_array_end(); ++it)
    {
        if(n.get_run_bounding_box(it,bb,mesh,run))
        {
//            if(get_stats)
//            {
//                cerr<<"run: "<<run.first<<" "<<run.second<<endl;
//                cerr<<"bbox: "<<bb<<endl;
//            }
//            int a; cin>>a;
            if(b.completely_contains(bb))
            {
                if(get_stats)
                    qS.box_completely_contains_bbox_num++;

//                if(get_stats)
//                    cerr<<"b completely_contains box: pushing -> "<<run.second-run.first<<endl;

                for(int t_id=run.first; t_id<=run.second; t_id++)
                {
                    if(get_stats)
                        qS.access_per_tetra[t_id]++;

                    if(!qS.checkTetra[t_id])
                    {
                        qS.checkTetra[t_id]=true;
                        qS.tetrahedra.push_back(t_id);

                        if(get_stats)
                        {
                            qS.tetra_compl_cont_bbox_num++;
                            qS.avoided_tetra_geom_tests_num++;
                        }
                    }
                }

            }
            else if(b.intersects(bb))
            {
                if(get_stats)
                    qS.box_intersect_bbox_num++;

//                if(get_stats)
//                    cerr<<"b intersects bbox: checking -> "<<run.second-run.first<<endl;

                for(int t_id=run.first; t_id<=run.second; t_id++)
                {
                    if(get_stats)
                        if(!qS.checkTetra[t_id])
                            qS.box_intersect_bbox_geom_tests_num++;
                    atomic_tetra_in_box_test(t_id,b,qS,mesh,get_stats);
                }
            }
            else if(get_stats) // bbox does not intesect the search box
            {
//                cerr<<"bbox and box do not intersect"<<endl;
                //if(get_stats)
                qS.box_no_intersect_bbox_num++;
                for(int t_id=run.first; t_id<=run.second; t_id++)
                {
                    /// == WARNING ==
                    /// computing the statistics like this (i.e., by flagging as visited also these tops)
                    /// affects the stat about the average geometric tests executed
                    ///
                    if(!qS.checkTetra[t_id] && !qS.avoid_to_check_tetra[t_id])
                    {
                        qS.avoid_to_check_tetra[t_id]=true;
                        //if(get_stats)
                        qS.avoided_tetra_geom_tests_num++;
                    }
                }
            }
        }
        else
        {
            atomic_tetra_in_box_test(*it,b,qS,mesh,get_stats);
        }
    }

}

template<class N> void Spatial_Queries::add_tetrahedra_to_box_query_result(N& n, QueryStatistics& qS, bool get_stats)
{
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& tet_id = itPair.first;
        if(get_stats)
            qS.access_per_tetra[*tet_id]++;

        if(!qS.checkTetra[*tet_id])
        {
            qS.checkTetra[*tet_id]=true;
            qS.tetrahedra.push_back(*tet_id);

            if(get_stats)
            {
                qS.tetra_compl_cont_leaf_num++;
                qS.avoided_tetra_geom_tests_num++;
            }
        }
    }
}

template<class N, class D> void Spatial_Queries::exec_line_query(N &n, Box &dom, int level, Box &b, QueryStatistics &qS, Mesh &mesh, D &division, bool get_stats)
{
    if(get_stats)
        qS.numNode++;

    if(!Geometry_Wrapper::line_in_box(b.get_min(),b.get_max(),dom))
    {
        return;
    }

    if (n.is_leaf())
    {
        if(get_stats)
            qS.numLeaf++;
        exec_line_query_leaf(n,b,qS,mesh,get_stats);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->exec_line_query(*n.get_son(i), son_dom, son_level, b, qS, mesh, division, get_stats);
        }
    }
}

template<class N> void Spatial_Queries::exec_line_query_leaf(N& n, Box &b, QueryStatistics& qS, Mesh& mesh, bool get_stats)
{
    Box bb;
    pair<int,int> run;

    for(vector<int>::iterator it=n.get_t_array_begin(); it!=n.get_t_array_end(); ++it)
    {
        if(n.get_run_bounding_box(it,bb,mesh,run))
        {
            if(Geometry_Wrapper::line_in_bounding_box(b.get_min(),b.get_max(),bb))
            {
                for(int t_id=run.first; t_id<=run.second; t_id++)
                    atomic_line_in_tetra_test(t_id,b,qS,mesh,get_stats);
            }
        }
        else
        {
            atomic_line_in_tetra_test(*it,b,qS,mesh,get_stats);
        }
    }
}

#endif // SPATIAL_QUERIES_H
