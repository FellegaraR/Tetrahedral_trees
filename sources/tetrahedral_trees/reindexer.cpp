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

#include "reindexer.h"

void Reindexer::update_mesh_vertices(Mesh &mesh)
{
    Vertex v = Vertex();
    vector<Vertex> newVertexOrder;
    newVertexOrder.assign(mesh.get_num_vertices(),v);

    for(int i=1;i<=mesh.get_num_vertices();i++)
    {
        if(coherent_indices[i-1] == -1)
        {
            cerr<<"[updateMesh_reorderVertices] INDENTIFIED ISOLATED VERTEX: "<<i<<" "<<coherent_indices[i-1]<<endl;
            int a; cin>>a;
        }
        newVertexOrder[coherent_indices[i-1]-1] = mesh.get_vertex(i);
    }
    mesh.reset_vertices();
    mesh.reserve_vertices_space(newVertexOrder.size());
    for(unsigned int j=0;j<newVertexOrder.size();j++)
        mesh.add_vertex(newVertexOrder.at(j));

    for(int i=1;i<=mesh.get_num_tetrahedra();i++)
    {
        Tetrahedron& t = mesh.get_tetrahedron(i);
        for(int j=0;j<t.vertices_num();j++)
            t.setTV(j,coherent_indices[t.TV(j)-1]);
    }
}

void Reindexer::update_mesh_tetrahedra(Mesh &mesh)
{
    Tetrahedron t = Tetrahedron();
    vector<Tetrahedron> newTopSimplexesOrder;
    newTopSimplexesOrder.assign(mesh.get_num_tetrahedra(),t);

    for(int i=1; i<=mesh.get_num_tetrahedra(); i++)
    {
        newTopSimplexesOrder[coherent_indices[i-1]-1] = mesh.get_tetrahedron(i);
    }

    mesh.reset_tetrahedra();
    mesh.reserve_tetrahedra_space(newTopSimplexesOrder.size());
    for(unsigned i=0; i<newTopSimplexesOrder.size(); i++)
        mesh.add_tetrahedron(newTopSimplexesOrder[i]);
}

void Reindexer::extract_leaves_tetra_association(Mesh &mesh)
{
    map<vector<pair<int,int> >,vector<int> > leaf_tetra_association;
    for(unsigned i=0; i<tetra_leaves_association.size(); i++)
    {
        int t_id = i + 1;
//        cout<<t_id<<" --> T: "<<mesh.get_tetrahedron(t_id)<<" L: ";
//        for(unsigned l=0; l<tetra_leaves_association[i].size(); l++)
//            cout<<"["<<tetra_leaves_association[i][l].first<<" "<<tetra_leaves_association[i][l].second<<"]"<<" ";
//        cout<<endl;
//        int a; cin>>a;
        leaf_tetra_association[tetra_leaves_association[i]].push_back(t_id);
    }

    for(map<vector<pair<int,int> >,vector<int> >::iterator iter=leaf_tetra_association.begin(); iter!=leaf_tetra_association.end(); ++iter)
    {
        const vector<int> &t_list = iter->second;
        for(unsigned t=0; t<t_list.size(); t++)
        {
            coherent_indices[t_list[t]-1] = indices_counter;
            indices_counter++;
        }
    }
}
