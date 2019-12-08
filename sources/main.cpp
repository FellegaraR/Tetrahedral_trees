#include "main_utility_functions.h"
template<class T> int main_template(T& tree, global_variables &variables);
int main_input_query_generation(global_variables &variables);

int main(int argc, char** argv)
{
    if(argc == 1)
    {
        print_help();
        return 0;
    }

    global_variables variables = global_variables();

    if (read_arguments(argc, argv,variables) == -1)
    {
        print_usage();
        return (EXIT_FAILURE);
    }

    if (variables.isTreeFile)
    {
        //nel caso di lettura da file recupero le info dal nome del file
        setParameters(variables);
    }

    //controllo che tutto sia inizializzato correttamente
    if (!variables.is_getInput)
    {
        if (!checkParameters(variables))
            return (EXIT_FAILURE);
    }

    if(variables.division_type == "ok")
    {
        if (variables.crit_type == "pr")
        {
            P_Tree<OK_Subdivision> tree = P_Tree<OK_Subdivision>(variables.vertices_per_leaf);
            return main_template(tree,variables);
        }
        else if (variables.crit_type == "pm")
        {
            PT_Tree<OK_Subdivision> tree = PT_Tree<OK_Subdivision>(variables.vertices_per_leaf, variables.tetrahedra_per_leaf);
            return main_template(tree,variables);
        }
        else if (variables.crit_type == "pm2")
        {
            T_Tree<OK_Subdivision> tree = T_Tree<OK_Subdivision>(variables.tetrahedra_per_leaf);
            return main_template(tree,variables);
        }
        else if (variables.crit_type == "pmr")
        {
            RT_Tree<OK_Subdivision> tree = RT_Tree<OK_Subdivision>(variables.tetrahedra_per_leaf);
            return main_template(tree,variables);
        }
        else
        {
            cerr << "Not a Valid Criterion Type: Use pr pm pm2 or pmr as criterion value" << endl;
            return -1;
        }
    }
    else if (variables.division_type == "kd")
    {
        if (variables.crit_type == "pr")
        {
            P_Tree<KD_Subdivision> tree = P_Tree<KD_Subdivision>(variables.vertices_per_leaf);
            return main_template(tree,variables);
        }
        else if (variables.crit_type == "pm")
        {
            PT_Tree<KD_Subdivision> tree = PT_Tree<KD_Subdivision>(variables.vertices_per_leaf, variables.tetrahedra_per_leaf);
            return main_template(tree,variables);
        }
        else if (variables.crit_type == "pm2")
        {
            T_Tree<KD_Subdivision> tree = T_Tree<KD_Subdivision>(variables.tetrahedra_per_leaf);
            return main_template(tree,variables);
        }
        else if (variables.crit_type == "pmr")
        {
            RT_Tree<KD_Subdivision> tree = RT_Tree<KD_Subdivision>(variables.tetrahedra_per_leaf);
            return main_template(tree,variables);
        }
        else
        {
            cerr << "Not a Valid Criterion Type: Use pr pm pm2 or pmr as criterion value" << endl;
            return -1;
        }
    }
    else if (variables.is_getInput)
    {
        return main_input_query_generation(variables);
    }
    else
    {
        cerr << "Not a Valid Division Type: Use kd or ok as division value" << endl;
        return -1;
    }

    return (EXIT_SUCCESS);
}

template<class T> int main_template(T& tree, global_variables &variables)
{
    Timer time;

    //Legge l'input
    if (!Reader::read_mesh(tree.get_mesh(), variables.mesh_path))
    {
        cout << "Error Loading .ts file. Execution Stopped." << endl;
        return -1;
    }

    stringstream base_info;
    base_info << variables.vertices_per_leaf << " " << variables.tetrahedra_per_leaf << " " << variables.crit_type << " ";

    if (variables.isTreeFile)
    {
        if (!Reader::read_tree(tree, tree.get_root(), variables.tree_path))
        {
            cerr << "Error Loading .tree file. Execution Stopped." << endl;
            return -1;
        }
    }
    else
    {
        //costruisco l'albero da zero
        stringstream tree_info;
        tree_info << base_info.str() << "Building ";
        time.start();
        tree.build_tree();
        time.stop();
        time.print_elapsed_time(tree_info.str());

        stringstream out;
        if (variables.crit_type == "pr")
            out << get_file_name(variables.mesh_path) << "_" << variables.division_type << "_" << variables.crit_type << "_v_" << variables.vertices_per_leaf << "_.tree";
        else if (variables.crit_type == "pm")
            out << get_file_name(variables.mesh_path) << "_" << variables.division_type << "_" << variables.crit_type << "_v_" << variables.vertices_per_leaf << "_t_" << variables.tetrahedra_per_leaf << "_.tree";
        else if (variables.crit_type == "pmr" || variables.crit_type == "pm2")
            out << get_file_name(variables.mesh_path) << "_" << variables.division_type << "_" << variables.crit_type << "_t_" << variables.tetrahedra_per_leaf << "_.tree";
        Writer::write_tree(out.str(), tree.get_root(), tree.get_decomposition());
    }


    if(variables.reindex)
    {
        time.start();
        Reindexer reindexer = Reindexer();
        reindexer.reindex_tree_and_mesh(tree);
        time.stop();
        time.print_elapsed_time("Index and Mesh Reindexing ");
    }

    Statistics stats;

    if (variables.is_index)
        stats.get_index_statistics(tree,variables.reindex);

    if (variables.query_type != NOTHING)
    {
        Spatial_Queries sq;
        Topological_Queries tq;

        cerr<<base_info.str()<<endl;
        if (variables.query_type == POINT)
            sq.exec_point_locations(tree,variables.query_path,stats);
        else if(variables.query_type == BOX)
            sq.exec_box_queries(tree,variables.query_path,stats);
        else if(variables.query_type == LINE)
        {
            //the face ordering is needed only by the line in tetra test
            Geometry_Wrapper::set_faces_ordering(tree.get_mesh());
            sq.exec_line_queries(tree,variables.query_path,stats);
        }
        else if(variables.query_type == WINDVT)
            tq.windowed_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_decomposition(),variables.query_path,variables.reindex);
        else if(variables.query_type == WINDDIST)
            tq.windowed_Distortion(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_decomposition(),variables.query_path,variables.reindex);
        else if(variables.query_type == WINDTT)
            tq.windowed_TT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_decomposition(),variables.query_path);
        else if(variables.query_type == LINETT)
        {
            //the face ordering is needed only by the line in tetra test
            Geometry_Wrapper::set_faces_ordering(tree.get_mesh());
            tq.linearized_TT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_decomposition(),variables.query_path);
        }
        else if(variables.query_type == BATCH)
        {
            tq.batched_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_decomposition(),variables.reindex);
            tq.batched_TT(tree.get_root(),tree.get_mesh(),tree.get_decomposition());
        }
    }

    return (EXIT_SUCCESS);
}

int main_input_query_generation(global_variables &variables)
{
    Mesh mesh;
    //Legge l'input
    if (!Reader::read_mesh(mesh, variables.mesh_path))
    {
        cout << "Error Loading .ts file. Execution Stopped." << endl;
        return -1;
    }

    if (variables.is_getInput)
    {
        if(variables.input_gen_type == "rand")
        {
            if(variables.query_type == POINT && variables.ratio == 0.0)
                Input_Generator::generate_random_point_inputs(mesh.get_domain(),variables.num_input_entries, get_file_name(variables.mesh_path));
            else if(variables.query_type == BOX && variables.ratio > 0.0)
                Input_Generator::generate_random_box_inputs(mesh.get_domain(),variables.ratio,variables.num_input_entries, get_file_name(variables.mesh_path));
            else if(variables.query_type == LINE && variables.ratio > 0.0)
                Input_Generator::generate_random_line_inputs(mesh.get_domain(),variables.ratio,variables.num_input_entries, get_file_name(variables.mesh_path));
        }
        else if(variables.input_gen_type == "near")
        {
            if(variables.query_type == POINT && variables.ratio == 0.0)
                Input_Generator::generate_near_point_inputs(mesh.get_domain(),variables.num_input_entries, mesh, get_file_name(variables.mesh_path));
            else if(variables.query_type == BOX && variables.ratio > 0.0)
                Input_Generator::generate_near_box_inputs(mesh.get_domain(),variables.ratio,variables.num_input_entries, mesh, get_file_name(variables.mesh_path));
            else if(variables.query_type == LINE && variables.ratio > 0.0)
                Input_Generator::generate_near_line_inputs(mesh.get_domain(),variables.ratio,variables.num_input_entries, mesh, get_file_name(variables.mesh_path));
        }
    }

    return EXIT_SUCCESS;
}
