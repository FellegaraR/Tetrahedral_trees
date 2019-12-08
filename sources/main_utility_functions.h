#ifndef MAIN_UTILITY_FUNCTIONS_H
#define MAIN_UTILITY_FUNCTIONS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <sstream>

#include "geometry/geometry_wrapper.h"
#include "io/reader.h"
#include "io/writer.h"
#include "queries/spatial_queries.h"
#include "queries/topological_queries.h"
#include "statistics/statistics.h"
#include "tetrahedral_trees/ok_subdivision.h"
#include "tetrahedral_trees/kd_subdivision.h"
#include "tetrahedral_trees/p_tree.h"
#include "tetrahedral_trees/pt_tree.h"
#include "tetrahedral_trees/t_tree.h"
#include "tetrahedral_trees/rt_tree.h"
#include "tetrahedral_trees/reindexer.h"
#include "utilities/input_generator.h"
#include "utilities/string_management.h"
#include "utilities/timer.h"

using namespace std;
using namespace string_management;

#define DEFAULT_X_PER_LEAF -1
#define DEFAULT "null"
enum TopoQueryType { POINT, LINE, BOX, WINDVT, WINDDIST, WINDTT, LINETT, BATCH, NOTHING };
#define BOLD  "\033[1m\033[33m" //for dark background shell
//#define BOLD "\033[1m\033[31m"  //for white background shell
#define RESET   "\033[0m"

/// GLOBAL VARIABLES
struct global_variables
{
    string mesh_path, query_path, exe_name, tree_path;
    string division_type;
    string crit_type;
    bool is_index, is_getInput, isTreeFile, reindex;
    int vertices_per_leaf;
    int tetrahedra_per_leaf;

    int num_input_entries;
    double ratio;
    TopoQueryType query_type;
    string input_gen_type;

    global_variables()
    {
        division_type = DEFAULT;
        crit_type = DEFAULT;
        vertices_per_leaf = DEFAULT_X_PER_LEAF;
        tetrahedra_per_leaf = DEFAULT_X_PER_LEAF;
        query_type = NOTHING;
        is_index = false;
        is_getInput = false;
        isTreeFile = false;
        reindex = false;

        num_input_entries = 0;
        input_gen_type = DEFAULT;
    }
};
///

int read_arguments(int argc, char** argv, global_variables &variables);
void print_usage();
void setParameters(global_variables &variables);
bool checkParameters(global_variables &variables);
void print_help();
void print_paragraph(string stringa, int cols);

void setParameters(global_variables &variables)
{
    vector<string> tokens;
    tokenize(variables.tree_path, tokens, "_");

    for (unsigned int i = 0; i < tokens.size(); i++)
    {
        if (tokens.at(i) == "pm" || tokens.at(i) == "pr" || tokens.at(i) == "cd" || tokens.at(i) == "pmr" || tokens.at(i) == "pm2")
            variables.crit_type = tokens.at(i);
        if (tokens.at(i) == "ok" || tokens.at(i) == "kd")
            variables.division_type = tokens.at(i);
        if (tokens.at(i) == "v")
            variables.vertices_per_leaf = atoi(tokens.at(i + 1).c_str());
        if (tokens.at(i) == "t")
            variables.tetrahedra_per_leaf = atoi(tokens.at(i + 1).c_str());
    }
}

bool checkParameters(global_variables &variables)
{
    if (variables.crit_type == DEFAULT || variables.division_type == DEFAULT)
    {
        cout << "Error initializing criterion or division type. Execution Stopped." << endl;
        print_usage();
        return false;
    }

    if (variables.crit_type == "pm" && (variables.vertices_per_leaf == DEFAULT_X_PER_LEAF || variables.tetrahedra_per_leaf == DEFAULT_X_PER_LEAF))
    {
        cout << "Error initializing vertices_per_leaf or tetrahedra_per_leaf. Execution Stopped." << endl;
        print_usage();
        return false;
    }
    if (variables.crit_type == "pr" && variables.vertices_per_leaf == DEFAULT_X_PER_LEAF)
    {
        cout << "Error initializing vertices_per_leaf. Execution Stopped." << endl;
        print_usage();
        return false;
    }
    if ((variables.crit_type == "pmr" || variables.crit_type == "pm2") && variables.tetrahedra_per_leaf == DEFAULT_X_PER_LEAF)
    {
        cout << "Error initializing tetrahedra_per_leaf. Execution Stopped." << endl;
        print_usage();
        return false;
    }
    return true;
}

int read_arguments(int argc, char** argv, global_variables &variables)
{
    string trash;
    variables.exe_name = argv[0];
    variables.exe_name = strip_path(variables.exe_name);

    for(int i=1; i<(argc); i++)
    {
        char* tag = argv[i];
        if(strcmp(tag, "-i") == 0)
        {
            variables.mesh_path = argv[i+1];
            i++;
        }
        else if(strcmp(tag, "-f") == 0)
        {
            variables.tree_path = argv[i+1];
            variables.isTreeFile = true;
            i++;
        }
        else if(strcmp(tag, "-d") == 0)
        {
            variables.division_type = argv[i+1];
            i++;
        }
        else if(strcmp(tag, "-c") == 0)
        {
            variables.crit_type = argv[i+1];
            i++;
        }
        else if(strcmp(tag, "-v") == 0)
        {
            variables.vertices_per_leaf = atoi(argv[i+1]);
            if (variables.vertices_per_leaf < 1) {
                cerr << "Error: the limit of vertices per leaf must be greater than 0" << endl;
                return -1;
            }
            i++;
        }
        else if(strcmp(tag, "-t") == 0)
        {
            variables.tetrahedra_per_leaf = atoi(argv[i+1]);
            if (variables.tetrahedra_per_leaf < 1) {
                cerr << "Error: the limit of tetrahedra per leaf must be greater than 0" << endl;
                return -1;
            }
            i++;
        }
        else if(strcmp(tag, "-s") == 0)
        {
            variables.is_index = true;
        }
        else if(strcmp(tag, "-r") == 0)
        {
            variables.reindex = true;
        }
        else if(strcmp(tag, "-q") == 0)
        {
            trash = argv[i+1];
            vector<string> tok;
            tokenize(trash,tok,"-");
            if(tok.size()==1 && tok[0]=="batch")
                variables.query_type = BATCH;
            else if(tok.size()<2)
                cerr<<"[-q argument] error when reading arguments"<<endl;
            else
            {
                if(tok[0] == "wvt")
                    variables.query_type = WINDVT;
                else if(tok[0] == "wdist")
                    variables.query_type = WINDDIST;
                else if(tok[0] == "wtt")
                    variables.query_type = WINDTT;
                else if(tok[0] == "ltt")
                    variables.query_type = LINETT;
                else if(tok[0] == "point")
                    variables.query_type = POINT;
                else if(tok[0] == "box")
                    variables.query_type = BOX;
                else if(tok[0] == "line")
                    variables.query_type = LINE;

                variables.query_path = tok[1];
            }
            i++;
        }
        else if(strcmp(tag, "-g") == 0)
        {
            trash = argv[i+1];
            vector<string> tok;
            tokenize(trash,tok,"-");
            if(tok.size()<4)
                cerr<<"[-g argument] error when reading arguments"<<endl;
            else
            {
                if(tok.at(0) == "point")
                    variables.query_type = POINT;
                else if(tok.at(0) == "line")
                    variables.query_type = LINE;
                else if(tok.at(0) == "box")
                    variables.query_type = BOX;

                variables.ratio = atof(tok[1].c_str());
                variables.num_input_entries = atoi(tok[2].c_str());
                variables.input_gen_type = tok[3];
                variables.is_getInput = true;
            }
            i++;
        }
    }

    return 0;
}

void print_usage()
{
    cerr<<"Wrong Usage. Run ./tetrahedraltrees for detailed instructions."<<endl;
}

void print_help()
{
    //annoying stuff to get the dimension of the output shell (!!! not sure it works on Mac,
    //everything based on the command tput. If it doesn't work the dimension si setted to 80 by default)
    FILE* fp;
    char path[1035];

    int cols;
    fp = popen("tput cols", "r");
    if(fp != NULL){
        fgets(path, sizeof(path)-1, fp);
        cols = atoi(path);
    }
    else{
        cols = 80;
    }

    printf(BOLD "\n  NAME:\n\n" RESET);
    printf("\tTetrahedral Trees library\n\n" RESET);

    printf(BOLD "  USAGE: \n\n" RESET);
    printf(BOLD "    ./tetrahedraltrees {<-v [kv] -t [kt] -c [crit] -d [div] | -f [tree_file]>\n" RESET);
    printf(BOLD "                       -q [op-file] -s -r} | {-g [query-ratio-quantity-type]}\n" RESET);
    printf(BOLD "                       -i [mesh_file]\n" RESET);

    printf(BOLD "    -v [kv]\n" RESET);
    print_paragraph("kv is the vertices threshold per leaf. This parameter is needed by P-Ttrees and PT-Ttrees.", cols);
    printf(BOLD "    -t [kt]\n" RESET);
    print_paragraph("kt is the tetrahedra threshold per leaf. This parameter is needed by RT-Ttrees, PT-Ttrees and T-Ttrees.", cols);
    printf(BOLD "    -c [crit]\n" RESET);
    print_paragraph("crit is the criterion type of the index. This can be P-Ttree (pr), RT-Ttree (pmr),  PT-Ttree (pm) or  T-Ttree (pm2).", cols);
    printf(BOLD "    -d [div]\n" RESET);
    print_paragraph("div is the division type of the index. This can be octree (ok) or kD-tree (kd).", cols);

    print_paragraph("NOTA: these arguments must be used in conjunction to create an index. "
                    "This operation generate as output a file containing the tetrahedral index.", cols);

    printf(BOLD "    -f [tree_file]\n" RESET);
    print_paragraph("reads an spatial index from an input file", cols);
    print_paragraph("tree_file contains a Tetrahedral tree index. This file has a fixed syntax of the name "
                    "that allows to recover the informations needed to get the tetrahedraltree index (i.e., kv, kt, division and critirion types)", cols);

    print_paragraph("NOTA: you can use -f argument [OR] {-v / -t / -c / -d} accordingly to the chosen criterion.", cols);

    printf(BOLD "    -q [op-file]\n" RESET);
    print_paragraph("executes a query op, picking the inputs from file", cols);
    print_paragraph("'op' can be: point - box - line - wvt - wdist - wtt - ltt \n"
                    "'point' stands for point location, 'box' for box query, 'line' for line query, "
                    "'wvt' for windowed VT query, 'wdist' windowed Distortion computation, "
                    "'wtt' for windowed TT query and 'ltt' for linearized TT query."
                    "'file' represent the path of the file that contains the inputs for the queries.", cols);

    printf(BOLD "    -g [query-ratio-quantity-type]\n" RESET);
    print_paragraph("generates a given number of input data for a specific query", cols);
    print_paragraph("query can be: point - box - line. "
                    "'ratio' is a number between 0 and 1, and and represents the percentage of the maximum side of the domain to pick. "
                    "'quantity' is a positive number that indicate the number of inputs to generate. "
                    "type can be: near - rand. 'near' stands for a point (picked randomly) that is near to the mesh, "
                    "while 'rand' stands for a point (picked randomly) that is inside the domain.", cols);
    print_paragraph("If 'query' is equal to point 'ratio' must be equal to 0, otherwise 'ratio' must be greater than 0.", cols);

    printf(BOLD "    -s\n" RESET);
    print_paragraph("computes the tetrahedral tree statistics.", cols);
    printf(BOLD "    -r\n" RESET);
    print_paragraph("activate the procedures to exploit the spatial coherence of the index and the mesh.", cols);
    printf(BOLD "    - i [mesh_file]\n" RESET);
    print_paragraph("reads the mesh_file containing the tetrahedral mesh.", cols);

    printf(BOLD "  EXAMPLE[1]: \n" RESET);
    printf("          ./tetrahedraltrees -v 20 -c pr -d ok -s -i mesh.ts\n");
    print_paragraph("reads the mesh [mesh.ts]. Then, builds a P-Ttree index with kv=20 and with subdivision octree. \n"
                    "Finally, it computes the index statistics (-s).", cols);

    printf(BOLD "  EXAMPLE[2]: \n" RESET);
    printf("          ./tetrahedraltrees -f tree_file -q wvt-boxfile -r -i mesh.ts\n");
    print_paragraph("reads the mesh [mesh.ts]. Then, reads the index from tree_file (obtaining the tree "
                    "parameters direcly from the file name) and spatially reordering it (-r). "
                    "Finally, it executes the windowed VT querys, using as query boxes those into 'boxfile'.", cols);

    printf(BOLD "  IMPLEMENTATION:\n" RESET);
    printf("          Author: Riccardo Fellegara\n");
    printf("          Group: G3 Geometry and Graphics Group\n");
    printf("          Man-page Last Update: May 2016\n\n");

    printf(BOLD "  DESCRIPTION: \n" RESET);
    print_paragraph("We address the problem of performing spatial queries on tetrahedral meshes. These latter "
                    "arise in several application domains including 3D GIS, scientific visualization, finite element "
                    "analysis. We have defined and implemented a family of spatial indexes, that we call tetrahedral "
                    "trees. Tetrahedral trees are based on a subdivision of a cubic domain containing the mesh "
                    "defined either by an octree or a 3D kD-tree. For each of them, we have four variants of "
                    "the spatial index, depending on four different subdivision criteria.", cols);
}

void print_paragraph(string stringa, int cols){
    if((int)stringa.size() < cols-20){
        printf("          %s\n\n",stringa.c_str());
    }
    else{
        float dim = (float)(stringa.size()/((float)cols-20));
        int dim_int = dim;
        if(dim > dim_int) dim_int++;
        for(int i=0; i<dim_int; i++){
            if(stringa.at((cols-20)*i) == ' ')
                printf("         %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
            else printf("          %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
        }
        printf("\n");
    }
}


#endif // MAIN_UTILITY_FUNCTIONS_H
