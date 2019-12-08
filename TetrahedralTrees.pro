#-------------------------------------------------
#
# Project created by QtCreator 2011-01-20T12:37:39
#
#-------------------------------------------------

#QT       += core
#QT       -= gui

TARGET = tetrahedral_trees
#CONFIG   += console
CONFIG   -= app_bundle
CONFIG -= qt

#TEMPLATE = app
LANGUAGE = C++

# Directories
DESTDIR = dist/
OBJECTS_DIR = build/

LIBS+= -lrt
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3 \
    -march=native

INCLUDEPATH += "sources"

SOURCES += \  
    sources/utilities/sorting.cpp \
    sources/geometry/geometry.cpp \
    sources/geometry/geometry_distortion.cpp \
    sources/geometry/geometry_wrapper.cpp \
    sources/io/reader.cpp \
    sources/io/writer.cpp \
    sources/queries/border_checker.cpp \
    sources/queries/topological_queries.cpp \
    sources/tetrahedral_trees/reindexer.cpp \
    sources/utilities/input_generator.cpp \
    sources/utilities/string_management.cpp \
    sources/utilities/timer.cpp \
    sources/main.cpp \
    sources/queries/spatial_queries.cpp \
    sources/statistics/statistics.cpp \
    sources/basic_types/tetrahedron.cpp \
    sources/tetrahedral_trees/kd_subdivision.cpp \
    sources/tetrahedral_trees/ok_subdivision.cpp \
    sources/tetrahedral_trees/node_t.cpp
    

HEADERS += \    
    sources/basic_types/box.h \
    sources/basic_types/mesh.h \
    sources/basic_types/point.h \
    sources/basic_types/tetrahedron.h \
    sources/basic_types/vertex.h \
    sources/geometry/geometry.h \
    sources/geometry/geometry_distortion.h \
    sources/geometry/geometry_wrapper.h \
    sources/io/reader.h \
    sources/io/writer.h \
    sources/queries/border_checker.h \
    sources/queries/topological_queries.h \
    sources/queries/topological_queries_windowed.h \
    sources/statistics/full_query_statistics.h \
    sources/statistics/index_statistics.h \
    sources/statistics/query_statistics.h \
    sources/statistics/statistics.h \
    sources/tetrahedral_trees/kd_subdivision.h \
    sources/tetrahedral_trees/node.h \
    sources/tetrahedral_trees/node_v.h \
    sources/tetrahedral_trees/ok_subdivision.h \
    sources/tetrahedral_trees/reindexer.h \
    sources/tetrahedral_trees/run_iterator.h \
    sources/tetrahedral_trees/subdivision.h \
    sources/utilities/input_generator.h \
    sources/utilities/sorting.h \
    sources/utilities/sorting_structure.h \
    sources/utilities/string_management.h \
    sources/utilities/timer.h \
    sources/tetrahedral_trees/tree.h \
    sources/tetrahedral_trees/pt_tree.h \
    sources/tetrahedral_trees/t_tree.h \
    sources/tetrahedral_trees/rt_tree.h \
    sources/tetrahedral_trees/p_tree.h \
    sources/main_utility_functions.h \
    sources/queries/spatial_queries.h \
    sources/tetrahedral_trees/node_t.h \
    sources/queries/topological_queries_batched.h \
    sources/basic_types/basic_types.h
    

