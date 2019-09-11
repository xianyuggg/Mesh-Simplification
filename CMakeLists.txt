cmake_minimum_required(VERSION 3.13)
project(MeshSImplification_refined)

set(CMAKE_CXX_STANDARD 11)

add_executable(MeshSImplification_refined main.cpp Mesh.cpp Mesh.h PairHeap.h FaceVertex.h Const.h)