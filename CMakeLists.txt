cmake_minimum_required(VERSION 3.31)
project(kur_such)

set(CMAKE_CXX_STANDARD 20)

add_executable(kur_such main.cpp
        Global.h
        BasisFunctions.cpp
        BasisFunctions.h
        LinearAlgebra.h
        LinearAlgebra.cpp
        FEM.h
        FEM.cpp
        IO.h
        IO.cpp)
