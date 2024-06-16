cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 17)

project(singular_matrix)

set(SOURCES
    main.cpp
    svd.cpp
    eigenvalues/Jacobi.cpp
    math/matrix.cpp
    tests/test.cpp
)
set(HEADERS
    svd.h
    eigenvalues/Jacobi.h
    math/math_utils.h
    math/matrix.h
    math/solving_equation.h
    tests/test.h
)
add_executable(singular_matrix ${SOURCES} ${HEADERS})

target_include_directories(singular_matrix PUBLIC eigenvalues PUBLIC math PUBLIC tests)

