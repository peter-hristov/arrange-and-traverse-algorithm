cmake_minimum_required(VERSION 3.18)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)

# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

set(PROJECT_NAME fv99)

project(${PROJECT_NAME})


message(STATUS "${PROJECT_NAME}")

#add_compile_options(-std=c++11)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenMP REQUIRED)

find_package(CGAL 6.0.1 REQUIRED)
find_package(VTK 9.4.1 REQUIRED)

set(SOURCES
    ./src/Data.cpp
    ./src/io.cpp
    ./src/TetMesh.cpp
    ./src/ReebSpace2.cpp
    ./src/Fiber.cpp
    ./src/Arrangement.cpp
    ./src/ReebSpace.cpp
    ./src/TracerVisualiserMain.cpp
    )

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_precompile_headers(${PROJECT_NAME} PRIVATE ./src/CGALTypedefs.h)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
    ${OPENGL_LIBRARIES}  # required so that include is "-framework OpenGL" on macos instad of "-lOpenGL"
    ${GLUT_LIBRARIES}
    ${VTK_LIBRARIES}
    CGAL::CGAL
    OpenMP::OpenMP_CXX
    )
