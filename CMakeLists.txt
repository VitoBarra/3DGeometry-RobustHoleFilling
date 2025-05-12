cmake_minimum_required(VERSION 3.28)
project(3DGeometry-RobustHoleFilling)

set(CMAKE_CXX_STANDARD 23)

#find_package(Eigen REQUIRED)

# Raccoglie tutti i file sorgente .cpp
file(GLOB_RECURSE SOURCE CONFIGURE_DEPENDS src/*.cpp libs/vcglib/wrap/ply/plylib.cpp)
# Raccoglie tutti i file header .hpp (o .h)
file(GLOB_RECURSE HEADER CONFIGURE_DEPENDS include/*.h)

# Usa i file raccolti nell'eseguibile
add_executable(${CMAKE_PROJECT_NAME} ${SOURCE} ${HEADER})

add_subdirectory("libs/vcglib")
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${VCGLib_SOURCE_DIR} )
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC vcglib )