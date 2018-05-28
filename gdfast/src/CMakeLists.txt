project(gdfast)
#SET(CMAKE_MODULE_PATH "/home/maarten/tmp/FindCUDA/CMake/cuda" ${CMAKE_MODULE_PATH})
#FIND_PACKAGE(CUDA)
#link_directories(${EXTRA_LIB_DIR} /usr/local/gsl/lib/)
#set(SOURCE_GDFAST jeans.cpp profile.cpp gdfast.cpp )
set(SOURCE parameters.cpp mem.cpp df.cpp profile.cpp profile_2d.cpp sparse_grid_dynamic.cpp sparse_grid.cpp coordinate_systems.cpp profile_axi.cpp profile3d.cpp transformation.cpp mesh2.cpp  poisson_fem2.cpp polynomial.cpp mesh.cpp functions.cpp galaxy.cpp aperture.cpp jeans.cpp gdfast.cpp orbit_integrator.cpp torus_spherical.cpp)
#ndarray/dtype.cpp ndarray/matrix.cpp  ndarray/numpy.cpp  ndarray/scalars.cpp  ndarray/ufunc.cpp ndarray/managers.cpp)
 

#  google_profiler.cpp orbit_integrator_cuda.cu)

set(SOURCE_SCHW gdfast_schw.cpp optimization_schw.cpp) # 

if ( ($ENV{USER} MATCHES "breddels"))
	message("username: $ENV{USER}")
	include_directories("$ENV{HOME}/stow/64/include")
	#include_directories("/opt/local/include/eigen2")
	link_directories("$ENV{HOME}/stow/64/lib")
	link_directories("$ENV{HOME}/stow/64/lib64")
endif ()
if ( ($ENV{USER} MATCHES "maartenbreddels") AND APPLE)
  message("username: $ENV{USER}")
  include_directories("$ENV{HOME}/stow/include")
  include_directories("/opt/local/include/eigen2")
  link_directories("$ENV{HOME}/stow/lib")
endif ()
if ( APPLE)
  message("apple: include directory /opt/local/include")
  include_directories("/opt/local/include")
  message("apple: library directory /opt/local/lib")
  link_directories("/opt/local/lib")
endif ()

find_package(PythonLibs 2.7 REQUIRED)
message("python: including directory ${PYTHON_INCLUDE_PATH}")
include_directories(${PYTHON_INCLUDE_PATH})
set(PYTHON_EXECUTABLE "python")
message("python binary: ${PYTHON_EXECUTABLE}")
execute_process(
    COMMAND
    ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
    OUTPUT_VARIABLE NUMPY_INCLUDE_PATH
    RESULT_VARIABLE NUMPY_ERR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

if(NUMPY_ERR)
  message("numpy not found")
else()
  message("python: including directory ${NUMPY_INCLUDE_PATH}")
  include_directories("${NUMPY_INCLUDE_PATH}")
endif()





add_library(gdfast SHARED ${SOURCE})
add_library(gdfast_schw SHARED ${SOURCE_SCHW})
target_link_libraries(gdfast_schw gsl gslcblas m nlopt)
target_link_libraries(gdfast_schw ${boost_python_lib})

if ( APPLE)
  target_link_libraries(gdfast python2.6) # only osx requires linking to python
  target_link_libraries(gdfast_schw python2.6) # only osx requires linking to python
endif()
find_library(boost_python_lib boost_python)
target_link_libraries(gdfast ${boost_python_lib} )
target_link_libraries(gdfast gsl gslcblas m nlopt) #profiler
#link_directories("/usr/lib64")
link_directories("${BREDHOME}/.local")
target_link_libraries(gdfast) #profiler




set_target_properties(gdfast PROPERTIES PREFIX "")
set_target_properties(gdfast PROPERTIES SUFFIX ".so")
set_target_properties(gdfast PROPERTIES LIBRARY_OUTPUT_DIRECTORY $ENV{HOME}/mab/python/mab/gd)

set_target_properties(gdfast_schw PROPERTIES PREFIX "")
set_target_properties(gdfast_schw PROPERTIES SUFFIX ".so")
set_target_properties(gdfast_schw PROPERTIES LIBRARY_OUTPUT_DIRECTORY $ENV{HOME}/mab/python/mab/gd)

#include_directories(/home/maarten/tmp/pyublas/src/cpp/)
#include_directories(${PYTHON_INCLUDE_PATH})
#include_directories(${interpol_SOURCE_DIR})
#include_directories(BEFORE ${EXTRA_INCLUDE_DIR} /usr/local/gsl/include)
#include_directories(BEFORE ${BOOST_INCLUDE_DIR})
#include_directories(BEFORE ${NUMPY_INCLUDE_DIR1})
#include_directories(BEFORE ${NUMPY_INCLUDE_DIR2})
#include_directories(BEFORE /Users/users/breddels/src/quadpp/src)
#include_directories(BEFORE /Users/users/breddels/local/lib/python2.6/site-packages/PyUblas-0.93.1-py2.6-linux-x86_64.egg/include)
#include_directories(BEFORE /Users/users/breddels/software/lib/python2.6/site-packages/PyUblas-0.93.1-py2.6-linux-x86_64.egg/include)



#set_target_properties(gdfast PROPERTIES CMAKE_CXX_FLAGS -g -Wall -mfpmath=sse,387 -msse3 -O0 -std=c++0xPREFIX "")

#add_library(gdfast_schw SHARED ${SOURCE_SCHW})
#target_link_libraries(gdfast_schw gsl gslcblas m nlopt)
#target_link_libraries(gdfast_schw ${boost_python_lib})

#target_link_libraries(gdfast gsl gslcblas m)
#SET(Boost_USE_STATIC_LIBS OFF)
#SET(Boost_USE_MULTITHREAD OFF)
#FIND_PACKAGE( Boost 1.34 EXACT COMPONENTS python)

#set_target_properties(gdfast_schw PROPERTIES PREFIX "")
#set_target_properties(gdfast_schw PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${BREDHOME}/python/mab/gd)
#set_target_properties(gdfast PROPERTIES LIBRARY_OUTPUT_DIRECTORY /home/maarten/python/mab/gd)


#add_executable(test-deal dealtest.cpp)
#target_link_libraries(test-deal base deal_II_3d deal_II_2d deal_II_1d)

add_executable(hyperquadtree hyperquadtree.cpp)
add_executable(cpptest cpptest.cpp)
add_executable(vectortest vectortest.cpp)
add_executable(sparse_grid2 sparse_grid2.cpp)
#add_executable(testpolynomials testpolynomials.cpp)
#add_executable(testautodiff2 testautodiff2.cpp)
#add_executable(torus torus.cpp)
#add_executable(torus_spherical torus_spherical.cpp)
#target_link_libraries(torus_spherical  ${boost_python_lib} )
#target_link_libraries(torus_spherical gsl gslcblas m nlopt)
#set_target_properties(torus_spherical PROPERTIES CMAKE_CXX_FLAGS -g -Wall -mfpmath=sse,387 -msse3 -O0 -std=c++0xPREFIX "")
#target_link_libraries(test-deal base deal_II_3d deal_II_2d deal_II_1d)
#set_target_properties(test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${go_SOURCE_DIR})
#set_target_properties(test PROPERTIES OUTPUT_NAME test2)
#CUDA_ADD_EXECUTABLE(test test-cuda2.cu )
