cmake_minimum_required(VERSION 2.6)

set(BREDHOME /Users/users/breddels)
if ( APPLE)
  message("apple: setting c++ compiler to /opt/local/bin/g++-mp-4.5")
  set(CMAKE_CXX_COMPILER "/opt/local/bin/g++-mp-4.5")
endif ()
#set(CMAKE_CXX_FLAGS "-Wall -O2 -std=c++0x -I/usr/lib64/qt-3.3/include" CACHE STRING	"bla" FORCE)
set(CMAKE_CXX_FLAGS "-Wall -O2 -std=c++0x -L/usr/lib64 -L/Users/users/breddels/.local/lib -I/Users/users/breddels/.local/lib/python2.7/site-packages/PyUblas-2011.1-py2.7-linux-x86_64.egg/include -I/Users/users/breddels/.local/include" CACHE STRING	"bla" FORCE)
set(EXTRA_LIB_DIR ${BREDHOME}/.local/lib CACHE PATH "extra directory to search for libraries")
set(EXTRA_INCLUDE_DIR ${BREDHOME}/stow/64/include CACHE PATH "extra directory to search for header files")
#set(BOOST_INCLUDE_DIR ${EXTRA_INCLUDE_DIR}/boost-1_35 CACHE PATH "boost include directory")
set(BOOST_PYTHON_LIBRARY boost_python CACHE PATH "boost::python library file")
set(NUMPY_INCLUDE_DIR1 ${BREDHOME}/stow/64/lib/python2.6/site-packages/numpy/core/include  CACHE PATH "extra include dir for numpy")
set(NUMPY_INCLUDE_DIR2 ${BREDHOME}/stow/64/lib/python2.6/site-packages/numpy/numarray/numpy CACHE PATH "extra include dir for numpy")

project(gdfast C)
add_subdirectory(src)



