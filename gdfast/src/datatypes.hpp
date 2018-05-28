#pragma once

#include <pyublas/numpy.hpp>

namespace gd {

typedef pyublas::numpy_vector<double> double_vector;
typedef pyublas::numpy_matrix<double> double_matrix;
typedef pyublas::numpy_vector<float> float_vector;
typedef pyublas::numpy_matrix<float> float_matrix;

typedef pyublas::numpy_vector<int> int_vector;
typedef pyublas::numpy_matrix<int> int_matrix;
}