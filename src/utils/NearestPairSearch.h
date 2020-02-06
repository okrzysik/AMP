#ifndef included_AMP_NearestPairSearch
#define included_AMP_NearestPairSearch

#include <iostream>
#include <stdlib.h>
#include <utility>


//! Function to compute the closest pair of points
/*!
 * This function will calculate the closest pair of points in a list
 * @param N         The number of points in the list
 * @param x         The coordinates of the verticies (NDIM x N)
 */
template<int NDIM, class TYPE>
inline std::pair<int, int> find_min_dist( const int N, const TYPE *x );


#include "AMP/utils/NearestPairSearch.hpp"

#endif
