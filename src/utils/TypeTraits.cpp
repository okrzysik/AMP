#include "AMP/utils/TypeTraits.h"

#include <array>
#include <complex>
#include <list>
#include <memory>
#include <set>
#include <string_view>
#include <type_traits>
#include <vector>


static_assert( !AMP::is_shared_ptr_v<double> );
static_assert( AMP::is_shared_ptr_v<std::shared_ptr<double>> );

static_assert( !AMP::is_vector_v<double> );
static_assert( AMP::is_vector_v<std::vector<double>> );

static_assert( !std::is_array_v<double> );
static_assert( std::is_array_v<double[4]> );

static_assert( !AMP::is_Array_v<double> );
static_assert( AMP::is_Array_v<AMP::Array<double>> );

static_assert( !AMP::is_container_v<double> );
static_assert( AMP::has_begin_v<std::vector<double>> );
static_assert( AMP::has_begin_v<std::array<double, 5>> );
static_assert( AMP::has_begin_v<std::list<double>> );
static_assert( AMP::has_begin_v<std::set<double>> );
