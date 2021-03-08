// This file contains helper functions and interfaces for reading/writing HDF5
#ifndef included_AMP_TypeTraits
#define included_AMP_TypeTraits

#include <memory>
#include <utility>
#include <vector>

#include "AMP/utils/ArraySize.h" // Forward declare Array


namespace AMP {


//! Checks whether T is a shared_ptr
template<typename T>
struct is_shared_ptr : std::false_type {
};
template<typename T>
struct is_shared_ptr<std::shared_ptr<T>> : std::true_type {
};


//! Checks whether T is a std::vector
template<typename T>
struct is_vector : std::false_type {
};
template<typename T>
struct is_vector<std::vector<T>> : std::true_type {
};


// Function to test if a type is a std::pair
template<typename>
struct is_pair : std::false_type {
};
template<typename T, typename U>
struct is_pair<std::pair<T, U>> : std::true_type {
};


//! Checks whether T is an AMP::Array
template<typename T>
struct is_Array : std::false_type {
};
template<typename T>
struct is_Array<Array<T>> : std::true_type {
};


//! Checks whether T has a size() function
template<typename T>
struct has_size {
private:
    typedef std::true_type yes;
    typedef std::false_type no;
    template<typename U>
    static auto test( int ) -> decltype( std::declval<U>().size() == 1, yes() );
    template<typename>
    static no test( ... );

public:
    static constexpr bool value = std::is_same<decltype( test<T>( 0 ) ), yes>::value;
};


} // namespace AMP


#endif
