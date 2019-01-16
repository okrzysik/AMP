#ifndef included_AMP_ArrayClass
#define included_AMP_ArrayClass


#include <array>
#include <cstring>
#include <functional>
#include <initializer_list>
#include <iosfwd>
#include <memory>
#include <vector>

#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"


namespace AMP {


#if defined( __CUDA_ARCH__ )
#include <cuda.h>
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif
#if defined( USING_GCC ) || defined( USING_CLANG )
#define ATTRIBUTE_INLINE __attribute__( ( always_inline ) )
#else
#define ATTRIBUTE_INLINE
#endif


#if ( defined( DEBUG ) || defined( _DEBUG ) ) && !defined( NDEBUG )
#define CHECK_ARRAY_LENGTH( i )                                      \
    do {                                                             \
        if ( i >= d_length )                                         \
            throw std::out_of_range( "Index exceeds array bounds" ); \
    } while ( 0 )
#else
#define CHECK_ARRAY_LENGTH( i ) \
    do {                        \
    } while ( 0 )
#endif


#define ArraySharedPtr AMP::shared_ptr


// Forward decleration
class FunctionTable;


//! Simple range class
template<class TYPE = size_t>
class Range final
{
public:
    //! Empty constructor
    constexpr Range() : i( 0 ), j( -1 ), k( 1 ) {}

    /*!
     * Create a range i:k:j (or i:j)
     * @param i_            Starting value
     * @param j_            Ending value
     * @param k_            Increment value
     */
    constexpr Range( TYPE i_, TYPE j_, TYPE k_ = 1 ) : i( i_ ), j( j_ ), k( k_ ) {}

    TYPE i, j, k;
};


//! Simple class to store the array dimensions
class ArraySize final
{
public:
    //! Empty constructor
    constexpr ArraySize();

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     */
    constexpr ArraySize( size_t N1 );

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     */
    constexpr ArraySize( size_t N1, size_t N2 );

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     */
    constexpr ArraySize( size_t N1, size_t N2, size_t N3 );

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     * @param N4            Number of elements in the fourth dimension
     */
    constexpr ArraySize( size_t N1, size_t N2, size_t N3, size_t N4 );

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     * @param N4            Number of elements in the fourth dimension
     * @param N5            Number of elements in the fifth dimension
     */
    constexpr ArraySize( size_t N1, size_t N2, size_t N3, size_t N4, size_t N5 );

    /*!
     * Create from initializer list
     * @param N             Size of the array
     */
    constexpr ArraySize( std::initializer_list<size_t> N );

    /*!
     * Create from raw pointer
     * @param ndim          Number of dimensions
     * @param dims          Dimensions
     */
    constexpr ArraySize( size_t ndim, const size_t *dims );

    /*!
     * Create from std::vector
     * @param N             Size of the array
     */
    ArraySize( const std::vector<size_t> &N );

    // Copy/assignment constructors
    constexpr ArraySize( const ArraySize &rhs ) = default;
    constexpr ArraySize( ArraySize &&rhs )      = default;
    constexpr ArraySize &operator=( const ArraySize &rhs ) = default;
    constexpr ArraySize &operator=( ArraySize &&rhs ) = default;

    /*!
     * Access the ith dimension
     * @param i             Index to access
     */
    constexpr inline size_t operator[]( size_t i ) const { return d_N[i]; }

    //! Return the number of dimensions
    constexpr HOST_DEVICE inline uint8_t ndim() const ATTRIBUTE_INLINE { return d_ndim; }

    //! Return the number of dimensions
    constexpr HOST_DEVICE inline size_t size() const ATTRIBUTE_INLINE { return d_ndim; }

    //! Return the total number of elements in the array
    constexpr HOST_DEVICE inline size_t length() const ATTRIBUTE_INLINE { return d_length; }

    //! Resize the dimension
    constexpr void resize( uint8_t dim, size_t N );

    /*!
     * Reshape the Array so that the number of dimensions is the
     *    max of ndim and the largest dim>1.
     * @param ndim          Desired number of dimensions
     */
    constexpr inline void setNdim( uint8_t ndim ) ATTRIBUTE_INLINE
    {
        d_ndim = std::max( ndim, d_ndim );
    }

    //! Returns an iterator to the beginning
    constexpr inline const size_t *begin() const ATTRIBUTE_INLINE { return d_N; }

    //! Returns an iterator to the end
    constexpr inline const size_t *end() const ATTRIBUTE_INLINE { return d_N + d_ndim; }

    // Check if two matrices are equal
    constexpr inline bool operator==( const ArraySize &rhs ) const ATTRIBUTE_INLINE
    {
        return d_ndim == rhs.d_ndim && memcmp( d_N, rhs.d_N, sizeof( d_N ) ) == 0;
    }

    // Check if two matrices are equal
    constexpr inline bool approxEqual( const ArraySize &rhs ) const ATTRIBUTE_INLINE
    {
        return ( length() == 0 && rhs.length() ) || memcmp( d_N, rhs.d_N, sizeof( d_N ) ) == 0;
    }

    //! Check if two matrices are not equal
    constexpr HOST_DEVICE inline bool operator!=( const ArraySize &rhs ) const ATTRIBUTE_INLINE
    {
        return d_ndim != rhs.d_ndim || memcmp( d_N, rhs.d_N, sizeof( d_N ) ) != 0;
    }

    //! Maximum supported dimension
    constexpr static uint8_t maxDim() ATTRIBUTE_INLINE { return 5u; }

    //! Get the index
    constexpr HOST_DEVICE inline size_t index( size_t i ) const ATTRIBUTE_INLINE
    {
        CHECK_ARRAY_LENGTH( i );
        return i;
    }

    //! Get the index
    constexpr HOST_DEVICE inline size_t index( size_t i1, size_t i2 ) const ATTRIBUTE_INLINE
    {
        size_t index = i1 + i2 * d_N[0];
        CHECK_ARRAY_LENGTH( index );
        return index;
    }

    //! Get the index
    constexpr HOST_DEVICE inline size_t
    index( size_t i1, size_t i2, size_t i3 ) const ATTRIBUTE_INLINE
    {
        size_t index = i1 + d_N[0] * ( i2 + d_N[1] * i3 );
        CHECK_ARRAY_LENGTH( index );
        return index;
    }

    //! Get the index
    constexpr HOST_DEVICE inline size_t
    index( size_t i1, size_t i2, size_t i3, size_t i4 ) const ATTRIBUTE_INLINE
    {
        size_t index = i1 + d_N[0] * ( i2 + d_N[1] * ( i3 + d_N[2] * i4 ) );
        CHECK_ARRAY_LENGTH( index );
        return index;
    }

    //! Get the index
    constexpr HOST_DEVICE inline size_t
    index( size_t i1, size_t i2, size_t i3, size_t i4, size_t i5 ) const ATTRIBUTE_INLINE
    {
        size_t index = i1 + d_N[0] * ( i2 + d_N[1] * ( i3 + d_N[2] * ( i4 + d_N[3] * i5 ) ) );
        CHECK_ARRAY_LENGTH( index );
        return index;
    }

private:
    uint8_t d_ndim;
    size_t d_length;
    size_t d_N[5];
};


/*!
 * Class Array is a multi-dimensional array class written by Mark Berrill
 */
template<class TYPE, class FUN = FunctionTable, class Allocator = std::allocator<TYPE>>
class Array final
{
public: // Constructors / assignment operators
    /*!
     * Create a new empty Array
     */
    Array();

    /*!
     * Create an Array with the given size
     * @param N             Size of the array
     */
    explicit Array( const ArraySize &N );

    /*!
     * Create a new 1D Array with the given number of elements
     * @param N             Number of elements in the array
     */
    explicit Array( size_t N );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     */
    explicit Array( size_t N_rows, size_t N_columns );

    /*!
     * Create a new 3D Array with the given number of rows and columns
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     */
    explicit Array( size_t N1, size_t N2, size_t N3 );

    /*!
     * Create a new 4D Array with the given number of rows and columns
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     * @param N4            Number of elements in the fourth dimension
     */
    explicit Array( size_t N1, size_t N2, size_t N3, size_t N4 );

    /*!
     * Create a new 4D Array with the given number of rows and columns
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     * @param N4            Number of elements in the fourth dimension
     * @param N5            Number of elements in the fifth dimension
     */
    explicit Array( size_t N1, size_t N2, size_t N3, size_t N4, size_t N5 );

    /*!
     * Create a multi-dimensional Array with the given number of elements
     * @param N             Number of elements in each dimension
     * @param data          Optional raw array to copy the src data
     */
    explicit Array( const std::vector<size_t> &N, const TYPE *data = NULL );

    /*!
     * Create a 1D Array with the range
     * @param range         Range of the data
     */
    explicit Array( const Range<TYPE> &range );

    /*!
     * Create a 1D Array using a string that mimic's MATLAB
     * @param range         Range of the data
     */
    explicit Array( std::string range );

    /*!
     * Create a 1D Array with the given initializer list
     * @param data          Input data
     */
    Array( std::initializer_list<TYPE> data );


    /*!
     * Copy constructor
     * @param rhs           Array to copy
     */
    Array( const Array &rhs );

    /*!
     * Move constructor
     * @param rhs           Array to copy
     */
    Array( Array &&rhs );

    /*!
     * Assignment operator
     * @param rhs           Array to copy
     */
    Array &operator=( const Array &rhs );

    /*!
     * Move assignment operator
     * @param rhs           Array to copy
     */
    Array &operator=( Array &&rhs );

    /*!
     * Assignment operator
     * @param rhs           std::vector to copy
     */
    Array &operator=( const std::vector<TYPE> &rhs );


public: // Views/copies/subset
    /*!
     * Create a 1D Array view to a raw block of data
     * @param N             Number of elements in the array
     * @param data          Pointer to the data
     */
    static ArraySharedPtr<Array> view( size_t N, ArraySharedPtr<TYPE> const &data );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     * @param data          Pointer to the data
     */
    static ArraySharedPtr<Array>
    view( size_t N_rows, size_t N_columns, ArraySharedPtr<TYPE> const &data );

    /*!
     * Create a new 3D Array view to a raw block of data
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     * @param data          Pointer to the data
     */
    static ArraySharedPtr<Array>
    view( size_t N1, size_t N2, size_t N3, ArraySharedPtr<TYPE> const &data );

    /*!
     * Create a multi-dimensional Array view to a raw block of data
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    static ArraySharedPtr<Array> view( const ArraySize &N, ArraySharedPtr<TYPE> const &data );


    /*!
     * Create a 1D Array view to a raw block of data
     * @param N             Number of elements in the array
     * @param data          Pointer to the data
     */
    static ArraySharedPtr<const Array> constView( size_t N,
                                                  ArraySharedPtr<const TYPE> const &data );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     * @param data          Pointer to the data
     */
    static ArraySharedPtr<const Array>
    constView( size_t N_rows, size_t N_columns, ArraySharedPtr<const TYPE> const &data );

    /*!
     * Create a new 3D Array view to a raw block of data
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     * @param data          Pointer to the data
     */
    static ArraySharedPtr<const Array>
    constView( size_t N1, size_t N2, size_t N3, ArraySharedPtr<const TYPE> const &data );

    /*!
     * Create a multi-dimensional Array view to a raw block of data
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    static ArraySharedPtr<const Array> constView( const ArraySize &N,
                                                  ArraySharedPtr<const TYPE> const &data );


    /*!
     * Make this object a view of the src
     * @param src           Source vector to take the view of
     */
    void view2( Array &src );

    /*!
     * Make this object a view of the data
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    void view2( const ArraySize &N, ArraySharedPtr<TYPE> const &data );

    /*!
     * Make this object a view of the raw data (expert use only).
     * Use view2( N, ArraySharedPtr(data,[](TYPE*){}) ) instead.
     *   Note: this interface is not recommended as it does not protect from
     *   the src data being deleted while still being used by the Array.
     *   Additionally for maximum performance it does not set the internal shared_ptr
     *   so functions like getPtr and resize will not work correctly.
     * @param ndim          Number of dimensions
     * @param dims          Number of elements in each dimension
     * @param data          Pointer to the data
     */
    void viewRaw( int ndim, const size_t *dims, TYPE *data );

    /*!
     * Make this object a view of the raw data (expert use only).
     * Use view2( N, ArraySharedPtr(data,[](TYPE*){}) ) instead.
     *   Note: this interface is not recommended as it does not protect from
     *   the src data being deleted while still being used by the Array.
     *   Additionally for maximum performance it does not set the internal shared_ptr
     *   so functions like getPtr and resize will not work correctly.
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    void viewRaw( const ArraySize &N, TYPE *data );

    /*!
     * Convert an array of one type to another.  This may or may not allocate new memory.
     * @param array         Input array
     */
    template<class TYPE2>
    static ArraySharedPtr<Array<TYPE2, FUN, Allocator>>
    convert( ArraySharedPtr<Array<TYPE, FUN, Allocator>> array );


    /*!
     * Convert an array of one type to another.  This may or may not allocate new memory.
     * @param array         Input array
     */
    template<class TYPE2>
    static ArraySharedPtr<const Array<TYPE2, FUN, Allocator>>
    convert( ArraySharedPtr<const Array<TYPE, FUN, Allocator>> array );


    /*!
     * Copy and convert data from another array to this array
     * @param array         Source array
     */
    template<class TYPE2>
    void copy( const Array<TYPE2, FUN, Allocator> &array );

    /*!
     * Copy and convert data from a raw vector to this array.
     *    Note: The current array must be allocated to the proper size first.
     * @param array         Source array
     */
    template<class TYPE2>
    void copy( const TYPE2 *array );

    /*!
     * Copy and convert data from this array to a raw vector.
     * @param array         Source array
     */
    template<class TYPE2>
    void copyTo( TYPE2 *array ) const;

    /*!
     * Copy and convert data from this array to a new array
     */
    template<class TYPE2>
    Array<TYPE2, FUN, Allocator> cloneTo() const;

    /*! swap the raw data pointers for the Arrays after checking for compatibility */
    void swap( Array &other );

    /*!
     * Fill the array with the given value
     * @param value         Value to fill
     */
    void fill( const TYPE &value );

    /*!
     * Scale the array by the given value
     * @param scale         Value to scale by
     */
    void scale( const TYPE &scale );

    /*!
     * Set the values of this array to pow(base, exp)
     * @param base        Base array
     * @param exp         Exponent value
     */
    void pow( const Array &base, const TYPE &exp );

    //! Destructor
    ~Array();


    //! Clear the data in the array
    void clear();


    //! Return the size of the Array
    inline int ndim() const { return d_size.ndim(); }


    //! Return the size of the Array
    inline ArraySize &size() { return d_size; }


    //! Return the size of the Array
    inline const ArraySize &size() const { return d_size; }


    //! Return the size of the Array
    inline size_t size( int d ) const { return d_size[d]; }


    //! Return the size of the Array
    inline size_t length() const { return d_size.length(); }


    //! Return true if the Array is empty
    inline bool empty() const { return d_size.length() == 0; }


    /*!
     * Resize the Array
     * @param N             NUmber of elements
     */
    void resize( size_t N );

    /*!
     * Resize the Array
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     */
    void resize( size_t N_rows, size_t N_columns );

    /*!
     * Resize the Array
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     */
    void resize( size_t N1, size_t N2, size_t N3 );

    /*!
     * Resize the Array
     * @param N             Number of elements in each dimension
     */
    void resize( const ArraySize &N );


    /*!
     * Resize the given dimension of the array
     * @param dim           The dimension to resize
     * @param N             Number of elements for the given dimension
     * @param value         Value to initialize new elements
     */
    void resizeDim( int dim, size_t N, const TYPE &value );


    /*!
     * Reshape the Array (total size of array will not change)
     * @param N             Number of elements in each dimension
     */
    void reshape( const ArraySize &N );


    /*!
     * Reshape the Array so that the number of dimensions is the
     *    max of ndim and the largest dim>1.
     * @param ndim          Desired number of dimensions
     */
    inline void setNdim( int ndim ) { d_size.setNdim( ndim ); }


    /*!
     * Subset the Array (total size of array will not change)
     * @param index         Index to subset (imin,imax,jmin,jmax,kmin,kmax,...)
     */
    template<class TYPE2 = TYPE>
    Array<TYPE2, FUN, Allocator> subset( const std::vector<size_t> &index ) const;


    /*!
     * Subset the Array (total size of array will not change)
     * @param index         Index to subset (ix:kx:jx,iy:ky:jy,...)
     */
    template<class TYPE2 = TYPE>
    Array<TYPE2, FUN, Allocator> subset( const std::vector<Range<size_t>> &index ) const;


    /*!
     * Copy data from an array into a subset of this array
     * @param index         Index of the subset (imin,imax,jmin,jmax,kmin,kmax,...)
     * @param subset        The subset array to copy from
     */
    template<class TYPE2>
    void copySubset( const std::vector<size_t> &index, const Array<TYPE2, FUN, Allocator> &subset );

    /*!
     * Copy data from an array into a subset of this array
     * @param index         Index of the subset
     * @param subset        The subset array to copy from
     */
    template<class TYPE2>
    void copySubset( const std::vector<Range<size_t>> &index,
                     const Array<TYPE2, FUN, Allocator> &subset );

    /*!
     * Add data from an array into a subset of this array
     * @param index         Index of the subset (imin,imax,jmin,jmax,kmin,kmax,...)
     * @param subset        The subset array to add from
     */
    void addSubset( const std::vector<size_t> &index, const Array<TYPE, FUN, Allocator> &subset );

    /*!
     * Add data from an array into a subset of this array
     * @param index         Index of the subset
     * @param subset        The subset array to add from
     */
    void addSubset( const std::vector<Range<size_t>> &index,
                    const Array<TYPE, FUN, Allocator> &subset );


public: // Accessors
    /*!
     * Access the desired element
     * @param i             The row index
     */
    HOST_DEVICE inline TYPE &operator()( size_t i ) ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i )];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     */
    HOST_DEVICE inline const TYPE &operator()( size_t i ) const ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i )];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    HOST_DEVICE inline TYPE &operator()( size_t i, size_t j ) ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i, j )];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    HOST_DEVICE inline const TYPE &operator()( size_t i, size_t j ) const ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i, j )];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     * @param k             The third index
     */
    HOST_DEVICE inline TYPE &operator()( size_t i, size_t j, size_t k ) ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i, j, k )];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     * @param k             The third index
     */
    HOST_DEVICE inline const TYPE &operator()( size_t i, size_t j, size_t k ) const ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i, j, k )];
    }

    /*!
     * Access the desired element
     * @param i1            The first index
     * @param i2            The second index
     * @param i3            The third index
     * @param i4            The fourth index
     */
    HOST_DEVICE inline TYPE &
    operator()( size_t i1, size_t i2, size_t i3, size_t i4 ) ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i1, i2, i3, i4 )];
    }

    /*!
     * Access the desired element
     * @param i1            The first index
     * @param i2            The second index
     * @param i3            The third index
     * @param i4            The fourth index
     */
    HOST_DEVICE inline const TYPE &
    operator()( size_t i1, size_t i2, size_t i3, size_t i4 ) const ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i1, i2, i3, i4 )];
    }

    /*!
     * Access the desired element
     * @param i1            The first index
     * @param i2            The second index
     * @param i3            The third index
     * @param i4            The fourth index
     * @param i5            The fifth index
     */
    HOST_DEVICE inline TYPE &
    operator()( size_t i1, size_t i2, size_t i3, size_t i4, size_t i5 ) ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i1, i2, i3, i4, i5 )];
    }

    /*!
     * Access the desired element
     * @param i1            The first index
     * @param i2            The second index
     * @param i3            The third index
     * @param i4            The fourth index
     * @param i5            The fifth index
     */
    HOST_DEVICE inline const TYPE &
    operator()( size_t i1, size_t i2, size_t i3, size_t i4, size_t i5 ) const ATTRIBUTE_INLINE
    {
        return d_data[d_size.index( i1, i2, i3, i4, i5 )];
    }

    /*!
     * Access the desired element as a raw pointer
     * @param i             The global index
     */
    HOST_DEVICE inline TYPE *ptr( size_t i ) ATTRIBUTE_INLINE
    {
        return i >= d_size.length() ? nullptr : &d_data[i];
    }

    /*!
     * Access the desired element as a raw pointer
     * @param i             The global index
     */
    HOST_DEVICE inline const TYPE *ptr( size_t i ) const ATTRIBUTE_INLINE
    {
        return i >= d_size.length() ? nullptr : &d_data[i];
    }

    //! Get iterator to beginning of data
    inline TYPE *begin() ATTRIBUTE_INLINE { return d_data; }

    //! Get iterator to beginning of data
    inline const TYPE *begin() const ATTRIBUTE_INLINE { return d_data; }

    //! Get iterator to beginning of data
    inline TYPE *end() ATTRIBUTE_INLINE { return d_data + d_size.length(); }

    //! Get iterator to beginning of data
    inline const TYPE *end() const ATTRIBUTE_INLINE { return d_data + d_size.length(); }

    //! Return the pointer to the raw data
    inline ArraySharedPtr<TYPE> getPtr() ATTRIBUTE_INLINE { return d_ptr; }

    //! Return the pointer to the raw data
    inline ArraySharedPtr<const TYPE> getPtr() const ATTRIBUTE_INLINE { return d_ptr; }

    //! Return the pointer to the raw data
    HOST_DEVICE inline TYPE *data() ATTRIBUTE_INLINE { return d_data; }

    //! Return the pointer to the raw data
    HOST_DEVICE inline const TYPE *data() const ATTRIBUTE_INLINE { return d_data; }


public: // Operator overloading
    //! Check if two matrices are equal
    // Equality means the dimensions and data have to be identical
    bool operator==( const Array &rhs ) const;

    //! Check if two matrices are not equal
    inline bool operator!=( const Array &rhs ) const { return !this->operator==( rhs ); }

    //! Add another array
    Array &operator+=( const Array &rhs );

    //! Subtract another array
    Array &operator-=( const Array &rhs );

    //! Add a scalar
    Array &operator+=( const TYPE &rhs );

    //! Subtract a scalar
    Array &operator-=( const TYPE &rhs );


public: // Math operations
    //! Concatenates the arrays along the dimension dim.
    static Array cat( const std::vector<Array> &x, int dim = 0 );

    //! Concatenates a given array with the current array
    void cat( const Array &x, int dim = 0 );

    //! Initialize the array with random values (defined from the function table)
    void rand();

    //! Return true if NaNs are present
    bool NaNs() const;

    //! Return the smallest value
    TYPE min() const;

    //! Return the largest value
    TYPE max() const;

    //! Return the sum of all elements
    TYPE sum() const;

    //! Return the mean of all elements
    TYPE mean() const;

    //! Return the min of all elements in a given direction
    Array min( int dir ) const;

    //! Return the max of all elements in a given direction
    Array max( int dir ) const;

    //! Return the sum of all elements in a given direction
    Array sum( int dir ) const;

    //! Return the smallest value
    TYPE min( const std::vector<size_t> &index ) const;

    //! Return the largest value
    TYPE max( const std::vector<size_t> &index ) const;

    //! Return the sum of all elements
    TYPE sum( const std::vector<size_t> &index ) const;

    //! Return the mean of all elements
    TYPE mean( const std::vector<size_t> &index ) const;

    //! Return the smallest value
    TYPE min( const std::vector<Range<size_t>> &index ) const;

    //! Return the largest value
    TYPE max( const std::vector<Range<size_t>> &index ) const;

    //! Return the sum of all elements
    TYPE sum( const std::vector<Range<size_t>> &index ) const;

    //! Return the mean of all elements
    TYPE mean( const std::vector<Range<size_t>> &index ) const;

    //! Find all elements that match the operator
    std::vector<size_t> find( const TYPE &value,
                              std::function<bool( const TYPE &, const TYPE & )> compare ) const;


    //! Print an array
    void
    print( std::ostream &os, const std::string &name = "A", const std::string &prefix = "" ) const;

    //! Multiply two arrays
    static Array multiply( const Array &a, const Array &b );

    //! Transpose an array
    Array reverseDim() const;

    //! Replicate an array a given number of times in each direction
    Array repmat( const std::vector<size_t> &N ) const;

    //! Coarsen an array using the given filter
    Array coarsen( const Array &filter ) const;

    //! Coarsen an array using the given filter
    Array coarsen( const std::vector<size_t> &ratio,
                   std::function<TYPE( const Array & )> filter ) const;

    /*!
     * Perform a element-wise operation y = f(x)
     * @param[in] fun           The function operation
     * @param[in] x             The input array
     */
    static Array transform( std::function<TYPE( const TYPE & )> fun, const Array &x );

    /*!
     * Perform a element-wise operation z = f(x,y)
     * @param[in] fun           The function operation
     * @param[in] x             The first array
     * @param[in] y             The second array
     */
    static Array transform( std::function<TYPE( const TYPE &, const TYPE & )> fun,
                            const Array &x,
                            const Array &y );

    /*!
     * axpby operation: this = alpha*x + beta*this
     * @param[in] alpha         alpha
     * @param[in] x             x
     * @param[in] beta          beta
     */
    void axpby( const TYPE &alpha, const Array &x, const TYPE &beta );

    /*!
     * Linear interpolation
     * @param[in] x             Position as a decimal index
     */
    inline TYPE interp( const std::vector<double> &x ) { return interp( x.data() ); }

    /*!
     * Linear interpolation
     * @param[in] x             Position as a decimal index
     */
    inline TYPE interp( const double *x );

    /**
     * \fn equals (Array & const rhs, TYPE tol )
     * \brief  Determine if two Arrays are equal using an absolute tolerance
     * \param[in] rhs Vector to compare to
     * \param[in] tol Tolerance of comparison
     * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
     */
    bool equals( const Array &rhs, TYPE tol = 0.000001 ) const;

private:
    ArraySize d_size;           // Size of each dimension
    TYPE *d_data;               // Raw pointer to data in array
    ArraySharedPtr<TYPE> d_ptr; // Shared pointer to data in array
    void allocate( const ArraySize &N );

public:
    template<class TYPE2, class FUN2, class Allocator2>
    inline bool sizeMatch( const Array<TYPE2, FUN2, Allocator2> &rhs ) const
    {
        return d_size == rhs.d_size;
    }

private:
    inline void checkSubsetIndex( const std::vector<Range<size_t>> &range ) const;
    inline std::vector<Range<size_t>> convert( const std::vector<size_t> &index ) const;
    static inline void getSubsetArrays( const std::vector<Range<size_t>> &range,
                                        std::array<size_t, 5> &first,
                                        std::array<size_t, 5> &last,
                                        std::array<size_t, 5> &inc,
                                        std::array<size_t, 5> &N );
};


} // namespace AMP
#include "AMP/utils/Array.hpp"
#include "AMP/utils/FunctionTable.h"

#endif
