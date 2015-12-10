#ifndef included_ArrayClass
#define included_ArrayClass

#include <iostream>
#include <memory>
#include <array>
#include <functional>
#include "shared_ptr.h"
#include "utils/Utilities.h"

namespace AMP {

#define ARRAY_NDIM_MAX 5        // Maximum number of dimensions supported


#define GET_ARRAY_INDEX(i1,i2,i3) i1+d_N[0]*(i2+d_N[1]*i3)
#if defined(DEBUG) || defined(_DEBUG)
    #define CHECK_ARRAY_INDEX(i1,i2,i3) \
        if ( GET_ARRAY_INDEX(i1,i2,i3)<0 || GET_ARRAY_INDEX(i1,i2,i3)>=d_length ) \
            AMP_ERROR("Index exceeds array bounds");
#else
    #define CHECK_ARRAY_INDEX(i1,i2,i3) 
#endif

#if defined(__CUDA_ARCH__)
    #include <cuda.h>
    #define HOST_DEVICE __host__ __device__
#else
    #define HOST_DEVICE
#endif


/*!
 * Class Array is a multi-dimensional array class written by Mark Berrill
 */
template<class TYPE>
class Array
{
public:

    /*!
     * Create a new empty Array
     */
    Array( );

    /*!
     * Create a new 1D Array with the given number of elements
     * @param N             Number of elements in the array
     */
    Array( size_t N );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     */
    Array( size_t N_rows, size_t N_columns );

    /*!
     * Create a new 3D Array with the given number of rows and columns
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     */
    Array( size_t N1, size_t N2, size_t N3 );

    /*!
     * Create a multi-dimensional Array with the given number of elements
     * @param N             Number of elements in each dimension
     * @param data          Optional raw array to copy the src data
     */
    Array( const std::vector<size_t>& N, const TYPE* data=NULL );

    /*!
     * Copy constructor
     * @param rhs           Array to copy
     */
    Array( const Array& rhs );

    /*!
     * Move constructor
     * @param rhs           Array to copy
     */
    Array( Array&& rhs );

    /*!
     * Assignment operator
     * @param rhs           Array to copy
     */
    Array& operator=( const Array& rhs );

    /*!
     * Move assignment operator
     * @param rhs           Array to copy
     */
    Array& operator=( Array&& rhs );
    
    /*!
     * Assignment operator
     * @param rhs           std::vector to copy
     */
    Array& operator=( const std::vector<TYPE>& rhs );
    
    
    /*!
     * Create a 1D Array view to a raw block of data
     * @param N             Number of elements in the array
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view( size_t N, std::shared_ptr<TYPE> const& data );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view( size_t N_rows, size_t N_columns, std::shared_ptr<TYPE> const& data );

    /*!
     * Create a new 3D Array view to a raw block of data
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view( size_t N1, size_t N2, size_t N3, std::shared_ptr<TYPE> const& data );

    /*!
     * Create a multi-dimensional Array view to a raw block of data
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view( const std::vector<size_t>& N, std::shared_ptr<TYPE> const& data );

   
    
    /*!
     * Create a 1D Array view to a raw block of data
     * @param N             Number of elements in the array
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView( size_t N, std::shared_ptr<const TYPE> const& data );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView( size_t N_rows, size_t N_columns, std::shared_ptr<const TYPE> const& data );

    /*!
     * Create a new 3D Array view to a raw block of data
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView( size_t N1, size_t N2, size_t N3, std::shared_ptr<const TYPE> const& data );

    /*!
     * Create a multi-dimensional Array view to a raw block of data
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView( const std::vector<size_t>& N, std::shared_ptr<const TYPE> const& data );


    /*!
     * Make this object a view of the src
     * @param src           Source vector to take the view of
     */
    void view2( Array& src );

    /*!
     * Make this object a view of the data
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    void view2( const std::vector<size_t>& N, std::shared_ptr<TYPE> const& data );

    /*!
     * Make this object a view of the raw data (expert use only).
     * Use view2( N, std::shared_ptr(data,[](TYPE*){}) ) instead.
     *   Note: this interface is not recommended as it does not protect from
     *   the src data being deleted while still being used by the Array.
     *   Additionally for maximum performance it does not set the internal shared_ptr
     *   so functions like getPtr and resize will not work correctly.
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    void viewRaw( const std::initializer_list<size_t>& N, TYPE* data );


    /*!
     * Convert an array of one type to another.  This may or may not allocate new memory.
     * @param array         Input array
     */
    template<class TYPE2>
    static std::shared_ptr<Array<TYPE2> > convert( std::shared_ptr<Array<TYPE> > array );


    /*!
     * Convert an array of one type to another.  This may or may not allocate new memory.
     * @param array         Input array
     */
    template<class TYPE2>
    static std::shared_ptr<const Array<TYPE2> > convert( std::shared_ptr<const Array<TYPE> > array );


    /*!
     * Copy and convert data from another array to this array
     * @param array         Source array
     */
    template<class TYPE2>
    void copy( const Array<TYPE2>& array );

    /*!
     * Copy and convert data from a raw vector to this array.
     *    Note: The current array must be allocated to the proper size first.
     * @param array         Source array
     */
    template<class TYPE2>
    void copy( const TYPE2* array );

    /*!
     * Copy and convert data from this array to a raw vector.
     * @param array         Source array
     */
    template<class TYPE2>
    void copyTo( TYPE2* array ) const;


    /*!
     * Fill the array with the given value
     * @param value         Value to fill
     */
    void fill( const TYPE& value );

    /*!
     * Scale the array by the given value
     * @param scale         Value to scale by
     */
    void scale( const TYPE& scale );


    //! Destructor
    ~Array( );


    //! Clear the data in the array
    void clear( );


    //! Return the size of the Array
    inline int ndim( ) const { return d_ndim; }


    //! Return the size of the Array
    inline std::vector<size_t> size( ) const { return std::vector<size_t>(d_N,d_N+d_ndim); }


    //! Return the size of the Array
    inline size_t size( int d ) const { return d_N[d]; }


    //! Return the size of the Array
    inline size_t length( ) const { return d_length; }


    //! Return true if the Array is empty
    inline bool empty( ) const { return d_length==0; }


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
    void resize( const std::vector<size_t>& N );

    /*!
     * Resize the given dimension of the array
     * @param dim           The dimension to resize
     * @param N             Number of elements for the given dimension
     * @param value         Value to initialize new elements
     */
    void resizeDim( int dim, size_t N, const TYPE& value );


    /*!
     * Reshape the Array (total size of array will not change)
     * @param N             Number of elements in each dimension
     */
    void reshape( const std::vector<size_t>& N );


    /*!
     * Subset the Array (total size of array will not change)
     * @param index         Index to subset (imin,imax,jmin,jmax,kmin,kmax,...)
     */
    Array<TYPE> subset( const std::vector<size_t>& index ) const;


    /*!
     * Copy data from a subset array
     * @param index         Index of the subset (imin,imax,jmin,jmax,kmin,kmax,...)
     * @param subset        The subset array to copy from
     */
    void copyFromSubset( const std::vector<size_t>& index, const Array<TYPE>& subset );


    /*!
     * Access the desired element
     * @param i             The row index
     */
    HOST_DEVICE inline TYPE& operator()( size_t i ) { CHECK_ARRAY_INDEX(i,0,0) return d_data[i]; }

    /*!
     * Access the desired element
     * @param i             The row index
     */
    HOST_DEVICE inline const TYPE& operator()( size_t i ) const { CHECK_ARRAY_INDEX(i,0,0) return d_data[i]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    HOST_DEVICE inline TYPE& operator()( size_t i, size_t j ) { CHECK_ARRAY_INDEX(i,j,0) return d_data[i+j*d_N[0]]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    HOST_DEVICE inline const TYPE& operator()( size_t i, size_t j ) const { CHECK_ARRAY_INDEX(i,j,0) return d_data[i+j*d_N[0]]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     * @param k             The third index
     */
    HOST_DEVICE inline TYPE& operator()( size_t i, size_t j, size_t k ) { CHECK_ARRAY_INDEX(i,j,k) return d_data[GET_ARRAY_INDEX(i,j,k)]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     * @param k             The third index
     */
    HOST_DEVICE inline const TYPE& operator()( size_t i, size_t j, size_t k ) const { CHECK_ARRAY_INDEX(i,j,k) return d_data[GET_ARRAY_INDEX(i,j,k)]; }


    //! Check if two matricies are equal
    bool operator==( const Array& rhs ) const;

    //! Check if two matricies are not equal
    inline bool operator!=( const Array& rhs ) const { return !this->operator==(rhs); }


    //! Return the pointer to the raw data
    inline std::shared_ptr<TYPE> getPtr( ) { return d_ptr; }

    //! Return the pointer to the raw data
    inline std::shared_ptr<const TYPE> getPtr( ) const { return d_ptr; }

    //! Return the pointer to the raw data
    HOST_DEVICE inline TYPE* get( ) { return d_data; }

    //! Return the pointer to the raw data
    HOST_DEVICE inline const TYPE* get( ) const { return d_data; }


    //! Return true if NaNs are present
    inline bool NaNs( ) const;

    //! Return the smallest value
    inline TYPE min( ) const;

    //! Return the largest value
    inline TYPE max( ) const;

    //! Return the sum of all elements
    inline TYPE sum( ) const;

    //! Return the min of all elements in a given direction
    Array<TYPE> min( int dir ) const;

    //! Return the max of all elements in a given direction
    Array<TYPE> max( int dir ) const;

    //! Return the sum of all elements in a given direction
    Array<TYPE> sum( int dir ) const;

    //! Find all elements that match the operator
    std::vector<size_t> find( const TYPE& value, 
        std::function<bool(const TYPE&,const TYPE&)> compare ) const; 

    //! Add another array
    Array& operator+=( const Array& rhs);

    //! Subtract another array
    Array& operator-=( const Array& rhs);

private:
    int d_ndim;                             // Number of dimensions in array
    size_t d_N[ARRAY_NDIM_MAX];             // Size of each dimension
    size_t d_length;                        // Total length of array
    TYPE *d_data;                           // Raw pointer to data in array
    std::shared_ptr<TYPE> d_ptr;            // Shared pointer to data in array
    void allocate( const std::vector<size_t>& N );
};

}

#include "Array.hpp"

#endif

