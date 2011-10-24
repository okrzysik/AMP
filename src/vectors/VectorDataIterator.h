#ifndef included_AMP_VectorIterators
#define included_AMP_VectorIterators

#include <iterator>

#include "utils/Utilities.h"

#include "VectorIndexer.h"

namespace AMP {
namespace LinearAlgebra {

  class MultiVector;
  class ConstVectorDataIterator;

  /**
    * \class VectorDataIterator
    * \brief  Iterator for local data in a vector
    *
    * \details Even though a vector may have non-contiguous storage of data, the
    * interface presents a contiguous block of memory to the user:  each element
    * in the vector is given an offset from 0 and the vector is packed.
    * This allows for a random access iterator on the data.
    *
    * Vector::begin() and Vector::end() return this class.  This class
    * uses the DataBlock interface in vectors to access data.  As a result,
    * for some non-AMP managed vectors, this class may not be the most efficient.
    */

  class VectorDataIterator
  {
    private:
      MultiVector *d_MultiVector;
      Vector   *d_Vec;
      double   *d_Block;
      size_t    d_CurBlock, d_CurOffset;
      VectorIndexer::shared_ptr   d_Indexer;

      size_t   dbSize () const;

      static void advance ( int , VectorDataIterator & );
      static void recede  ( int , VectorDataIterator & );
      static int subtract ( const VectorDataIterator & , const VectorDataIterator & );

    public:
      /** \brief  Convenince typedef for testing
        */
      typedef  Vector        vector_type;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  int            difference_type;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  double         value_type;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  double *       pointer;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  double &       reference;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  std::random_access_iterator_tag  iterator_category;

      /** \brief  Default constructor
        * \details  Creates an iterator that points nowhere in particular
        */
      VectorDataIterator ();

      /** \brief Copy constructor
        * \param[in] rhs iterator to copy
        * \details  Copies an iterator
        */
      VectorDataIterator ( const VectorDataIterator &rhs );

      /** \brief Constructor from a vector
        * \param[in] p  a (non-reference counted) pointer to the vector being iterated over
        * \param[in] block the contiguous block in the vector that holds the data being pointed to
        * \param[in] offset the offset of the data being pointed to.
        * \param[in] ndx An optional indexer for handling iterators on sparse vectors.
        * \details Vector::begin() and Vector::end() call this function to instantiate a new iterator.
        * This method should not be invoked outside of this use.
        */
      VectorDataIterator ( Vector *p , int block , int offset , VectorIndexer::shared_ptr ndx = VectorIndexer::shared_ptr () );

      /** \brief Constructor used to advance an iterator to the next vector in a multivector
        * \param[in]  i  The iterator to advance
        * \param[in]  p  The multivector pointed to
        * \param[in]  blockNum  The next block to use in the multivector
        */
      VectorDataIterator ( VectorDataIterator i , MultiVector *p , size_t blockNum );

      /** \brief Dereference the iterator
        * \return Value pointed to by the iterator
        * \warning Setting values in the vector with the iterator requires firing of dataChanged event.
        * \see PetscVector
        * \see DataChangeListener
        * \details  This returns a reference to the data pointed to by the iterator
        */
      double   & operator * ();

      /** \brief Test for equality
        * \return True iff iterators point to the same place on the same vector
        * \details Returns true iff rhs points to the exact same offset of the same vector.
        * It is possible for two iterators to point to the same spot in memory and return false.  This
        * is due to the abstraction of a contiguous block presented by Vector.
        */
      bool       operator == ( const VectorDataIterator &rhs ) const;

      /** \brief Test for inequality
        * \returns True iff !(*this == rhs)
        * \details Returns !(*this == rhs)
        */
      bool       operator != ( const VectorDataIterator &rhs ) const;

      /** \brief Increment the iterator
        * \returns a reference to this iterator
        */
      VectorDataIterator  &operator++ ();

      /** \brief Increment the iterator
        * \returns a copy of this iterator
        */
      VectorDataIterator   operator++ ( int  );

      /** \brief Decrement the iterator
        * \returns a reference to this iterator
        */
      VectorDataIterator  &operator-- ();

      /** \brief Decrement the iterator
        * \returns a copy of this iterator
        */
      VectorDataIterator   operator-- ( int  );

      /** \brief Add a constant to this iterator
        * \param[in] i Offset to move forward
        * \returns A new iterator i steps forward
        * \details  Equivalent to
        * \code
          VectorDataIterator t = *this;
          for ( int k = 0 ; k != i ; k++ )
            ++t;
          \endcode
        * except it is \f$O(1)\f$ in most cases.
        */
      VectorDataIterator   operator + ( int i );

      /** \brief Add a constant to this iterator
        * \param[in] i Offset to move forward
        * \returns A new iterator i steps forward
        * \details  Equivalent to
        * \code
          for ( int k = 0 ; k != i ; k++ )
            operator ++ ();
          \endcode
        * except it is \f$O(1)\f$ in most cases.
        */
      VectorDataIterator  &operator += ( int i );

      /** \brief Subtract a constant to this iterator
        * \param[in] i Offset to move backward
        * \returns A new iterator i steps backward
        * \details  Equivalent to
        * \code
          VectorDataIterator t = *this;
          for ( int k = 0 ; k != i ; k++ )
            --t;
          \endcode
        * except it is \f$O(1)\f$ in most cases.
        */
      VectorDataIterator   operator - ( int i );

      /** \brief Subtract a constant to this iterator
        * \param [in] i Offset to move backward
        * \returns A new iterator i steps backward
        * \details  Equivalent to
        * \code
          for ( int k = 0 ; k != i ; k++ )
            operator -- ();
          \endcode
        * except it is \f$O(1)\f$ in most cases.
        */
      VectorDataIterator  &operator -= ( int i );

      /** \brief Compute distance between two iterators
        * \param  rhs The other iterator to measure distance to
        * \returns Number of entities between iterators.  If rhs == this,
        * the answer is 0. Else, if rhs is after this, then the answer is
        * positive, otherwise it is negative.
        * \details If rhs is after this, this is equivalent to
        * \code
          int answer = 0;
          VectorDataIterator  t = *this;
          while ( t != rhs )
          {
            t++;
            answer++;
          }
          return answer;
          \endcode
        * except it is \f$O(1)\f$ in most cases.  If rhs is before this,
        * the equivalent code is the same except t and answer are
        * decremented.
        */
      int                  operator - ( const VectorDataIterator &rhs ) const;


      /** \brief Return data a distance from the current iterator
        * \param i Offset from this iterator
        * \return The data pointed to at this + i
        * \details  Equivalent to
        * \code
          VectorDataIterator t = *this + i;
          return *t;
          \endcode
        */
      double              &operator [] ( int i );

      friend class ConstVectorDataIterator;
  };

  /**
    * \class ConstVectorDataIterator
    * \brief  Iterator for local data in a vector
    *
    * \details Even though a vector may have non-contiguous storage of data, the
    * interface presents a contiguous block of memory to the user:  each element
    * in the vector is given an offset from 0 and the vector is packed.
    * This allows for a random access iterator on the data.
    *
    * Vector::begin() and Vector::end() return this class.  This class
    * uses the DataBlock interface in vectors to access data.  As a result,
    * for some non-AMP managed vectors, this class may not be the most efficient.
    */

  class ConstVectorDataIterator
  {
    private:
      MultiVector * d_MultiVector;
      const Vector   *d_Vec;
      const double   *d_Block;
      size_t    d_CurBlock, d_CurOffset;
      VectorIndexer::shared_ptr  d_Indexer;

      size_t   dbSize () const;

      static void advance ( int , ConstVectorDataIterator & );
      static void recede  ( int , ConstVectorDataIterator & );
      static int subtract ( const ConstVectorDataIterator & , const ConstVectorDataIterator & );

    public:
      /** \brief  Convenince typedef for testing
        */
      typedef  const Vector   vector_type;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  int            difference_type;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  const double         value_type;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  const double *       pointer;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  const double &       reference;
      /** \brief  Required typedef for iterator_traits
        */
      typedef  std::random_access_iterator_tag  iterator_category;


      /** \brief  Default constructor
        * \details  Creates an iterator that points nowhere in particular
        */
      ConstVectorDataIterator ();

      /** \brief Copy constructor
        * \param rhs iterator to copy
        * \details  Copies an iterator
        */
      ConstVectorDataIterator ( const ConstVectorDataIterator &rhs );

      /** \brief Copy constructor
        * \param rhs iterator to copy
        * \details  Copies an iterator
        */
      ConstVectorDataIterator ( const VectorDataIterator &rhs );

      /** \brief Constructor from a vector
        * \param p  a (non-reference counted) pointer to the vector being iterated over
        $ \param block  The block of the vector being indexed
        * \param offset the offset of the data being pointed to.
        * \param ndx a vector indexer if p is a sparse vector
        * \details Vector::begin() and Vector::end() call this function to instantiate a new iterator.
        * This method should not be invoked outside of this use.
        */
      ConstVectorDataIterator ( const Vector *p , int block , int offset , VectorIndexer::shared_ptr ndx = VectorIndexer::shared_ptr () );

      /** \brief Constructor used to advance an iterator to the next vector in a multivector
        * \param[in]  i  The iterator to advance
        * \param[in]  p  The multivector pointed to
        * \param[in]  blockNum  The next block to use in the multivector
        */
      ConstVectorDataIterator ( ConstVectorDataIterator i , const MultiVector * const p , size_t blockNum );

      //ConstVectorDataIterator &operator = ( const ConstVectorDataIterator &rhs );

      /** \brief Dereference the iterator
        */
      const double   &operator * ();

      /** \brief Test for equality
        * \details Returns true iff rhs points to the exact same offset of the same vector.
        * It is possible for two iterators to point to the same spot in memory and return false.  This
        * is due to the abstraction of a contiguous block presented by Vector.
        */
      bool       operator == ( const ConstVectorDataIterator &rhs ) const;

      /** \brief Test for inequality
        * \details Returns !(*this == rhs)
        */
      bool       operator != ( const ConstVectorDataIterator &rhs ) const;


      /** \brief Increment the iterator
        * \returns a reference to this iterator
        */
      ConstVectorDataIterator  &operator++ ();

      /** \brief Increment the iterator
        * \returns a copy of this iterator
        */
      ConstVectorDataIterator   operator++ ( int  );

      /** \brief Decrement the iterator
        * \returns a reference to this iterator
        */
      ConstVectorDataIterator  &operator-- ();

      /** \brief Decrement the iterator
        * \returns a copy of this iterator
        */
      ConstVectorDataIterator   operator-- ( int  );

      /** \brief Add a constant to this iterator
        * \param i Offset to move forward
        * \returns A new iterator i steps forward
        * \details  Equivalent to
        * \code
          ConstVectorDataIterator t = *this;
          for ( int k = 0 ; k != i ; k++ )
            ++t;
          \endcode
        * except it is \f$O(1)\f$ in most cases.
        */
      ConstVectorDataIterator   operator + ( int i );

      /** \brief Move this iterator forward
        * \param i Offset to move
        * \returns A reference to this iterator
        * \details  Equivalent to
        * \code
          for ( int k = 0 ; k != i ; k++ )
            opertor ++ ();
          \endcode
        * except it is \f$O(1)\f$ in most cases.
        */
      ConstVectorDataIterator  &operator += ( int i );

      /** \brief Subtract a constant from this iterator
        */
      ConstVectorDataIterator   operator - ( int );

      /** \brief Subtract a constant from this iterator
        */
      ConstVectorDataIterator  &operator -= ( int );

      /** \brief Compute distance between two iterators
        * \param  rhs The other iterator to measure distance to
        * \returns Number of entities between iterators.  If rhs == this,
        * the answer is 0. Else, if rhs is after this, then the answer is
        * positive, otherwise it is negative.
        * \details If rhs is after this, this is equivalent to
        * \code
          int answer = 0;
          ConstVectorDataIterator  t = *this;
          while ( t != rhs )
          {
            t++;
            answer++;
          }
          return answer;
          \endcode
        * except it is \f$O(1)\f$ in most cases.  If rhs is before this,
        * the equivalent code is the same except t and answer are
        * decremented.
        */
      int               operator - ( const ConstVectorDataIterator &rhs ) const;

      /** \brief Return data a distance from the current iterator
        * \param i Offset from this iterator
        * \return The data pointed to at this + i
        * \details  Equivalent to
        * \code
          ConstVectorDataIterator t = *this + i;
          return *t;
          \endcode
        */
      const double             &operator [] ( int i );
  };

}
}

#include "VectorDataIterator.inline.h"

#endif
