#ifndef included_AMP_ArrayVector
#define included_AMP_ArrayVector

#include "SimpleVector.h"
#include "utils/Array.h"

namespace AMP {
namespace LinearAlgebra {

/** \brief A core-local vector
  * \details This is a Vector that implements the Vector interface for a std::vector<double>.
  */
template <typename T>
class ArrayVector : public SimpleVector<T>
{
private:
    AMP::Array<T> d_array;

    ArrayVector();
    ArrayVector( const ArrayVector & );

public:
    /** \brief    Create a ArrayVector
      * \details  This is the factory method for the ArrayVector.  It returns the shared pointer
      * to be used in the code
      * \param    localSize  The number of elements in the vector on this processor
      * \param    var The variable associated with the new vector
      */
    static Vector::shared_ptr create( const std::vector<size_t> &localSize,
                                      Variable::shared_ptr var );

    /** \brief    Create a ArrayVector
      * \details  This is the factory method for the ArrayVector.  It returns the shared pointer
      * to be used in the code
      * \param    localSize  The number of elements in the vector on this processor
      * \param    var The variable associated with the new vector
      * \param    comm The variable associated with the new vector
      */
    static Vector::shared_ptr
    create( const std::vector<size_t> &localSize, Variable::shared_ptr var, AMP_MPI comm );

    /** \brief    Create a ArrayVector
      * \details  This is the factory method for the ArrayVector.  It returns the shared pointer
      * to be used in the code that spans a comm and contains ghost values.
      * \param    var The variable associated with the new vector
      * \param    DOFs The DOFManager
      * \param    commlist The communication list
      */
    static Vector::shared_ptr create( Variable::shared_ptr var,
                                      AMP::Discretization::DOFManager::shared_ptr DOFs,
                                      AMP::LinearAlgebra::CommunicationList::shared_ptr commlist );

    /** \brief  Destructor
      */
    virtual ~ArrayVector() {}

    virtual std::string type() const override { return "ArrayVector"; }

    using Vector::cloneVector;
    virtual Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const override;
    using Vector::copyVector;
    virtual void copyVector( Vector::const_shared_ptr src_vec ) override;
    virtual void swapVectors( Vector &other ) override;
    virtual void aliasVector( Vector &other ) override;

    //! return a non-const reference to the internal data container
    Array<T> &getArray( void ) { return d_array; }

    //! return a const reference to the internal data container
    const Array<T> &getArray( void ) const { return d_array; }
};
}
}

#include "ArrayVector.hpp"

#endif
