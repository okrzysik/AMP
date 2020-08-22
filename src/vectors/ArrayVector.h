#ifndef included_AMP_ArrayVector
#define included_AMP_ArrayVector

#include <string>

#include "AMP/utils/Array.h"
#include "AMP/utils/FunctionTable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/ArrayVectorData.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"
#include "AMP/vectors/operations/VectorOperationsDefault.hpp"


namespace AMP {

namespace LinearAlgebra {

/** \brief A core-local vector
 * \details This is a Vector that implements the Vector interface for a std::vector<double>.
 */
template<typename T, typename FUN = FunctionTable, typename Allocator = std::allocator<T>>
class ArrayVector : public Vector
{
private:

    ArrayVector();
    explicit ArrayVector( const ArrayVector & );

public:
    explicit ArrayVector( std::shared_ptr<ArrayVectorData<T, FUN, Allocator>> data );

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

    std::string type() const override { return "ArrayVector"; }

    using Vector::cloneVector;
    Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector &other ) override;
    /**
     * \brief This method is used to implement the assemble interface
     * of PETSc.
     * \details  This method is empty except for instantiations of NativePetscVector
     */
    void assemble() override { AMP_ERROR( "Not implemented" ); }

    //! resize the ArrayVector and reset the internal data structures
    void resize( const std::vector<size_t> &localDims );

};

} // namespace LinearAlgebra
} // namespace AMP

#include "ArrayVector.hpp"

#endif
