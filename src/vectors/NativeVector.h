#ifndef included_AMP_NativeVector
#define included_AMP_NativeVector

#include "Vector.h"

namespace AMP {
namespace LinearAlgebra {


//! A typedef for NativeVectorParameters
typedef VectorParameters NativeVectorParameters;


/** \brief A light wrapper around a vector from a TPL
 * \details  If the memory management of a Vector is controlled by a TPL, then
 * it is a NativeVector.
 */

class NativeVector : public Vector
{
public:
    //! Convenience typedef
    typedef AMP::shared_ptr<NativeVectorParameters> parameters_ptr;
    //! brief Convenience typedef
    typedef NativeVectorParameters parameters;


    /** \brief Construct a NativeVector from a set of parameters
     * \param[in] params  The parameters describing the NativeVector
     */
    explicit NativeVector( parameters_ptr params );

    //! Virtual Destuctor
    virtual ~NativeVector() {}

    /** \brief Get a ManagedVector copy of this vector (copies data)
     * \param[in] comm The communicator to build the vector on
     * \return A ManagedVector copy of this vector.
     */
    virtual Vector::shared_ptr getManagedVectorCopy( AMP_MPI comm ) = 0;

    /** \brief Get a ManagedVector duplicate of this vector (does not copy data)
     * \param[in] comm The communicator to build the vector on
     * \return A ManagedVector duplicate (data not copied) of this vector.
     */
    virtual Vector::shared_ptr getManagedVectorDuplicate( AMP_MPI comm ) = 0;

protected:
    //! Empty constructor.
    NativeVector();
};


} // namespace LinearAlgebra
} // namespace AMP

#include "NativeVector.inline.h"

#endif
