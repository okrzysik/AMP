#ifndef included_AMP_NativeThyraVector
#define included_AMP_NativeThyraVector

#include "AMP/vectors/Vector.h"

namespace AMP {
namespace LinearAlgebra {


/** \class NativeThyraVector
 * \brief An AMP Vector that uses Thyra for parallel data management, linear algebra,
 * etc.
 * \details  This is an AMP wrapper to Thyra.  This is different from ManagedThyraVector
 * in that this class does not replace calls to Vec*.  Rather, it wraps these calls.
 * This class is used when Thyra is chosen as the default linear algebra engine.
 *
 * This class is not to be used directly, just through base class interfaces.
 * \see ThyraVector
 * \see ManagedThyraVector
 */
class NativeThyraVector : public Vector
{
public:
    /** \brief Construct a wrapper for a Thyra Vec from a set of parameters
     * \param[in] params The parameters describing the Vec
     */
    explicit NativeThyraVector( VectorParameters::shared_ptr params );
    explicit NativeThyraVector( std::shared_ptr<VectorData> data );

    //! Destructor
    virtual ~NativeThyraVector();

    //! Overloaded functions
    std::string type() const override { return "Native Thyra Vector"; }
    using Vector::cloneVector;
    Vector::shared_ptr cloneVector( const Variable::shared_ptr ) const override;
    void swapVectors( Vector &other ) override;

protected:
    //! Empty constructor.
    NativeThyraVector();
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
