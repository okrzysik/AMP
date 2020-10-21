#ifndef included_AMP_ManagedEpetraVector
#define included_AMP_ManagedEpetraVector

#include "EpetraVector.h"
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

namespace AMP {
namespace LinearAlgebra {

/** \class ManagedEpetraVector
 * \brief  Vector capable of returning an Epetra_Vector from a ManagedVector
 * \details  One popular package of numerical algorithms is Trilinos.  Many
 * of the pieces of Trilinos rely on a basic vector called an Epetra vector.
 * This class will provide an Epetra_Vector view on a vector.  This should
 * not be used explicitly.  Rather, the EpetraVector interface provides
 * the getEpetra_Vector interface.  In this case, the class will populate
 * an Epetra_Vector with a view of an array that the ManagedVector has.
 *
 * This class is the class returned by EpetraVector::view() and
 * EpetraVector::constView().
 *
 * \see EpetraVector
 */
class ManagedEpetraVector : public Vector, public EpetraVector
{
public:
    /** \brief Create a view of a vector
     * \param[in] alias  Vector to view
     */
    explicit ManagedEpetraVector( Vector::shared_ptr alias );

    virtual ~ManagedEpetraVector();

    // These methods are adequately documented in a base class
    std::string type() const override
    {
        return "Managed Epetra Vector" + d_VectorData->VectorDataName();
    }

    std::unique_ptr<Vector> rawClone( const Variable::shared_ptr var ) const override;
    void swapVectors( Vector &other ) override;

    using Vector::copyVector;
    void copyVector( Vector::const_shared_ptr vec ) override;

    Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name ) override;
    Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const override;

    Epetra_Vector &getEpetra_Vector() override;
    const Epetra_Vector &getEpetra_Vector() const override;

    std::shared_ptr<Vector> getManagedVec() override { return shared_from_this(); }
    std::shared_ptr<const Vector> getManagedVec() const override { return shared_from_this(); }
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
