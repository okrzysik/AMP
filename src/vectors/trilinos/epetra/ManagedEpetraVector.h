#ifndef included_AMP_ManagedEpetraVector
#define included_AMP_ManagedEpetraVector

#include "EpetraVector.h"
#include "EpetraVectorEngine.h"
#include "vectors/ManagedVector.h"
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


class ManagedEpetraVector : public ManagedVector, public EpetraVector
{
public:
    /** \brief Create a ManagedEpetraVector from a set of parameters
      * \param[in] params  A VectorParameters class used to construct this vector
      */
    explicit ManagedEpetraVector( VectorParameters::shared_ptr params );

    /** \brief Create a view of a vector
      * \param[in] alias  Vector to view
      */
    explicit ManagedEpetraVector( Vector::shared_ptr alias );


    // These methods are adequately documented in a base class
    virtual std::string type() const override;

    using Vector::cloneVector;
    virtual Vector::shared_ptr cloneVector( const Variable::shared_ptr var ) const override;
    virtual void copy( const VectorOperations &vec ) override;

    Epetra_Vector &getEpetra_Vector() override;
    const Epetra_Vector &getEpetra_Vector() const override;
    virtual void assemble() override;

protected:
    virtual ManagedVector *getNewRawPtr() const override;
};
}
}

#include "ManagedEpetraVector.inline.h"

#endif
