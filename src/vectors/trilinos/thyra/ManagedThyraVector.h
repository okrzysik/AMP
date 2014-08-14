#ifndef included_AMP_ManagedThyraVector
#define included_AMP_ManagedThyraVector

// AMP includes
#include "vectors/ManagedVector.h"
#include "vectors/trilinos/thyra/ThyraVector.h"



namespace AMP {
namespace LinearAlgebra {


//! ManagedThyraVectorParameters 
typedef ManagedVectorParameters ManagedThyraVectorParameters;


/** \class ManagedThyraVector
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
class ManagedThyraVector : public ManagedVector, public ThyraVector
{
public:

    /** \brief Create a ManagedThyraVector from a set of parameters
      * \param[in] params  A VectorParameters class used to construct this vector
      */
    ManagedThyraVector( VectorParameters::shared_ptr  params );

    /** \brief Create a view of a vector
      * \param[in] alias  Vector to view
      */
    ManagedThyraVector( Vector::shared_ptr  alias );

    //! Destructor
    virtual ~ManagedThyraVector();

    // These methods are adequately documented in a base class
    virtual std::string type() const;

    using Vector::cloneVector;
    virtual Vector::shared_ptr  cloneVector( const Variable::shared_ptr var ) const;
    virtual void copyVector( Vector::const_shared_ptr vec );
    virtual void assemble() {}
    ManagedVector* getNewRawPtr() const;

protected:


};

}
}

#endif
