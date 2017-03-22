#ifndef included_AMP_ThyraVector
#define included_AMP_ThyraVector

// AMP includes
#include "vectors/DataChangePassThrough.h"
#include "vectors/Vector.h"

// Thyra includes
DISABLE_WARNINGS
#include "Thyra_VectorDefaultBase_decl.hpp"
ENABLE_WARNINGS


namespace AMP {
namespace LinearAlgebra {


/**
  *  \class  ThyraVector
  *  \brief  ThyraVector is a bridge between AMP::LinearAlgebra::Vector and
  *  the Thyra::VectorDefaultBase data structure.
  *
  *  A ThyraVector has a Thyra::VectorDefaultBase data structure.  Given an
  *  AMP::LinearAlgebra::Vector, this class can create a Thyra view without
  *  copying the data.  As such, this class serves three purposes:
  *  -# Provides a Thyra Vector for derived classes to use, fill, manage, etc.
  *  -# Provides an interface for accessing this Thyra Vector independent of derived classes
  *  -# Provides a static method for creating a Thyra view of an AMP Vector.
  *
  *  This allows the Castable class to be used to verify correctness of code.  For instance,
  *  given a Vector shared pointer, it is possible to get the Thyra Vector safely thusly
  \code
     Vector::shared_ptr  vector;
     vector->castTo<ThyraVector>().getVec();
  \endcode
  *  The castTo ensures that the Thyra Vector exists.  If the Vec does not exist, the castTo will
  *  throw an error.  If, on the other hand, you have an arbitrary AMP Vector that may
  *  or may not have a Thyra Vector associated with it, you can use the static members
  *  ThyraVector::view or ThyraVector::constView
  *  to create the Thyra Vector if it doesn't already exist, give or take some edge cases.
  */
class ThyraVector : public DataChangePassThrough
{
public:
    //!  Destructor
    virtual ~ThyraVector();

    /**
      *  \brief  Obtain Thyra Vector for use in Thyra routines
      *
      *  This function is used to get a Thyra vector.  The following idiom
      *  should be used since it fails gracefully.  In this function,
      *  a view may be created before the Vec is extracted
      */
    virtual Teuchos::RCP<Thyra::VectorBase<double>> getVec();

    /**
      *  \brief  Obtain Thyra Vector for use in Thyra routines
      *
      *  This function is used to get a Thyra vector.  The following idiom
      *  should be used since it fails gracefully.  In this function,
      *  a view may be created before the Vec is extracted
      */
    virtual Teuchos::RCP<const Thyra::VectorBase<double>> getVec() const;

    /**
      *  \brief  If needed, create a Thyra wrapper for AmpVector.  Otherwise, return AmpVector.
      *  \param  AmpVector  a shared pointer to a Vector
      *
      *  \details The function attempts to return a view with the least amount of work.
      *  IT WILL NEVER COPY DATA.
      *  - If AmpVector is already a ThyraVector, it is returned.
      *  - Else, if AmpVector is a ManagedVector, it is wrapped and returned
      *  - Else, if AmpVector can be used as a VectorEngine, a new ManagedVector
      *    is created and returned.
      *  Otherwise, this function will throw an error.
      */
    static Vector::shared_ptr view( Vector::shared_ptr AmpVector );

    /**
      *  \brief  If needed, create a const Thyra wrapper for AmpVector.  Otherwise, return
     * AmpVector.
      *  \param  AmpVector  a shared pointer to a Vector
      *
      *  \details The function attempts to return a view with the least amount of work.
      *  IT WILL NEVER COPY DATA.
      *  - If AmpVector is already a ThyraVector, it is returned.
      *  - Else, if AmpVector is a ManagedVector, it is wrapped and returned
      *  - Else, if AmpVector can be used as a VectorEngine, a new ManagedVector
      *    is created and returned.
      *  Otherwise, this function will throw an error.
      */
    static Vector::const_shared_ptr constView( Vector::const_shared_ptr AmpVector );


    //! Return an AMP Vector from the Thyra::VectorBase
    static AMP::LinearAlgebra::Vector::shared_ptr view( Thyra::VectorBase<double> *vec );

    //! Return an AMP Vector from the Thyra::VectorBase
    static AMP::LinearAlgebra::Vector::const_shared_ptr
    constView( const Thyra::VectorBase<double> *vec );

protected:
    /**
      *  \brief  Thyra Vector holding data in the vector
      *
      *  Whether created with VecCreate (called Native) or
      *  a view of an AMP:Vector (called Managed), this pointer
      *  is what is used when calling the Thyra Vector interface
      */
    Teuchos::RCP<Thyra::VectorBase<double>> d_thyraVec;

    /**
      *  \brief  Swap the underlying PETSc Vec with another
      *  AMP::LinearAlgebra::Vector.
      */
    void swapThyraVec( ThyraVector &rhs ) { std::swap( d_thyraVec, rhs.d_thyraVec ); }

    /**
      *  \brief  Construct a PetscVector
      *
      *  This can only be called by a derived class or the static function below.  There is
      *  no need to create this vector directly since it is virtual.
      */
    ThyraVector();
};
}
}


#endif
