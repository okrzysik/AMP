#ifndef included_AMP_PetscVector
#define included_AMP_PetscVector

#include "vectors/DataChangeListener.h"
#include "vectors/Vector.h"
#include "vectors/petsc/PetscHelpers.h"

#include "petscvec.h"


namespace AMP {
namespace LinearAlgebra {

/**  \class  PetscVectorParameters
 *  \brief  Parameters class to construct a PETSc Vector
 *
 *  Since the constructor of PetscVector is protected, this class is of little
 *  use to any but Backplane developers.  Also, since PetscVector
 *  is almost trivial, so is the parameters class.
 */
class PetscVectorParameters
{
};


/**
 *  \class  PetscVector
 *  \brief  PetscVector is a bridge between AMP::LinearAlgebra::Vector and
 *  the PETSc Vec data structure.
 *
 *  A PetscVector has a Vec data structure.  Given an
 *  AMP::LinearAlgebra::Vector, this class can create a PETSc view without
 *  copying the data.  As such, this class serves three purposes:
 *  -# Provides a PETSc Vec for derived classes to use, fill, manage, etc.
 *  -# Provides an interface for accessing this PETSc Vec independent of derived classes
 *  -# Provides a static method for creating a PETSc view of an AMP Vector.
 *
 */

class PetscVector : public DataChangeListener
{
private:
    PetscRandom *d_PetscRandom;

protected:
    /**
     *  \brief  PETSc Vec holding data in the vector
     *
     *  Whether created with VecCreate (called Native) or
     *  a view of an AMP:Vector (called Managed), this pointer
     *  is what is used when calling the PETSc Vec interface
     */
    Vec d_petscVec;

    /**
     *  \brief Retrieve a valide PETSc random context
     *  \param  comm  The communicator to create the context around.
     *
     *  If PetscRandomCreate has not been called, this will
     *  call it.
     */
    PetscRandom &getPetscRandom( AMP_MPI comm = AMP_MPI( AMP_COMM_NULL ) );

    /**
     *  \brief  Swap the underlying PETSc Vec with another
     *  AMP::LinearAlgebra::Vector.
     */
    void swapPetscVec( PetscVector &rhs ) { std::swap( d_petscVec, rhs.d_petscVec ); }

    /**
     *  \brief  Construct a PetscVector
     *
     *  This can only be called by a derived class or the static function below.  There is
     *  no need to create this vector directly since it is virtual.
     */
    PetscVector();

public:
    /**
     *  \brief  Destructor
     */
    virtual ~PetscVector();

    /**
      *  \brief  Obtain PETSc Vec for use in PETSc routines
      *
      *  This function is used to get a PETSc vector.  The following idiom
      *  should be used since it fails gracefully.  In this function,
      *  a view may be created before the Vec is extracted
      *\code
      double  DoPETScMax ( Vector::shared_ptr  &in )
      {
        double   ans;
        // Create a PETSc Vec if necessary
        Vector::shared_ptr in_petsc_view = PetscVector::view( in );

        // Extract the Vec
        Vec  in_vec = dynamic_pointer_cast<PetscVector>(in_petsc_view)->getVec();

        // Perform a PETSc operation
        VecMax ( in_vec , &abs );
        return ans;
      }
      \endcode
      */
    virtual Vec &getVec();

    /**
      *  \brief  Obtain PETSc Vec for use in PETSc routines
      *
      *  This function is used to get a PETSc vector.  The following idiom
      *  should be used since it fails gracefully.  In this function,
      *  a view may be created before the Vec is extracted
      *\code
      double  DoPETScMax ( Vector::shared_ptr  &in )
      {
        double   ans;
        // Create a PETSc Vec if necessary
        Vector::shared_ptr in_petsc_view = PetscVector::view( in );

        // Extract the Vec
        Vec  in_vec = dynamic_pointer_cast<PetscVector>(in_petsc_view)->getVec();

        // Perform a PETSc operation
        VecMax ( in_vec , &abs );
        return ans;
      }
      \endcode
      */
    virtual Vec getVec() const;

    /**
     *  \brief  If needed, create a PETSc wrapper for AmpVector.  Otherwise, return AmpVector.
     *  \param  AmpVector  a shared pointer to a Vector
     *
     *  \details The function attempts to return a view with the least amount of work.
     *  IT WILL NEVER COPY DATA.
     *  - If AmpVector is already a PetscVector, it is returned.
     *  - Else, if AmpVector is a ManagedVector, it is wrapped and returned
     *  - Else, if AmpVector can be used as a VectorEngine, a new ManagedVector
     *    is created and returned.
     *  Otherwise, this function will throw an error.
     */
    static Vector::shared_ptr view( Vector::shared_ptr AmpVector );

    /**
     *  \brief  If needed, create a const PETSc wrapper for AmpVector.  Otherwise, return
     * AmpVector.
     *  \param  AmpVector  a shared pointer to a Vector
     *
     *  \details The function attempts to return a view with the least amount of work.
     *  IT WILL NEVER COPY DATA.
     *  - If AmpVector is already a PetscVector, it is returned.
     *  - Else, if AmpVector is a ManagedVector, it is wrapped and returned
     *  - Else, if AmpVector can be used as a VectorEngine, a new ManagedVector
     *    is created and returned.
     *  Otherwise, this function will throw an error.
     */
    static Vector::const_shared_ptr constView( Vector::const_shared_ptr AmpVector );


    /**
     *  \brief  Check if petsc is holding a view that might prevent us from deleting the vector
     *  \details This function checks if petsc might be holding a view of the vector
     *    that would prevent us from deleting the vector.  This function returns false
     *    if we can safely delete the vector.
     */
    virtual bool petscHoldsView() const = 0;


    virtual void dataChanged();
};
} // namespace LinearAlgebra
} // namespace AMP

#include "PetscVector.inline.h"

#endif
