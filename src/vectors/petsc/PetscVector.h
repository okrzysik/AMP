#ifndef included_AMP_PetscVector
#define included_AMP_PetscVector

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/petsc/PetscHelpers.h"

#include "petscvec.h"


namespace AMP::LinearAlgebra {


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
class PetscVector final
{
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
        Vec  in_vec = std::dynamic_pointer_cast<PetscVector>(in_petsc_view)->getVec();

        // Perform a PETSc operation
        VecMax ( in_vec , &abs );
        return ans;
      }
      \endcode
      */
    inline Vec &getVec() { return d_Vec; }

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
        Vec  in_vec = std::dynamic_pointer_cast<PetscVector>(in_petsc_view)->getVec();

        // Perform a PETSc operation
        VecMax ( in_vec , &abs );
        return ans;
      }
      \endcode
      */
    inline const Vec &getVec() const { return d_Vec; }

    /**
     *  \brief  If needed, create a PETSc wrapper for AmpVector.  Otherwise, return AmpVector.
     *  \details The function attempts to return a view with the least amount of work.
     *     It will never copy data.  If the vector cannot be wrapped it wll return an error.
     *  \param  AmpVector  a shared pointer to a Vector
     */
    static std::shared_ptr<PetscVector> view( Vector::shared_ptr AmpVector );

    /**
     *  \brief  If needed, create a PETSc wrapper for AmpVector.  Otherwise, return AmpVector.
     *  \details The function attempts to return a view with the least amount of work.
     *     It will never copy data.  If the vector cannot be wrapped it wll return an error.
     *  \param  AmpVector  a shared pointer to a Vector
     */
    static std::shared_ptr<const PetscVector> constView( Vector::const_shared_ptr AmpVector );


public:
    inline Vec &getNativeVec() { return d_Vec; }
    inline const Vec &getNativeVec() const { return d_Vec; }
    inline std::shared_ptr<Vector> getManagedVec() { return d_vector; }
    inline std::shared_ptr<const Vector> getManagedVec() const { return d_vector; }


protected:
    //! Empty constructor
    PetscVector();

    //! Default constructor
    explicit PetscVector( std::shared_ptr<Vector> vec );

protected:
    Vec d_Vec;
    std::shared_ptr<Vector> d_vector;
};


} // namespace AMP::LinearAlgebra

#endif
