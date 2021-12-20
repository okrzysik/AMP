#ifndef included_AMP_PetscMatrix
#define included_AMP_PetscMatrix

#include "AMP/matrices/Matrix.h"
#include "AMP/vectors/petsc/PetscHelpers.h"

#include "petscmat.h"


namespace AMP {
namespace LinearAlgebra {


/**
 *  \class  PetscMatrix
 *  \brief  PetscMatrix is a bridge between AMP::LinearAlgebra::Matrix and
 *  the PETSc Vec data structure.
 *
 *  A PetscMatrix has a Mat data structure.  Given an
 *  AMP::LinearAlgebra::Matrix, this class can create a PETSc view without
 *  copying the data.  As such, this class serves three purposes:
 *  -# Provides a PETSc Mat for derived classes to use, fill, manage, etc.
 *  -# Provides an interface for accessing this PETSc Mat independent of derived classes
 *  -# Provides a static method for creating a PETSc view of an AMP Matrix.
 *
 */
class PetscMatrix final
{
public:
    /**
     *  \brief  Destructor
     */
    virtual ~PetscMatrix();

    /**
     *  \brief  Obtain PETSc Mat for use in PETSc routines
     *
     *  This function is used to get a PETSc matrix.  The following idiom
     *  should be used since it fails gracefully.  In this function,
     *  a view may be created before the Mat is extracted
     */
    inline Mat &getMat() { return d_Mat; }

    /**
     *  \brief  Obtain PETSc Mat for use in PETSc routines
     *
     *  This function is used to get a PETSc matrix.  The following idiom
     *  should be used since it fails gracefully.  In this function,
     *  a view may be created before the Mat is extracted
     */
    inline const Mat &getMat() const { return d_Mat; }

    /**
     *  \brief  If needed, create a PETSc wrapper for AmpMatrix.  Otherwise, return AmpMatrix.
     *  \details The function attempts to return a view with the least amount of work.
     *     It will never copy data.  If the matrix cannot be wrapped it wll return an error.
     *  \param  AmpMatrix  a shared pointer to a Matrix
     */
    static std::shared_ptr<PetscMatrix> view( std::shared_ptr<Matrix> AmpMatrix );

    /**
     *  \brief  If needed, create a PETSc wrapper for AmpMatrix.  Otherwise, return AmpMatrix.
     *  \details The function attempts to return a view with the least amount of work.
     *     It will never copy data.  If the matrix cannot be wrapped it wll return an error.
     *  \param  AmpMatrix  a shared pointer to a Matrix
     */
    static std::shared_ptr<const PetscMatrix> constView( std::shared_ptr<const Matrix> AmpMatrix );


public:
    inline Mat &getNativeMat() { return d_Mat; }
    inline const Mat &getNativeMat() const { return d_Mat; }
    inline std::shared_ptr<Matrix> getManagedMat() { return d_matrix; }
    inline std::shared_ptr<const Matrix> getManagedMat() const { return d_matrix; }


protected:
    //! Empty constructor
    PetscMatrix();

    //! Default constructor
    explicit PetscMatrix( std::shared_ptr<Matrix> mat );

protected:
    Mat d_Mat;
    std::shared_ptr<Matrix> d_matrix;
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
