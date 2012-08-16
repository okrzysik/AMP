
#ifndef included_AMP_Petsc_Matrix
#define included_AMP_Petsc_Matrix

//AMP files
#include "matrices/petsc/PetscMatrix.h"

namespace AMP {
namespace LinearAlgebra {

  /** \class NativePetscMatrix
    * \brief  This is a thin wrapper around PETSc Mat
    * \details  As opposed to ManagedPetscMatrix, this is a
    * thin wrapper around a PETSc Mat.
    */
  class NativePetscMatrix : public PetscMatrix 
  {
    protected:
      /** \brief Unused default constructor
        */
      NativePetscMatrix();

      virtual void multiply ( shared_ptr other_op , shared_ptr &result );

    public:

      /** \brief  Construct a matrix from a PETSc Mat.
        * \param[in] m  The Mat to wrap
        * \param[in] dele  Let this class deallocate the Mat
        */
      NativePetscMatrix(Mat m , bool dele = false );
      

      /** \brief Destructor
        */
      virtual ~NativePetscMatrix();

      /** \brief Create a NativePetscMatrix with the same non-zero
        * structure
        * \param[in] m  The matrix to duplicate
        * \return A new matrix with the same non-zero structure
        */
      static Matrix::shared_ptr   duplicateMat ( Mat m );

      /** \brief Copy data from a PETSc Mat
        * \param[in] m  The matrix with the data
        */
      void                        copyFromMat ( Mat m );

      virtual void mult ( Vector::const_shared_ptr in , Vector::shared_ptr out );
      virtual void multTranspose ( Vector::const_shared_ptr in , Vector::shared_ptr out );

      virtual shared_ptr  cloneMatrix () const;


      virtual Vector::shared_ptr getRightVector ();
      virtual Vector::shared_ptr getLeftVector ();
      virtual Discretization::DOFManager::shared_ptr  getRightDOFManager ();
      virtual Discretization::DOFManager::shared_ptr  getLeftDOFManager ();

      virtual size_t  numGlobalRows ();
      virtual size_t  numGlobalColumns ();

      virtual void  scale ( double alpha );
      virtual void  axpy ( double alpha , const Matrix &x );

      virtual void  addValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values );

      virtual void  setValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values );

      virtual void getRowByGlobalID ( int row , std::vector<unsigned int> &cols , std::vector<double> &values ) const;

      virtual void setScalar ( double );
      virtual void setDiagonal ( const Vector::shared_ptr &in );

      virtual void makeConsistent ();
      virtual Vector::shared_ptr  extractDiagonal ( Vector::shared_ptr p = Vector::shared_ptr() );
      virtual double L1Norm () const;


    private:

      /** \brief  Communicator for the matrix
        */
      AMP_MPI d_comm;


  };

}
}

#include "NativePetscMatrix.inline.h"

#endif



