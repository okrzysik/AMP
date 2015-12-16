#ifndef included_AMP_ManagedPetscVector
#define included_AMP_ManagedPetscVector

#include "vectors/ManagedVector.h"
#include "vectors/petsc/PetscVector.h"
#include "vectors/ExternalVectorDeleter.h"


extern "C"{
#include "petsc.h"
#include "petscvec.h"
}



namespace AMP {
namespace LinearAlgebra {


/** \typedef ManagedPetscVectorParameters
* \brief   Requirements for constructing a ManagedPetscVector
*/
typedef ManagedVectorParameters ManagedPetscVectorParameters;


/** \class ManagedPetscVector
 * \brief A class that provides a PETSc vector interfaced to a ManagedVector.
 * \details  This class provides a PETSc Vec specially configured to call 
 * through ManagedVector.
 *
 * In general, the end user will not have to use this class.  This class is 
 * returned by PetscVector::view() and PetscVector::constView(); 
 *
 * \see PetscVector
 */
class ManagedPetscVector: public ManagedVector, public PetscVector
{
private:
      bool     d_bMadeWithPetscDuplicate;

protected:
      /** \brief  Convenience typedef fpr a ManagedVector
        */
      typedef  ManagedVector   ParentVector;

      /** \brief Populate PETSc data structures with functions that call
        * back into the Vector interface
        */
      void initPetsc();


public:


      /** \brief Construct a new ManagedPetscVector given a set of parameters
        * \param[in] params  The parameters describing the new vector
        */
      explicit ManagedPetscVector( VectorParameters::shared_ptr params );

      /** \brief Construct a view of another vector
        * \param[in] alias The vector to view
        */
      explicit ManagedPetscVector ( Vector::shared_ptr alias );

      /** \brief Method to create a duplicate of this vector for VecDuplicate
        * \return Raw pointer to a new vector.  This does not copy data
        */
      ManagedPetscVector  *petscDuplicate ();

      /** \brief Identifies whether this vector was created through the 
        * VecDuplicate interface
        * \return true if constructed with VecDuplicate.  False otherwise
        */
      bool          constructedWithPetscDuplicate ();

      /** \brief Destructor
        */
      virtual ~ManagedPetscVector();

      /** \brief Create an exact clone of this vector.
        * \return A raw pointer to a clone of this vector
        */
      ManagedPetscVector  *rawClone () const;


      /** \brief Copy data from a PETSc Vec to an AMP Vector
        * \param[out] dest  Vector to copy to
        * \param[in]  src   Vec to copy from
        */
      static void  copyFromPetscVec ( Vector &dest , Vec src );

      /** \brief Create data from a PETSc Vec, but do not copy data
        * \param[in]  src   Vec to copy from
        * \param[in]  comm  The AMP_MPI to create the AMP Vector on.
        * \return A new AMP vector with identical data distribution to the 
        * PETSc Vec.
        */
      static Vector::shared_ptr  createFromPetscVec ( Vec src , AMP_MPI &comm );
      

      // These are adequately documented in a base class.
public:
      virtual void  swapVectors ( Vector &other );
      using Vector::cloneVector;
      virtual Vector::shared_ptr  cloneVector ( const Variable::shared_ptr p ) const;

      virtual std::string type() const;
      void  assemble();
 
      virtual bool petscHoldsView() const;

protected:
      virtual ManagedVector *getNewRawPtr () const;

};


}
}

#include "ManagedPetscVector.inline.h"

#endif

