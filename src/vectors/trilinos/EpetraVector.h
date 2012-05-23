#ifndef included_AMP_EpetraVector
#define included_AMP_EpetraVector

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "vectors/Vector.h"

namespace AMP {
namespace LinearAlgebra {

  /**  \class  EpetraVectorParameters
    *  \brief  Parameters class to construct an EpetraVector
    */
  class EpetraVectorParameters : public VectorParameters {};


  /**  \class EpetraVector
    *  \brief A class that manages an Epetra_Vector
    *
    *  \see PetscVector
    *  \see SundialsVector
    *
    *  \details An EpetraVector presents an Epetra_Vector class.  Given an
    *  AMP::LinearAlgebra::Vector, this class can create an Epetra view without
    *  copying the data.  As such, this class serves three purposes:
    *  -# Provides an Epetra_Vector for derived classes to use, fill, manage, etc.
    *  -# Provides an interface for accessing this Epetra_Vector independent of base or derived classes
    *  -# Provides a static method for creating an Epetra_Vector view of an AMP Vector.
    *
    *  This allows the Castable class to be used to verify correctness of code.  For instance,
    *  given a Vector shared pointer, it is possible to get the Epetra_Vector safely thusly
    \code
       Vector::shared_ptr  vector;
       vector->castTo<EpetraVector>().getEpetra_Vector();
    \endcode
    *  The castTo ensures that the Epetra_Vector exists.  An Epetra_Vector expects data
    *  in a contiguous block.  If the vector does not contain a contiuous block of data, then
    *  an Epetra_Vector cannot be created.
    */
  class EpetraVector 
  {
    protected:
      /**  \brief Constructor
        */
      EpetraVector ();

    public:
      /**  \brief Destructor
        */
      virtual ~EpetraVector();

      /**
        *  \brief  Obtain Epetra_Vector for use in Trilinos routines
        *
        *  \details This function is used to get a Epetra vector.  The 
        *  following idiom should be used since it fails gracefully.  In 
        *  this function, a view may be created before the Vec is extracted
        *  \see view()
        *  \returns Epetra_Vector wrapper for this vector
        *\code
        double  DoEpetraMax ( Vector::shared_ptr  &in )
        {
          double   ans;

         // Create an Epetra_Vector, if necessary
          Vector::shared_ptr  in_epetra_view = EpetraVector::view ( in );  

         // Extract the Epetra_Vector
          Epetra_Vector &in_vec = in_epetra_view->castTo<EpetraVector>().getEpetra_Vector ();    

         // Perform an Epetra_Vector operation
          in_vec.MaxValue ( &abs );
          return ans;
        }
        \endcode
        */
      virtual       Epetra_Vector &getEpetra_Vector ()        = 0;

      /**
        *  \brief  Obtain Epetra_Vector for use in Trilinos routines
        *
        *  \see view()
        *  \returns Epetra_Vector wrapper for this vector
        *  \details This function is used to get a Epetra vector.  The 
        *  following idiom should be used since it fails gracefully.  In 
        *  this function, a view may be created before the Epetra_Vector is extracted
        *\code
        double  DoEpetraMax ( Vector::shared_ptr  &in )
        {
          double   ans;

         // Create an Epetra_Vector, if necessary
          Vector::shared_ptr  in_epetra_view = EpetraVector::view ( in );  

         // Extract the Epetra_Vector
          Epetra_Vector &in_vec = in_epetra_view->castTo<EpetraVector>().getEpetra_Vector ();    

         // Perform an Epetra_Vector operation
          in_vec.MaxValue ( &abs );
          return ans;
        }
        \endcode
        */
      virtual const Epetra_Vector &getEpetra_Vector () const  = 0;

      /**
        *  \brief  Obtain a view of a vector with an Epetra_Vector wrapper
        *  \param[in] vec  The vector to get an Epetra_Vector view of.
        *  \return A Vector::shared_ptr guaranteed to have an Epetra_Vector
        *   wrapper available through the getEpetra_Vector() interface.
        *  \see getEpetra_Vector()
        *  \details  If the vector has an Epetra_Vector wrapper already
        *  created, it is returned.  Otherwise, it will try to create an
        *  Epetra_Vector wrapper around the Vector.  If it fails, an
        *  exception is thrown.
        */
      static Vector::shared_ptr  view ( Vector::shared_ptr vec );

      /**
        *  \brief  Obtain a view of a vector with an Epetra_Vector wrapper
        *  \param[in] vec The vector to get an Epetra_Vector view of.
        *  \return A Vector::shared_ptr guaranteed to have an Epetra_Vector
        *   wrapper available through the getEpetra_Vector() interface.
        *  \see getEpetra_Vector()
        *  \details  If the vector has an Epetra_Vector wrapper already
        *  created, it is returned.  Otherwise, it will try to create an
        *  Epetra_Vector wrapper around the Vector.  If it fails, an
        *  exception is thrown.
        */
      static const Vector::shared_ptr  constView ( const Vector::shared_ptr vec );

  };

}
}


#endif
