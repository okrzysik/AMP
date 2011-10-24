#include "ManagedEpetraVector.h"
#include "utils/Utilities.h"


namespace AMP {
namespace LinearAlgebra {



  ManagedEpetraVector::ManagedEpetraVector ( VectorParameters::shared_ptr params ) : ManagedVector ( params ) ,
                                                         EpetraVector ()
  {
  }


  ManagedEpetraVector::ManagedEpetraVector ( shared_ptr alias ) : ManagedVector ( alias ) ,
                                                         EpetraVector ()
  {
  }

  void ManagedEpetraVector::copyVector(const Vector &vec)
  {
    // there must be a more sensible way of doing this but I can't find the documentation - BP
    if ( vec.isA<ManagedEpetraVector>() )
    {
      double scale = 1.0;
      getEpetra_Vector().Scale(scale, vec.castTo<EpetraVector>().getEpetra_Vector());
    }
    else
    {
      Vector::copyVector ( vec );
    }
  }
  

}
}

