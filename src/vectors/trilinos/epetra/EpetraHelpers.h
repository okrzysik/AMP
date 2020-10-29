// This file contains so definitions and wrapper functions for Epetra
#ifndef AMP_EpetraHelpers
#define AMP_EpetraHelpers


#include <memory>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>


namespace AMP::LinearAlgebra {

class Vector;


/********************************************************
 * Get an Epetra vector from an AMP vector               *
 ********************************************************/
std::shared_ptr<Epetra_Vector> getEpetra( std::shared_ptr<Vector> vec );


} // namespace AMP::LinearAlgebra

#endif
