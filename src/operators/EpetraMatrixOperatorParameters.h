#ifndef included_EpetraMatrixOperatorParameters_h
#define included_EpetraMatrixOperatorParameters_h

#include "matrices/trilinos/EpetraMatrix.h"

namespace AMP {
namespace Operator {
  class EpetraMatrixOperatorParameters : public OperatorParameters
  {
    public:
      Epetra_CrsMatrix    *d_Matrix;

      EpetraMatrixOperatorParameters ( const boost::shared_ptr<Database> &db ) : OperatorParameters ( db ) {}
  };

}
}

#endif
