#ifndef included_EpetraMatrixOperatorParameters_h
#define included_EpetraMatrixOperatorParameters_h

#include "AMP/matrices/trilinos/EpetraMatrix.h"

namespace AMP {
namespace Operator {
class EpetraMatrixOperatorParameters : public OperatorParameters
{
public:
    Epetra_CrsMatrix *d_Matrix;

    explicit EpetraMatrixOperatorParameters( const AMP::shared_ptr<Database> &db )
        : OperatorParameters( db ), d_Matrix( nullptr )
    {
    }
};
} // namespace Operator
} // namespace AMP

#endif
