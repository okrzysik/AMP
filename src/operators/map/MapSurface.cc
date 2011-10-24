#include "MapSurface.h"
#include "utils/Utilities.h"
#include "ampmesh/MeshUtils.h"


namespace AMP {
namespace Operator {


void MapSurface :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a , const double b )
    {
          AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );

          inpVec = u->subsetVectorForVariable(d_inpVariable);

          boost::shared_ptr<AMP::LinearAlgebra::Vector> nullVec;

          mapMaster->apply(nullVec, inpVec, nullVec, 1.0, 0.0 );
          mapTarget->apply(nullVec, gap1DVec, nullVec, 1.0, 0.0 );
    }


}
}


