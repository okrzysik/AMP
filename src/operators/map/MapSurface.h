#ifndef included_AMP_MapSurface
#define included_AMP_MapSurface

#include "boost/shared_ptr.hpp"

#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "operators/map/Map3Dto1D.h"
#include "operators/map/Map1Dto3D.h"
#include "vectors/Vector.h"
#include "vectors/Variable.h"
#include "vectors/SimpleVector.h"

#include <string>

#ifdef DEBUG_CHECK_ASSERTIONS
#include <cassert>
#endif


namespace AMP {
namespace Operator {
 
class MapSurface : public MapOperator
{
public :

    MapSurface(const boost::shared_ptr<OperatorParameters> & params);
    virtual ~MapSurface() { }

    virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

    boost::shared_ptr<AMP::LinearAlgebra::Vector> getBoundaryVector(const AMP::LinearAlgebra::Vector::shared_ptr &u) {
        return (u->subsetVectorForVariable(d_outVariable)) ;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
        return d_inpVariable;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_outVariable;
    }

    void setVector(AMP::LinearAlgebra::Vector::shared_ptr scratchVec)
    {
        outVec = scratchVec;
        mapTarget->setVector(outVec);
    }

    /*
    boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& )
    {
        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));

        boost::shared_ptr<MapOperatorParameters> outParams(new MapOperatorParameters(tmp_db));

        outParams->d_BoundaryId = (mapMasterParams->d_db)->getInteger("BoundaryId");

        return outParams;

    } */

protected :

    boost::shared_ptr<AMP::LinearAlgebra::Vector> gap1DVec; 
    boost::shared_ptr<AMP::LinearAlgebra::Variable> gapVariable;

    boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;
    boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

    AMP::LinearAlgebra::Vector::shared_ptr inpVec; 
    AMP::LinearAlgebra::Vector::shared_ptr outVec;

private :

    boost::shared_ptr<Map3Dto1D> mapMaster;
    boost::shared_ptr<MapOperatorParameters> mapMasterParams;
    boost::shared_ptr<Map1Dto3D> mapTarget;
    boost::shared_ptr<MapOperatorParameters> mapTargetParams;
};


}
}

#endif
