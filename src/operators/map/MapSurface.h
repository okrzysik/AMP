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

      MapSurface(const boost::shared_ptr<OperatorParameters> & params) 
        : MapOperator (params)
      {

        boost::shared_ptr<MapOperatorParameters> myparams = 
          boost::dynamic_pointer_cast<MapOperatorParameters>(params);

        //Construct for Map3Dto1D, Map1Dto3D
        boost::shared_ptr<AMP::Database> mapMaster_db = myparams->d_db->getDatabase("Map3Dto1D");
        mapMasterParams.reset(new AMP::Operator::MapOperatorParameters(mapMaster_db));
        mapMasterParams->d_MeshAdapter    = myparams->d_MeshAdapter;
        mapMaster.reset(new Map3Dto1D (mapMasterParams));

        boost::shared_ptr<AMP::Database> mapTarget_db = myparams->d_db->getDatabase("Map1Dto3D");
        mapTargetParams.reset(new AMP::Operator::MapOperatorParameters(mapTarget_db));
        mapTargetParams->d_MapAdapter = myparams->d_MapAdapter;
        mapTarget.reset(new Map1Dto3D (mapTargetParams));

        mapMaster->setZLocations(mapTarget->getZLocations());

        std::string inpVar = mapMaster_db->getString("InputVariable"); 
        d_inpVariable.reset(new AMP::Mesh::NodalScalarVariable(inpVar,myparams->d_MeshAdapter));

        std::string outVar = mapTarget_db->getString("OutputVariable"); 
        d_outVariable.reset(new AMP::Mesh::NodalScalarVariable(outVar,myparams->d_MapAdapter));

        std::string gapVar = mapMaster_db->getString("OutputVariable"); 
        gapVariable.reset(new AMP::LinearAlgebra::Variable(gapVar));
        gap1DVec = AMP::LinearAlgebra::SimpleVector::create( mapTarget->getNumZlocations(), gapVariable );

        mapMaster->setVector(gap1DVec);

      }

      virtual ~MapSurface() { }

      virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
              AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

      boost::shared_ptr<AMP::LinearAlgebra::Vector> getBoundaryVector(const AMP::LinearAlgebra::Vector::shared_ptr &u) {
        return (u->subsetVectorForVariable(d_outVariable)) ;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
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

      }
*/
    protected :

      boost::shared_ptr<AMP::LinearAlgebra::Vector> gap1DVec; 
      boost::shared_ptr<AMP::LinearAlgebra::Variable> gapVariable;

      boost::shared_ptr<AMP::Mesh::NodalScalarVariable> d_inpVariable;
      boost::shared_ptr<AMP::Mesh::NodalScalarVariable> d_outVariable;

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
