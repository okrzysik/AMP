#include "MapSurface.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Operator {


MapSurface::MapSurface(const boost::shared_ptr<OperatorParameters> & params):
    MapOperator (params)
{

    boost::shared_ptr<MapOperatorParameters> myparams = 
        boost::dynamic_pointer_cast<MapOperatorParameters>(params);

    // Construct for Map3Dto1D, Map1Dto3D
    boost::shared_ptr<AMP::Database> mapMaster_db = myparams->d_db->getDatabase("Map3Dto1D");
    mapMasterParams.reset(new AMP::Operator::MapOperatorParameters(mapMaster_db));
    mapMasterParams->d_Mesh = myparams->d_Mesh;
    mapMasterParams->d_MapMesh = myparams->d_Mesh;
    mapMasterParams->d_MapComm = myparams->d_Mesh->getComm();
    mapMaster.reset(new Map3Dto1D (mapMasterParams));

    boost::shared_ptr<AMP::Database> mapTarget_db = myparams->d_db->getDatabase("Map1Dto3D");
    mapTargetParams.reset(new AMP::Operator::MapOperatorParameters(mapTarget_db));
    mapTargetParams->d_Mesh = myparams->d_MapMesh;
    mapTargetParams->d_MapMesh = myparams->d_MapMesh;
    mapTargetParams->d_MapComm = myparams->d_MapMesh->getComm();
    mapTarget.reset(new Map1Dto3D (mapTargetParams));

    mapMaster->setZLocations(mapTarget->getZLocations());

    std::string inpVar = mapMaster_db->getString("InputVariable"); 
    d_inpVariable.reset(new AMP::LinearAlgebra::Variable(inpVar));

    std::string outVar = mapTarget_db->getString("OutputVariable"); 
    d_outVariable.reset(new AMP::LinearAlgebra::Variable(outVar));

    std::string gapVar = mapMaster_db->getString("OutputVariable"); 
    gapVariable.reset(new AMP::LinearAlgebra::Variable(gapVar));
    gap1DVec = AMP::LinearAlgebra::SimpleVector::create( mapTarget->getNumZlocations(), gapVariable );

    mapMaster->setVector(gap1DVec);

}



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


