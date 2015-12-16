#include "MapSurface.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Operator {


MapSurface::MapSurface( const AMP::shared_ptr<OperatorParameters> &params ) : MapOperator( params )
{

    AMP::shared_ptr<MapOperatorParameters> myparams =
        AMP::dynamic_pointer_cast<MapOperatorParameters>( params );

    // Construct for Map3Dto1D, Map1Dto3D
    AMP::shared_ptr<AMP::Database> mapMaster_db = myparams->d_db->getDatabase( "Map3Dto1D" );
    mapMasterParams.reset( new AMP::Operator::MapOperatorParameters( mapMaster_db ) );
    mapMasterParams->d_Mesh    = myparams->d_Mesh;
    mapMasterParams->d_MapMesh = myparams->d_Mesh;
    mapMasterParams->d_MapComm = myparams->d_Mesh->getComm();
    mapMaster.reset( new Map3Dto1D( mapMasterParams ) );

    AMP::shared_ptr<AMP::Database> mapTarget_db = myparams->d_db->getDatabase( "Map1Dto3D" );
    mapTargetParams.reset( new AMP::Operator::MapOperatorParameters( mapTarget_db ) );
    mapTargetParams->d_Mesh    = myparams->d_MapMesh;
    mapTargetParams->d_MapMesh = myparams->d_MapMesh;
    mapTargetParams->d_MapComm = myparams->d_MapMesh->getComm();
    mapTarget.reset( new Map1Dto3D( mapTargetParams ) );

    mapMaster->setZLocations( mapTarget->getZLocations() );

    std::string inpVar = mapMaster_db->getString( "InputVariable" );
    d_inpVariable.reset( new AMP::LinearAlgebra::Variable( inpVar ) );

    std::string outVar = mapTarget_db->getString( "OutputVariable" );
    d_outVariable.reset( new AMP::LinearAlgebra::Variable( outVar ) );

    std::string gapVar = mapMaster_db->getString( "OutputVariable" );
    gapVariable.reset( new AMP::LinearAlgebra::Variable( gapVar ) );
    gap1DVec = AMP::LinearAlgebra::SimpleVector<double>::create( mapTarget->getNumZlocations(),
                                                                 gapVariable );

    mapMaster->setVector( gap1DVec );
}


void MapSurface::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr )
{
    AMP_INSIST( ( ( u.get() ) != NULL ), "NULL Solution Vector" );

    inpVec = subsetInputVector( u );

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> nullVec;

    mapMaster->apply( inpVec, nullVec );
    mapTarget->apply( gap1DVec, nullVec );
}
}
}
