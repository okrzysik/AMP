#include "MapOperator.h"
#include "AMP/utils/Utilities.h"

namespace AMP {
namespace Operator {


void MapOperator::reset( const AMP::shared_ptr<OperatorParameters> &params )
{
    AMP::shared_ptr<MapOperatorParameters> myparams =
        AMP::dynamic_pointer_cast<MapOperatorParameters>( params );

    AMP_INSIST( myparams.get() != nullptr, "NULL parameter" );
    AMP_INSIST( myparams->d_db.get() != nullptr, "NULL database" );
    AMP_INSIST( !myparams->d_MapComm.isNull(), "NULL communicator" );

    AMP_INSIST( ( myparams->d_db )->keyExists( "BoundaryId" ), "Key ''tflow_id'' is missing!" );
    d_boundaryId = ( myparams->d_db )->getScalar<int>( "BoundaryId" );

    d_MapComm = myparams->d_MapComm;
    d_MapMesh = myparams->d_MapMesh;
}
} // namespace Operator
} // namespace AMP
