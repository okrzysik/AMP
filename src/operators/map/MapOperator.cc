#include "MapOperator.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {


void MapOperator :: reset(const AMP::shared_ptr<OperatorParameters>& params)
{
    AMP::shared_ptr<MapOperatorParameters> myparams =
    AMP::dynamic_pointer_cast<MapOperatorParameters>(params);

    AMP_INSIST( myparams.get()!=NULL, "NULL parameter" );
    AMP_INSIST( myparams->d_db.get()!=NULL, "NULL database" );
    AMP_INSIST( !myparams->d_MapComm.isNull(), "NULL communicator" );

    AMP_INSIST( (myparams->d_db)->keyExists("BoundaryId"), "Key ''tflow_id'' is missing!" );
    d_boundaryId = (myparams->d_db)->getInteger("BoundaryId");

    d_MapComm = myparams->d_MapComm;
    d_MapMesh = myparams->d_MapMesh;
}


}
}

