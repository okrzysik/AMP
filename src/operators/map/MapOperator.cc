#include "MapOperator.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {

  void MapOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
  {
      boost::shared_ptr<MapOperatorParameters> myparams =
          boost::dynamic_pointer_cast<MapOperatorParameters>(params);

      AMP_INSIST( ((myparams.get()) != NULL), "NULL parameter" );
      AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );

      AMP_INSIST( (myparams->d_db)->keyExists("BoundaryId"), "Key ''tflow_id'' is missing!" );
      d_boundaryId = (myparams->d_db)->getInteger("BoundaryId");
  }

}
}

