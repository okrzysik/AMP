#ifndef included_AMP_MapOperator
#define included_AMP_MapOperator


#include "utils/InputDatabase.h"

#include "utils/shared_ptr.h"

#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "operators/map/MapOperatorParameters.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"

#include <string>

#ifdef DEBUG_CHECK_ASSERTIONS

#endif


namespace AMP {
namespace Operator {


/**
 * Class MapOperator is the base class for various mapping alogorithms. This
 * class stores a pointer to the mapAdapter to which the solution has to be
 * mapped from Operator's meshAdapter.
 */

class MapOperator : public Operator
{
public:
    /**
      Constructor calls the reset member which reads the information about the boundary.
      */
    explicit MapOperator( const AMP::shared_ptr<OperatorParameters> &params ) : Operator( params )
    {
        AMP::shared_ptr<MapOperatorParameters> myparams =
            AMP::dynamic_pointer_cast<MapOperatorParameters>( params );
        d_boundaryId = 0;
        reset( myparams );
    }

    /**
      Destructor
      */
    virtual ~MapOperator() {}

    virtual void reset( const AMP::shared_ptr<OperatorParameters> &params ) override;

    virtual AMP::LinearAlgebra::Variable::shared_ptr createInputVariable( const std::string &, int )
    {
        // Implemented in derived classes
        AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
        return emptyPointer;
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable( const std::string &,
                                                                           int )
    {
        // Implemented in derived classes
        AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
        return emptyPointer;
    }

    virtual void setInputVariableName( const std::string &, int )
    {
        // Implemented in derived classes.
    }

    virtual void setOutputVariableName( const std::string &, int )
    {
        // Implemented in derived classes.
    }

protected:
    unsigned int d_boundaryId;

    AMP::Mesh::Mesh::shared_ptr d_MapMesh;

    // Communicator for the Map
    AMP_MPI d_MapComm;

private:
};
} // namespace Operator
} // namespace AMP

#endif
