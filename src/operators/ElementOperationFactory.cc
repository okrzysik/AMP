
#include "ElementOperationFactory.h"
#include "utils/Utilities.h"

#ifdef USE_EXT_LIBMESH
#include "operators/diffusion/DiffusionElement.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionNonlinearElement.h"
#include "operators/libmesh/MassLinearElement.h"
#include "operators/libmesh/SourceNonlinearElement.h"
#include "operators/mechanics/MechanicsElement.h"
#include "operators/mechanics/MechanicsLinearElement.h"
#include "operators/mechanics/MechanicsLinearUpdatedLagrangianElement.h"
#include "operators/mechanics/MechanicsNonlinearElement.h"
#include "operators/mechanics/MechanicsNonlinearUpdatedLagrangianElement.h"
#endif


#define resetElementOperation( NAME )                                  \
    do {                                                               \
        if ( name == #NAME ) retElementOp.reset( new NAME( params ) ); \
    } while ( 0 )


namespace AMP {
namespace Operator {


AMP::shared_ptr<ElementOperation>
ElementOperationFactory::createElementOperation( AMP::shared_ptr<Database> elementOperationDb )
{
    AMP::shared_ptr<ElementOperation> retElementOp;
    AMP::shared_ptr<ElementOperationParameters> params;

    AMP_INSIST( elementOperationDb.get() != NULL,
                "ElementOperationFactory::createElementOperation:: NULL Database object input" );

    std::string name = elementOperationDb->getString( "name" );
    NULL_USE( name );

    params.reset( new ElementOperationParameters( elementOperationDb ) );

#ifdef USE_EXT_LIBMESH
    resetElementOperation( MechanicsLinearElement );
    resetElementOperation( MechanicsNonlinearElement );
    resetElementOperation( MechanicsLinearUpdatedLagrangianElement );
    resetElementOperation( MechanicsNonlinearUpdatedLagrangianElement );
    resetElementOperation( DiffusionLinearElement );
    resetElementOperation( DiffusionNonlinearElement );
    resetElementOperation( MassLinearElement );
    resetElementOperation( SourceNonlinearElement );
#endif

    return retElementOp;
}


} // namespace Operator
} // namespace AMP
