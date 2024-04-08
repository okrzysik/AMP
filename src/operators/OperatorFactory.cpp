#include "AMP/operators/OperatorFactory.h"
#include "AMP/AMP_TPLs.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/CoupledOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/Operator.h"
#include "AMP/operators/boundary/ColumnBoundaryOperator.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"

#ifdef AMP_USE_LIBMESH
    #include "AMP/operators/boundary/libmesh/RobinMatrixCorrection.h"
    #include "AMP/operators/boundary/libmesh/RobinVectorCorrection.h"
    #include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
    #include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
    #include "AMP/operators/subchannel/FlowFrapconJacobian.h"
    #include "AMP/operators/subchannel/FlowFrapconOperator.h"
    #include "AMP/operators/subchannel/SubchannelFourEqLinearOperator.h"
    #include "AMP/operators/subchannel/SubchannelTwoEqLinearOperator.h"
#endif


namespace AMP::Operator {


// Macro to register an operator
#define REGISTER_OPERATOR( NAME )                                         \
    static struct NAME##_INIT {                                           \
        NAME##_INIT()                                                     \
        {                                                                 \
            auto fun = []( std::shared_ptr<OperatorParameters> params ) { \
                return std::make_unique<NAME>( params );                  \
            };                                                            \
            OperatorFactory::registerFactory( #NAME, fun );               \
        }                                                                 \
    } NAME##_INIT


// Create the operator
std::unique_ptr<Operator> OperatorFactory::create( std::shared_ptr<OperatorParameters> params )
{
    AMP_ASSERT( params );
    auto db = params->d_db;
    AMP_ASSERT( db );
    auto objectName = db->getString( "name" );
    return create( objectName, params );
}


// Register operators
REGISTER_OPERATOR( CoupledOperator );
REGISTER_OPERATOR( ColumnOperator );
REGISTER_OPERATOR( LinearBVPOperator );
REGISTER_OPERATOR( ColumnBoundaryOperator );
REGISTER_OPERATOR( DirichletMatrixCorrection );
REGISTER_OPERATOR( DirichletVectorCorrection );
#ifdef AMP_USE_LIBMESH
REGISTER_OPERATOR( DiffusionLinearFEOperator );
REGISTER_OPERATOR( MechanicsLinearFEOperator );
REGISTER_OPERATOR( RobinMatrixCorrection );
REGISTER_OPERATOR( RobinVectorCorrection );
REGISTER_OPERATOR( FlowFrapconJacobian );
REGISTER_OPERATOR( FlowFrapconOperator );
REGISTER_OPERATOR( SubchannelTwoEqLinearOperator );
REGISTER_OPERATOR( SubchannelFourEqLinearOperator );
#endif

} // namespace AMP::Operator
