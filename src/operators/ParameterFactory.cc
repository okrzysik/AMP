#include "operators/ParameterFactory.h"

#ifdef USE_EXT_LIBMESH
#include "operators/NeutronicsRhsParameters.h"
#include "operators/boundary/DirichletMatrixCorrectionParameters.h"
#include "operators/mechanics/MechanicsLinearFEOperatorParameters.h"
#include "operators/mechanics/MechanicsNonlinearFEOperatorParameters.h"
#endif


#define resetParameters( NAME )                                      \
    do {                                                             \
        if ( name == #NAME )                                         \
            retParameters.reset( new NAME##Parameters( input_db ) ); \
    } while ( 0 )


namespace AMP {
namespace Operator {


AMP::shared_ptr<OperatorParameters>
ParameterFactory::createParameter( AMP::shared_ptr<AMP::Database> input_db,
                                   AMP::Mesh::Mesh::shared_ptr mesh )
{
    AMP::shared_ptr<OperatorParameters> retParameters;
    std::string name;

    AMP_INSIST( input_db.get() != nullptr,
                "ParameterFactory::createParameter:: NULL Database object input" );
    AMP_INSIST( input_db->keyExists( "name" ),
                "ParameterFactory::createParameter:: key 'name' must be a part of database " );

    name = input_db->getString( "name" );
    NULL_USE( name );

#ifdef USE_EXT_LIBMESH
    resetParameters( DirichletMatrixCorrection );
    resetParameters( MechanicsLinearFEOperator );
    resetParameters( MechanicsNonlinearFEOperator );
    resetParameters( NeutronicsRhs );
    AMP_ASSERT( retParameters != nullptr );
#endif
    retParameters->d_Mesh = mesh;

    return retParameters;
}


} // namespace Operator
} // namespace AMP
